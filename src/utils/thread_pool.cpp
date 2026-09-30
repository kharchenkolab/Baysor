#include "baysor/utils/thread_pool.h"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <exception>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#ifdef EIGEN_GEMM_THREADPOOL
#include <Eigen/Core>
#include <Eigen/ThreadPool>
#endif

namespace baysor {

namespace {

// ---------------------------------------------------------------------------
// Thread-local state
// ---------------------------------------------------------------------------

// Worker index of the current thread: [0, n_workers) on pool workers, and
// n_workers - 1 on the calling thread (which participates as the last worker).
// Non-pool threads that submit work reuse index n_workers - 1 as well; only
// one parallel region runs at a time, so indices never collide.
thread_local int t_worker_index = 0;
thread_local int t_region_depth = 0;

// ---------------------------------------------------------------------------
// Job: one parallel region
// ---------------------------------------------------------------------------

struct Job {
    std::function<void(std::int64_t, std::int64_t, int)> fn;
    std::int64_t begin = 0;
    std::int64_t end = 0;
    std::int64_t chunk = 1;
    Scheduling sched = Scheduling::Dynamic;
    std::int64_t n_chunks = 0;
    std::vector<std::int64_t> static_starts;  // Static scheduling: start index per chunk
    std::vector<std::int64_t> static_lens;    // Static scheduling: length per chunk

    std::atomic<std::int64_t> next{0};
    std::atomic<int> active{0};               // participants not finished with this job
    std::atomic<bool> failed{false};
    std::exception_ptr error;
    std::mutex error_mutex;

    std::atomic<bool> done{false};
    std::mutex done_mutex;
    std::condition_variable done_cv;

    void record_error() {
        std::lock_guard<std::mutex> lk(error_mutex);
        if (!error) {
            error = std::current_exception();
        }
        failed.store(true, std::memory_order_release);
    }
};

// ---------------------------------------------------------------------------
// Global pool
// ---------------------------------------------------------------------------

class ThreadPool {
public:
    explicit ThreadPool(int n_workers) : n_workers_(std::max(1, n_workers)) {
        // n_workers_ threads total; the submitting thread participates as the
        // last worker, so only n_workers_ - 1 background threads are spawned.
        int bg = n_workers_ - 1;
        workers_.reserve(static_cast<size_t>(bg));
        for (int i = 0; i < bg; ++i) {
            workers_.emplace_back([this, i] { worker_loop(i); });
        }
    }

    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lk(mutex_);
            stop_ = true;
            ++seq_;
        }
        cv_.notify_all();
        for (auto& w : workers_) {
            if (w.joinable()) w.join();
        }
    }

    int n_workers() const { return n_workers_; }

    void run(Job& job) {
        job.active.store(n_workers_, std::memory_order_relaxed);
        {
            std::lock_guard<std::mutex> lk(mutex_);
            job_ = &job;
            ++seq_;
        }
        cv_.notify_all();

        // The submitting thread participates as the last worker. Remember its
        // previous worker index so nested serial runs inside this region use
        // the same scratch-buffer slot as the chunk bodies themselves.
        int my_id = n_workers_ - 1;
        int old_id = t_worker_index;
        t_worker_index = my_id;
        process_job(job, my_id);
        t_worker_index = old_id;

        finish_participation(job);
        wait_for_done(job);
        if (job.error) {
            std::rethrow_exception(job.error);
        }
    }

    static void finish_participation(Job& job) {
        if (job.active.fetch_sub(1, std::memory_order_acq_rel) == 1) {
            job.done.store(true, std::memory_order_release);
            std::lock_guard<std::mutex> lk(job.done_mutex);
            job.done_cv.notify_one();
        }
    }

private:
    void worker_loop(int id) {
        t_worker_index = id;
        std::uint64_t seen = 0;
        for (;;) {
            Job* job = nullptr;
            {
                std::unique_lock<std::mutex> lk(mutex_);
                cv_.wait(lk, [&] { return stop_ || seq_ != seen; });
                if (stop_) return;
                seen = seq_;
                job = job_;
            }
            process_job(*job, id);
            finish_participation(*job);
        }
    }

    // Claim and execute chunks until the job is exhausted or has failed.
    static void process_job(Job& job, int worker_id) {
        for (;;) {
            if (job.failed.load(std::memory_order_acquire)) break;
            std::int64_t idx = job.next.fetch_add(1, std::memory_order_relaxed);
            if (idx >= job.n_chunks) break;

            std::int64_t b, e;
            if (job.sched == Scheduling::Static) {
                b = job.static_starts[static_cast<size_t>(idx)];
                e = b + job.static_lens[static_cast<size_t>(idx)];
            } else {
                b = job.begin + idx * job.chunk;
                e = std::min(job.end, b + job.chunk);
            }

            ++t_region_depth;
            try {
                job.fn(b, e, worker_id);
            } catch (...) {
                job.record_error();
            }
            --t_region_depth;
        }
    }

    void wait_for_done(Job& job) {
        // Short bounded spin first (helps small regions), then block on the CV
        // so an idle pool burns no CPU.
        for (int s = 0; s < 4096; ++s) {
            if (job.done.load(std::memory_order_acquire)) return;
            std::atomic_signal_fence(std::memory_order_seq_cst);
        }
        std::unique_lock<std::mutex> lk(job.done_mutex);
        job.done_cv.wait(lk, [&] { return job.done.load(std::memory_order_acquire); });
    }

    int n_workers_;
    std::vector<std::thread> workers_;
    std::mutex mutex_;
    std::condition_variable cv_;
    Job* job_ = nullptr;
    std::uint64_t seq_ = 0;
    bool stop_ = false;
};

std::unique_ptr<ThreadPool>& pool_storage() {
    static std::unique_ptr<ThreadPool> pool;
    return pool;
}

std::mutex& region_mutex() {
    static std::mutex m;
    return m;
}

int& configured_threads() {
    static int n = 0;  // 0 = not configured yet
    return n;
}

int effective_threads() {
    int n = configured_threads();
    if (n <= 0) {
        n = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
    }
    return std::max(1, n);
}

#ifdef EIGEN_GEMM_THREADPOOL
// Eigen's own GEMM pool (Eigen >= 3.4.90). Sized to `--threads` minus one,
// mirroring the Baysor pool (background threads plus the calling thread).
// Baysor never runs Eigen GEMM inside a pool region and never runs a pool
// region inside a GEMM, so the two pools are never active at the same time
// and the total runnable thread count stays bounded by `--threads`.
//
// Eigen::setGemmThreadPool never destroys the previously registered pool (its
// parallelizer keeps dereferencing the registered pointer even when
// single-threaded), so pools are kept alive here and reused by size. A process
// that configures a handful of distinct thread counts parks at most one idle
// set of pool threads per distinct size; idle threads block on a condition
// variable and burn no CPU.
struct EigenPoolRegistry {
    std::vector<std::unique_ptr<Eigen::ThreadPool>> owned;
    Eigen::ThreadPool* registered = nullptr;
    int registered_threads = 0;
};

EigenPoolRegistry& eigen_pool_registry() {
    static EigenPoolRegistry reg;
    return reg;
}

void configure_eigen_gemm(int n_threads) {
    auto& reg = eigen_pool_registry();
    if (n_threads >= 3) {
        if (reg.registered != nullptr && reg.registered_threads == n_threads - 1) {
            return;  // already configured for this size
        }
        reg.owned.push_back(std::make_unique<Eigen::ThreadPool>(n_threads - 1));
        Eigen::ThreadPool* pool = reg.owned.back().get();
        Eigen::setGemmThreadPool(pool);
        reg.registered = pool;
        reg.registered_threads = n_threads - 1;
    } else {
        // 1-2 threads: run GEMM inline on the calling thread (bitwise identical
        // to the serial path). setNbThreads(1) makes the parallelizer take the
        // serial code path; the registered pool (if any) is left in place but
        // never scheduled on.
        reg.registered_threads = 0;
        Eigen::setNbThreads(1);
    }
}
#endif

} // namespace

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

int thread_pool_size() {
    return effective_threads();
}

void set_thread_pool_size(int n_threads) {
    if (n_threads <= 0) n_threads = 1;
    std::lock_guard<std::mutex> lk(region_mutex());
    configured_threads() = n_threads;
    // Rebuild the pool. Safe because no parallel region can be active while
    // holding region_mutex.
    pool_storage() = std::make_unique<ThreadPool>(n_threads);
#ifdef EIGEN_GEMM_THREADPOOL
    configure_eigen_gemm(n_threads);
#endif
}

bool inside_parallel_region() {
    return t_region_depth > 0;
}

int current_worker_index() {
    return t_worker_index;
}

void run_parallel_chunks(std::int64_t begin, std::int64_t end, std::int64_t chunk,
                         Scheduling sched,
                         const std::function<void(std::int64_t, std::int64_t, int)>& fn) {
    if (end <= begin) return;

    if (chunk <= 0) chunk = 1;
    std::int64_t n = end - begin;
    std::int64_t n_chunks = (sched == Scheduling::Static)
        ? 1
        : (n + chunk - 1) / chunk;

    // Serial fallback: nested regions and single-threaded runs execute inline
    // on the calling thread, in index order.
    if (t_region_depth > 0 || effective_threads() <= 1) {
        int worker = std::min(t_worker_index, effective_threads() - 1);
        if (sched == Scheduling::Static) {
            ++t_region_depth;
            try {
                fn(begin, end, worker);
            } catch (...) {
                --t_region_depth;
                throw;
            }
            --t_region_depth;
        } else {
            for (std::int64_t c = 0; c < n_chunks; ++c) {
                std::int64_t b = begin + c * chunk;
                std::int64_t e = std::min(end, b + chunk);
                ++t_region_depth;
                try {
                    fn(b, e, worker);
                } catch (...) {
                    --t_region_depth;
                    throw;
                }
                --t_region_depth;
            }
        }
        return;
    }

    // Parallel region.
    std::lock_guard<std::mutex> region_lk(region_mutex());

    ThreadPool* pool = pool_storage().get();
    if (pool == nullptr) {
        pool_storage() = std::make_unique<ThreadPool>(effective_threads());
        pool = pool_storage().get();
    }

    Job job;
    job.fn = fn;
    job.begin = begin;
    job.end = end;
    job.chunk = chunk;
    job.sched = sched;

    int n_workers = pool->n_workers();
    if (sched == Scheduling::Static) {
        // Contiguous blocks, one per worker, as evenly as possible (OpenMP
        // `schedule(static)` semantics).
        std::int64_t blocks = std::min<std::int64_t>(n_workers, n);
        std::int64_t base = n / blocks;
        std::int64_t rem = n % blocks;
        job.n_chunks = blocks;
        job.static_starts.reserve(static_cast<size_t>(blocks));
        job.static_lens.reserve(static_cast<size_t>(blocks));
        std::int64_t pos = begin;
        for (std::int64_t b = 0; b < blocks; ++b) {
            std::int64_t len = base + (b < rem ? 1 : 0);
            job.static_starts.push_back(pos);
            job.static_lens.push_back(len);
            pos += len;
        }
    } else {
        job.n_chunks = n_chunks;
    }

    pool->run(job);
}

} // namespace baysor
