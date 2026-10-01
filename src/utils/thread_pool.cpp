#include "baysor/utils/thread_pool.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#ifdef __linux__
#include <climits>
#include <linux/futex.h>
#include <sys/syscall.h>
#include <unistd.h>
#endif

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

// Chunk length: `chunk` (at least 1), or one chunk per worker for Static.
std::int64_t chunk_length(std::int64_t n, std::int64_t chunk, Scheduling sched, int n_workers) {
    return std::max<std::int64_t>(sched == Scheduling::Static ? (n + n_workers - 1) / n_workers : chunk, 1);
}

// Spin budget before blocking at job hand-off, job completion and region
// barriers: it absorbs the wake-up latency of back-to-back parallel loops.
constexpr int kSpinUs = 20;

// `us`, unless BAYSOR_POOL_SPIN_US overrides it: the profiling suite sets 0
// under Valgrind, where the spin (bounded by a wall-clock deadline) would
// make instruction counts depend on host timing.
int spin_budget(int us) {
    static const int env_us = [] {
        const char* env = std::getenv("BAYSOR_POOL_SPIN_US");
        try {
            return env != nullptr ? std::max(std::stoi(env), -1) : -1;
        } catch (const std::logic_error&) {  // not a number, or out of range
            return -1;
        }
    }();
    return env_us >= 0 ? env_us : us;
}

// Busy-waits up to `us` microseconds for `ready()`; returns its last value.
template <class Pred>
bool spin_until(Pred ready, int us) {
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::microseconds(us);
    while (!ready()) {
        if (std::chrono::steady_clock::now() >= deadline) return false;
#if defined(__x86_64__) || defined(__i386__)
        __builtin_ia32_pause();
#elif defined(__aarch64__) || defined(__arm__)
        __asm__ __volatile__("yield" ::: "memory");
#else
        std::this_thread::yield();
#endif
    }
    return true;
}

#ifdef __linux__
// Futex wait/wake on a 32-bit atomic (private to this process). Sleeping
// barrier waiters are released with one FUTEX_WAKE instead of a
// condition-variable broadcast that makes them contend for one mutex.
void futex_wait(std::atomic<std::uint32_t>* word, std::uint32_t expected) {
    syscall(SYS_futex, reinterpret_cast<std::uint32_t*>(word), FUTEX_WAIT_PRIVATE,
            expected, nullptr, nullptr, 0);
}

void futex_wake_all(std::atomic<std::uint32_t>* word) {
    syscall(SYS_futex, reinterpret_cast<std::uint32_t*>(word), FUTEX_WAKE_PRIVATE,
            INT_MAX, nullptr, nullptr, 0);
}
#endif

// ---------------------------------------------------------------------------
// Job: one parallel region
// ---------------------------------------------------------------------------

struct Job {
    std::function<void(std::int64_t, std::int64_t, int)> fn;
    std::int64_t begin = 0;
    std::int64_t end = 0;
    std::int64_t chunk = 1;
    std::int64_t n_chunks = 0;
    // Persistent region: every participant runs `(*region)(worker_index)`
    // exactly once instead of claiming chunks.
    const std::function<void(int)>* region = nullptr;

    std::atomic<std::int64_t> next{0};
    std::atomic<int> active{0};               // participants not finished with this job
    std::atomic<bool> failed{false};
    std::exception_ptr error;
    std::mutex error_mutex;

    std::atomic<bool> done{false};

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
            workers_.push_back(std::make_unique<Worker>());
        }
        for (int i = 0; i < bg; ++i) {
            workers_[static_cast<size_t>(i)]->thread =
                std::thread([this, i] { worker_loop(*workers_[static_cast<size_t>(i)], i); });
        }
    }

    ~ThreadPool() {
        for (auto& w : workers_) {
            {
                std::lock_guard<std::mutex> lk(w->mutex);
                w->stop.store(true, std::memory_order_release);
                ++w->gen;
            }
            w->cv.notify_one();
        }
        for (auto& w : workers_) {
            if (w->thread.joinable()) w->thread.join();
        }
    }

    int n_workers() const { return n_workers_; }

    void run(Job& job) {
        // Wake only as many workers as there is work for (the submitting
        // thread participates as the last worker), and only the workers we
        // wake take part in the job: no thundering herd, and tiny regions run
        // on the calling thread alone.
        int bg = n_workers_ - 1;
        int wake = (job.region != nullptr)
            ? bg
            : static_cast<int>(std::min<std::int64_t>(bg, job.n_chunks));
        job.active.store(wake + 1, std::memory_order_relaxed);

        for (int k = 0; k < wake; ++k) {
            Worker& w = *workers_[static_cast<size_t>((wake_cursor_ + k) % workers_.size())];
            std::lock_guard<std::mutex> lk(w.mutex);
            w.job = &job;
            ++w.gen;
            w.cv.notify_one();
        }
        if (!workers_.empty()) {
            wake_cursor_ = (wake_cursor_ + wake) % static_cast<int>(workers_.size());
        }

        // The submitting thread participates as the last worker. Remember its
        // previous worker index so nested serial runs inside this region use
        // the same scratch-buffer slot as the chunk bodies themselves.
        int my_id = n_workers_ - 1;
        int old_id = t_worker_index;
        t_worker_index = my_id;
        process_job(job, my_id);
        t_worker_index = old_id;

        finish_participation(*this, job);
        wait_for_done(job);
        if (job.error) {
            std::rethrow_exception(job.error);
        }
    }

private:
    struct Worker {
        std::thread thread;
        std::mutex mutex;
        std::condition_variable cv;
        std::atomic<std::uint64_t> gen{0};  // bumped when a job (or stop) is assigned
        Job* job = nullptr;
        std::atomic<bool> stop{false};
    };

    static void finish_participation(ThreadPool& pool, Job& job) {
        if (job.active.fetch_sub(1, std::memory_order_acq_rel) == 1) {
            // The completion signal lives on the pool, not on the job: the
            // submitter may destroy its stack-local Job as soon as `done` is
            // set, so the last participant must not touch the job afterwards
            // (not even its mutex/condition_variable).
            job.done.store(true, std::memory_order_release);
            std::lock_guard<std::mutex> lk(pool.done_mutex_);
            pool.done_cv_.notify_one();
        }
    }

    void worker_loop(Worker& w, int id) {
        t_worker_index = id;
        // Start from 0, never from the current generation: the pool may be
        // destroyed (or a job assigned) before this thread first runs, and
        // initializing `seen` from the already-bumped generation would lose
        // the stop/job signal.
        std::uint64_t seen = 0;
        for (;;) {
            std::uint64_t g = w.gen.load(std::memory_order_acquire);
            if (g != seen) {
                seen = g;
                if (w.stop.load(std::memory_order_acquire)) return;
                Job* job = w.job;
                process_job(*job, id);
                finish_participation(*this, *job);
                continue;
            }

            const auto woken = [&] { return w.gen.load(std::memory_order_acquire) != seen; };
            if (spin_until(woken, spin_budget(kSpinUs))) continue;
            std::unique_lock<std::mutex> lk(w.mutex);
            w.cv.wait(lk, woken);
        }
    }

    // Claim and execute chunks until the job is exhausted or has failed.
    static void process_job(Job& job, int worker_id) {
        if (job.region != nullptr) {
            ++t_region_depth;
            try {
                (*job.region)(worker_id);
            } catch (...) {
                job.record_error();
            }
            --t_region_depth;
            return;
        }
        for (;;) {
            if (job.failed.load(std::memory_order_acquire)) break;
            std::int64_t idx = job.next.fetch_add(1, std::memory_order_relaxed);
            if (idx >= job.n_chunks) break;

            const std::int64_t b = job.begin + idx * job.chunk;
            const std::int64_t e = std::min(job.end, b + job.chunk);
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
        const auto done = [&] { return job.done.load(std::memory_order_acquire); };
        if (spin_until(done, spin_budget(kSpinUs))) return;
        std::unique_lock<std::mutex> lk(done_mutex_);
        done_cv_.wait(lk, done);
    }

    int n_workers_;
    std::vector<std::unique_ptr<Worker>> workers_;
    int wake_cursor_ = 0;
    std::mutex done_mutex_;
    std::condition_variable done_cv_;
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
        n = default_thread_count();
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
// Persistent regions: shared state
// ---------------------------------------------------------------------------

namespace detail {

struct RegionShared {
    int n_workers = 1;
    int spin_us = 0;  // barrier spin budget

    // Dynamic loops: one monotonically increasing chunk counter for the whole
    // region. Every participant claims chunks of a loop until it overshoots,
    // which it does exactly once per loop, so after a loop with C chunks the
    // counter has advanced by C + n_workers on every participant's view.
    alignas(64) std::atomic<std::int64_t> next_chunk{0};

    // Centralized barrier: arrival count plus generation number; waiters
    // spin briefly on the generation, then sleep (futex on Linux, else a
    // condition variable).
    alignas(64) std::atomic<int> arrived{0};
    alignas(64) std::atomic<std::uint32_t> generation{0};
    std::atomic<int> sleepers{0};
#ifndef __linux__
    std::mutex mutex;
    std::condition_variable cv;
#endif

    std::atomic<bool> failed{false};
    std::mutex error_mutex;
    std::exception_ptr error;

    void record_error() {
        std::lock_guard<std::mutex> lk(error_mutex);
        if (!error) {
            error = std::current_exception();
        }
        failed.store(true, std::memory_order_release);
    }
};

} // namespace detail

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

int default_thread_count() {
#ifdef __linux__
    // Prefer physical cores over logical CPUs (hardware_concurrency counts
    // SMT siblings; on this workload Hyper-Threading adds little and doubles
    // wake-up costs). Count unique (package, core) pairs in sysfs.
    try {
        std::set<std::string> cores;
        for (int cpu = 0; cpu < 4096; ++cpu) {
            std::ifstream pkg("/sys/devices/system/cpu/cpu" + std::to_string(cpu) +
                              "/topology/physical_package_id");
            if (!pkg.is_open()) break;
            std::string package_id, core_id;
            std::getline(pkg, package_id);
            std::ifstream core("/sys/devices/system/cpu/cpu" + std::to_string(cpu) +
                               "/topology/core_id");
            if (!std::getline(core, core_id)) break;
            cores.insert(package_id + ":" + core_id);
        }
        if (!cores.empty()) {
            return static_cast<int>(cores.size());
        }
    } catch (...) {
        // fall through to hardware_concurrency
    }
#endif
    return static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
}

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

    // Serial fallback: nested regions and single-threaded runs execute inline
    // on the calling thread, in index order.
    if (t_region_depth > 0 || effective_threads() <= 1) {
        chunk = chunk_length(end - begin, chunk, sched, 1);
        int worker = std::min(t_worker_index, effective_threads() - 1);
        for (std::int64_t b = begin; b < end; b += chunk) {
            ++t_region_depth;
            try {
                fn(b, std::min(end, b + chunk), worker);
            } catch (...) {
                --t_region_depth;
                throw;
            }
            --t_region_depth;
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
    job.chunk = chunk_length(end - begin, chunk, sched, pool->n_workers());
    job.n_chunks = (end - begin + job.chunk - 1) / job.chunk;
    pool->run(job);
}

// ---------------------------------------------------------------------------
// Persistent regions
// ---------------------------------------------------------------------------

bool ParallelRegion::cancelled() const {
    return shared_ != nullptr && shared_->failed.load(std::memory_order_acquire);
}

void ParallelRegion::run_guarded(const std::function<void()>& fn) {
    if (shared_ == nullptr) {
        fn();  // serial region: exceptions propagate directly
        return;
    }
    try {
        fn();
    } catch (...) {
        shared_->record_error();
    }
}

void ParallelRegion::for_chunks(std::int64_t begin, std::int64_t end, std::int64_t chunk,
                                Scheduling sched,
                                const std::function<void(std::int64_t, std::int64_t, int)>& fn) {
    chunk = chunk_length(end - begin, chunk, sched, n_workers_);

    if (shared_ == nullptr) {
        // Serial region: the run_parallel_chunks serial path.
        for (std::int64_t b = begin; b < end; b += chunk) {
            fn(b, std::min(end, b + chunk), worker_);
        }
        return;
    }

    if (end > begin) {
        const std::int64_t n_chunks = (end - begin + chunk - 1) / chunk;
        for (;;) {
            const std::int64_t idx =
                shared_->next_chunk.fetch_add(1, std::memory_order_relaxed) - chunk_base_;
            if (idx >= n_chunks) break;
            if (cancelled()) continue;  // keep claiming so the counter stays consistent
            const std::int64_t b = begin + idx * chunk;
            const std::int64_t e = std::min(end, b + chunk);
            run_guarded([&]() { fn(b, e, worker_); });
        }
        chunk_base_ += n_chunks + n_workers_;
    }
    barrier();
}

void ParallelRegion::barrier() {
    if (shared_ == nullptr) return;
    detail::RegionShared& sh = *shared_;

    const std::uint32_t gen = sh.generation.load(std::memory_order_acquire);
    if (sh.arrived.fetch_add(1, std::memory_order_acq_rel) == n_workers_ - 1) {
        // Last to arrive: reset the count (nobody can arrive at the next
        // barrier before observing the new generation) and release everyone.
        // seq_cst on generation/sleepers: a waiter that registers as a
        // sleeper either sees the new generation or is seen here.
        sh.arrived.store(0, std::memory_order_relaxed);
#ifdef __linux__
        sh.generation.store(gen + 1, std::memory_order_seq_cst);
        if (sh.sleepers.load(std::memory_order_seq_cst) > 0) futex_wake_all(&sh.generation);
#else
        bool wake;
        {
            std::lock_guard<std::mutex> lk(sh.mutex);
            sh.generation.store(gen + 1, std::memory_order_seq_cst);
            wake = sh.sleepers.load(std::memory_order_relaxed) > 0;
        }
        if (wake) sh.cv.notify_all();
#endif
        return;
    }

    const auto released = [&] { return sh.generation.load(std::memory_order_acquire) != gen; };
    if (spin_until(released, sh.spin_us)) return;
#ifdef __linux__
    sh.sleepers.fetch_add(1, std::memory_order_seq_cst);
    while (sh.generation.load(std::memory_order_seq_cst) == gen) {
        futex_wait(&sh.generation, gen);
    }
    sh.sleepers.fetch_sub(1, std::memory_order_relaxed);
#else
    std::unique_lock<std::mutex> lk(sh.mutex);
    sh.sleepers.fetch_add(1, std::memory_order_relaxed);
    sh.cv.wait(lk, released);
    sh.sleepers.fetch_sub(1, std::memory_order_relaxed);
#endif
}

void parallel_region(const std::function<void(ParallelRegion&)>& body) {
    if (t_region_depth > 0 || effective_threads() <= 1) {
        ParallelRegion region(nullptr, 0, 1);
        ++t_region_depth;
        try {
            body(region);
        } catch (...) {
            --t_region_depth;
            throw;
        }
        --t_region_depth;
        return;
    }

    std::lock_guard<std::mutex> region_lk(region_mutex());

    ThreadPool* pool = pool_storage().get();
    if (pool == nullptr) {
        pool_storage() = std::make_unique<ThreadPool>(effective_threads());
        pool = pool_storage().get();
    }

    detail::RegionShared shared;
    shared.n_workers = pool->n_workers();
    // No barrier spinning with more workers than physical cores: spinning
    // waiters take CPU time from the preempted participants they wait for.
    static const int physical_cores = default_thread_count();
    shared.spin_us = spin_budget(shared.n_workers > physical_cores ? 0 : kSpinUs);
    const std::function<void(int)> participant = [&](int worker) {
        ParallelRegion region(&shared, worker, shared.n_workers);
        body(region);
    };

    Job job;
    job.region = &participant;
    job.n_chunks = shared.n_workers;
    pool->run(job);  // rethrows an exception escaping a participant's body

    if (shared.error) {
        std::rethrow_exception(shared.error);
    }
}

} // namespace baysor
