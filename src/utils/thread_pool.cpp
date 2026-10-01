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

// Worker index of the current thread inside a parallel loop or region, -1
// outside. Pool threads keep their index for life; the submitting thread is
// the last worker. Calls nested inside a loop or region run serially.
thread_local int t_worker = -1;

bool runs_serially() {
    return t_worker >= 0 || thread_pool_size() <= 1;
}

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
// Sleeping barrier waiters are released with one FUTEX_WAKE instead of a
// condition-variable broadcast that makes them contend for one mutex.
void futex_wait(std::atomic<std::uint32_t>& word, std::uint32_t expected) {
    syscall(SYS_futex, reinterpret_cast<std::uint32_t*>(&word), FUTEX_WAIT_PRIVATE,
            expected, nullptr, nullptr, 0);
}

void futex_wake_all(std::atomic<std::uint32_t>& word) {
    syscall(SYS_futex, reinterpret_cast<std::uint32_t*>(&word), FUTEX_WAKE_PRIVATE,
            INT_MAX, nullptr, nullptr, 0);
}
#endif

// The first exception thrown in a parallel loop or region.
struct FirstError {
    std::atomic<bool> failed{false};
    std::mutex mutex;
    std::exception_ptr error;

    void record() {
        std::lock_guard<std::mutex> lk(mutex);
        if (!error) error = std::current_exception();
        failed.store(true, std::memory_order_release);
    }

    void rethrow() const {
        if (error) std::rethrow_exception(error);
    }
};

// n_workers - 1 background threads; the thread that submits a job takes part
// as worker n_workers - 1.
class ThreadPool {
public:
    explicit ThreadPool(int n_workers) : workers_(static_cast<size_t>(n_workers - 1)) {
        for (int i = 0; i < n_workers - 1; ++i) {
            workers_[i].thread = std::thread([this, i] { worker_loop(i); });
        }
    }

    ~ThreadPool() {
        stop_ = true;
        for (auto& w : workers_) wake(w);
        for (auto& w : workers_) w.thread.join();
    }

    int n_workers() const { return static_cast<int>(workers_.size()) + 1; }

    // Runs `body(worker)` on the calling thread and on `n_wake` background
    // workers, and returns when all have finished; `body` must not throw.
    // Waking only the workers there is work for avoids a thundering herd;
    // the rotation spreads small jobs over the workers.
    void run(int n_wake, const std::function<void(int)>& body) {
        body_ = &body;
        active_.store(n_wake, std::memory_order_relaxed);
        for (int k = 0; k < n_wake; ++k) wake(workers_[(cursor_ + k) % workers_.size()]);
        cursor_ = (cursor_ + n_wake) % workers_.size();

        t_worker = n_workers() - 1;
        body(t_worker);
        t_worker = -1;

        const auto done = [this] { return active_.load(std::memory_order_acquire) == 0; };
        if (!spin_until(done, spin_budget(kSpinUs))) {
            std::unique_lock<std::mutex> lk(done_mutex_);
            done_cv_.wait(lk, done);
        }
    }

private:
    struct alignas(64) Worker {
        std::thread thread;
        std::mutex mutex;
        std::condition_variable cv;
        std::atomic<std::uint64_t> gen{0};  // bumped once per job, and to stop
    };

    static void wake(Worker& w) {
        std::lock_guard<std::mutex> lk(w.mutex);
        ++w.gen;
        w.cv.notify_one();
    }

    void worker_loop(int id) {
        t_worker = id;
        Worker& w = workers_[id];
        for (std::uint64_t jobs = 0;; ++jobs) {
            const auto woken = [&] { return w.gen.load(std::memory_order_acquire) != jobs; };
            if (!spin_until(woken, spin_budget(kSpinUs))) {
                std::unique_lock<std::mutex> lk(w.mutex);
                w.cv.wait(lk, woken);
            }
            if (stop_) return;
            (*body_)(id);
            if (active_.fetch_sub(1, std::memory_order_acq_rel) == 1) {
                std::lock_guard<std::mutex> lk(done_mutex_);
                done_cv_.notify_one();
            }
        }
    }

    std::vector<Worker> workers_;
    size_t cursor_ = 0;
    // Set before waking workers, which read them after waking.
    const std::function<void(int)>* body_ = nullptr;
    bool stop_ = false;
    std::atomic<int> active_{0};  // woken workers still running `body_`
    std::mutex done_mutex_;
    std::condition_variable done_cv_;
};

std::mutex g_run_mutex;  // one job on the pool at a time; guards g_pool
std::unique_ptr<ThreadPool> g_pool;

int& configured_threads() {
    static int n = default_thread_count();
    return n;
}

ThreadPool& pool() {
    if (!g_pool) g_pool = std::make_unique<ThreadPool>(thread_pool_size());
    return *g_pool;
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

namespace detail {

struct RegionShared {
    int n_workers = 1;
    int spin_us = 0;  // barrier spin budget

    // One chunk counter for all dynamic loops of the region (see for_chunks).
    alignas(64) std::atomic<std::int64_t> next_chunk{0};

    // Centralized barrier: waiters spin briefly on the generation, then sleep
    // (futex on Linux, else a condition variable).
    alignas(64) std::atomic<int> arrived{0};
    alignas(64) std::atomic<std::uint32_t> generation{0};
    std::atomic<int> sleepers{0};
#ifndef __linux__
    std::mutex mutex;
    std::condition_variable cv;
#endif

    FirstError error;
};

} // namespace detail

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
    return configured_threads();
}

void set_thread_pool_size(int n_threads) {
    std::lock_guard<std::mutex> lk(g_run_mutex);
    configured_threads() = std::max(1, n_threads);
    g_pool.reset();  // recreated with the new size on first use
#ifdef EIGEN_GEMM_THREADPOOL
    configure_eigen_gemm(configured_threads());
#endif
}

void run_parallel_chunks(std::int64_t begin, std::int64_t end, std::int64_t chunk,
                         Scheduling sched,
                         const std::function<void(std::int64_t, std::int64_t, int)>& fn) {
    if (end <= begin) return;

    if (runs_serially()) {
        chunk = chunk_length(end - begin, chunk, sched, 1);
        const int worker = std::max(t_worker, 0);
        for (std::int64_t b = begin; b < end; b += chunk) fn(b, std::min(end, b + chunk), worker);
        return;
    }

    std::lock_guard<std::mutex> lk(g_run_mutex);
    ThreadPool& p = pool();
    chunk = chunk_length(end - begin, chunk, sched, p.n_workers());
    const std::int64_t n_chunks = (end - begin + chunk - 1) / chunk;
    std::atomic<std::int64_t> next{0};
    FirstError error;
    p.run(static_cast<int>(std::min<std::int64_t>(p.n_workers() - 1, n_chunks)), [&](int worker) {
        for (std::int64_t c; !error.failed.load(std::memory_order_acquire) &&
                             (c = next.fetch_add(1, std::memory_order_relaxed)) < n_chunks;) {
            const std::int64_t b = begin + c * chunk;
            try {
                fn(b, std::min(end, b + chunk), worker);
            } catch (...) {
                error.record();
            }
        }
    });
    error.rethrow();
}

// ---------------------------------------------------------------------------
// Persistent regions
// ---------------------------------------------------------------------------

bool ParallelRegion::cancelled() const {
    return shared_ != nullptr && shared_->error.failed.load(std::memory_order_acquire);
}

void ParallelRegion::run_guarded(const std::function<void()>& fn) {
    if (shared_ == nullptr) {
        fn();  // serial region: exceptions propagate directly
        return;
    }
    try {
        fn();
    } catch (...) {
        shared_->error.record();
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
        if (sh.sleepers.load(std::memory_order_seq_cst) > 0) futex_wake_all(sh.generation);
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
        futex_wait(sh.generation, gen);
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
    if (runs_serially()) {
        ParallelRegion region(nullptr, 0, 1);
        body(region);
        return;
    }

    std::lock_guard<std::mutex> lk(g_run_mutex);
    ThreadPool& p = pool();
    detail::RegionShared shared;
    shared.n_workers = p.n_workers();
    // No barrier spinning with more workers than physical cores: spinning
    // waiters take CPU time from the preempted participants they wait for.
    static const int physical_cores = default_thread_count();
    shared.spin_us = spin_budget(shared.n_workers > physical_cores ? 0 : kSpinUs);
    p.run(shared.n_workers - 1, [&](int worker) {
        ParallelRegion region(&shared, worker, shared.n_workers);
        try {
            body(region);
        } catch (...) {
            shared.error.record();
        }
    });
    shared.error.rethrow();
}

} // namespace baysor
