#pragma once

// Baysor's persistent thread pool. It runs all internal parallel loops and,
// through the subpar/umappp hooks at the end of this file, those of the
// FetchContent dependencies.
//  - One global pool sized by `set_thread_pool_size` (`--threads`); idle
//    workers block, so an idle pool burns no CPU.
//  - With 1 thread, and for calls nested in a parallel loop or region,
//    everything runs inline on the calling thread in index order.
//  - Dynamic chunk boundaries depend only on (begin, end, chunk), never on the
//    thread count or scheduling: the chunk-keyed E-step RNG relies on this.
//  - Exceptions thrown by loop bodies are rethrown on the calling thread.

#include <algorithm>
#include <cstdint>
#include <functional>
#include <type_traits>
#include <vector>

namespace baysor {

enum class Scheduling {
    Static,   ///< one contiguous chunk per worker; `chunk` is ignored
    Dynamic   ///< chunks of `chunk` indices handed out in index order
};

/// Number of threads of the global pool (>= 1).
int thread_pool_size();

/// Number of physical CPU cores when detectable (SMT siblings add little here
/// and double wake-up costs), otherwise std::thread::hardware_concurrency().
int default_thread_count();

/// Sizes the global pool; values < 1 mean 1. Must not run concurrently with
/// parallel work; in practice called once at start-up from `--threads`.
void set_thread_pool_size(int n_threads);

/// Runs `fn(chunk_begin, chunk_end, worker)` on disjoint chunks covering
/// [begin, end). `worker` is in [0, thread_pool_size()) and indexes per-worker
/// scratch buffers; nested calls keep the worker of the enclosing chunk.
void run_parallel_chunks(std::int64_t begin, std::int64_t end, std::int64_t chunk,
                         Scheduling sched,
                         const std::function<void(std::int64_t, std::int64_t, int)>& fn);

namespace detail {

// Chunk body running a per-index `fn(index)` or `fn(index, worker)`.
template <class Fn>
auto per_index(Fn& fn) {
    return [&fn](std::int64_t b, std::int64_t e, int worker) {
        for (std::int64_t i = b; i < e; ++i) {
            if constexpr (std::is_invocable_v<Fn&, std::int64_t, int>) {
                fn(i, worker);
            } else {
                fn(i);
            }
        }
    };
}

} // namespace detail

/// Per-index loop over [begin, end) in dynamic chunks of `chunk` indices;
/// `fn` takes (index) or (index, worker).
template <class Fn>
void parallel_for(std::int64_t begin, std::int64_t end, std::int64_t chunk, Fn&& fn) {
    run_parallel_chunks(begin, end, chunk, Scheduling::Dynamic, detail::per_index(fn));
}

/// parallel_for with one contiguous chunk per worker.
template <class Fn>
void parallel_for_static(std::int64_t begin, std::int64_t end, Fn&& fn) {
    run_parallel_chunks(begin, end, 0, Scheduling::Static, detail::per_index(fn));
}

/// Deterministic reduction over [begin, end). `fn(chunk_begin, chunk_end, acc)`
/// accumulates into one of at most 64 buckets initialized to `init`, which are
/// then folded in index order: `out = combine(out, bucket)`, starting from
/// `out = init`. The buckets depend on the range size and `bucket_size` but
/// never on the thread count; with 1 thread there is a single bucket.
template <class T, class Fn, class Combine>
T parallel_reduce(std::int64_t begin, std::int64_t end, std::int64_t bucket_size,
                  T init, Fn&& fn, Combine&& combine) {
    if (end <= begin) return init;
    const std::int64_t n = end - begin;
    const std::int64_t n_buckets = (bucket_size <= 0 || thread_pool_size() <= 1)
        ? 1 : std::min<std::int64_t>((n + bucket_size - 1) / bucket_size, 64);
    const std::int64_t stride = (n + n_buckets - 1) / n_buckets;

    std::vector<T> accs(static_cast<size_t>(n_buckets), init);
    run_parallel_chunks(begin, end, stride, Scheduling::Dynamic,
        [&](std::int64_t b, std::int64_t e, int) {
            fn(b, e, accs[static_cast<size_t>((b - begin) / stride)]);
        });
    for (auto& acc : accs) init = combine(init, acc);
    return init;
}

// ---------------------------------------------------------------------------
// Persistent parallel regions
// ---------------------------------------------------------------------------
// `parallel_region(body)` wakes the pool once and runs `body(region)` on every
// worker (OpenMP `#pragma omp parallel`), so a sequence of work-shared loops
// with short serial glue costs one wake-up instead of one per loop.
//  - Every participant must call the work-sharing constructs (`for_chunks`,
//    `for_each`, `barrier`, `single`) in the same order with the same
//    arguments; control flow that depends on shared data must read it only
//    after a barrier that orders its last write.
//  - An exception in a loop body or `single` block cancels the region: the
//    remaining constructs are skipped, `cancelled()` is true on every
//    participant after the next barrier, and the first exception is rethrown
//    when the region ends.
//  - With 1 thread, or nested in a parallel loop or region, the body runs once
//    inline with `n_workers() == 1`.

namespace detail { struct RegionShared; }

class ParallelRegion {
public:
    ParallelRegion(detail::RegionShared* shared, int worker, int n_workers)
        : shared_(shared), worker_(worker), n_workers_(n_workers) {}
    ParallelRegion(const ParallelRegion&) = delete;
    ParallelRegion& operator=(const ParallelRegion&) = delete;

    /// Index of this participant in [0, n_workers()). The master, which runs
    /// `single` blocks, is the calling thread: n_workers() - 1.
    int worker_index() const { return worker_; }
    int n_workers() const { return n_workers_; }
    bool is_master() const { return worker_ == n_workers_ - 1; }

    /// Work-shared loop with the chunks of run_parallel_chunks, then a barrier.
    /// `worker` is worker_index().
    void for_chunks(std::int64_t begin, std::int64_t end, std::int64_t chunk,
                    Scheduling sched,
                    const std::function<void(std::int64_t, std::int64_t, int)>& fn);

    /// Per-index for_chunks in dynamic chunks; `fn` takes (index) or (index, worker).
    template <class Fn>
    void for_each(std::int64_t begin, std::int64_t end, std::int64_t chunk, Fn&& fn) {
        for_chunks(begin, end, chunk, Scheduling::Dynamic, detail::per_index(fn));
    }

    /// Waits until every participant has arrived.
    void barrier();

    /// Runs `fn()` on the master only, then a barrier.
    template <class Fn>
    void single(Fn&& fn) {
        if (is_master() && !cancelled()) {
            try {
                fn();
            } catch (...) {
                record_error();
            }
        }
        barrier();
    }

    /// True once a construct has failed.
    bool cancelled() const;

private:
    void record_error();

    detail::RegionShared* shared_;
    int worker_;
    int n_workers_;
    std::int64_t chunk_base_ = 0;  // this participant's view of the shared chunk counter
};

/// Runs `body` once on every pool worker (see above) and returns when all
/// participants have finished; rethrows the first exception.
void parallel_region(const std::function<void(ParallelRegion&)>& body);

// ---------------------------------------------------------------------------
// Custom-parallelization hooks of the FetchContent dependencies
// ---------------------------------------------------------------------------
// src/processing/data_processing/umap_wrappers.cpp defines
// SUBPAR_CUSTOM_PARALLELIZE_RANGE / _SIMPLE and UMAPPP_CUSTOM_PARALLEL as these
// functions, so umappp, knncolle, irlba and CppKmeans run on the Baysor pool.

/// `subpar::parallelize_range` with subpar's split: worker w of
/// W = min(num_workers, num_tasks) gets num_tasks / W tasks, plus one while
/// w < num_tasks % W, so per-worker data indexed by `w` behaves as in subpar.
template <class Task_, class Run_>
void subpar_parallelize_range(int num_workers, Task_ num_tasks, Run_ run_task_range) {
    if (num_tasks <= Task_(0)) return;
    if (num_workers <= 1 || num_tasks == Task_(1)) {
        run_task_range(0, Task_(0), num_tasks);
        return;
    }
    const auto n = static_cast<std::int64_t>(num_tasks);
    const std::int64_t n_workers = std::min<std::int64_t>(num_workers, n);
    const std::int64_t per_worker = n / n_workers;
    const std::int64_t rem = n % n_workers;
    parallel_for(0, n_workers, 1, [&](std::int64_t w) {
        run_task_range(static_cast<int>(w), static_cast<Task_>(w * per_worker + std::min(w, rem)),
                       static_cast<Task_>(per_worker + (w < rem ? 1 : 0)));
    });
}

/// `subpar::parallelize_simple`: one task per worker.
template <class Task_, class Run_>
void subpar_parallelize_simple(Task_ num_tasks, Run_ run_task) {
    parallel_for(0, static_cast<std::int64_t>(num_tasks), 1, [&](std::int64_t t) {
        run_task(static_cast<Task_>(t));
    });
}

/// umappp's `UMAPPP_CUSTOM_PARALLEL`: `run(first, last)` on contiguous ranges.
template <class Task_, class Run_>
void umappp_parallel_range(Task_ num_tasks, Run_ run, int num_threads) {
    subpar_parallelize_range(num_threads, num_tasks,
                             [&](int, Task_ start, Task_ len) { run(start, start + len); });
}

} // namespace baysor
