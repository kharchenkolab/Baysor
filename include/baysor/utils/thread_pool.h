#pragma once

// Baysor's own persistent thread pool, replacing OpenMP for all internal
// parallel loops (and, via the subpar custom-parallel hooks, for the
// FetchContent dependencies that would otherwise spawn threads per call).
//
// Design notes:
//  - One global pool, sized once at start-up from `--threads` (see
//    `set_thread_pool_size`). Workers are persistent and block on a condition
//    variable when idle, so an idle pool burns no CPU.
//  - `parallel_for` supports static and dynamic chunking. Dynamic chunking
//    hands out fixed-size chunks through an atomic counter; chunk boundaries
//    therefore do not depend on the number of threads or on scheduling, which
//    is what makes the chunk-keyed RNG scheme in the E-step deterministic.
//  - `run_parallel_chunks` additionally reports the executing worker index so
//    bodies can use per-worker scratch buffers.
//  - A thread-local "inside a parallel region" flag makes nested parallel calls
//    run serially on the current thread (replacing `omp_in_parallel`).
//  - Exceptions thrown in worker tasks are rethrown on the calling thread.
//  - With 1 thread no worker threads exist at all: everything runs inline on
//    the calling thread in index order, preserving the exact serial code path
//    (single-threaded results stay bitwise identical to the OpenMP build with
//    OMP_NUM_THREADS=1).
//  - `parallel_reduce` accumulates into a fixed number of buckets that depends
//    on the data size but not on the thread count, and merges the buckets in
//    index order. This keeps floating-point reductions deterministic and
//    (for >1 thread) independent of the thread count. With 1 thread the whole
//    range is a single bucket, i.e. a plain sequential accumulation.

#include <algorithm>
#include <cstdint>
#include <functional>
#include <type_traits>
#include <vector>

namespace baysor {

enum class Scheduling {
    Static,   ///< contiguous blocks, one per worker (OpenMP `schedule(static)`)
    Dynamic   ///< fixed-size chunks handed out via an atomic counter (OpenMP `schedule(dynamic, chunk)`)
};

/// Number of worker threads of the global pool (>= 1). This is the value the
/// pool was configured with, even on nested/serial paths.
int thread_pool_size();

/// (Re)configure the global pool to `n_threads` workers. Values <= 1 configure
/// the serial fallback (no worker threads). Must be called before parallel
/// work; in practice once at start-up from `--threads`.
void set_thread_pool_size(int n_threads);

/// True while executing inside a parallel region (on a pool worker or in a
/// serial nested run). Nested parallel calls run serially.
bool inside_parallel_region();

/// Index in [0, thread_pool_size()) of the current worker thread, for indexing
/// per-worker scratch buffers. Returns 0 on non-pool threads (the calling
/// thread, which only ever runs when the pool is serial, i.e. with 1 thread).
int current_worker_index();

/// Core primitive: run `fn(chunk_begin, chunk_end, worker_index)` on disjoint
/// half-open index ranges covering [begin, end).
///  - Dynamic: ranges are chunks of `chunk` indices handed out via an atomic
///    counter; chunk boundaries are fixed given (begin, end, chunk).
///  - Static: `chunk` is ignored and [begin, end) is split into contiguous
///    blocks (one per worker) as evenly as possible.
/// With 1 thread (or when nested inside another parallel region) everything
/// runs inline on the calling thread in index order.
void run_parallel_chunks(std::int64_t begin, std::int64_t end, std::int64_t chunk,
                         Scheduling sched,
                         const std::function<void(std::int64_t, std::int64_t, int)>& fn);

namespace detail {

// Detects bodies taking (index) vs (index, worker_index).
template <class Fn, class = void>
struct takes_worker : std::false_type {};
template <class Fn>
struct takes_worker<Fn, std::void_t<decltype(std::declval<Fn&>()(std::int64_t(0), 0))>>
    : std::true_type {};

} // namespace detail

/// Per-index parallel loop over [begin, end). `fn` takes either (index) or
/// (index, worker_index). `chunk` is the dynamic chunk size (see
/// run_parallel_chunks).
template <class Fn>
void parallel_for(std::int64_t begin, std::int64_t end, std::int64_t chunk, Fn&& fn,
                  Scheduling sched = Scheduling::Dynamic) {
    auto body = [f = std::forward<Fn>(fn)](std::int64_t b, std::int64_t e, int worker) mutable {
        if constexpr (detail::takes_worker<Fn>::value) {
            for (std::int64_t i = b; i < e; ++i) f(i, worker);
        } else {
            for (std::int64_t i = b; i < e; ++i) f(i);
        }
    };
    run_parallel_chunks(begin, end, chunk, sched, body);
}

/// Per-index serial-or-parallel loop for `schedule(static)` sites.
template <class Fn>
void parallel_for_static(std::int64_t begin, std::int64_t end, Fn&& fn) {
    parallel_for(begin, end, /*chunk=*/0, std::forward<Fn>(fn), Scheduling::Static);
}

/// Deterministic reduction over [begin, end). The range is split into a fixed
/// number of buckets (depending on the range size and `bucket_size`, but never
/// on the thread count; with 1 thread there is exactly one bucket covering the
/// whole range). `fn(bucket_begin, bucket_end, T& acc)` accumulates into a
/// private bucket initialized to `init`, and the buckets are merged in index
/// order with `combine(acc, bucket)`. The final fold is
/// `out = init; out = combine(out, bucket_i)` in bucket order, so with one
/// bucket this is exactly the sequential accumulation.
template <class T, class Fn, class Combine>
T parallel_reduce(std::int64_t begin, std::int64_t end, std::int64_t bucket_size,
                  T init, Fn&& fn, Combine&& combine) {
    if (end <= begin) return init;

    constexpr std::int64_t kMaxBuckets = 64;
    std::int64_t n = end - begin;
    std::int64_t n_buckets;
    if (bucket_size <= 0 || thread_pool_size() <= 1) {
        n_buckets = 1;
    } else {
        n_buckets = (n + bucket_size - 1) / bucket_size;
        if (n_buckets > kMaxBuckets) n_buckets = kMaxBuckets;
        if (n_buckets < 1) n_buckets = 1;
    }
    std::int64_t stride = (n + n_buckets - 1) / n_buckets;

    std::vector<T> accs(static_cast<size_t>(n_buckets), init);
    run_parallel_chunks(begin, end, stride, Scheduling::Dynamic,
        [&](std::int64_t b, std::int64_t e, int) {
            std::int64_t bucket = (b - begin) / stride;
            fn(b, e, accs[static_cast<size_t>(bucket)]);
        });

    T out = init;
    for (auto& acc : accs) {
        out = combine(out, acc);
    }
    return out;
}

} // namespace baysor

// ---------------------------------------------------------------------------
// Custom-parallelization hooks for the FetchContent dependencies
// ---------------------------------------------------------------------------
// Defined as SUBPAR_CUSTOM_PARALLELIZE_RANGE / SUBPAR_CUSTOM_PARALLELIZE_SIMPLE
// / UMAPPP_CUSTOM_PARALLEL (see src/processing/data_processing/umap_wrappers.cpp,
// which must include this header before the dependency headers) so that
// umappp, knncolle, irlba and CppKmeans run their parallel loops on the Baysor
// pool instead of OpenMP or per-call std::thread spawns.

namespace baysor {

/// Drop-in replacement for `subpar::parallelize_range`. Keeps subpar's default
/// contiguous split (worker w gets `num_tasks / num_workers` tasks plus one
/// extra while w < `num_tasks % num_workers`), so per-worker data indexed by
/// `w` behaves exactly as with subpar's own scheme. `run_task_range(w, start,
/// length)` is invoked once per worker slot.
template <class Task_, class Run_>
void subpar_parallelize_range(int num_workers, Task_ num_tasks, Run_ run_task_range) {
    if (num_tasks <= Task_(0)) return;
    if (num_workers <= 1 || num_tasks == Task_(1)) {
        run_task_range(0, Task_(0), num_tasks);
        return;
    }
    std::int64_t W = std::min<std::int64_t>(num_workers, static_cast<std::int64_t>(num_tasks));
    std::int64_t tpu = static_cast<std::int64_t>(num_tasks) / W;
    std::int64_t rem = static_cast<std::int64_t>(num_tasks) % W;
    parallel_for(0, W, 1, [&](std::int64_t w) {
        Task_ start = static_cast<Task_>(w * tpu + std::min<std::int64_t>(w, rem));
        Task_ len = static_cast<Task_>(tpu + (w < rem ? 1 : 0));
        run_task_range(static_cast<int>(w), start, len);
    });
}

/// Drop-in replacement for `subpar::parallelize_simple` (1:1 task/worker).
template <class Task_, class Run_>
void subpar_parallelize_simple(Task_ num_tasks, Run_ run_task) {
    if (num_tasks <= Task_(0)) return;
    parallel_for(0, static_cast<std::int64_t>(num_tasks), 1, [&](std::int64_t w) {
        run_task(static_cast<Task_>(w));
    });
}

/// Drop-in replacement for umappp's `UMAPPP_CUSTOM_PARALLEL`: split
/// [0, num_tasks) into contiguous ranges and call `run(first, last)` on each.
template <class Task_, class Run_>
void umappp_parallel_range(Task_ num_tasks, Run_ run, int num_threads) {
    if (num_tasks <= Task_(0)) return;
    std::int64_t W = std::min<std::int64_t>(
        std::max(1, num_threads), static_cast<std::int64_t>(num_tasks));
    std::int64_t tpu = static_cast<std::int64_t>(num_tasks) / W;
    std::int64_t rem = static_cast<std::int64_t>(num_tasks) % W;
    parallel_for(0, W, 1, [&](std::int64_t w) {
        std::int64_t start = w * tpu + std::min<std::int64_t>(w, rem);
        std::int64_t len = tpu + (w < rem ? 1 : 0);
        run(static_cast<Task_>(start), static_cast<Task_>(start + len));
    });
}

} // namespace baysor
