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
//  - Parallel calls nested inside a parallel loop or region run serially on
//    the current thread.
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
    Static,   ///< one contiguous chunk per worker; `chunk` is ignored
    Dynamic   ///< chunks of `chunk` indices handed out in index order
};

/// Number of worker threads of the global pool (>= 1). This is the value the
/// pool was configured with, even on nested/serial paths.
int thread_pool_size();

/// Default thread count when nothing is configured explicitly: the number of
/// physical CPU cores when detectable (SMT siblings add little here and double
/// wake-up costs), otherwise std::thread::hardware_concurrency().
int default_thread_count();

/// (Re)configure the global pool to `n_threads` workers. Values <= 1 configure
/// the serial fallback (no worker threads). Must be called before parallel
/// work; in practice once at start-up from `--threads`.
void set_thread_pool_size(int n_threads);

/// Core primitive: run `fn(chunk_begin, chunk_end, worker_index)` on disjoint
/// half-open index ranges covering [begin, end). `worker_index` is in
/// [0, thread_pool_size()) and indexes per-worker scratch buffers; nested
/// calls keep the worker of the enclosing chunk.
///  - Dynamic: ranges are chunks of `chunk` indices handed out via an atomic
///    counter; chunk boundaries are fixed given (begin, end, chunk).
///  - Static: `chunk` is ignored and [begin, end) is split into one
///    contiguous chunk per worker.
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

// ---------------------------------------------------------------------------
// Persistent parallel regions
// ---------------------------------------------------------------------------
// `parallel_region(body)` wakes the pool once and runs `body(region)` on every
// worker (OpenMP `#pragma omp parallel`). Inside, the participants execute
// work-shared loops, barriers and single blocks, so a sequence of parallel
// phases with short serial glue costs one wake-up instead of one per loop.
//
// Rules (as in OpenMP):
//  - Every participant must call the work-sharing constructs (`for_chunks`,
//    `for_each`, `barrier`, `single`) in the same order with the same
//    arguments; control flow that depends on shared data must read it only
//    after a barrier that orders its last write.
//  - Exceptions thrown inside `for_chunks` bodies and `single` blocks are
//    recorded, the remaining constructs are skipped (`cancelled()` becomes
//    true after the next barrier), and the first exception is rethrown on the
//    calling thread when the region ends. The body itself must not throw
//    outside these constructs.
//  - With 1 thread, or when called inside another parallel region, the body
//    runs once on the calling thread with `n_workers() == 1`; loops then run
//    their chunks inline in index order, exactly like run_parallel_chunks.
//  - Nested parallel_for / run_parallel_chunks calls inside the region run
//    serially on the calling participant.

namespace detail { struct RegionShared; }

class ParallelRegion {
public:
    /// Index of this participant in [0, n_workers()). The calling thread is
    /// n_workers() - 1 (it also runs `single` blocks).
    int worker_index() const { return worker_; }
    int n_workers() const { return n_workers_; }
    bool is_master() const { return worker_ == n_workers_ - 1; }

    /// Work-shared loop with the chunks of run_parallel_chunks, then a barrier.
    /// `fn(chunk_begin, chunk_end, worker_index)`.
    void for_chunks(std::int64_t begin, std::int64_t end, std::int64_t chunk,
                    Scheduling sched,
                    const std::function<void(std::int64_t, std::int64_t, int)>& fn);

    /// Per-index work-shared loop; `fn` takes (index) or (index, worker_index).
    template <class Fn>
    void for_each(std::int64_t begin, std::int64_t end, std::int64_t chunk, Fn&& fn,
                  Scheduling sched = Scheduling::Dynamic) {
        for_chunks(begin, end, chunk, sched,
            [&fn](std::int64_t b, std::int64_t e, int worker) {
                if constexpr (detail::takes_worker<Fn>::value) {
                    for (std::int64_t i = b; i < e; ++i) fn(i, worker);
                } else {
                    for (std::int64_t i = b; i < e; ++i) fn(i);
                }
            });
    }

    /// All participants wait until every participant has arrived.
    void barrier();

    /// Run `fn()` on the master participant only, then barrier.
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

    /// True once a construct has failed (stable between barriers).
    bool cancelled() const;

    ParallelRegion(detail::RegionShared* shared, int worker, int n_workers)
        : shared_(shared), worker_(worker), n_workers_(n_workers) {}
    ParallelRegion(const ParallelRegion&) = delete;
    ParallelRegion& operator=(const ParallelRegion&) = delete;

private:
    void record_error();

    detail::RegionShared* shared_;
    int worker_;
    int n_workers_;
    std::int64_t chunk_base_ = 0;   // this participant's view of the shared chunk counter
};

/// Run `body` once on every pool worker (see above). Blocks until all
/// participants have returned; rethrows the first recorded exception.
void parallel_region(const std::function<void(ParallelRegion&)>& body);

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
