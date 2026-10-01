// Unit tests for Baysor's persistent thread pool (baysor/utils/thread_pool.h).
#include <gtest/gtest.h>

#include "baysor/utils/thread_pool.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <functional>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using baysor::ParallelRegion;
using baysor::Scheduling;
using Chunk = std::pair<std::int64_t, std::int64_t>;
using ChunkFn = std::function<void(std::int64_t, std::int64_t, int)>;

// Sets the global pool size for the lifetime of the guard.
struct PoolSizeGuard {
    explicit PoolSizeGuard(int n) { baysor::set_thread_pool_size(n); }
    ~PoolSizeGuard() { baysor::set_thread_pool_size(old); }
    int old = baysor::thread_pool_size();
};

std::vector<int> range(int begin, int end) {
    std::vector<int> v(end - begin);
    std::iota(v.begin(), v.end(), begin);
    return v;
}

// [begin, end) cut into chunks of `chunk` indices.
std::vector<Chunk> fixed_chunks(std::int64_t begin, std::int64_t end, std::int64_t chunk) {
    std::vector<Chunk> chunks;
    for (std::int64_t b = begin; b < end; b += chunk) chunks.emplace_back(b, std::min(end, b + chunk));
    return chunks;
}

// The sorted chunks that `loop` passes to its chunk body, whose worker index
// must be in range.
std::vector<Chunk> chunks_of(const std::function<void(const ChunkFn&)>& loop) {
    std::mutex m;
    std::vector<Chunk> chunks;
    loop([&](std::int64_t b, std::int64_t e, int worker) {
        EXPECT_GE(worker, 0);
        EXPECT_LT(worker, baysor::thread_pool_size());
        std::lock_guard<std::mutex> lk(m);
        chunks.emplace_back(b, e);
    });
    std::sort(chunks.begin(), chunks.end());
    return chunks;
}

std::vector<Chunk> region_chunks(std::int64_t begin, std::int64_t end, std::int64_t chunk,
                                 Scheduling sched) {
    return chunks_of([&](const ChunkFn& fn) {
        baysor::parallel_region([&](ParallelRegion& r) { r.for_chunks(begin, end, chunk, sched, fn); });
    });
}

std::vector<Chunk> pool_chunks(std::int64_t begin, std::int64_t end, std::int64_t chunk,
                               Scheduling sched) {
    return chunks_of([&](const ChunkFn& fn) { baysor::run_parallel_chunks(begin, end, chunk, sched, fn); });
}

// Nested loops, reductions and regions must run inline (taking the pool for
// them would deadlock), in index order, on the worker of the enclosing chunk.
void expect_nested_calls_inline(int outer_worker) {
    std::vector<int> order;
    baysor::parallel_for(0, 10, 3, [&](int j, int w) {
        EXPECT_EQ(w, outer_worker);
        order.push_back(j);
    });
    EXPECT_EQ(order, range(0, 10));
    EXPECT_EQ(baysor::parallel_reduce<std::int64_t>(0, 10, 4, 0,
                  [](std::int64_t b, std::int64_t e, std::int64_t& acc) { acc += e - b; },
                  std::plus<std::int64_t>()),
              10);
    baysor::parallel_region([&](ParallelRegion& inner) {
        EXPECT_EQ(inner.n_workers(), 1);
        inner.for_each(10, 15, 2, [&](int j) { order.push_back(j); });
    });
    EXPECT_EQ(order, range(0, 15));
}

} // namespace

TEST(ThreadPool, SerialPoolRunsInlineInIndexOrder) {
    for (int n_threads : {0, 1}) {
        PoolSizeGuard guard(n_threads);
        EXPECT_EQ(baysor::thread_pool_size(), 1);
        std::vector<int> order;
        baysor::parallel_for(0, 100, 16, [&](int i, int w) {
            EXPECT_EQ(w, 0);
            order.push_back(i);
        });
        baysor::parallel_for_static(100, 150, [&](int i) { order.push_back(i); });
        EXPECT_EQ(order, range(0, 150));
    }
}

TEST(ThreadPool, ChunksPartitionTheRangeIndependentlyOfThreadCount) {
    for (int n_threads : {1, 2, 3, 8}) {
        PoolSizeGuard guard(n_threads);
        SCOPED_TRACE("threads " + std::to_string(n_threads));
        // Dynamic: fixed chunks; the range is deliberately not chunk-aligned.
        EXPECT_EQ(pool_chunks(3, 10'003, 64, Scheduling::Dynamic), fixed_chunks(3, 10'003, 64));
        EXPECT_EQ(region_chunks(3, 10'003, 64, Scheduling::Dynamic), fixed_chunks(3, 10'003, 64));
        // Static: one contiguous chunk per worker.
        const auto per_worker = fixed_chunks(0, 1001, (1001 + n_threads - 1) / n_threads);
        EXPECT_EQ(pool_chunks(0, 1001, 0, Scheduling::Static), per_worker);
        EXPECT_EQ(region_chunks(0, 1001, 0, Scheduling::Static), per_worker);
        // Empty ranges never call the body.
        EXPECT_TRUE(pool_chunks(5, 5, 4, Scheduling::Dynamic).empty());
        EXPECT_TRUE(pool_chunks(5, 5, 0, Scheduling::Static).empty());
    }
}

TEST(ThreadPool, ParallelWorkersActuallyParticipate) {
    PoolSizeGuard guard(4);
    std::atomic<int> active{0};
    std::atomic<int> peak_active{0};
    std::atomic<bool> first{true};

    // The first body to start waits (bounded) until a second body overlaps it,
    // which proves that the loop is executed by more than one thread.
    baysor::parallel_for(0, 8, 1, [&](int) {
        int a = active.fetch_add(1) + 1;
        int prev = peak_active.load();
        while (a > prev && !peak_active.compare_exchange_weak(prev, a)) {}
        if (first.exchange(false)) {
            auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(2);
            while (peak_active.load() < 2 && std::chrono::steady_clock::now() < deadline) {
                std::this_thread::yield();
            }
        }
        active.fetch_sub(1);
    });
    EXPECT_GE(peak_active.load(), 2);
}

TEST(ThreadPool, NestedCallsRunInlineOnTheEnclosingWorker) {
    PoolSizeGuard guard(4);
    std::atomic<int> outer{0};
    baysor::parallel_for(0, 32, 4, [&](int, int w) {
        expect_nested_calls_inline(w);
        ++outer;
    });
    baysor::parallel_region([&](ParallelRegion& r) {
        r.for_each(0, 8, 1, [&](int, int w) {
            expect_nested_calls_inline(w);
            ++outer;
        });
    });
    EXPECT_EQ(outer.load(), 40);
}

TEST(ThreadPool, ExceptionsPropagateToCaller) {
    for (int n_threads : {1, 3}) {
        PoolSizeGuard guard(n_threads);
        EXPECT_THROW(baysor::parallel_for(0, 1000, 16, [](int i) {
                         if (i == 500) throw std::runtime_error("dynamic");
                     }),
                     std::runtime_error);
        EXPECT_THROW(baysor::parallel_for_static(0, 1000, [](int i) {
                         if (i == 17) throw std::runtime_error("static");
                     }),
                     std::runtime_error);
        EXPECT_THROW(baysor::parallel_reduce<double>(0, 1000, 16, 0.0,
                         [](std::int64_t, std::int64_t, double&) { throw std::runtime_error("reduce"); },
                         std::plus<double>()),
                     std::runtime_error);
        EXPECT_THROW(baysor::parallel_for(0, 8, 1, [](int) {
                         baysor::parallel_for(0, 4, 1, [](int j) {
                             if (j == 2) throw std::runtime_error("nested");
                         });
                     }),
                     std::runtime_error);

        // The pool stays usable.
        std::atomic<int> count{0};
        baysor::parallel_for(0, 100, 8, [&](int) { count.fetch_add(1); });
        EXPECT_EQ(count.load(), 100);
    }
}

TEST(ThreadPool, ParallelReduceFoldsFixedBucketsInIndexOrder) {
    // Each bucket collects the chunks it accumulated; the fold concatenates.
    using Chunks = std::vector<Chunk>;
    const auto reduce = [] {
        return baysor::parallel_reduce<Chunks>(0, 10'000, 16, Chunks{},
            [](std::int64_t b, std::int64_t e, Chunks& acc) { acc.emplace_back(b, e); },
            [](Chunks a, const Chunks& b) {
                a.insert(a.end(), b.begin(), b.end());
                return a;
            });
    };
    {
        PoolSizeGuard guard(1);
        EXPECT_EQ(reduce(), (Chunks{{0, 10'000}}));
    }
    for (int n_threads : {2, 3, 8}) {
        PoolSizeGuard guard(n_threads);
        // 625 buckets of 16 capped at 64 buckets of ceil(10000 / 64) = 157.
        EXPECT_EQ(reduce(), fixed_chunks(0, 10'000, 157)) << "threads " << n_threads;
    }
}

TEST(ThreadPool, ManySmallLoopsDoNotDeadlock) {
    PoolSizeGuard guard(3);
    for (int r = 0; r < 500; ++r) {
        std::atomic<int> count{0};
        baysor::parallel_for(0, 3, 1, [&](int) { count.fetch_add(1); });
        EXPECT_EQ(count.load(), 3);
    }
}

TEST(ThreadPool, SubparHooksKeepSubparSplit) {
    PoolSizeGuard guard(4);
    std::mutex m;
    std::vector<std::tuple<int, int, int>> ranges;
    baysor::subpar_parallelize_range(3, 10, [&](int w, int start, int len) {
        std::lock_guard<std::mutex> lk(m);
        ranges.emplace_back(w, start, len);
    });
    std::sort(ranges.begin(), ranges.end());
    EXPECT_EQ(ranges, (std::vector<std::tuple<int, int, int>>{{0, 0, 4}, {1, 4, 3}, {2, 7, 3}}));

    EXPECT_EQ(chunks_of([](const ChunkFn& fn) {
                  baysor::umappp_parallel_range(std::size_t(10),
                      [&](std::size_t first, std::size_t last) { fn(first, last, 0); }, 3);
              }),
              (std::vector<Chunk>{{0, 4}, {4, 7}, {7, 10}}));
}

// ============================================================================
// Persistent regions (parallel_region / ParallelRegion)
// ============================================================================

TEST(ThreadPoolRegion, OneThreadRunsBodyOnceInlineInIndexOrder) {
    PoolSizeGuard guard(1);
    int calls = 0;
    std::vector<int> order;
    baysor::parallel_region([&](ParallelRegion& r) {
        ++calls;
        EXPECT_EQ(r.n_workers(), 1);
        EXPECT_EQ(r.worker_index(), 0);
        EXPECT_TRUE(r.is_master());
        r.for_each(0, 50, 7, [&](int i) { order.push_back(i); });
        r.for_chunks(50, 60, 0, Scheduling::Static, [&](std::int64_t b, std::int64_t e, int) {
            for (std::int64_t i = b; i < e; ++i) order.push_back(static_cast<int>(i));
        });
        r.barrier();
        r.single([&]() { order.push_back(60); });
    });
    EXPECT_EQ(calls, 1);
    EXPECT_EQ(order, range(0, 61));
}

TEST(ThreadPoolRegion, EveryWorkerRunsTheBodyOnce) {
    for (int n_threads : {2, 3, 8}) {
        PoolSizeGuard guard(n_threads);
        std::vector<std::atomic<int>> seen(n_threads);
        std::atomic<int> masters{0};
        baysor::parallel_region([&](ParallelRegion& r) {
            ASSERT_EQ(r.n_workers(), n_threads);
            ASSERT_GE(r.worker_index(), 0);
            ASSERT_LT(r.worker_index(), n_threads);
            seen[r.worker_index()].fetch_add(1);
            if (r.is_master()) masters.fetch_add(1);
        });
        for (int w = 0; w < n_threads; ++w) {
            EXPECT_EQ(seen[w].load(), 1) << "threads " << n_threads << " worker " << w;
        }
        EXPECT_EQ(masters.load(), 1);
    }
}

TEST(ThreadPoolRegion, SequencesOfLoopsCoverEachIndexOnceWithBarriersBetween) {
    PoolSizeGuard guard(4);
    constexpr int n = 3001;
    std::vector<int> a(n, 0), b(n, 0);
    std::atomic<int> violations{0};
    std::vector<long> per_worker(4, 0);
    long total = -1;
    baysor::parallel_region([&](ParallelRegion& r) {
        for (int rep = 0; rep < 50; ++rep) {
            // Dynamic loop with an odd chunk, then a loop that reads what the
            // previous loop wrote at other indices (needs the implicit barrier).
            r.for_each(0, n, 13, [&](int i) { a[i] = rep * n + i; });
            r.for_each(0, n, 29, [&](int i) {
                int j = n - 1 - i;
                if (a[j] != rep * n + j) violations.fetch_add(1);
                b[i] += 1;
            });
            // Empty loops must not desynchronise the participants.
            r.for_each(5, 5, 4, [&](int) { violations.fetch_add(1); });
        }
        r.for_chunks(0, n, 0, Scheduling::Static, [&](std::int64_t lo, std::int64_t hi, int w) {
            EXPECT_EQ(w, r.worker_index());
            for (std::int64_t i = lo; i < hi; ++i) per_worker[w] += b[i];
        });
        r.single([&]() { total = per_worker[0] + per_worker[1] + per_worker[2] + per_worker[3]; });
        // After single's barrier every participant sees the result.
        if (total != 50L * n) violations.fetch_add(1);
    });
    EXPECT_EQ(violations.load(), 0);
    EXPECT_EQ(total, 50L * n);
    for (int i = 0; i < n; ++i) ASSERT_EQ(b[i], 50) << i;
}

TEST(ThreadPoolRegion, ExceptionsCancelTheRegionAndPropagate) {
    for (int n_threads : {1, 4}) {
        PoolSizeGuard guard(n_threads);
        std::atomic<int> after{0};
        EXPECT_THROW(
            baysor::parallel_region([&](ParallelRegion& r) {
                r.for_each(0, 100, 4, [&](int i) {
                    if (i == 37) throw std::runtime_error("boom");
                });
                // cancelled() is consistent across participants after the
                // loop's barrier.
                if (!r.cancelled()) after.fetch_add(1);
                r.single([&]() { after.fetch_add(100); });  // skipped once cancelled
            }),
            std::runtime_error);
        EXPECT_EQ(after.load(), 0) << "threads " << n_threads;

        // A throwing single block is reported as well.
        EXPECT_THROW(
            baysor::parallel_region([&](ParallelRegion& r) {
                r.single([&]() { throw std::logic_error("single"); });
                r.for_each(0, 10, 1, [&](int) { after.fetch_add(1); });
            }),
            std::logic_error);
        EXPECT_EQ(after.load(), 0) << "threads " << n_threads;

        // The pool stays usable.
        std::atomic<int> count{0};
        baysor::parallel_region([&](ParallelRegion& r) {
            r.for_each(0, 64, 4, [&](int) { count.fetch_add(1); });
        });
        EXPECT_EQ(count.load(), 64);
    }
}

TEST(ThreadPoolRegion, ManyRegionsWithLongSerialSectionsDoNotDeadlock) {
    PoolSizeGuard guard(4);
    long sum = 0;
    for (int rep = 0; rep < 200; ++rep) {
        std::vector<int> v(64, 0);
        baysor::parallel_region([&](ParallelRegion& r) {
            r.for_each(0, 64, 1, [&](int i) { v[i] = i; });
            r.single([&]() {
                // Longer than the spin budget every few reps: waiters block.
                if (rep % 20 == 0) std::this_thread::sleep_for(std::chrono::milliseconds(1));
                for (int x : v) sum += x;
            });
            r.barrier();
        });
    }
    EXPECT_EQ(sum, 200L * (63 * 64 / 2));
}

TEST(ThreadPoolRegion, OversubscribedRegionsSleepAtBarriersAndStayCorrect) {
    // More workers than physical cores: barrier waiters do not spin and
    // always go to sleep (futex on Linux), the path most sensitive to lost
    // wake-ups.
    const int n_threads = 2 * baysor::default_thread_count() + 1;
    PoolSizeGuard guard(n_threads);
    std::vector<long> per_worker(n_threads, 0);
    long expected = 0;
    for (int rep = 0; rep < 100; ++rep) {
        baysor::parallel_region([&](ParallelRegion& r) {
            for (int phase = 0; phase < 10; ++phase) {
                r.for_each(0, 97, 3, [&](int i, int w) { per_worker[w] += i; });
                r.barrier();
            }
        });
        expected += 10L * (96 * 97 / 2);
    }
    long total = 0;
    for (long v : per_worker) total += v;
    EXPECT_EQ(total, expected);
}
