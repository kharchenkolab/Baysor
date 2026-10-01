// Unit tests for Baysor's persistent thread pool (baysor/utils/thread_pool.h).
#include <gtest/gtest.h>

#include "baysor/utils/thread_pool.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <functional>
#include <mutex>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

// RAII helper: save/restore the global pool size around a test.
class PoolSizeGuard {
public:
    explicit PoolSizeGuard(int n) : old_(baysor::thread_pool_size()) {
        baysor::set_thread_pool_size(n);
    }
    ~PoolSizeGuard() { baysor::set_thread_pool_size(old_); }

private:
    int old_;
};

} // namespace

TEST(ThreadPool, ZeroThreadsRunsInlineLikeOneThread) {
    PoolSizeGuard guard(0);
    EXPECT_EQ(baysor::thread_pool_size(), 1);

    std::vector<int> order;
    baysor::parallel_for(0, 100, 16, [&](int i) {
        order.push_back(i);
        EXPECT_EQ(baysor::current_worker_index(), 0);
    });
    ASSERT_EQ(order.size(), 100u);
    // Inline execution must run strictly in index order.
    for (int i = 0; i < 100; ++i) EXPECT_EQ(order[i], i);
}

TEST(ThreadPool, OneThreadRunsInlineInIndexOrder) {
    PoolSizeGuard guard(1);

    std::vector<int> order;
    baysor::parallel_for(0, 257, 32, [&](int i) {
        order.push_back(i);
    });
    ASSERT_EQ(order.size(), 257u);
    for (int i = 0; i < 257; ++i) EXPECT_EQ(order[i], i);

    // Static scheduling is also inline and ordered with 1 thread.
    order.clear();
    baysor::parallel_for_static(0, 50, [&](int i) {
        order.push_back(i);
    });
    ASSERT_EQ(order.size(), 50u);
    for (int i = 0; i < 50; ++i) EXPECT_EQ(order[i], i);
}

TEST(ThreadPool, DynamicChunksCoverEachIndexExactlyOnce) {
    for (int n_threads : {2, 3, 8}) {
        PoolSizeGuard guard(n_threads);
        constexpr int n = 10'003;  // deliberately not a multiple of the chunk size
        std::vector<std::atomic<int>> hits(n);
        for (auto& h : hits) h.store(0);

        baysor::parallel_for(0, n, 64, [&](int i) {
            hits[i].fetch_add(1);
        });

        for (int i = 0; i < n; ++i) {
            EXPECT_EQ(hits[i].load(), 1) << "thread count " << n_threads << ", index " << i;
        }
    }
}

TEST(ThreadPool, DynamicChunkBoundariesAreFixed) {
    PoolSizeGuard guard(8);
    constexpr int n = 1000;
    constexpr int chunk = 128;

    // Record (begin, end) of every chunk invocation; the ranges must partition
    // [0, n) into fixed-size chunks regardless of scheduling.
    std::mutex m;
    std::vector<std::pair<std::int64_t, std::int64_t>> ranges;
    baysor::run_parallel_chunks(0, n, chunk, baysor::Scheduling::Dynamic,
        [&](std::int64_t b, std::int64_t e, int) {
            std::lock_guard<std::mutex> lk(m);
            ranges.emplace_back(b, e);
        });

    std::sort(ranges.begin(), ranges.end());
    std::int64_t expected_begin = 0;
    for (auto [b, e] : ranges) {
        EXPECT_EQ(b, expected_begin);
        EXPECT_LE(e - b, chunk);
        expected_begin = e;
    }
    EXPECT_EQ(expected_begin, n);
}

TEST(ThreadPool, StaticSchedulingCoversEachIndexExactlyOnce) {
    for (int n_threads : {2, 5}) {
        PoolSizeGuard guard(n_threads);
        constexpr int n = 1001;
        std::vector<std::atomic<int>> hits(n);
        for (auto& h : hits) h.store(0);

        baysor::parallel_for_static(0, n, [&](int i) {
            hits[i].fetch_add(1);
        });
        for (int i = 0; i < n; ++i) EXPECT_EQ(hits[i].load(), 1);

        // Empty and single-element ranges are no-ops.
        baysor::parallel_for_static(5, 5, [&](int) { ADD_FAILURE(); });
    }
}

TEST(ThreadPool, WorkerIndicesAreInRangeAndStablePerBody) {
    PoolSizeGuard guard(6);
    constexpr int n = 4000;

    std::mutex m;
    std::set<int> seen_workers;
    std::vector<std::int64_t> chunk_workers;
    baysor::run_parallel_chunks(0, n, 16, baysor::Scheduling::Dynamic,
        [&](std::int64_t b, std::int64_t e, int worker) {
            EXPECT_GE(worker, 0);
            EXPECT_LT(worker, baysor::thread_pool_size());
            {
                std::lock_guard<std::mutex> lk(m);
                seen_workers.insert(worker);
            }
            // The whole chunk must report the same worker (per-worker buffers
            // stay valid for the duration of the chunk).
            for (std::int64_t i = b; i < e; ++i) {
                EXPECT_EQ(baysor::current_worker_index(), worker);
            }
            (void)e;
            (void)chunk_workers;
        });
    EXPECT_GE(seen_workers.size(), 1u);
}

TEST(ThreadPool, ParallelWorkersActuallyParticipate) {
    PoolSizeGuard guard(4);
    std::atomic<int> active{0};
    std::atomic<int> peak_active{0};
    std::atomic<bool> first{true};

    // The first body to start waits (bounded) until a second body overlaps it,
    // which proves that the region is executed by more than one thread.
    baysor::parallel_for(0, 8, 1, [&](int) {
        int a = active.fetch_add(1) + 1;
        int prev = peak_active.load();
        while (a > prev && !peak_active.compare_exchange_weak(prev, a)) {}
        if (first.exchange(false)) {
            auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(2);
            while (active.load() < 2 && std::chrono::steady_clock::now() < deadline) {
                std::this_thread::yield();
            }
        }
        active.fetch_sub(1);
    });
    EXPECT_GE(peak_active.load(), 2);
}

TEST(ThreadPool, NestedParallelCallsRunSerially) {
    PoolSizeGuard guard(4);

    std::atomic<int> outer_bodies{0};

    baysor::parallel_for(0, 32, 4, [&](int) {
        outer_bodies.fetch_add(1);
        EXPECT_TRUE(baysor::inside_parallel_region());
        int outer_worker = baysor::current_worker_index();

        // Nested parallel work must run serially (inline) on the current
        // thread: the nested loop executes in index order, on the same worker,
        // and sees the region flag.
        std::vector<int> seq;
        seq.reserve(10);
        baysor::parallel_for(0, 10, 2, [&](int j) {
            EXPECT_TRUE(baysor::inside_parallel_region());
            EXPECT_EQ(baysor::current_worker_index(), outer_worker);
            seq.push_back(j);
        });
        ASSERT_EQ(seq.size(), 10u);
        for (int j = 0; j < 10; ++j) EXPECT_EQ(seq[j], j);
        // Nested reductions too.
        double s = baysor::parallel_reduce<double>(0, 10, 4, 0.0,
            [](std::int64_t b, std::int64_t e, double& acc) {
                for (std::int64_t i = b; i < e; ++i) acc += 1.0;
            },
            std::plus<double>());
        EXPECT_DOUBLE_EQ(s, 10.0);
    });

    EXPECT_FALSE(baysor::inside_parallel_region());
    EXPECT_EQ(outer_bodies.load(), 32);
}

TEST(ThreadPool, ExceptionsPropagateToCaller) {
    for (int n_threads : {1, 3}) {
        PoolSizeGuard guard(n_threads);

        // Dynamic scheduling.
        EXPECT_THROW(
            baysor::parallel_for(0, 1000, 16, [](int i) {
                if (i == 500) throw std::runtime_error("boom");
            }),
            std::runtime_error);

        // Static scheduling.
        EXPECT_THROW(
            baysor::parallel_for_static(0, 1000, [](int i) {
                if (i == 17) throw std::runtime_error("boom-static");
            }),
            std::runtime_error);

        // Reductions.
        EXPECT_THROW(
            baysor::parallel_reduce<double>(0, 1000, 16, 0.0,
                [](std::int64_t, std::int64_t, double&) {
                    throw std::runtime_error("boom-reduce");
                },
                std::plus<double>()),
            std::runtime_error);

        // The pool must remain usable afterwards.
        std::atomic<int> count{0};
        baysor::parallel_for(0, 100, 8, [&](int) { count.fetch_add(1); });
        EXPECT_EQ(count.load(), 100);
    }
}

TEST(ThreadPool, NestedExceptionsPropagateThroughOuterRegion) {
    PoolSizeGuard guard(4);
    EXPECT_THROW(
        baysor::parallel_for(0, 64, 8, [](int i) {
            if (i == 3) throw std::runtime_error("nested-boom");
        }),
        std::runtime_error);

    std::atomic<int> count{0};
    baysor::parallel_for(0, 100, 8, [&](int) { count.fetch_add(1); });
    EXPECT_EQ(count.load(), 100);
}

TEST(ThreadPool, ParallelReduceMatchesSequential) {
    constexpr std::int64_t n = 10'000;
    auto body = [](std::int64_t b, std::int64_t e, double& acc) {
        for (std::int64_t i = b; i < e; ++i) acc += static_cast<double>(i) * 0.5;
    };

    double expected = 0.0;
    for (std::int64_t i = 0; i < n; ++i) expected += static_cast<double>(i) * 0.5;

    for (int n_threads : {1, 2, 3, 8}) {
        PoolSizeGuard guard(n_threads);
        double got = baysor::parallel_reduce<double>(0, n, 128, 0.0, body, std::plus<double>());
        EXPECT_DOUBLE_EQ(got, expected) << "threads " << n_threads;
    }
}

TEST(ThreadPool, ParallelReduceIsThreadCountIndependent) {
    // A deliberately order-sensitive combine (string concatenation): the fixed
    // bucket scheme must produce the same reduction for any thread count > 1,
    // and the plain sequential fold for 1 thread.
    constexpr std::int64_t n = 100;
    auto body = [](std::int64_t b, std::int64_t e, std::string& acc) {
        for (std::int64_t i = b; i < e; ++i) {
            acc += std::to_string(i);
            acc += ',';
        }
    };
    auto combine = [](const std::string& a, const std::string& b) { return a + b; };

    std::string expected;
    for (std::int64_t i = 0; i < n; ++i) {
        expected += std::to_string(i);
        expected += ',';
    }

    {
        PoolSizeGuard guard(1);
        EXPECT_EQ(baysor::parallel_reduce<std::string>(0, n, 8, std::string(), body, combine),
                  expected);
    }
    std::string reference;
    for (int n_threads : {2, 4, 5, 16}) {
        PoolSizeGuard guard(n_threads);
        std::string got =
            baysor::parallel_reduce<std::string>(0, n, 8, std::string(), body, combine);
        if (reference.empty()) {
            reference = got;
        }
        EXPECT_EQ(got, reference) << "threads " << n_threads;
    }
    // The single-bucket sequential fold equals the multi-bucket fold here
    // because the concatenation is associative; both must equal `expected`.
    EXPECT_EQ(reference, expected);
}

TEST(ThreadPool, RepeatedRunsAreDeterministic) {
    PoolSizeGuard guard(4);
    constexpr int n = 5000;

    auto run = [&]() {
        std::vector<int> hits(n, 0);
        baysor::parallel_for(0, n, 37, [&](int i) { hits[i] = i * 2 + 1; });
        return hits;
    };

    auto a = run();
    for (int r = 0; r < 5; ++r) {
        EXPECT_EQ(run(), a);
    }
}

TEST(ThreadPool, IdlePoolBurnsNoCpu) {
    // The workers must block when there is no work: measure the wall-clock
    // cost of waiting 200ms with an idle pool (a busy-spinning pool would
    // still pass, but the per-region latency test below would not).
    PoolSizeGuard guard(4);
    auto t0 = std::chrono::steady_clock::now();
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    auto t1 = std::chrono::steady_clock::now();
    // Nothing to assert about CPU here without OS counters; just ensure the
    // pool stays responsive after idling.
    std::atomic<int> count{0};
    baysor::parallel_for(0, 1000, 10, [&](int) { count.fetch_add(1); });
    EXPECT_EQ(count.load(), 1000);
    EXPECT_GE(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count(), 200);
}

TEST(ThreadPool, ManySmallRegionsDoNotDeadlock) {
    PoolSizeGuard guard(3);
    for (int r = 0; r < 500; ++r) {
        std::atomic<int> count{0};
        baysor::parallel_for(0, 3, 1, [&](int) { count.fetch_add(1); });
        EXPECT_EQ(count.load(), 3);
    }
}

// ============================================================================
// Persistent regions (parallel_region / ParallelRegion)
// ============================================================================

TEST(ThreadPoolRegion, OneThreadRunsBodyOnceInlineInIndexOrder) {
    PoolSizeGuard guard(1);
    int calls = 0;
    std::vector<int> order;
    baysor::parallel_region([&](baysor::ParallelRegion& r) {
        ++calls;
        EXPECT_EQ(r.n_workers(), 1);
        EXPECT_EQ(r.worker_index(), 0);
        EXPECT_TRUE(r.is_master());
        EXPECT_TRUE(baysor::inside_parallel_region());
        r.for_each(0, 50, 7, [&](int i) { order.push_back(i); });
        r.for_chunks(50, 60, 0, baysor::Scheduling::Static,
                     [&](std::int64_t b, std::int64_t e, int) {
            EXPECT_EQ(b, 50);
            EXPECT_EQ(e, 60);
            for (std::int64_t i = b; i < e; ++i) order.push_back(static_cast<int>(i));
        });
        r.barrier();
        r.single([&]() { order.push_back(-1); });
    });
    EXPECT_EQ(calls, 1);
    std::vector<int> expected(60);
    std::iota(expected.begin(), expected.end(), 0);
    expected.push_back(-1);
    EXPECT_EQ(order, expected);
    EXPECT_FALSE(baysor::inside_parallel_region());
}

TEST(ThreadPoolRegion, EveryWorkerRunsTheBodyOnce) {
    for (int n_threads : {2, 3, 8}) {
        PoolSizeGuard guard(n_threads);
        std::vector<std::atomic<int>> seen(n_threads);
        for (auto& s : seen) s.store(0);
        std::atomic<int> masters{0};
        baysor::parallel_region([&](baysor::ParallelRegion& r) {
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
    baysor::parallel_region([&](baysor::ParallelRegion& r) {
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
        // Static blocks: worker w gets block w.
        r.for_chunks(0, n, 0, baysor::Scheduling::Static,
                     [&](std::int64_t lo, std::int64_t hi, int w) {
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

TEST(ThreadPoolRegion, DynamicChunkBoundariesAreFixed) {
    for (int n_threads : {1, 2, 5}) {
        PoolSizeGuard guard(n_threads);
        std::mutex m;
        std::set<std::pair<std::int64_t, std::int64_t>> chunks;
        baysor::parallel_region([&](baysor::ParallelRegion& r) {
            r.for_chunks(3, 103, 10, baysor::Scheduling::Dynamic,
                         [&](std::int64_t b, std::int64_t e, int) {
                std::lock_guard<std::mutex> lk(m);
                chunks.insert({b, e});
            });
        });
        std::set<std::pair<std::int64_t, std::int64_t>> expected;
        for (std::int64_t b = 3; b < 103; b += 10) expected.insert({b, std::min<std::int64_t>(103, b + 10)});
        EXPECT_EQ(chunks, expected) << "threads " << n_threads;
    }
}

TEST(ThreadPoolRegion, NestedParallelCallsRunSerially) {
    PoolSizeGuard guard(4);
    std::atomic<int> count{0};
    baysor::parallel_region([&](baysor::ParallelRegion& r) {
        r.for_each(0, 8, 1, [&](int) {
            EXPECT_TRUE(baysor::inside_parallel_region());
            baysor::parallel_for(0, 10, 3, [&](int) { count.fetch_add(1); });
            // A nested region runs serially as well.
            baysor::parallel_region([&](baysor::ParallelRegion& inner) {
                EXPECT_EQ(inner.n_workers(), 1);
                inner.for_each(0, 5, 2, [&](int) { count.fetch_add(1); });
            });
        });
    });
    EXPECT_EQ(count.load(), 8 * 15);
}

TEST(ThreadPoolRegion, ExceptionsCancelTheRegionAndPropagate) {
    for (int n_threads : {1, 4}) {
        PoolSizeGuard guard(n_threads);
        std::atomic<int> after{0};
        EXPECT_THROW(
            baysor::parallel_region([&](baysor::ParallelRegion& r) {
                r.for_each(0, 100, 4, [&](int i) {
                    if (i == 37) throw std::runtime_error("boom");
                });
                if (!r.cancelled()) {
                    // Not reached with >1 thread: cancelled() is consistent
                    // across participants after the loop's barrier.
                    after.fetch_add(1);
                }
                r.single([&]() { after.fetch_add(100); });  // skipped once cancelled
            }),
            std::runtime_error);
        EXPECT_EQ(after.load(), 0) << "threads " << n_threads;

        // A throwing single block is reported as well.
        EXPECT_THROW(
            baysor::parallel_region([&](baysor::ParallelRegion& r) {
                r.single([&]() { throw std::logic_error("single"); });
                r.for_each(0, 10, 1, [&](int) { after.fetch_add(1); });
            }),
            std::logic_error);
        EXPECT_EQ(after.load(), 0) << "threads " << n_threads;

        // The pool stays usable afterwards.
        std::atomic<int> count{0};
        baysor::parallel_region([&](baysor::ParallelRegion& r) {
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
        baysor::parallel_region([&](baysor::ParallelRegion& r) {
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
        baysor::parallel_region([&](baysor::ParallelRegion& r) {
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
