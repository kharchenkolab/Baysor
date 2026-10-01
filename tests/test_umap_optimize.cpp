// Baysor's UMAP layout optimiser (umap_optimize.h) against umappp's.

#include "baysor/processing/data_processing/umap_optimize.h"
#include "baysor/processing/data_processing/umap_wrappers.h"

#include <gtest/gtest.h>

#include <cstring>
#include <random>
#include <vector>

namespace {

using Graph = umappp::NeighborList<int, double>;

// A random weighted graph shaped like a symmetrised UMAP kNN graph: every
// vertex has 5-20 neighbours with similarities in (0, 1].
Graph random_graph(int n, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> deg(5, 20), node(0, n - 1);
    std::uniform_real_distribution<double> w(0.01, 1.0);
    Graph g(n);
    for (int i = 0; i < n; ++i) {
        const int k = deg(rng);
        for (int j = 0; j < k; ++j) {
            int t = node(rng);
            if (t == i) t = (t + 1) % n;
            g[i].emplace_back(t, w(rng));
        }
    }
    return g;
}

std::vector<double> random_embedding(int n, int ndim, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> u(-10.0, 10.0);
    std::vector<double> e(static_cast<size_t>(n) * ndim);
    for (auto& v : e) v = u(rng);
    return e;
}

bool bitwise_equal(const std::vector<double>& a, const std::vector<double>& b) {
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

// Runs umappp's serial optimiser and Baysor's on the same input.
void expect_same_as_umappp(int ndim, double a, double b, int epochs, int epoch_limit) {
    const int n = 300;
    const Graph g = random_graph(n, 11);
    const std::vector<double> init = random_embedding(n, ndim, 12);

    auto setup_ref = umappp::internal::similarities_to_epochs<int, double>(g, epochs, 5.0);
    auto setup_new = setup_ref;
    std::vector<double> emb_ref = init, emb_new = init;
    std::mt19937_64 rng_ref(42), rng_new(42);

    umappp::internal::optimize_layout<int, double>(
        ndim, emb_ref.data(), setup_ref, a, b, 1.0, 1.0, rng_ref, epoch_limit);
    baysor::umap_detail::optimize_layout_dispatch(
        static_cast<size_t>(ndim), emb_new.data(), setup_new, a, b, 1.0, 1.0, rng_new, epoch_limit);

    EXPECT_FALSE(bitwise_equal(emb_ref, init)) << "the optimiser did not move the embedding";
    EXPECT_TRUE(bitwise_equal(emb_ref, emb_new)) << "ndim " << ndim;
    EXPECT_EQ(setup_ref.current_epoch, setup_new.current_epoch);
    EXPECT_TRUE(bitwise_equal(setup_ref.epoch_of_next_sample, setup_new.epoch_of_next_sample));
    EXPECT_TRUE(bitwise_equal(setup_ref.epoch_of_next_negative_sample, setup_new.epoch_of_next_negative_sample));
    // Both consumed the same random numbers.
    EXPECT_EQ(rng_ref(), rng_new());
}

} // namespace

TEST(UmapOptimize, FixedDimensionMatchesUmapppBitwise) {
    // a, b of find_ab(2.0, 0.1) (the NCV colours) and of find_ab(1.0, 0.1).
    for (int ndim : {2, 3}) {
        expect_same_as_umappp(ndim, 0.544663, 0.842052, 50, 0);
        expect_same_as_umappp(ndim, 1.576943, 0.895061, 50, 0);
    }
}

TEST(UmapOptimize, RuntimeDimensionMatchesUmapppBitwise) {
    expect_same_as_umappp(1, 0.544663, 0.842052, 30, 0);
    expect_same_as_umappp(4, 0.544663, 0.842052, 30, 0);
}

TEST(UmapOptimize, EpochLimitMatchesUmapppBitwise) {
    expect_same_as_umappp(3, 0.544663, 0.842052, 50, 20);
}
