// Baysor's UMAP layout optimiser (umap_optimize.h) against umappp's.

#include "baysor/processing/data_processing/umap_optimize.h"
#include "baysor/processing/data_processing/umap_wrappers.h"

#include <gtest/gtest.h>

#include <cmath>
#include <cstring>
#include <limits>
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
    baysor::umap_detail::optimize_layout_dispatch</*FastPow_=*/false>(
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

TEST(UmapOptimize, FastPowMatchesStdPow) {
    // d2 is floored at DBL_EPSILON by quick_squared_distance; cover that up to
    // far beyond any embedding distance, for the b values of the UMAP kernels.
    std::mt19937_64 rng(5);
    std::uniform_real_distribution<double> log_x(std::log(std::numeric_limits<double>::epsilon()), 30.0);
    double max_rel = 0;
    for (double b : {0.5, 0.79, 0.842052, 0.895061, 1.0, 1.2, 1.9}) {
        for (int i = 0; i < 200000; ++i) {
            const double x = std::exp(log_x(rng));
            const double ref = std::pow(x, b);
            max_rel = std::max(max_rel, std::abs(baysor::umap_detail::fast_pow(x, b) - ref) / ref);
        }
        EXPECT_NEAR(baysor::umap_detail::fast_pow(1.0, b), 1.0, 3e-16);
    }
    EXPECT_LT(max_rel, 1e-13);
}

TEST(UmapOptimize, FastPowTablesMatchLibm) {
    const auto& t = baysor::umap_detail::fast_pow_tables;
    for (int i = 0; i < 256; ++i) {
        const double c = 1.0 + (i + 0.5) / 256.0;
        EXPECT_EQ(t.inv_c[i], 1.0 / c);
        EXPECT_NEAR(t.ln_c[i], std::log(c), 4.5e-16 * std::log(c)) << i;  // <= 2 ulp
    }
    for (int j = 0; j < 256; ++j) {
        EXPECT_NEAR(t.exp2_j256[j], std::exp2(j / 256.0), 4.5e-16 * std::exp2(j / 256.0)) << j;  // <= 2 ulp
    }
}

TEST(UmapOptimize, FastPowFallsBackOutsideItsDomain) {
    using baysor::umap_detail::fast_pow;
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double sub = std::numeric_limits<double>::denorm_min() * 12345;
    EXPECT_EQ(fast_pow(0.0, 0.8), std::pow(0.0, 0.8));
    EXPECT_EQ(fast_pow(sub, 0.8), std::pow(sub, 0.8));
    EXPECT_EQ(fast_pow(inf, 0.8), inf);
    EXPECT_TRUE(std::isnan(fast_pow(nan, 0.8)));
    EXPECT_TRUE(std::isnan(fast_pow(-2.0, 0.8)));
    EXPECT_TRUE(std::isnan(fast_pow(2.0, nan)));
    // |b ln x| beyond the exp() range of the approximation.
    EXPECT_EQ(fast_pow(1e300, 1.5), std::pow(1e300, 1.5));
    EXPECT_EQ(fast_pow(1e-300, 1.5), std::pow(1e-300, 1.5));
    // Still inside the range: approximated. The error grows with |b ln x|
    // (rounding of the double y = b ln x): ~1e-13 here, ~1e-15 for UMAP's d2.
    EXPECT_NEAR(fast_pow(1e300, 1.0) / 1e300, 1.0, 1e-12);
}

TEST(UmapOptimize, FastPowLayoutStaysCloseToExactLayout) {
    // fast_pow only changes the last bits of the gradient coefficients; over a
    // short run the layout stays close to the one with std::pow.
    const int n = 300, ndim = 3;
    const Graph g = random_graph(n, 21);
    const std::vector<double> init = random_embedding(n, ndim, 22);
    auto setup_exact = umappp::internal::similarities_to_epochs<int, double>(g, 20, 5.0);
    auto setup_fast = setup_exact;
    std::vector<double> emb_exact = init, emb_fast = init;
    std::mt19937_64 rng_exact(3), rng_fast(3);
    baysor::umap_detail::optimize_layout_dispatch<false>(
        ndim, emb_exact.data(), setup_exact, 0.544663, 0.842052, 1.0, 1.0, rng_exact, 0);
    baysor::umap_detail::optimize_layout_dispatch<true>(
        ndim, emb_fast.data(), setup_fast, 0.544663, 0.842052, 1.0, 1.0, rng_fast, 0);
    double max_diff = 0;
    for (size_t i = 0; i < emb_exact.size(); ++i) {
        ASSERT_TRUE(std::isfinite(emb_fast[i]));
        max_diff = std::max(max_diff, std::abs(emb_exact[i] - emb_fast[i]));
    }
    EXPECT_LT(max_diff, 1e-3);
}
