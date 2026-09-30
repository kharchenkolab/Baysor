// Determinism tests for the multi-threaded E-step RNG scheme: streams are
// keyed by (iteration, fixed-size chunk index), so results must not depend on
// scheduling or on the number of threads, and must change with the iteration
// salt. With 1 thread the global stream is used in index order instead.
#include <gtest/gtest.h>

#include "baysor/processing/bmm_algorithm/bmm_algorithm.h"
#include "baysor/processing/distributions/categorical_smoothed.h"
#include "baysor/processing/distributions/mv_normal.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/processing/models/bmm_data.h"
#include "baysor/processing/models/component.h"
#include "baysor/utils/general.h"
#include "baysor/utils/thread_pool.h"

#include <cmath>
#include <cstdint>
#include <vector>

namespace {

class PoolSizeGuard {
public:
    explicit PoolSizeGuard(int n) : old_(baysor::thread_pool_size()) {
        baysor::set_thread_pool_size(n);
    }
    ~PoolSizeGuard() { baysor::set_thread_pool_size(old_); }

private:
    int old_;
};

baysor::Component<2> make_component(
    double center_x, double center_y, double variance, int n_samples, int guid
) {
    Eigen::Vector2d center;
    center << center_x, center_y;
    Eigen::Matrix2d sigma = Eigen::Matrix2d::Identity() * variance;
    baysor::MvNormal<2> pos_params(center, sigma);
    baysor::CategoricalSmoothed comp_params(1, 1.0);
    comp_params.set_dense_counts({1.0f});
    baysor::Component<2> comp(pos_params, comp_params, std::nullopt, guid);
    comp.prior_probability = 1.0;
    comp.n_samples = n_samples;
    comp.confidence = 1.0;
    return comp;
}

// A chain of molecules whose two components have exactly identical densities,
// so every E-step assignment draw is a fair coin flip and the RNG scheme is
// observable. Deliberately larger than the E-step chunk size (1024) so that
// multiple per-chunk streams are used.
baysor::BmmData<2> make_symmetric_chain(int n) {
    baysor::BmmData<2> data;
    data.position_data.resize(2, n);
    for (int i = 0; i < n; ++i) {
        data.position_data(0, i) = 0.1 * i;
        data.position_data(1, i) = 0.01 * ((i * 37) % 13);
    }
    data.composition_data.assign(n, 0);
    data.confidence.assign(n, 1.0);

    std::vector<int> src, dst;
    std::vector<double> wt;
    // Two edge sets (i, i+1) and (i, i+2): every molecule sees neighbors from
    // both classes with equal weight, so every E-step draw is a fair coin flip.
    for (int i = 0; i + 1 < n; ++i) {
        src.push_back(i);
        dst.push_back(i + 1);
        wt.push_back(1.0);
    }
    for (int i = 0; i + 2 < n; ++i) {
        src.push_back(i);
        dst.push_back(i + 2);
        wt.push_back(1.0);
    }
    data.adj_list = baysor::AdjList::from_edge_list(
        src.data(), dst.data(), wt.data(), static_cast<int>(src.size()), n);

    // Two wide, identical Gaussians: their pdfs are effectively constant and
    // equal along the chain.
    data.components.push_back(make_component(0.0, 0.0, 1e8, n / 2, 1));
    data.components.push_back(make_component(0.0, 0.0, 1e8, n / 2, 2));

    data.assignment.resize(n);
    for (int i = 0; i < n; ++i) data.assignment[i] = (i % 2) + 1;
    data.max_component_guid = 2;
    data.cluster_penalty_mult = 0.25;
    data.use_gene_smoothing = true;
    data.mrf_strength = 0.1;
    data.real_edge_weight = 1.0;
    return data;
}

// Three compact blobs with within-blob chain edges; positions are generated
// deterministically (no RNG) so fixtures are reproducible.
baysor::BmmData<2> make_blob_data(int n_per_blob) {
    const int n = 3 * n_per_blob;
    baysor::BmmData<2> data;
    data.position_data.resize(2, n);
    data.assignment.resize(n);
    for (int b = 0; b < 3; ++b) {
        double cx = 10.0 * b;
        for (int j = 0; j < n_per_blob; ++j) {
            int i = b * n_per_blob + j;
            double u = static_cast<double>((i * 2654435761u) % 1000u) / 1000.0 - 0.5;
            double v = static_cast<double>((i * 40503u) % 1000u) / 1000.0 - 0.5;
            data.position_data(0, i) = cx + 2.0 * u;
            data.position_data(1, i) = 2.0 * v;
            data.assignment[i] = b + 1;
        }
    }
    data.composition_data.assign(n, 0);
    data.confidence.assign(n, 1.0);

    std::vector<int> src, dst;
    std::vector<double> wt;
    for (int b = 0; b < 3; ++b) {
        for (int j = 0; j + 1 < n_per_blob; ++j) {
            src.push_back(b * n_per_blob + j);
            dst.push_back(b * n_per_blob + j + 1);
            wt.push_back(1.0);
        }
    }
    data.adj_list = baysor::AdjList::from_edge_list(
        src.data(), dst.data(), wt.data(), static_cast<int>(src.size()), n);

    for (int b = 0; b < 3; ++b) {
        data.components.push_back(make_component(10.0 * b, 0.0, 1.0, n_per_blob, b + 1));
    }
    data.max_component_guid = 3;
    data.cluster_penalty_mult = 0.25;
    data.use_gene_smoothing = true;
    data.mrf_strength = 0.1;
    data.real_edge_weight = 1.0;
    return data;
}

} // namespace

TEST(RngDeterminism, EstepIsThreadCountIndependent) {
    constexpr int n = 3000;
    std::vector<int> reference;
    for (int n_threads : {2, 3, 5, 8}) {
        PoolSizeGuard guard(n_threads);
        auto data = make_symmetric_chain(n);
        baysor::expect_dirichlet_spatial(data, /*stochastic=*/true, /*rng_salt=*/7);
        if (reference.empty()) {
            reference = data.assignment;
        }
        EXPECT_EQ(data.assignment, reference) << "threads " << n_threads;
    }

    // The fixture must actually exercise sampling: the assignment should not
    // be trivially unchanged everywhere.
    auto data = make_symmetric_chain(n);
    EXPECT_NE(data.assignment, reference);
}

TEST(RngDeterminism, EstepIsDeterministicAcrossRepeatedRuns) {
    PoolSizeGuard guard(4);
    constexpr int n = 2500;
    std::vector<int> reference;
    for (int r = 0; r < 5; ++r) {
        auto data = make_symmetric_chain(n);
        baysor::expect_dirichlet_spatial(data, /*stochastic=*/true, /*rng_salt=*/3);
        if (reference.empty()) reference = data.assignment;
        EXPECT_EQ(data.assignment, reference);
    }
}

TEST(RngDeterminism, EstepStreamsDifferBetweenIterations) {
    PoolSizeGuard guard(4);
    constexpr int n = 3000;
    auto data_a = make_symmetric_chain(n);
    auto data_b = make_symmetric_chain(n);
    baysor::expect_dirichlet_spatial(data_a, /*stochastic=*/true, /*rng_salt=*/1);
    baysor::expect_dirichlet_spatial(data_b, /*stochastic=*/true, /*rng_salt=*/2);
    // Different iterations must draw from different streams.
    EXPECT_NE(data_a.assignment, data_b.assignment);
}

TEST(RngDeterminism, SingleThreadedEstepStaysReproducible) {
    PoolSizeGuard guard(1);
    constexpr int n = 1500;
    // With 1 thread the E-step draws from the global stream in index order
    // (unchanged from the OpenMP build), so reproducibility holds for the same
    // global seed.
    std::vector<int> reference;
    for (int r = 0; r < 3; ++r) {
        baysor::reset_global_xoshiro_rng(1);
        auto data = make_symmetric_chain(n);
        baysor::expect_dirichlet_spatial(data, /*stochastic=*/true, /*rng_salt=*/0);
        if (reference.empty()) reference = data.assignment;
        EXPECT_EQ(data.assignment, reference);
    }
}

TEST(RngDeterminism, BmmIsThreadCountIndependent) {
    std::vector<int> reference_assignment;
    std::size_t reference_components = 0;
    for (int n_threads : {2, 4, 8}) {
        PoolSizeGuard guard(n_threads);
        auto data = make_blob_data(800);
        baysor::bmm(data,
                    /*min_molecules_drop=*/2,
                    /*n_iters=*/6,
                    /*assignment_history_depth=*/0,
                    /*verbose=*/false,
                    /*component_split_step=*/3,
                    /*refine=*/false);
        if (reference_assignment.empty()) {
            reference_assignment = data.assignment;
            reference_components = data.components.size();
        }
        EXPECT_EQ(data.assignment, reference_assignment) << "threads " << n_threads;
        EXPECT_EQ(data.components.size(), reference_components) << "threads " << n_threads;
    }
}
