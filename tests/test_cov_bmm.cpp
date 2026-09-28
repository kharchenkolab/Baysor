// Coverage tests for the BMM segmentation core (COV-2):
//   src/processing/bmm_algorithm/{bmm_algorithm,tracing,history_analysis,
//                                 compartment_segmentation}.cpp
#include <gtest/gtest.h>

#include "baysor/processing/bmm_algorithm/bmm_algorithm.h"
#include "baysor/processing/bmm_algorithm/compartment_segmentation.h"
#include "baysor/processing/bmm_algorithm/history_analysis.h"
#include "baysor/processing/bmm_algorithm/tracing.h"

#include "test_cov_helpers.h"

#include <Eigen/Dense>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace {

using baysor::AdjList;
using baysor::BmmData;
using baysor::CategoricalSmoothed;
using baysor::Component;
using baysor::MvNormal;
using baysor::ShapePrior;

// Three disconnected pairs of molecules, assigned to three components
// (sizes 3, 3, 2) — mirroring the stable fixture used elsewhere.
BmmData<2> cov2_make_disconnected_data() {
    BmmData<2> data;
    data.position_data.resize(2, 8);
    data.position_data <<
        0.0,  0.1,  0.2,   10.0, 10.1, 10.2,   20.0, 20.1,
        0.0,  0.0,  0.1,    0.0,  0.0,  0.1,    0.0,  0.0;
    data.composition_data = {0, 0, 0, 1, 1, 1, 0, 0};
    data.confidence.assign(8, 1.0);

    const int edge_src[] = {0, 1, 3, 4, 6};
    const int edge_dst[] = {1, 2, 4, 5, 7};
    const double edge_wt[] = {1.0, 1.0, 1.0, 1.0, 1.0};
    data.adj_list = AdjList::from_edge_list(edge_src, edge_dst, edge_wt, 5, 8);

    ShapePrior<2> prior;
    prior.std_values << 0.25, 0.25;
    prior.std_value_stds << 0.05, 0.05;
    prior.n_samples = 3;

    const Eigen::Matrix2d sigma = Eigen::Matrix2d::Identity() * 0.05;
    const Eigen::Vector2d centers[] = {
        (Eigen::Vector2d() << 0.1, 0.03).finished(),
        (Eigen::Vector2d() << 10.1, 0.03).finished(),
        (Eigen::Vector2d() << 20.05, 0.0).finished()
    };

    for (int ci = 0; ci < 3; ++ci) {
        MvNormal<2> pos_params(centers[ci], sigma);
        CategoricalSmoothed comp_params(2, 1.0);
        comp_params.set_uniform_counts(1.0f);
        data.components.emplace_back(pos_params, comp_params, prior, ci + 1);
    }

    data.assignment = {1, 1, 1, 2, 2, 2, 3, 3};
    data.max_component_guid = 3;
    data.noise_position_density = 1e-6;
    data.noise_density = 1e-6;
    data.prior_seg_confidence = 0.2;
    data.cluster_penalty_mult = 0.25;
    data.use_gene_smoothing = true;
    data.min_nuclei_frac = 0.1;
    data.mrf_strength = 0.1;
    data.real_edge_weight = 1.0;
    return data;
}

// Molecule 0 sits between component 1 (main segment 1) and component 2
// (main segment 2). With prior_seg_confidence > 0 the under-segmentation
// branch of adjust_densities_by_prior_segmentation penalizes component 2 for
// molecule 0, flipping the winner.
BmmData<2> cov2_make_undersegmentation_data(bool with_prior_segmentation) {
    BmmData<2> data;
    data.position_data.resize(2, 3);
    data.position_data <<
        0.0, 0.0,  0.1,
        0.0, 0.05, 0.0;
    data.composition_data = {0, 0, 0};
    data.confidence.assign(3, 1.0);

    const int edge_src[] = {0, 0};
    const int edge_dst[] = {1, 2};
    const double edge_wt[] = {1.0, 1.0};
    data.adj_list = AdjList::from_edge_list(edge_src, edge_dst, edge_wt, 2, 3);

    // Component 1 is far from molecule 0, component 2 sits right on it:
    // without the prior adjustment molecule 0 clearly prefers component 2.
    Eigen::Vector2d mu1, mu2;
    mu1 << 0.25, 0.18;
    mu2 << 0.10, 0.0;
    const Eigen::Matrix2d sigma = Eigen::Matrix2d::Identity() * 0.04;
    CategoricalSmoothed comp_params(1, 1.0);
    comp_params.set_dense_counts({1.0f});
    data.components.emplace_back(MvNormal<2>(mu1, sigma), comp_params, std::nullopt, 1);
    data.components.emplace_back(MvNormal<2>(mu2, sigma), comp_params, std::nullopt, 2);

    data.assignment = {1, 1, 2};
    data.max_component_guid = 2;
    data.mrf_strength = 0.1;
    data.use_gene_smoothing = true;
    data.prior_seg_confidence = 0.5;

    if (with_prior_segmentation) {
        data.segment_per_molecule = {1, 1, 2};
        data.n_molecules_per_segment = {2, 1};
        data.update_n_mols_per_segment();
    }
    return data;
}

// Two well-separated groups of molecules in 3D with prior segmentation:
// exercises every <3> instantiation of the segmentation loop.
BmmData<3> cov2_make_3d_data() {
    BmmData<3> data;
    constexpr int n = 12;
    data.position_data.resize(3, n);
    for (int i = 0; i < n; ++i) {
        if (i < 8) {
            data.position_data(0, i) = 0.06 * i - 0.2;
            data.position_data(1, i) = 0.04 * (i % 3);
            data.position_data(2, i) = 0.02 * i;
        } else {
            data.position_data(0, i) = 10.0 + 0.06 * (i - 8);
            data.position_data(1, i) = 0.04 * ((i - 8) % 3);
            data.position_data(2, i) = 0.02 * (i - 8);
        }
    }
    data.composition_data = {0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0};
    data.confidence.assign(n, 1.0);

    std::vector<int> src, dst;
    std::vector<double> wts;
    for (int i = 0; i < 7; ++i) {
        src.push_back(i);
        dst.push_back(i + 1);
        wts.push_back(1.0);
    }
    for (int i = 8; i < 11; ++i) {
        src.push_back(i);
        dst.push_back(i + 1);
        wts.push_back(1.0);
    }
    data.adj_list = AdjList::from_edge_list(
        src.data(), dst.data(), wts.data(),
        static_cast<int>(src.size()), n);

    ShapePrior<3> prior;
    prior.std_values << 0.3, 0.3, 0.3;
    prior.std_value_stds << 0.05, 0.05, 0.05;
    prior.n_samples = 10;

    const Eigen::Matrix3d sigma = Eigen::Matrix3d::Identity() * 0.02;
    Eigen::Vector3d center_a, center_b;
    center_a << 0.0, 0.0, 0.0;
    center_b << 10.0, 0.0, 0.0;
    CategoricalSmoothed comp_params(2, 1.0);
    comp_params.set_dense_counts({1.0f, 1.0f});
    data.components.emplace_back(MvNormal<3>(center_a, sigma), comp_params, prior, 1);
    data.components.emplace_back(MvNormal<3>(center_b, sigma), comp_params, prior, 2);

    data.assignment = {1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2};
    data.max_component_guid = 2;
    data.noise_position_density = 1e-6;
    data.prior_seg_confidence = 0.3;
    data.cluster_penalty_mult = 0.25;
    data.use_gene_smoothing = true;
    data.min_nuclei_frac = 0.1;
    data.mrf_strength = 0.1;
    data.real_edge_weight = 1.0;
    data.segment_per_molecule = {1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2};
    data.n_molecules_per_segment = {8, 4};
    return data;
}

Component<2> cov2_make_component_2d(int guid) {
    Eigen::Vector2d mu = Eigen::Vector2d::Zero();
    const Eigen::Matrix2d sigma = Eigen::Matrix2d::Identity();
    CategoricalSmoothed comp_params(1, 1.0);
    comp_params.set_dense_counts({1.0f});
    return Component<2>(MvNormal<2>(mu, sigma), comp_params, std::nullopt, guid);
}

Component<3> cov2_make_component_3d(int guid) {
    Eigen::Vector3d mu = Eigen::Vector3d::Zero();
    const Eigen::Matrix3d sigma = Eigen::Matrix3d::Identity();
    CategoricalSmoothed comp_params(1, 1.0);
    comp_params.set_dense_counts({1.0f});
    return Component<3>(MvNormal<3>(mu, sigma), comp_params, std::nullopt, guid);
}

} // namespace

// ============================================================================
// maximize / M-step
// ============================================================================

TEST(Cov2Bmm, MaximizeTracksClusterPerCellMode) {
    auto data = cov2_make_disconnected_data();
    data.cluster_per_molecule = {1, 1, 1, 2, 2, 2, 1, 1};
    data.cluster_per_cell.clear();

    baysor::maximize(data, /*freeze_composition=*/false, /*freeze_position=*/false);

    // Mode of cluster_per_molecule per component: comp1 {1,1,1} -> 1,
    // comp2 {2,2,2} -> 2, comp3 {1,1} -> 1.
    EXPECT_EQ(data.cluster_per_cell, (std::vector<int>{1, 2, 1}));
}

TEST(Cov2Bmm, NonFiniteNoiseDensityIsWarnedAndResetToZero) {
    // The test name promises the warning; capture it instead of only pinning
    // the reset value.
    auto sink = std::make_shared<baysor_test::CapturingSink>();
    baysor_test::LoggerGuard logger(sink);

    // NaN path: infinite position density times an empty component list.
    auto nan_data = BmmData<2>{};
    nan_data.noise_position_density = std::numeric_limits<double>::infinity();
    nan_data.noise_density = 42.0;  // must be overwritten by the M-step
    baysor::maximize(nan_data, /*freeze_composition=*/false, /*freeze_position=*/false);
    EXPECT_DOUBLE_EQ(nan_data.noise_density, 0.0);  // NaN was warned about and reset
    EXPECT_NE(sink->data().find("Unexpected noise density"), std::string::npos)
        << sink->data();
    sink->clear();

    // Inf path: infinite position density times a positive composition density.
    auto inf_data = BmmData<2>{};
    inf_data.noise_position_density = std::numeric_limits<double>::infinity();
    inf_data.noise_density = 42.0;
    inf_data.position_data.resize(2, 2);
    inf_data.position_data << 0.0, 0.1, 0.0, 0.0;
    inf_data.composition_data = {0, 0};
    inf_data.confidence = {1.0, 1.0};
    inf_data.assignment = {1, 1};
    Eigen::Vector2d mu = Eigen::Vector2d::Zero();
    CategoricalSmoothed comp_params(2, 1.0);
    comp_params.set_dense_counts({1.0f, 1.0f});
    inf_data.components.emplace_back(MvNormal<2>(mu, Eigen::Matrix2d::Identity()),
                                     comp_params, std::nullopt, 1);

    baysor::maximize(inf_data, /*freeze_composition=*/false, /*freeze_position=*/false);
    EXPECT_DOUBLE_EQ(inf_data.noise_density, 0.0);
    EXPECT_NE(sink->data().find("Unexpected noise density"), std::string::npos)
        << sink->data();
}

// ============================================================================
// Prior segmentation adjustment in the E-step
// ============================================================================

TEST(Cov2Bmm, UnderSegmentationBranchPenalizesForeignMainSegment) {
    auto no_prior = cov2_make_undersegmentation_data(/*with_prior_segmentation=*/false);
    auto with_prior = cov2_make_undersegmentation_data(/*with_prior_segmentation=*/true);

    baysor::expect_dirichlet_spatial(no_prior, /*stochastic=*/false);
    baysor::expect_dirichlet_spatial(with_prior, /*stochastic=*/false);

    // Without prior segmentation molecule 0 follows the nearest density
    // (component 2); with prior segmentation component 2's main segment (2)
    // differs from molecule 0's segment (1), so the under-segmentation
    // penalty flips molecule 0 to component 1.
    EXPECT_EQ(no_prior.assignment[0], 2);
    EXPECT_EQ(with_prior.assignment[0], 1);
    EXPECT_NE(no_prior.assignment[0], with_prior.assignment[0]);
}

// ============================================================================
// drop_unused_components
// ============================================================================

TEST(Cov2Bmm, DropUnusedComponentsResizesSegmentAndClusterBookkeeping) {
    auto data = cov2_make_disconnected_data();
    data.main_segment_per_cell = {1, 2, 3};
    data.cluster_per_cell = {1, 1, 2};

    // Component 3 has only 2 molecules; drop threshold 3 removes it.
    baysor::drop_unused_components(data, /*min_n_samples=*/3);

    EXPECT_EQ(data.n_components(), 2);
    EXPECT_EQ(data.assignment, (std::vector<int>{1, 1, 1, 2, 2, 2, 0, 0}));
    EXPECT_EQ(data.main_segment_per_cell, (std::vector<int>{1, 2}));
    EXPECT_EQ(data.cluster_per_cell, (std::vector<int>{1, 1}));

    // Nothing to drop: bookkeeping arrays stay untouched.
    auto stable = cov2_make_disconnected_data();
    stable.main_segment_per_cell = {1, 2, 3};
    stable.cluster_per_cell = {1, 1, 2};
    baysor::drop_unused_components(stable, /*min_n_samples=*/2);
    EXPECT_EQ(stable.n_components(), 3);
    EXPECT_EQ(stable.main_segment_per_cell, (std::vector<int>{1, 2, 3}));
    EXPECT_EQ(stable.cluster_per_cell, (std::vector<int>{1, 1, 2}));
}

// ============================================================================
// estimate_assignment_by_history
// ============================================================================

TEST(Cov2Bmm, EstimateAssignmentByHistoryMajorityVote) {
    auto data = cov2_make_disconnected_data();
    // Drop component 3 (guid 3) from the current component list so history
    // entries with guid 3 no longer match any current component.
    data.components.pop_back();
    data.assignment.resize(5, 1);
    data.position_data.conservativeResize(Eigen::NoChange, 5);
    data.composition_data.resize(5);
    data.confidence.resize(5);
    data.components[0].guid = 10;
    data.components[1].guid = 20;

    data.assignment_history = {
        {10, 10,  0, 99, 99},
        {10, 20,  0, 99, 99},
        {10, 10,  0, 20, 99},
    };

    auto [reassign, match_frac] = baysor::estimate_assignment_by_history(data);

    ASSERT_EQ(reassign.size(), 5u);
    ASSERT_EQ(match_frac.size(), 5u);
    EXPECT_EQ(reassign[0], 1);              // guid 10 in all three entries
    EXPECT_NEAR(match_frac[0], 1.0, 1e-12);
    EXPECT_EQ(reassign[1], 1);              // 10 wins 2:1 over 20
    EXPECT_NEAR(match_frac[1], 2.0 / 3.0, 1e-12);
    EXPECT_EQ(reassign[2], 0);              // noise (guid 0) majority
    EXPECT_NEAR(match_frac[2], 1.0, 1e-12);
    EXPECT_EQ(reassign[3], 2);              // dropped guid 99 ignored, 20 valid
    EXPECT_NEAR(match_frac[3], 1.0, 1e-12);
    EXPECT_EQ(reassign[4], 0);              // no valid history entry at all
    EXPECT_NEAR(match_frac[4], 0.0, 1e-12);
}

TEST(Cov2Bmm, EstimateAssignmentByHistoryFallsBackWithoutHistory) {
    auto data = cov2_make_disconnected_data();
    ASSERT_TRUE(data.assignment_history.empty());

    auto [reassign, match_frac] = baysor::estimate_assignment_by_history(data);

    EXPECT_EQ(reassign, data.assignment);
    ASSERT_EQ(match_frac.size(), data.assignment.size());
    for (double c : match_frac) {
        EXPECT_DOUBLE_EQ(c, 0.5);
    }
}

// ============================================================================
// Full bmm() loop with verbose diagnostics, tracing and refinement
// ============================================================================

TEST(Cov2Bmm, VerboseLoopWithHistoryAndRefineKeepsStableComponents) {
    baysor_test::GlobalRngGuard rng_guard;  // bmm() advances the global RNG
    auto data = cov2_make_disconnected_data();

    baysor::bmm(data,
                /*min_molecules_drop=*/2,
                /*n_iters=*/3,
                /*assignment_history_depth=*/3,
                /*verbose=*/true,
                /*component_split_step=*/3,
                /*refine=*/true,
                /*freeze_composition=*/false,
                /*freeze_position=*/false,
                /*freeze_components=*/false,
                /*tol=*/0.0,
                /*min_molecules_display=*/5);  // differs from drop threshold

    EXPECT_EQ(data.n_components(), 3);
    EXPECT_EQ(data.assignment, (std::vector<int>{1, 1, 1, 2, 2, 2, 3, 3}));
    ASSERT_EQ(data.assignment_history.size(), 3u);
    EXPECT_EQ(data.n_components_trace.size(), 4u);  // initial + 3 iterations
    // Refinement populated per-molecule confidence from the history vote.
    ASSERT_EQ(data.assignment_confidence.size(), 8u);
    for (double c : data.assignment_confidence) {
        EXPECT_DOUBLE_EQ(c, 1.0);
    }
}

TEST(Cov2Bmm, VerboseLoopWithDropOneOmitsRedundantDiagnosticColumns) {
    baysor_test::GlobalRngGuard rng_guard;  // bmm() advances the global RNG
    auto sink = std::make_shared<baysor_test::CapturingSink>();
    baysor_test::LoggerGuard logger(sink);
    auto data = cov2_make_disconnected_data();

    baysor::bmm(data,
                /*min_molecules_drop=*/1,
                /*n_iters=*/1,
                /*assignment_history_depth=*/0,
                /*verbose=*/true,
                /*component_split_step=*/3,
                /*refine=*/false,
                /*freeze_composition=*/false,
                /*freeze_position=*/false,
                /*freeze_components=*/false,
                /*tol=*/0.0,
                /*min_molecules_display=*/0);  // equals drop threshold

    EXPECT_EQ(data.n_components(), 3);
    EXPECT_TRUE(data.assignment_history.empty());
    EXPECT_EQ(data.n_components_trace.size(), 2u);  // initial + 1 iteration

    // Verbose diagnostics were captured: the per-iteration line must exist
    // with its base columns, while the >=drop / >=disp columns are omitted
    // because both thresholds coincide (min_molecules_drop == disp_thresh == 1).
    const std::string log = sink->data();
    EXPECT_NE(log.find("Iter"), std::string::npos) << log;
    EXPECT_NE(log.find("noise="), std::string::npos) << log;
    EXPECT_NE(log.find("total="), std::string::npos) << log;
    EXPECT_EQ(log.find("_segs="), std::string::npos) << log;
    EXPECT_EQ(log.find("_cells="), std::string::npos) << log;
}

TEST(Cov2Bmm, ThreeDimensionalLoopWithSegmentsAndRefineStaysStable) {
    baysor_test::GlobalRngGuard rng_guard;  // bmm() advances the global RNG
    auto data = cov2_make_3d_data();

    baysor::bmm(data,
                /*min_molecules_drop=*/2,
                /*n_iters=*/3,
                /*assignment_history_depth=*/2,
                /*verbose=*/true,
                /*component_split_step=*/3,
                /*refine=*/true,
                /*freeze_composition=*/false,
                /*freeze_position=*/false,
                /*freeze_components=*/false,
                /*tol=*/0.0,
                /*min_molecules_display=*/3);

    // Both spatial groups are internally connected and far apart, so the
    // segmentation must stay at two components with the original assignment.
    EXPECT_EQ(data.n_components(), 2);
    EXPECT_EQ(data.assignment,
              (std::vector<int>{1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2}));
    EXPECT_EQ(data.assignment_history.size(), 2u);  // trimmed to depth
    EXPECT_EQ(data.n_components_trace.size(), 4u);
    ASSERT_EQ(data.assignment_confidence.size(), 12u);
    for (double c : data.assignment_confidence) {
        EXPECT_DOUBLE_EQ(c, 1.0);
    }
    // History stores global GUIDs only (plus 0 for noise).
    for (const auto& row : data.assignment_history) {
        ASSERT_EQ(row.size(), 12u);
        for (int guid : row) {
            EXPECT_TRUE(guid == 0 || guid == 1 || guid == 2);
        }
    }
    // Prior-segmentation bookkeeping was refreshed during the loop.
    ASSERT_EQ(data.main_segment_per_cell.size(), 2u);
    EXPECT_EQ(data.main_segment_per_cell[0], 1);
    EXPECT_EQ(data.main_segment_per_cell[1], 2);
}

// ============================================================================
// tracing
// ============================================================================

TEST(Cov2Trace, AssignmentHistoryTrimsOldestEntriesBeyondDepth) {
    auto data = cov2_make_disconnected_data();
    data.assignment = {1, 2, 3, 0, 1, 2, 3, 1};

    baysor::trace_assignment_history(data, /*assignment_history_depth=*/2);
    baysor::trace_assignment_history(data, /*assignment_history_depth=*/2);
    EXPECT_EQ(data.assignment_history.size(), 2u);

    baysor::trace_assignment_history(data, /*assignment_history_depth=*/1);
    // Both stale entries are erased, only the newest snapshot remains.
    ASSERT_EQ(data.assignment_history.size(), 1u);
    EXPECT_EQ(data.assignment_history.back(),
              (std::vector<int>{1, 2, 3, 0, 1, 2, 3, 1}));

    // Non-positive depth records nothing.
    baysor::trace_assignment_history(data, /*assignment_history_depth=*/0);
    EXPECT_EQ(data.assignment_history.size(), 1u);
}

TEST(Cov2Trace, EstimateComponentLifespanHandlesUnbrokenAndBrokenStreaks) {
    EXPECT_TRUE(baysor::estimate_component_lifespan({}).empty());

    // Unbroken streaks: guid 1 present in all three snapshots, guid 3 only in
    // the last two (it disappears going backward).
    const std::vector<std::vector<int>> unbroken = {
        {1, 1, 2},
        {1, 1, 3},
        {1, 1, 3, 0},
    };
    auto life = baysor::estimate_component_lifespan(unbroken);
    ASSERT_EQ(life.size(), 2u);
    EXPECT_EQ(life.at(1), 3);
    EXPECT_EQ(life.at(3), 2);

    // Broken streaks: guid 2 vanishes in the oldest snapshot while guid 1 is
    // absent from the middle one, so neither survives the full history.
    const std::vector<std::vector<int>> broken = {
        {1, 0},
        {2, 0},
        {1, 1, 2, 0},
    };
    auto life_broken = baysor::estimate_component_lifespan(broken);
    ASSERT_EQ(life_broken.size(), 2u);
    EXPECT_EQ(life_broken.at(1), 1);
    EXPECT_EQ(life_broken.at(2), 2);
}

// ============================================================================
// history_analysis (unimplemented stub)
// ============================================================================

TEST(Cov2Hist, ReassignMoleculesWithHistoryReturnsEmptyForBothInstantiations) {
    auto data2 = cov2_make_disconnected_data();
    auto [assign2, conf2] = baysor::reassign_molecules_with_history(
        data2, /*n_stable_iters=*/5);
    EXPECT_TRUE(assign2.empty());
    EXPECT_TRUE(conf2.empty());

    auto data3 = cov2_make_3d_data();
    auto [assign3, conf3] = baysor::reassign_molecules_with_history(
        data3, /*n_stable_iters=*/5, /*min_molecules_per_cell=*/3,
        /*max_iters=*/50, /*outlier_confidence_threshold=*/0.5);
    EXPECT_TRUE(assign3.empty());
    EXPECT_TRUE(conf3.empty());
}

// ============================================================================
// compartment_segmentation (unimplemented stubs)
// ============================================================================

TEST(Cov2Comp, CompartmentStubsReturnEmptyAndLeaveGraphUntouched) {
    Eigen::MatrixXd pos(2, 4);
    pos << 0.0, 1.0, 2.0, 3.0,
           0.0, 0.0, 1.0, 1.0;
    const std::vector<int> genes = {1, 2, 1, 2};
    const std::vector<std::string> gene_names = {"GeneA", "GeneB"};
    const std::vector<std::string> nuclei_genes = {"GeneA"};
    const std::vector<std::string> cyto_genes = {"GeneB"};

    auto [probs, locked] = baysor::init_nuclei_cyto_compartments(
        pos, genes, gene_names, nuclei_genes, cyto_genes, /*scale=*/5.0);
    EXPECT_EQ(probs.rows(), 0);
    EXPECT_EQ(probs.cols(), 0);
    EXPECT_TRUE(locked.empty());

    const int edge_src[] = {0, 1};
    const int edge_dst[] = {1, 2};
    const double edge_wt[] = {1.0, 1.0};
    auto adj = AdjList::from_edge_list(edge_src, edge_dst, edge_wt, 2, 4);
    const auto adj_before = adj;
    const std::vector<double> confidence(4, 1.0);
    Eigen::MatrixXd assignment_probs = Eigen::MatrixXd::Zero(3, 4);
    const std::vector<bool> is_locked(4, false);

    auto result = baysor::segment_molecule_compartments(
        assignment_probs, is_locked, adj, confidence);
    EXPECT_TRUE(result.assignment.empty());
    EXPECT_EQ(result.assignment_probs.rows(), 0);
    EXPECT_EQ(result.assignment_probs.cols(), 0);
    EXPECT_TRUE(result.diffs.empty());
    EXPECT_TRUE(result.change_fracs.empty());

    const std::vector<double> nuclei_probs = {0.9, 0.8, 0.1, 0.2};
    const std::vector<double> cyto_probs = {0.1, 0.2, 0.9, 0.8};
    baysor::adjust_mrf_with_compartments(adj, nuclei_probs, cyto_probs);
    EXPECT_EQ(adj.indptr, adj_before.indptr);
    EXPECT_EQ(adj.indices, adj_before.indices);
    EXPECT_EQ(adj.weights, adj_before.weights);
}
