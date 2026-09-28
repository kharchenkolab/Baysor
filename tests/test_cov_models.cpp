// Coverage tests for models and distributions (COV-2):
//   src/processing/models/{adj_list,bmm_data,component}.cpp
//   src/processing/distributions/{mv_normal,categorical_smoothed}.cpp
//   plus the inline accessors in include/baysor/processing/{models,distributions}
#include <gtest/gtest.h>

#include "baysor/processing/distributions/categorical_smoothed.h"
#include "baysor/utils/general.h"
#include "baysor/processing/distributions/mv_normal.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/processing/models/bmm_data.h"
#include "baysor/processing/models/component.h"

#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using baysor::AdjList;
using baysor::BmmData;
using baysor::CategoricalSmoothed;
using baysor::Component;
using baysor::MvNormal;

Component<2> cov2_make_component2d(int guid, double prior_probability = 1.0,
                                   double confidence = 1.0) {
    Eigen::Vector2d mu = Eigen::Vector2d::Zero();
    CategoricalSmoothed comp_params(2, 1.0);
    comp_params.set_dense_counts({1.0f, 1.0f});
    Component<2> comp(MvNormal<2>(mu, Eigen::Matrix2d::Identity()), comp_params,
                      std::nullopt, guid);
    comp.prior_probability = prior_probability;
    comp.confidence = confidence;
    return comp;
}

Component<3> cov2_make_component3d(int guid) {
    Eigen::Vector3d mu = Eigen::Vector3d::Zero();
    CategoricalSmoothed comp_params(2, 1.0);
    comp_params.set_dense_counts({1.0f, 1.0f});
    Component<3> comp(MvNormal<3>(mu, Eigen::Matrix3d::Identity()), comp_params,
                      std::nullopt, guid);
    return comp;
}

} // namespace

// ============================================================================
// AdjList
// ============================================================================

TEST(Cov2Models, FromEdgeListEmptyEdgeListProducesZeroIndptr) {
    // n_edges == 0: all vertices are isolated.
    auto adj = AdjList::from_edge_list(nullptr, nullptr, nullptr, /*n_edges=*/0,
                                        /*n_verts=*/4);
    const std::vector<int32_t> expected(5, 0);
    EXPECT_EQ(adj.indptr, expected);
    EXPECT_EQ(adj.n_molecules(), 4);
    EXPECT_EQ(adj.nnz(), 0);
    for (int i = 0; i < 4; ++i) {
        EXPECT_EQ(adj.neighbor_count(i), 0);
    }

    // n_verts == 0: the edge list is ignored entirely.
    const int src[] = {0, 1};
    const int dst[] = {1, 0};
    const double wts[] = {1.0, 1.0};
    auto empty = AdjList::from_edge_list(src, dst, wts, /*n_edges=*/2,
                                          /*n_verts=*/0);
    EXPECT_EQ(empty.indptr, (std::vector<int32_t>{0}));
    EXPECT_EQ(empty.n_molecules(), 0);
}

TEST(Cov2Models, BmmDataAccessorsReportDimensions) {
    BmmData<2> data2;
    EXPECT_EQ(data2.n_molecules(), 0);
    EXPECT_EQ(data2.n_components(), 0);
    EXPECT_EQ(data2.n_genes(), 0);

    data2.position_data.resize(2, 5);
    CategoricalSmoothed params2(7, 1.0);
    data2.components.emplace_back(
        MvNormal<2>(Eigen::Vector2d::Zero(), Eigen::Matrix2d::Identity()),
        params2, std::nullopt, 1);
    EXPECT_EQ(data2.n_molecules(), 5);
    EXPECT_EQ(data2.n_components(), 1);
    EXPECT_EQ(data2.n_genes(), 7);

    BmmData<3> data3;
    EXPECT_EQ(data3.n_molecules(), 0);
    EXPECT_EQ(data3.n_components(), 0);
    EXPECT_EQ(data3.n_genes(), 0);

    data3.position_data.resize(3, 4);
    CategoricalSmoothed params3(3, 1.0);
    data3.components.emplace_back(
        MvNormal<3>(Eigen::Vector3d::Zero(), Eigen::Matrix3d::Identity()),
        params3, std::nullopt, 1);
    data3.components.emplace_back(
        MvNormal<3>(Eigen::Vector3d::Zero(), Eigen::Matrix3d::Identity()),
        params3, std::nullopt, 2);
    EXPECT_EQ(data3.n_molecules(), 4);
    EXPECT_EQ(data3.n_components(), 2);
    EXPECT_EQ(data3.n_genes(), 3);
}

TEST(Cov2Models, UpdateNMolsPerSegmentBreaksTiesByLargerSegment) {
    BmmData<2> data;
    data.position_data.resize(2, 3);
    data.position_data.setZero();
    data.composition_data = {0, 0, 0};
    data.confidence = {1.0, 1.0, 1.0};
    data.assignment = {1, 1, 1};
    data.max_component_guid = 1;
    data.components.push_back(cov2_make_component2d(1));

    // Segment 1 is fully inside the cell (1/1) and so is segment 2 (2/2):
    // equal fractions, but the larger segment must win the tie-break.
    data.segment_per_molecule = {1, 2, 2};
    data.n_molecules_per_segment = {1, 2};
    data.update_n_mols_per_segment();

    const auto& seg_map = data.components[0].n_molecules_per_segment;
    ASSERT_EQ(seg_map.size(), 2u);
    EXPECT_EQ(seg_map.at(1), 1);
    EXPECT_EQ(seg_map.at(2), 2);

    ASSERT_EQ(data.main_segment_per_cell.size(), 1u);
    EXPECT_EQ(data.main_segment_per_cell[0], 2);
}

// ============================================================================
// Component
// ============================================================================

TEST(Cov2Models, PdfPositionOnlyMatchesPdfWithMissingGene) {
    // 2D: shift the mean to (1, 2) so the evaluated density is not the mode.
    auto comp2 = cov2_make_component2d(/*guid=*/1, /*prior_probability=*/0.5);
    Eigen::Vector2d mu;
    mu << 1.0, 2.0;
    const Eigen::Matrix2d sigma = Eigen::Matrix2d::Identity() * 0.5;
    comp2.position_params = MvNormal<2>(mu, sigma);

    const double x2[2] = {1.5, 2.5};
    // MvNormal's divider matches the Julia parity quirk: 0.5*log((2*pi)^3*det)
    // for every dimensionality, hence (2*pi)^1.5 * sqrt(det) in the pdf.
    const double expected2 =
        0.5 * std::exp(-0.5) /
        (std::pow(2.0 * baysor::kPi, 1.5) * std::sqrt(0.25));  // prior * N(mu, 0.5 I)
    EXPECT_NEAR(comp2.pdf_position_only(x2), expected2, 1e-12);
    EXPECT_DOUBLE_EQ(comp2.pdf_position_only(x2), comp2.pdf(x2, /*gene=*/-1));

    // 3D at the mode of the standard normal.
    auto comp3 = cov2_make_component3d(/*guid=*/2);
    comp3.prior_probability = 0.25;
    const double x3[3] = {0.0, 0.0, 0.0};
    const double expected3 = 0.25 * std::pow(2.0 * baysor::kPi, -1.5);
    EXPECT_NEAR(comp3.pdf_position_only(x3), expected3, 1e-12);
    EXPECT_DOUBLE_EQ(comp3.pdf_position_only(x3), comp3.pdf(x3, /*gene=*/-1));
}

TEST(Cov2Models, Component2DContiguousMaximizeUsesIntegerQuantile) {
    auto comp = cov2_make_component2d(/*guid=*/1);

    double pos[] = {0.0, 0.0,   1.0, 2.0,   2.0, 4.0};
    const int genes[] = {0, 1, 0};
    double nuclei[] = {0.1, 0.5, 0.9};

    comp.maximize(pos, /*stride=*/2, genes, /*n_points=*/3,
                  nuclei, /*min_nuclei_frac=*/0.5,
                  /*freeze_position=*/false, /*freeze_composition=*/false);

    EXPECT_EQ(comp.n_samples, 3);
    // Position parameters use the nuclei probabilities as per-point weights
    // (floored at 0.01): (0*0.1 + 1*0.5 + 2*0.9) / 1.5 and
    // (0*0.1 + 2*0.5 + 4*0.9) / 1.5.
    EXPECT_NEAR(comp.position_params.mu(0), 2.3 / 1.5, 1e-12);
    EXPECT_NEAR(comp.position_params.mu(1), 4.6 / 1.5, 1e-12);
    // 1 - 0.5 = 0.5 quantile over three points lands exactly on the middle
    // value (integer index), which is what quantile_vec returns.
    EXPECT_DOUBLE_EQ(comp.confidence, 0.5);
    const auto dense = comp.composition_params.dense_counts();
    ASSERT_EQ(dense.size(), 2u);
    EXPECT_NEAR(dense[0], 2.0f, 1e-6);
    EXPECT_NEAR(dense[1], 1.0f, 1e-6);
}

TEST(Cov2Models, Component3DContiguousMaximizeKeepsConfidenceWithoutNuclei) {
    auto comp = cov2_make_component3d(/*guid=*/3);
    comp.confidence = 0.7;

    double pos[] = {0.0, 0.0, 0.0,
                    1.0, 1.0, 1.0,
                    2.0, 2.0, 2.0,
                    3.0, 3.0, 3.0};
    const int genes[] = {0, 0, 1, 1};

    comp.maximize(pos, /*stride=*/3, genes, /*n_points=*/4,
                  /*nuclei_probs=*/nullptr, /*min_nuclei_frac=*/0.1,
                  /*freeze_position=*/false, /*freeze_composition=*/false);

    EXPECT_EQ(comp.n_samples, 4);
    EXPECT_NEAR(comp.position_params.mu(0), 1.5, 1e-12);
    EXPECT_NEAR(comp.position_params.mu(1), 1.5, 1e-12);
    EXPECT_NEAR(comp.position_params.mu(2), 1.5, 1e-12);
    // No nuclei probabilities: confidence is left untouched.
    EXPECT_DOUBLE_EQ(comp.confidence, 0.7);
    const auto dense = comp.composition_params.dense_counts();
    ASSERT_EQ(dense.size(), 2u);
    EXPECT_NEAR(dense[0], 2.0f, 1e-6);
    EXPECT_NEAR(dense[1], 2.0f, 1e-6);
    // Degenerate (collinear) sample covariance is repaired to be usable.
    EXPECT_TRUE(std::isfinite(comp.position_params.sigma(0, 0)));
    EXPECT_GT(comp.position_params.sigma.determinant(), 0.0);
}

// ============================================================================
// MvNormal
// ============================================================================

TEST(Cov2Dist, DefaultConstructorsMatchJuliaParityNormaliser) {
    MvNormal<2> d2;
    EXPECT_TRUE(d2.mu.isZero());
    EXPECT_TRUE(d2.sigma.isIdentity());
    EXPECT_TRUE(d2.sigma_inv.isIdentity());

    // Julia parity quirk: norm_pdf_divider hardcodes (2*pi)^3 even in 2D,
    // so the default-constructed divider is 1.5*log(2*pi), not log(2*pi).
    EXPECT_NEAR(d2.pdf_divider, 1.5 * std::log(2.0 * baysor::kPi), 1e-12);
    const double x2[2] = {0.0, 0.0};
    EXPECT_NEAR(d2.pdf(x2), std::pow(2.0 * baysor::kPi, -1.5), 1e-12);
    EXPECT_NEAR(d2.log_pdf(x2), -1.5 * std::log(2.0 * baysor::kPi), 1e-12);

    MvNormal<3> d3;
    EXPECT_TRUE(d3.mu.isZero());
    EXPECT_TRUE(d3.sigma.isIdentity());
    EXPECT_TRUE(d3.sigma_inv.isIdentity());
    EXPECT_NEAR(d3.pdf_divider, 1.5 * std::log(2.0 * baysor::kPi), 1e-12);

    const double x3[3] = {0.0, 0.0, 0.0};
    EXPECT_NEAR(d3.pdf(x3), std::pow(2.0 * baysor::kPi, -1.5), 1e-12);
    EXPECT_NEAR(d3.log_pdf(x3), -1.5 * std::log(2.0 * baysor::kPi), 1e-12);

    // Off-mode density decays with the squared Mahalanobis distance.
    const double y2[2] = {1.0, 0.0};
    EXPECT_NEAR(d2.pdf(y2), std::exp(-0.5) * std::pow(2.0 * baysor::kPi, -1.5), 1e-12);
}

TEST(Cov2Dist, MaximizeWithTooFewPointsKeepsCovariance) {
    Eigen::Vector2d mu;
    mu << 5.0, 5.0;
    Eigen::Matrix2d sigma = Eigen::Vector2d(2.0, 3.0).asDiagonal();
    MvNormal<2> mv(mu, sigma);

    // Two points in 2D: not enough to estimate a covariance.
    double pts[] = {0.0, 0.0,   2.0, 2.0};
    mv.maximize(pts, /*n_points=*/2, /*stride=*/2);

    EXPECT_NEAR(mv.mu(0), 1.0, 1e-12);
    EXPECT_NEAR(mv.mu(1), 1.0, 1e-12);
    EXPECT_DOUBLE_EQ(mv.sigma(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(mv.sigma(1, 1), 3.0);
    EXPECT_DOUBLE_EQ(mv.sigma(0, 1), 0.0);
    // The cache still matches the retained covariance.
    EXPECT_DOUBLE_EQ(mv.sigma_inv(0, 0), 0.5);
    EXPECT_DOUBLE_EQ(mv.sigma_inv(1, 1), 1.0 / 3.0);
}

TEST(Cov2Dist, MaximizeWithoutShapePriorKeepsSampleCovariance) {
    MvNormal<2> mv;
    double pts[] = {0.0, 0.0,   1.0, 0.0,   0.0, 1.0,   1.0, 1.0};
    mv.maximize(pts, /*n_points=*/4, /*stride=*/2,
                /*center_probs=*/nullptr, /*shape_prior=*/nullptr,
                /*n_samples=*/-1);

    EXPECT_NEAR(mv.mu(0), 0.5, 1e-12);
    EXPECT_NEAR(mv.mu(1), 0.5, 1e-12);
    EXPECT_NEAR(mv.sigma(0, 0), 0.25, 1e-12);
    EXPECT_NEAR(mv.sigma(1, 1), 0.25, 1e-12);
    EXPECT_NEAR(mv.sigma(0, 1), 0.0, 1e-12);
    EXPECT_NEAR(mv.sigma_inv(0, 0), 4.0, 1e-12);
}

TEST(Cov2Dist, MaximizeIndexedWithoutShapePriorMatchesContiguousResult) {
    Eigen::MatrixXd pm(2, 4);
    pm << 0.0, 1.0, 0.0, 1.0,
          0.0, 0.0, 1.0, 1.0;
    // Process the points in reverse order; the moments are order-invariant.
    const int ids[] = {3, 2, 1, 0};

    MvNormal<2> indexed;
    indexed.maximize_indexed(pm, ids, /*n_points=*/4,
                             /*center_probs=*/nullptr, /*shape_prior=*/nullptr,
                             /*n_samples=*/-1);

    MvNormal<2> contiguous;
    double pts[] = {0.0, 0.0,   1.0, 0.0,   0.0, 1.0,   1.0, 1.0};
    contiguous.maximize(pts, /*n_points=*/4, /*stride=*/2);

    EXPECT_NEAR(indexed.mu(0), contiguous.mu(0), 1e-12);
    EXPECT_NEAR(indexed.mu(1), contiguous.mu(1), 1e-12);
    EXPECT_NEAR(indexed.sigma(0, 0), contiguous.sigma(0, 0), 1e-12);
    EXPECT_NEAR(indexed.sigma(1, 1), contiguous.sigma(1, 1), 1e-12);
    EXPECT_NEAR(indexed.sigma(0, 1), contiguous.sigma(0, 1), 1e-12);
    EXPECT_NEAR(indexed.mu(0), 0.5, 1e-12);
    EXPECT_NEAR(indexed.sigma(0, 0), 0.25, 1e-12);
}

TEST(Cov2Dist, MaximizeDegeneratePointsYieldsPositiveDefiniteCovariance) {
    MvNormal<2> mv;
    // All four points are identical: the sample covariance is exactly zero
    // and must be nudged to positive definiteness.
    double pts[] = {3.0, 3.0,   3.0, 3.0,   3.0, 3.0,   3.0, 3.0};
    mv.maximize(pts, /*n_points=*/4, /*stride=*/2);

    EXPECT_NEAR(mv.mu(0), 3.0, 1e-12);
    EXPECT_NEAR(mv.mu(1), 3.0, 1e-12);
    EXPECT_GT(mv.sigma(0, 0), 0.0);
    EXPECT_GT(mv.sigma(1, 1), 0.0);
    EXPECT_DOUBLE_EQ(mv.sigma(0, 1), 0.0);
    EXPECT_GT(mv.sigma.determinant(), 0.0);
    Eigen::LLT<Eigen::Matrix2d> llt(mv.sigma);
    EXPECT_EQ(llt.info(), Eigen::Success);
    EXPECT_TRUE(std::isfinite(mv.sigma_inv(0, 0)));
}

TEST(Cov2Dist, AdjustCovMatrixRepairsIndefiniteMatrix) {
    Eigen::Matrix2d sigma;
    sigma << 1.0, 2.0,
             2.0, 1.0;  // eigenvalues 3 and -1: not positive definite

    baysor::adjust_cov_matrix<2>(sigma);

    // Only the diagonal may change; the off-diagonal stays put while the
    // diagonal grows until LLT succeeds with a positive determinant.
    EXPECT_DOUBLE_EQ(sigma(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(sigma(1, 0), 2.0);
    EXPECT_GT(sigma(0, 0), 1.0);
    EXPECT_GT(sigma(1, 1), 1.0);
    EXPECT_GT(sigma.determinant(), 0.0);
    Eigen::LLT<Eigen::Matrix2d> llt(sigma);
    EXPECT_EQ(llt.info(), Eigen::Success);
}

// ============================================================================
// CategoricalSmoothed
// ============================================================================

TEST(Cov2Dist, SmoothedCategoricalUninformedPriorIsUniform) {
    CategoricalSmoothed fresh(5, 1.0);
    EXPECT_EQ(fresh.size(), 5);
    EXPECT_DOUBLE_EQ(fresh.pdf(0), 0.2);
    EXPECT_DOUBLE_EQ(fresh.pdf(4), 0.2);
    EXPECT_DOUBLE_EQ(fresh.pdf(-1), 1.0);  // missing gene is uninformative

    CategoricalSmoothed empty_panel(0, 1.0);
    EXPECT_EQ(empty_panel.size(), 0);
    EXPECT_DOUBLE_EQ(empty_panel.pdf(0), 1.0);
}

TEST(Cov2Dist, SmoothedCategoricalPdfWithoutSmoothingUsesRawFractions) {
    CategoricalSmoothed params(3, 1.0);
    params.set_dense_counts({2.0f, 0.0f, 1.0f});

    EXPECT_DOUBLE_EQ(params.pdf(0, /*use_smoothing=*/false), 2.0 / 3.0);
    EXPECT_DOUBLE_EQ(params.pdf(2, /*use_smoothing=*/false), 1.0 / 3.0);
    EXPECT_DOUBLE_EQ(params.pdf(1, /*use_smoothing=*/false), 0.0);

    // Laplace smoothing lifts only zero counts up to `smooth`.
    EXPECT_DOUBLE_EQ(params.pdf(1, /*use_smoothing=*/true), 1.0 / 4.0);
    EXPECT_DOUBLE_EQ(params.pdf(0, /*use_smoothing=*/true), 2.0 / 4.0);
}

TEST(Cov2Dist, SmoothedCategoricalSetUniformZeroClearsCounts) {
    CategoricalSmoothed params(4, 1.0);
    params.set_dense_counts({1.0f, 1.0f, 1.0f, 1.0f});
    ASSERT_GT(params.sum_counts, 0.0);

    params.set_uniform_counts(0.0f);
    EXPECT_DOUBLE_EQ(params.sum_counts, 0.0);
    EXPECT_EQ(params.n_genes, 0);
    EXPECT_DOUBLE_EQ(params.base_count, 0.0);
    EXPECT_DOUBLE_EQ(params.pdf(0), 0.25);  // back to the uniform prior
}

TEST(Cov2Dist, SmoothedCategoricalSetDenseCountsRejectsMismatchedSize) {
    CategoricalSmoothed params(3, 1.0);
    try {
        params.set_dense_counts({1.0f, 2.0f});
        FAIL() << "expected std::runtime_error for size mismatch";
    } catch (const std::runtime_error& e) {
        EXPECT_STREQ(e.what(),
                     "CategoricalSmoothed::set_dense_counts size mismatch");
    }
    // State is untouched by the rejected call.
    EXPECT_DOUBLE_EQ(params.sum_counts, 0.0);
    EXPECT_EQ(params.n_genes, 0);
}

TEST(Cov2Dist, SmoothedCategoricalResetClearsState) {
    CategoricalSmoothed params(3, 1.0);
    params.set_dense_counts({1.0f, 2.0f, 3.0f});
    params.set_uniform_counts(0.5f);
    ASSERT_EQ(params.n_genes, 3);

    params.reset();
    EXPECT_DOUBLE_EQ(params.sum_counts, 0.0);
    EXPECT_EQ(params.n_genes, 0);
    EXPECT_DOUBLE_EQ(params.base_count, 0.0);
    EXPECT_TRUE(params.gene_ids.empty());
    EXPECT_TRUE(params.counts.empty());
    EXPECT_DOUBLE_EQ(params.pdf(2), 1.0 / 3.0);
}
