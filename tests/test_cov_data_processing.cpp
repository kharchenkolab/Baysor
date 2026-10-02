// Coverage tests for the data-processing modules (COV-3): initialization,
// umap_wrappers, boundary_estimation, noise_estimation, triangulation,
// neighborhood_composition, utils, convex_hull.

#include <gtest/gtest.h>

#include "baysor/data_loading/data.h"
#include "baysor/processing/data_processing/boundary_estimation.h"
#include "baysor/processing/data_processing/boundary_estimation_internal.h"
#include "baysor/processing/data_processing/initialization.h"
#include "baysor/processing/data_processing/neighborhood_composition.h"
#include "baysor/processing/data_processing/noise_estimation.h"
#include "baysor/processing/data_processing/triangulation.h"
#include "baysor/processing/data_processing/umap_wrappers.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/processing/models/bmm_data.h"
#include "baysor/processing/utils/convex_hull.h"
#include "baysor/processing/utils/utils.h"

#include "test_cov_helpers.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace {

using baysor::AdjList;
using baysor::AdjacencyType;

// ---------------------------------------------------------------------------
// Small local helpers
// ---------------------------------------------------------------------------

AdjList cov3_chain_adj(int n_molecules) {
    std::vector<int> src;
    std::vector<int> dst;
    std::vector<double> wts;
    for (int i = 0; i < n_molecules - 1; ++i) {
        src.push_back(i);
        dst.push_back(i + 1);
        wts.push_back(1.0);
    }
    return AdjList::from_edge_list(src.data(), dst.data(), wts.data(),
                                   static_cast<int>(src.size()), n_molecules);
}

// Shoelace signed area (positive = counter-clockwise).
double cov3_signed_area(const Eigen::MatrixXd& poly) {
    double s = 0.0;
    const int n = static_cast<int>(poly.cols());
    for (int i = 0; i < n; ++i) {
        const int j = (i + 1) % n;
        s += poly(0, i) * poly(1, j) - poly(1, i) * poly(0, j);
    }
    return s / 2.0;
}

// Canonical (sorted) list of edges for comparing edge sets.
std::vector<std::pair<int, int>> cov3_sorted_edges(
    const std::vector<std::pair<int, int>>& edges
) {
    std::vector<std::pair<int, int>> out;
    out.reserve(edges.size());
    for (const auto& e : edges) {
        out.emplace_back(std::min(e.first, e.second), std::max(e.first, e.second));
    }
    std::sort(out.begin(), out.end());
    return out;
}

// Undirected edge set of an adjacency result.
std::set<std::pair<int, int>> edge_set(const baysor::AdjacencyResult& r) {
    std::set<std::pair<int, int>> out;
    for (size_t i = 0; i < r.edge_src.size(); ++i) out.insert(std::minmax(r.edge_src[i], r.edge_dst[i]));
    return out;
}

} // namespace

// ============================================================================
// convex_hull.cpp
// ============================================================================

TEST(Cov3Data_ConvexHull, EmptyInputReturnsEmptyMatrix) {
    Eigen::MatrixXd pts(2, 0);
    const Eigen::MatrixXd hull = baysor::convex_hull(pts);
    EXPECT_EQ(hull.rows(), 0);
    EXPECT_EQ(hull.cols(), 0);
}

TEST(Cov3Data_ConvexHull, FewerThanThreePointsReturnedAsIs) {
    Eigen::MatrixXd one(2, 1);
    one << 3.0, 4.0;
    const Eigen::MatrixXd hull1 = baysor::convex_hull(one);
    ASSERT_EQ(hull1.rows(), 2);
    ASSERT_EQ(hull1.cols(), 1);
    EXPECT_DOUBLE_EQ(hull1(0, 0), 3.0);
    EXPECT_DOUBLE_EQ(hull1(1, 0), 4.0);

    Eigen::MatrixXd two(2, 2);
    two << 0.0, 1.0,
           0.0, 1.0;
    const Eigen::MatrixXd hull2 = baysor::convex_hull(two);
    ASSERT_EQ(hull2.rows(), 2);
    ASSERT_EQ(hull2.cols(), 2);
    EXPECT_TRUE((hull2 - two).isZero(0.0));
}

TEST(Cov3Data_ConvexHull, SquareHullKeepsExactlyTheCornersClockwise) {
    // A square with two interior points, and a square with every corner duplicated.
    Eigen::MatrixXd interior(2, 6), duplicated(2, 8);
    interior << 0.0, 4.0, 4.0, 0.0, 1.0, 2.0,
                0.0, 0.0, 4.0, 4.0, 1.0, 3.0;
    duplicated << 0.0, 4.0, 4.0, 0.0, 0.0, 4.0, 4.0, 0.0,
                  0.0, 0.0, 4.0, 4.0, 0.0, 0.0, 4.0, 4.0;
    const std::vector<std::pair<double, double>> corners = {
        {0.0, 0.0}, {0.0, 4.0}, {4.0, 0.0}, {4.0, 4.0}};

    for (const Eigen::MatrixXd& pts : {interior, duplicated}) {
        const Eigen::MatrixXd hull = baysor::convex_hull(pts);
        ASSERT_EQ(hull.rows(), 2);
        ASSERT_EQ(hull.cols(), 4);
        std::vector<std::pair<double, double>> got;
        for (int i = 0; i < hull.cols(); ++i) got.emplace_back(hull(0, i), hull(1, i));
        std::sort(got.begin(), got.end());
        EXPECT_EQ(got, corners);
        EXPECT_NEAR(cov3_signed_area(hull), -16.0, 1e-12);  // clockwise
        EXPECT_NEAR(baysor::polygon_area(hull), 16.0, 1e-12);
    }
}

TEST(Cov3Data_ConvexHull, CollinearPointsReduceToEndpoints) {
    Eigen::MatrixXd pts(2, 3);
    pts << 0.0, 1.0, 2.0,
           0.0, 0.0, 0.0;

    const Eigen::MatrixXd hull = baysor::convex_hull(pts);
    ASSERT_EQ(hull.rows(), 2);
    // The middle point is popped (cross == 0), leaving the two endpoints.
    ASSERT_EQ(hull.cols(), 2);
    // Degenerate polygon: zero area.
    EXPECT_DOUBLE_EQ(baysor::polygon_area(hull), 0.0);

    std::vector<double> xs = {hull(0, 0), hull(0, 1)};
    std::sort(xs.begin(), xs.end());
    EXPECT_DOUBLE_EQ(xs[0], 0.0);
    EXPECT_DOUBLE_EQ(xs[1], 2.0);
}

TEST(Cov3Data_ConvexHull, PolygonAreaShoelaceBasics) {
    Eigen::MatrixXd tri(2, 3);
    tri << 0.0, 1.0, 0.0,
           0.0, 0.0, 1.0;
    EXPECT_NEAR(baysor::polygon_area(tri), 0.5, 1e-12);

    Eigen::MatrixXd rect(2, 4);
    rect << 0.0, 3.0, 3.0, 0.0,
            0.0, 0.0, 2.0, 2.0;
    EXPECT_NEAR(baysor::polygon_area(rect), 6.0, 1e-12);

    // Orientation must not matter (absolute value).
    Eigen::MatrixXd rev(2, 4);
    rev << 0.0, 3.0, 3.0, 0.0,
           2.0, 2.0, 0.0, 0.0;
    EXPECT_NEAR(baysor::polygon_area(rev), 6.0, 1e-12);

    // Fewer than 3 vertices: area is 0 by definition.
    Eigen::MatrixXd seg(2, 2);
    seg << 0.0, 1.0,
           0.0, 0.0;
    EXPECT_DOUBLE_EQ(baysor::polygon_area(seg), 0.0);
    Eigen::MatrixXd empty(2, 0);
    EXPECT_DOUBLE_EQ(baysor::polygon_area(empty), 0.0);
}

// ============================================================================
// umap_wrappers.cpp
// ============================================================================

TEST(Cov3Data_Umap, PrecomputedEmptyMatrixReturnsEmptyEmbedding) {
    Eigen::MatrixXd dist(0, 0);
    const Eigen::MatrixXd emb = baysor::umap_embed_precomputed(dist, 3, 5, 10, 1);
    EXPECT_EQ(emb.rows(), 3);
    EXPECT_EQ(emb.cols(), 0);
}

TEST(Cov3Data_Umap, PrecomputedRunsDeterministicallyAndKeepsSize) {
    // Two tight clusters in 2D -> symmetric distance matrix.
    constexpr int n = 40;
    Eigen::MatrixXd pts(2, n);
    for (int i = 0; i < n; ++i) {
        const double base = (i < n / 2) ? 0.0 : 50.0;
        pts(0, i) = base + static_cast<double>(i % 5);
        pts(1, i) = base + static_cast<double>((i / 5) % 5);
    }
    Eigen::MatrixXd dist(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            dist(i, j) = (pts.col(i) - pts.col(j)).norm();
        }
    }

    const Eigen::MatrixXd emb1 = baysor::umap_embed_precomputed(dist, 2, 5, 30, 7);
    const Eigen::MatrixXd emb2 = baysor::umap_embed_precomputed(dist, 2, 5, 30, 7);

    ASSERT_EQ(emb1.rows(), 2);
    ASSERT_EQ(emb1.cols(), n);
    ASSERT_EQ(emb2.rows(), 2);
    ASSERT_EQ(emb2.cols(), n);
    EXPECT_TRUE(emb1.allFinite());
    // Same seed must reproduce the same embedding bit-for-bit.
    EXPECT_TRUE((emb1 - emb2).isZero(0.0));
    // The embedding must not be degenerate (all points identical).
    EXPECT_GT(emb1.row(0).maxCoeff() - emb1.row(0).minCoeff(), 1e-6);
}

// ============================================================================
// utils.cpp
// ============================================================================

TEST(Cov3Data_Utils, KnnParallelReturnsEmptyOnEmptyInputs) {
    Eigen::MatrixXd empty_tree(2, 0);
    Eigen::MatrixXd queries(2, 3);
    queries.setZero();
    auto r1 = baysor::knn_parallel(empty_tree, queries, 3, true);
    EXPECT_TRUE(r1.indices.empty());
    EXPECT_TRUE(r1.distances.empty());
    EXPECT_EQ(r1.n, 0);
    EXPECT_EQ(r1.k, 0);

    Eigen::MatrixXd tree(2, 3);
    tree.setZero();
    Eigen::MatrixXd empty_query(2, 0);
    auto r2 = baysor::knn_parallel(tree, empty_query, 2, true);
    EXPECT_TRUE(r2.indices.empty());
    EXPECT_EQ(r2.n, 0);

    // k <= 0 is rejected as well.
    auto r3 = baysor::knn_parallel(tree, queries, 0, true);
    EXPECT_TRUE(r3.indices.empty());
    EXPECT_EQ(r3.n, 0);
}

TEST(Cov3Data_Utils, KnnParallelRowsAreContiguousAndKIsClamped) {
    Eigen::MatrixXd pts(2, 20);
    for (int i = 0; i < 20; ++i) {
        pts.col(i) << static_cast<double>(i % 5), static_cast<double>(i / 5);
    }
    auto result = baysor::knn_parallel(pts, pts, 4, true);
    ASSERT_EQ(result.n, 20);
    ASSERT_EQ(result.k, 4);
    ASSERT_EQ(result.indices.size(), 80u);
    ASSERT_EQ(result.distances.size(), 80u);
    EXPECT_EQ(result.idx_row(19), result.indices.data() + 76);
    EXPECT_EQ(result.dist_row(19), result.distances.data() + 76);
    EXPECT_EQ(baysor::knn_parallel(pts, pts, 50, true).k, 20);
}

TEST(Cov3Data_Utils, KnnKthDistancesMatchesFullKnnResult) {
    // Deterministic pseudo-random points with duplicated columns (distance
    // ties): the streamed kth-distance result must be bitwise equal to the
    // row of the full knn_parallel result for every k/kth combination.
    Eigen::MatrixXd pts(2, 300);
    std::uint64_t s = 987654321;
    auto next = [&s]() {
        s = s * 6364136223846793005ULL + 1442655040888963407ULL;
        return static_cast<double>((s >> 16) & 0xFFFFFF) / 16777216.0;
    };
    for (int i = 0; i < 300; ++i) {
        double v = next() * 50.0;
        pts.col(i) << v, next() * 50.0;
    }
    for (int dup = 0; dup < 30; ++dup) {
        pts.col(dup * 10) = pts.col(dup * 10 + 1);  // exact duplicate pair
    }

    for (int k : {1, 2, 7, 40, 300, 400}) {
        auto full = baysor::knn_parallel(pts, pts, k, /*sorted=*/true);
        const int k_eff = std::min(k, 300);
        for (int kth : {0, 1, k_eff / 2, k_eff - 1, k_eff, 2 * k_eff}) {
            auto streamed = baysor::knn_kth_distances(pts, k, kth);
            ASSERT_EQ(static_cast<int>(streamed.size()), 300) << "k=" << k;
            const int kth_eff = std::min(std::max(kth, 0), k_eff - 1);
            for (int i = 0; i < 300; ++i) {
                ASSERT_EQ(streamed[i], full.dist_row(i)[kth_eff])
                    << "k=" << k << " kth=" << kth << " i=" << i;
            }
        }
    }

    // Degenerate inputs mirror knn_parallel's empty result.
    Eigen::MatrixXd no_points(2, 0);
    EXPECT_TRUE(baysor::knn_kth_distances(no_points, 3, 1).empty());
    EXPECT_TRUE(baysor::knn_kth_distances(pts, 0, 0).empty());
}

// ============================================================================
// triangulation.cpp
// ============================================================================

TEST(Cov3Data_Triangulation, NormalizePointsJittersDuplicates) {
    // Draws from the global RNG to jitter duplicates: restore it afterwards.
    baysor_test::GlobalRngGuard rng_guard;
    Eigen::MatrixXd pts(2, 4);
    pts.col(0) << 5.0, 5.0;
    pts.col(1) << 5.0, 5.0;   // exact duplicate of col 0
    pts.col(2) << 6.0, 5.0;
    pts.col(3) << 5.0, 6.0;

    const Eigen::MatrixXd out = baysor::normalize_points(pts);
    ASSERT_EQ(out.rows(), 2);
    ASSERT_EQ(out.cols(), 4);

    // Without jitter both duplicates would map to exactly the same value;
    // normalize_points perturbs them by < 1e-4 to break KNN self-loops.
    const double dx = out(0, 0) - out(0, 1);
    const double dy = out(1, 0) - out(1, 1);
    const double delta = std::sqrt(dx * dx + dy * dy);
    EXPECT_GT(delta, 0.0);
    EXPECT_LT(delta, 1e-4);

    // Non-duplicate points keep their relative placement inside [1.01, ~1.9].
    EXPECT_GT(out.minCoeff(), 1.0);
    EXPECT_LT(out.maxCoeff(), 2.0);
    // Point 2 stays to the right of the duplicated pair.
    EXPECT_GT(out(0, 2), out(0, 0));
}

TEST(Cov3Data_Triangulation, FilterLongEdgesDropsEverythingWhenMadIsZero) {
    // A perfect square: 4 equal sides + 1 diagonal => 5 (odd) edges.
    // All side log-distances equal the median, so with MAD == 0 the strict
    // `< median` filter drops every edge.
    Eigen::MatrixXd pts(2, 4);
    pts << 0.0, 1.0, 0.0, 1.0,
           0.0, 0.0, 1.0, 1.0;

    auto unfiltered = baysor::adjacency_list(pts, /*filter=*/false);
    EXPECT_EQ(static_cast<int>(unfiltered.edge_src.size()), 5);

    auto filtered = baysor::adjacency_list(pts, /*filter=*/true);
    EXPECT_EQ(static_cast<int>(filtered.edge_src.size()), 0);
    EXPECT_TRUE(filtered.edge_dists.empty());
}

TEST(Cov3Data_Triangulation, FilterKeepsASubsetOfEdgesForIrregularPoints) {
    // 5 irregular points: an odd number of Delaunay edges (7), non-zero MAD.
    Eigen::MatrixXd pts(2, 5);
    pts << 0.0, 1.0, 0.2, 1.1, 0.5,
           0.0, 0.1, 1.0, 0.9, 1.7;

    auto unfiltered = baysor::adjacency_list(pts, /*filter=*/false);
    auto filtered = baysor::adjacency_list(pts, /*filter=*/true);

    const int n_un = static_cast<int>(unfiltered.edge_src.size());
    const int n_fi = static_cast<int>(filtered.edge_src.size());
    EXPECT_EQ(n_un, 7);
    EXPECT_LE(n_fi, n_un);
    EXPECT_GT(n_fi, 0);

    // Every kept edge must come from the unfiltered set with the same length.
    for (int i = 0; i < n_fi; ++i) {
        const int a = filtered.edge_src[i];
        const int b = filtered.edge_dst[i];
        bool found = false;
        for (int j = 0; j < n_un; ++j) {
            const int ua = unfiltered.edge_src[j];
            const int ub = unfiltered.edge_dst[j];
            if (std::min(a, b) == std::min(ua, ub) &&
                std::max(a, b) == std::max(ua, ub)) {
                found = true;
                EXPECT_NEAR(filtered.edge_dists[i], unfiltered.edge_dists[j], 1e-12);
                break;
            }
        }
        EXPECT_TRUE(found) << "edge " << i << " not present in unfiltered result";
    }
}

TEST(Cov3Data_Triangulation, TriangulationTypeIsCoercedToKnnIn3D) {
    Eigen::MatrixXd pts(3, 6);
    pts << 0.0, 1.0, 0.0, 0.0, 1.0, 1.0,
           0.0, 0.0, 1.0, 0.0, 1.0, 0.0,
           0.0, 0.0, 0.0, 1.0, 0.0, 1.0;

    auto knn = baysor::adjacency_list(pts, /*filter=*/false, 2.0, 3, AdjacencyType::Knn);
    auto coerced = baysor::adjacency_list(pts, /*filter=*/false, 2.0, 3,
                                          AdjacencyType::Triangulation);

    // In 3D the Triangulation request is silently downgraded to KNN.
    EXPECT_EQ(coerced.edge_src, knn.edge_src);
    EXPECT_EQ(coerced.edge_dst, knn.edge_dst);
    ASSERT_EQ(coerced.edge_dists.size(), knn.edge_dists.size());
    for (size_t i = 0; i < knn.edge_dists.size(); ++i) {
        EXPECT_DOUBLE_EQ(coerced.edge_dists[i], knn.edge_dists[i]);
    }
}

// The Triangulation+KNN merge is the sorted union of both sources.
TEST(Cov3Data_Triangulation, BothTypeIsSortedUnionOfKnnAndTriangulationEdges) {
    Eigen::MatrixXd pts(2, 6);
    pts << 0.0, 1.0, 2.0, 0.0, 1.0, 2.0,
           0.0, 0.2, 0.1, 1.0, 1.2, 0.9;

    auto knn = baysor::adjacency_list(pts, /*filter=*/false, 2.0, 2, AdjacencyType::Knn);
    auto tri = baysor::adjacency_list(pts, /*filter=*/false, 2.0, 2, AdjacencyType::Triangulation);
    auto both = baysor::adjacency_list(pts, /*filter=*/false, 2.0, 2, AdjacencyType::Both);

    EXPECT_GT(knn.edge_src.size(), 0u);
    EXPECT_GT(tri.edge_src.size(), 0u);

    std::vector<std::pair<int, int>> merged;
    for (size_t i = 0; i < both.edge_src.size(); ++i) merged.push_back(std::minmax(both.edge_src[i], both.edge_dst[i]));
    EXPECT_EQ(std::adjacent_find(merged.begin(), merged.end(), std::greater_equal<>()), merged.end())
        << "merged edges not strictly increasing";

    auto union_set = edge_set(knn);
    const auto tri_set = edge_set(tri);
    union_set.insert(tri_set.begin(), tri_set.end());
    EXPECT_EQ(edge_set(both), union_set);
}

// ============================================================================
// noise_estimation.cpp
// ============================================================================

TEST(Cov3Data_Noise, FitHardAssignsExtremeOutlierToNoise) {
    // One distance that clearly exceeds q90 + 3 * init_std.
    const std::vector<double> edge_lengths = {1.0, 1.0, 1.0, 1.0, 1.0, 1000.0};
    auto adj = cov3_chain_adj(6);

    auto result = baysor::fit_noise_probabilities(edge_lengths, adj, nullptr,
                                                  /*max_iters=*/200, /*tol=*/0.005,
                                                  /*verbose=*/false);

    ASSERT_EQ(result.assignment.size(), 6u);
    EXPECT_EQ(result.assignment[5], 2);           // outlier -> noise
    for (int i = 0; i < 5; ++i) EXPECT_EQ(result.assignment[i], 1);
    ASSERT_EQ(result.assignment_probs.rows(), 6);
    // Outlier row is hard-reassigned to noise probability 1 (post-normalization).
    EXPECT_GT(result.assignment_probs(5, 1), 0.99);
    EXPECT_LT(result.assignment_probs(5, 0), 0.01);
    // Signal component must keep the lower mean.
    EXPECT_LT(result.signal_mu, result.noise_mu);
    EXPECT_TRUE(std::isfinite(result.signal_sigma));
    EXPECT_TRUE(std::isfinite(result.noise_sigma));
}

TEST(Cov3Data_Noise, DegenerateEqualInputsStillYieldOrderedSignalAndNoise) {
    // All distances identical (and negative): both initial Gaussians collapse
    // to zero width, the zero-weight M-step fallback sets the first component
    // mean to 0. The final guard must still order signal_mu <= noise_mu.
    const std::vector<double> edge_lengths(6, -5.0);
    auto adj = cov3_chain_adj(6);

    auto result = baysor::fit_noise_probabilities(edge_lengths, adj, nullptr,
                                                  /*max_iters=*/100, /*tol=*/0.005,
                                                  /*verbose=*/false);

    EXPECT_LE(result.signal_mu, result.noise_mu);
    EXPECT_TRUE(std::isfinite(result.signal_mu));
    EXPECT_TRUE(std::isfinite(result.noise_mu));
    EXPECT_TRUE(std::isfinite(result.signal_sigma));
    EXPECT_TRUE(std::isfinite(result.noise_sigma));
    ASSERT_FALSE(result.diffs.empty());
    EXPECT_LE(result.diffs.back(), 0.005);        // the EM loop reported convergence
    ASSERT_EQ(result.assignment.size(), 6u);
    EXPECT_EQ(result.assignment[0], 1);
    // Both component densities underflow to zero for this input, so the row
    // normalization is skipped; do not pin the resulting raw values (they are
    // an implementation detail), only require finite, in-range entries.
    for (int i = 0; i < 6; ++i) {
        EXPECT_TRUE(std::isfinite(result.assignment_probs(i, 0)));
        EXPECT_TRUE(std::isfinite(result.assignment_probs(i, 1)));
        EXPECT_GE(result.assignment_probs(i, 0), 0.0);
        EXPECT_GE(result.assignment_probs(i, 1), 0.0);
    }
}

TEST(Cov3Data_Noise, EstimateConfidenceDetailsUsesDefaultNnId) {
    // 20 molecules on a 5x4 grid, one gene.
    baysor::MoleculeData data;
    for (int i = 0; i < 20; ++i) {
        data.x.push_back(static_cast<double>(i % 5));
        data.y.push_back(static_cast<double>(i / 5));
        data.gene.push_back(1);
    }
    data.gene_names = {"GeneA"};

    auto details = baysor::estimate_confidence_details(data, /*nn_id=*/0);

    // nn_id <= 0 falls back to max(n_genes / 10, 10).
    EXPECT_EQ(details.nn_id, 10);
    EXPECT_EQ(details.edge_lengths.size(), 20u);
    ASSERT_EQ(details.fit_result.assignment_probs.rows(), 20);
    EXPECT_EQ(details.fit_result.assignment_probs.cols(), 2);
}

TEST(Cov3Data_Noise, EstimateConfidenceDetailsFloorsAtPriorConfidence) {
    // First 10 molecules carry a prior segment; the last 10 do not.
    baysor::MoleculeData data;
    for (int i = 0; i < 20; ++i) {
        data.x.push_back(static_cast<double>(i % 5));
        data.y.push_back(static_cast<double>(i / 5));
        data.gene.push_back(1 + (i % 2));
    }
    data.gene_names = {"GeneA", "GeneB"};
    data.prior_segmentation.assign(20, 0);
    for (int i = 0; i < 10; ++i) data.prior_segmentation[i] = 1;

    const double prior_confidence = 0.5;
    auto details = baysor::estimate_confidence_details(data, /*nn_id=*/6, prior_confidence);

    ASSERT_EQ(details.nn_id, 6);
    ASSERT_EQ(details.fit_result.assignment_probs.rows(), 20);
    // Molecules inside the prior segment get a signal-probability floor of
    // prior_confidence^2 = 0.25; unassigned molecules get no floor.
    const double floor = prior_confidence * prior_confidence;
    for (int i = 0; i < 10; ++i) {
        EXPECT_GE(details.fit_result.assignment_probs(i, 0), floor - 1e-12)
            << "molecule " << i;
    }
    for (int i = 10; i < 20; ++i) {
        // Unassigned molecules get no floor; assert the real value: the
        // per-molecule row must be a proper distribution over the two
        // components (a range check [0, 1] alone would pass trivially).
        EXPECT_NEAR(details.fit_result.assignment_probs(i, 0) +
                        details.fit_result.assignment_probs(i, 1),
                    1.0, 1e-9) << "molecule " << i;
    }
}

// ============================================================================
// initialization.cpp
// ============================================================================

TEST(Cov3Data_Initialization, BuildMoleculeGraphWithSingleMoleculeHasNoEdges) {
    baysor::MoleculeData data;
    data.x = {1.0};
    data.y = {2.0};
    data.gene = {1};
    data.gene_names = {"A"};

    auto adj = baysor::build_molecule_graph(data);
    EXPECT_EQ(adj.n_molecules(), 1);
    EXPECT_EQ(adj.nnz(), 0);
    ASSERT_EQ(adj.indptr.size(), 2u);
    EXPECT_EQ(adj.indptr[0], 0);
    EXPECT_EQ(adj.indptr[1], 0);
}

TEST(Cov3Data_Initialization, CellCenters2DUsesIsotropicCovarianceWithConfidences) {
    // 6 molecules near the origin, 4 near (10,10); last one is low-confidence.
    Eigen::MatrixXd pos(2, 10);
    const double ax[6] = {1.0, -1.0, 0.0, 0.0, 0.5, -0.5};
    const double ay[6] = {0.0, 0.0, 1.0, -1.0, 0.0, 0.0};
    for (int i = 0; i < 6; ++i) {
        pos(0, i) = ax[i];
        pos(1, i) = ay[i];
    }
    for (int i = 6; i < 10; ++i) {
        pos(0, i) = 10.0 + (i - 6);
        pos(1, i) = 10.0 + (i - 6);
    }
    std::vector<double> confidences(10, 0.9);
    confidences[9] = 0.1;  // excluded from center selection (< 0.25)

    const double scale = 0.4;
    auto init = baysor::cell_centers_uniformly<2>(pos, 3, &confidences, scale);

    ASSERT_EQ(init.centers.rows(), 3);
    ASSERT_EQ(init.centers.cols(), 2);
    EXPECT_EQ(init.covs.size(), 3u);
    ASSERT_EQ(init.assignment.size(), 10u);

    // scale > 0 => every covariance is exactly scale^2 * I.
    const Eigen::Matrix2d expected_cov = Eigen::Matrix2d::Identity() * scale * scale;
    for (const auto& cov : init.covs) {
        EXPECT_TRUE(cov.isApprox(expected_cov, 1e-15));
    }

    // Centers: one low-sum molecule from group A, another A molecule, and the
    // max-sum molecule of group B (index 8 — index 9 has low confidence).
    EXPECT_EQ(init.assignment[8], 3);
    for (int i = 6; i < 10; ++i) EXPECT_EQ(init.assignment[i], 3);
    for (int i = 0; i < 6; ++i) EXPECT_GE(init.assignment[i], 1);
    for (int i = 0; i < 6; ++i) EXPECT_LE(init.assignment[i], 3);
}

TEST(Cov3Data_Initialization, CellCenters2DComputesSampleCovariancesWithFallback) {
    // Group A: 2 points near the origin (invalid for N=2 -> identity
    // fallback). Group B: 3 points near (10,10) with a non-trivial sample
    // covariance. The coordinate sums interleave so that the evenly-spaced
    // center selection picks one center per group: centers are index 0
    // (sum 0) and index 2 (sum 20).
    Eigen::MatrixXd pos(2, 5);
    pos.col(0) << 0.0, 0.0;
    pos.col(1) << 0.4, 0.0;
    pos.col(2) << 10.0, 10.0;
    pos.col(3) << 10.4, 10.0;
    pos.col(4) << 10.0, 10.4;

    auto init = baysor::cell_centers_uniformly<2>(pos, 2, /*confidences=*/nullptr,
                                                  /*scale=*/-1.0);

    ASSERT_EQ(init.centers.rows(), 2);
    EXPECT_EQ(init.covs.size(), 2u);
    ASSERT_EQ(init.assignment.size(), 5u);
    EXPECT_EQ(init.assignment[0], 1);
    EXPECT_EQ(init.assignment[1], 1);
    EXPECT_EQ(init.assignment[2], 2);
    EXPECT_EQ(init.assignment[3], 2);
    EXPECT_EQ(init.assignment[4], 2);

    // Invalid cluster (2 points <= N): fallback from the median eigenvalue
    // of the valid cluster {0.01778, 0.05333}, clamped to >= 1 => identity.
    const Eigen::Matrix2d& cov_a = init.covs[0];
    EXPECT_TRUE(cov_a.isApprox(Eigen::Matrix2d::Identity(), 1e-15));

    // Valid cluster: exact sample covariance of group B.
    // diag(8/225), off-diag -4/225.
    const Eigen::Matrix2d& cov_b = init.covs[1];
    EXPECT_NEAR(cov_b(0, 0), 8.0 / 225.0, 1e-12);
    EXPECT_NEAR(cov_b(1, 1), 8.0 / 225.0, 1e-12);
    EXPECT_NEAR(cov_b(0, 1), -4.0 / 225.0, 1e-12);
    EXPECT_NEAR(cov_b(1, 0), -4.0 / 225.0, 1e-12);
}

TEST(Cov3Data_Initialization, CellCenters3DValidAndInvalidClusters) {
    // A: 6 axis points around origin, B: 6 axis points (radius 2) at x=20,
    // C: 2 points at x=40 (invalid for N=3 -> fallback).
    Eigen::MatrixXd pos(3, 14);
    int col = 0;
    const double rad_a = 1.0, rad_b = 2.0;
    for (double s : {1.0, -1.0}) {
        pos.col(col++) << s * rad_a, 0.0, 0.0;
        pos.col(col++) << 0.0, s * rad_a, 0.0;
        pos.col(col++) << 0.0, 0.0, s * rad_a;
    }
    for (double s : {1.0, -1.0}) {
        pos.col(col++) << 20.0 + s * rad_b, 0.0, 0.0;
        pos.col(col++) << 20.0, s * rad_b, 0.0;
        pos.col(col++) << 20.0, 0.0, s * rad_b;
    }
    pos.col(col++) << 40.0, 0.0, 0.0;
    pos.col(col++) << 40.5, 0.5, 0.5;

    auto init = baysor::cell_centers_uniformly<3>(pos, 3, /*confidences=*/nullptr,
                                                  /*scale=*/-1.0);

    ASSERT_EQ(init.centers.rows(), 3);
    EXPECT_EQ(init.covs.size(), 3u);
    ASSERT_EQ(init.assignment.size(), 14u);
    for (int i = 0; i < 6; ++i) EXPECT_EQ(init.assignment[i], 1);
    for (int i = 6; i < 12; ++i) EXPECT_EQ(init.assignment[i], 2);
    EXPECT_EQ(init.assignment[12], 3);
    EXPECT_EQ(init.assignment[13], 3);

    // Sample covariances of A and B are (2/6)I and (8/6)I.
    EXPECT_TRUE(init.covs[0].isApprox(Eigen::Matrix3d::Identity() * (2.0 / 6.0), 1e-12));
    EXPECT_TRUE(init.covs[1].isApprox(Eigen::Matrix3d::Identity() * (8.0 / 6.0), 1e-12));

    // Fallback for the 2-point cluster: median of all eigenvalues
    // ({1/3 x3, 4/3 x3} sorted, middle = 4/3), clamped to >= 1 => (4/3) I.
    EXPECT_TRUE(init.covs[2].isApprox(Eigen::Matrix3d::Identity() * (4.0 / 3.0), 1e-12));
}

TEST(Cov3Data_Initialization, CellCenters3DAllClustersTooSmallUseIdentityFallback) {
    // Two clusters of 3 points each: 3 <= N=3, so no cluster is valid and
    // the median-eigenvalue list stays empty => Ones fallback (identity).
    Eigen::MatrixXd pos(3, 6);
    pos.col(0) << 1.0, 0.0, 0.0;
    pos.col(1) << -1.0, 0.0, 0.0;
    pos.col(2) << 0.0, 1.0, 0.0;
    pos.col(3) << 20.0, 0.0, 0.0;
    pos.col(4) << 21.0, 0.0, 0.0;
    pos.col(5) << 20.0, 1.0, 0.0;

    auto init = baysor::cell_centers_uniformly<3>(pos, 2, nullptr, -1.0);

    ASSERT_EQ(init.covs.size(), 2u);
    EXPECT_TRUE(init.covs[0].isApprox(Eigen::Matrix3d::Identity(), 1e-15));
    EXPECT_TRUE(init.covs[1].isApprox(Eigen::Matrix3d::Identity(), 1e-15));
    ASSERT_EQ(init.assignment.size(), 6u);
    EXPECT_EQ(init.assignment[0], 1);
    EXPECT_EQ(init.assignment[3], 2);
}

TEST(Cov3Data_Initialization, CellCentersDegenerateWhenAllConfidencesAreLow) {
    // n_clusters <= 1 now raises Julia's "n must be > 1" error, so request 2
    // centers: no molecule clears the 0.25 confidence threshold, the center
    // selection degenerates to an empty set and every molecule is assigned
    // to cluster 1.
    Eigen::MatrixXd pos(2, 4);
    pos << 0.0, 1.0, 0.0, 1.0,
           0.0, 0.0, 1.0, 1.0;
    std::vector<double> low_conf(4, 0.1);

    auto init = baysor::cell_centers_uniformly<2>(pos, /*n_clusters=*/2, &low_conf, 1.0);

    EXPECT_EQ(init.centers.rows(), 0);
    EXPECT_EQ(init.centers.cols(), 2);
    EXPECT_TRUE(init.covs.empty());
    ASSERT_EQ(init.assignment.size(), 4u);
    for (int a : init.assignment) EXPECT_EQ(a, 1);
}

TEST(Cov3Data_Initialization, InitializeBmmData2DWithoutPriorAndDefaultConfidence) {
    baysor::MoleculeData data;
    data.x = {0.0, 0.5, 1.0, 3.0, 3.5, 4.0};
    data.y = {0.0, 0.5, 0.0, 0.0, 0.5, 1.0};
    data.gene = {1, 2, 0, 1, 2, 1};     // gene 0 must map to -1
    data.gene_names = {"A", "B"};
    data.cluster = {1, 1, 1, 2, 2, 2};
    data.nuclei_probs = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    // confidence intentionally empty -> defaults to 0.95 everywhere

    auto adj = baysor::build_molecule_graph(data);
    auto bm = baysor::initialize_bmm_data<2>(
        data, adj, /*n_cells_init=*/10, /*scale=*/0.25, /*scale_std=*/"25%",
        /*prior_seg_confidence=*/0.7, /*min_molecules_per_cell=*/5, /*verbose=*/true);

    EXPECT_EQ(bm.n_molecules(), 6);
    EXPECT_EQ(bm.n_components(), 6);           // n_cells_init clamped to n_mols
    EXPECT_EQ(bm.max_component_guid, 6);
    ASSERT_EQ(bm.position_data.rows(), 2);
    ASSERT_EQ(bm.position_data.cols(), 6);
    EXPECT_DOUBLE_EQ(bm.position_data(0, 5), 4.0);
    EXPECT_EQ(bm.composition_data, (std::vector<int>{0, 1, -1, 0, 1, 0}));
    EXPECT_EQ(bm.confidence, std::vector<double>(6, 0.95));
    EXPECT_EQ(bm.cluster_per_molecule, (std::vector<int>{1, 1, 1, 2, 2, 2}));
    EXPECT_EQ(bm.nuclei_prob_per_molecule,
              (std::vector<double>{0.1, 0.2, 0.3, 0.4, 0.5, 0.6}));
    EXPECT_TRUE(bm.segment_per_molecule.empty());
    EXPECT_TRUE(bm.n_molecules_per_segment.empty());

    ASSERT_EQ(bm.assignment.size(), 6u);
    std::vector<int> sorted_assign = bm.assignment;
    std::sort(sorted_assign.begin(), sorted_assign.end());
    EXPECT_EQ(sorted_assign, (std::vector<int>{1, 2, 3, 4, 5, 6}));

    EXPECT_EQ(bm.adj_list.n_molecules(), adj.n_molecules());
    EXPECT_GT(bm.adj_list.nnz(), 0);

    // noise_position_density = exp(-0.5 * N * 9) / (2*pi*scale^2)^(N/2)
    const double expected_npd =
        std::exp(-0.5 * 2 * 9.0) / std::pow(2.0 * baysor::kPi * 0.25 * 0.25, 2 / 2.0);
    EXPECT_NEAR(bm.noise_position_density, expected_npd, 1e-15);
    EXPECT_DOUBLE_EQ(bm.noise_density, 0.0);
    EXPECT_DOUBLE_EQ(bm.prior_seg_confidence, 0.7);
    EXPECT_DOUBLE_EQ(bm.mrf_strength, 0.1);
    EXPECT_TRUE(bm.use_gene_smoothing);
    EXPECT_DOUBLE_EQ(bm.min_nuclei_frac, 0.1);
    EXPECT_DOUBLE_EQ(bm.cluster_penalty_mult, 0.25);
    EXPECT_DOUBLE_EQ(bm.real_edge_weight, 1.0);

    // One molecule per initial cell, uniform gene prior, isotropic
    // covariance, and a shape prior built from scale / scale_std.
    ASSERT_EQ(bm.components.size(), 6u);
    for (int ci = 0; ci < 6; ++ci) {
        const auto& comp = bm.components[ci];
        EXPECT_EQ(comp.guid, ci + 1);
        EXPECT_EQ(comp.n_samples, 1);
        EXPECT_EQ(comp.composition_params.size(), 2);
        const auto counts = comp.composition_params.dense_counts();
        ASSERT_EQ(counts.size(), 2u);
        for (float c : counts) EXPECT_FLOAT_EQ(c, 1.0f);
        ASSERT_TRUE(comp.shape_prior.has_value());
        EXPECT_DOUBLE_EQ(comp.shape_prior->std_values(0), 0.25);
        EXPECT_DOUBLE_EQ(comp.shape_prior->std_value_stds(0), 0.0625);  // 25% of 0.25
        EXPECT_EQ(comp.shape_prior->n_samples, 5);
        const Eigen::Matrix2d expected_sigma = Eigen::Matrix2d::Identity() * 0.0625;
        EXPECT_TRUE(comp.position_params.sigma.isApprox(expected_sigma, 1e-15));
    }
}

TEST(Cov3Data_Initialization, InitializeBmmData3DWithPriorSegmentationAndConfidence) {
    baysor::MoleculeData data;
    data.x = {0.0, 1.0, 2.0, 10.0, 11.0, 12.0, 20.0, 21.0, 22.0};
    data.y = {0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0};
    data.z = {0.0, 0.5, 1.0, 2.0, 2.5, 3.0, 4.0, 4.5, 5.0};
    data.gene = {1, 2, 3, 1, 2, 3, 1, 2, 3};
    data.gene_names = {"A", "B", "C"};
    data.confidence = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};
    data.prior_segmentation = {1, 1, 1, 2, 2, 2, 3, 3, 0};

    auto adj = baysor::build_molecule_graph(data);
    auto bm = baysor::initialize_bmm_data<3>(
        data, adj, /*n_cells_init=*/3, /*scale=*/0.5, /*scale_std=*/"0.3",
        /*prior_seg_confidence=*/0.4, /*min_molecules_per_cell=*/3, /*verbose=*/false);

    EXPECT_EQ(bm.n_molecules(), 9);
    EXPECT_EQ(bm.n_components(), 3);
    EXPECT_EQ(bm.max_component_guid, 3);
    ASSERT_EQ(bm.position_data.rows(), 3);

    // Provided confidence is copied verbatim.
    EXPECT_EQ(bm.confidence, data.confidence);
    EXPECT_EQ(bm.composition_data, (std::vector<int>{0, 1, 2, 0, 1, 2, 0, 1, 2}));

    // Prior segmentation bookkeeping.
    EXPECT_EQ(bm.segment_per_molecule, data.prior_segmentation);
    ASSERT_EQ(bm.n_molecules_per_segment.size(), 3u);
    EXPECT_EQ(bm.n_molecules_per_segment[0], 3);
    EXPECT_EQ(bm.n_molecules_per_segment[1], 3);
    EXPECT_EQ(bm.n_molecules_per_segment[2], 2);
    ASSERT_EQ(bm.main_segment_per_cell.size(), 3u);
    for (int ms : bm.main_segment_per_cell) {
        EXPECT_GE(ms, 0);
        EXPECT_LE(ms, 3);
    }
    // update_n_mols_per_segment(): per-component segment counts sum up to the
    // number of that component's molecules that carry a prior segment.
    ASSERT_EQ(bm.assignment.size(), 9u);
    std::vector<int> expected_per_component(3, 0);
    for (int i = 0; i < 9; ++i) {
        if (bm.assignment[i] > 0 && data.prior_segmentation[i] > 0) {
            expected_per_component[bm.assignment[i] - 1]++;
        }
    }
    for (int ci = 0; ci < 3; ++ci) {
        int total = 0;
        for (const auto& [seg, n] : bm.components[ci].n_molecules_per_segment) total += n;
        EXPECT_EQ(total, expected_per_component[ci]) << "component " << ci;
    }

    // scale_std = "0.3" (absolute), scale = 0.5.
    for (const auto& comp : bm.components) {
        ASSERT_TRUE(comp.shape_prior.has_value());
        EXPECT_DOUBLE_EQ(comp.shape_prior->std_values(0), 0.5);
        EXPECT_DOUBLE_EQ(comp.shape_prior->std_value_stds(0), 0.3);
        EXPECT_EQ(comp.shape_prior->n_samples, 3);
        const Eigen::Matrix3d expected_cov =
            Eigen::Matrix3d::Identity() * 0.5 * 0.5;
        EXPECT_TRUE(comp.position_params.sigma.isApprox(expected_cov, 1e-15));
        EXPECT_EQ(comp.composition_params.size(), 3);
        EXPECT_EQ(comp.guid > 0, true);
        const auto counts = comp.composition_params.dense_counts();
        ASSERT_EQ(counts.size(), 3u);
        for (float c : counts) EXPECT_FLOAT_EQ(c, 1.0f);
    }
}

// ============================================================================
// boundary_estimation.cpp
// ============================================================================

TEST(Cov3Data_Boundary, DirectCall2DUsesDefaultCellNames) {
    // Two well-separated 2x2 squares, no cell_names pointer.
    Eigen::MatrixXd pos(2, 8);
    pos << 0.0, 2.0, 2.0, 0.0, 10.0, 12.0, 12.0, 10.0,
           0.0, 0.0, 2.0, 2.0,  0.0,  0.0,  2.0,  2.0;
    std::vector<int> labels = {1, 1, 1, 1, 2, 2, 2, 2};

    auto polys = baysor::boundary_polygons(pos, labels);

    ASSERT_EQ(polys.size(), 2u);
    ASSERT_EQ(polys.count("1"), 1u);
    ASSERT_EQ(polys.count("2"), 1u);
    ASSERT_EQ(polys.at("1").cols(), 4);
    ASSERT_EQ(polys.at("2").cols(), 4);
    EXPECT_NEAR(baysor::polygon_area(polys.at("1")), 4.0, 1e-9);
    EXPECT_NEAR(baysor::polygon_area(polys.at("2")), 4.0, 1e-9);
}

TEST(Cov3Data_Boundary, UnsupportedDimensionReturnsEmpty) {
    Eigen::MatrixXd pos(1, 5);
    pos.setZero();
    std::vector<int> labels(5, 1);
    auto polys = baysor::boundary_polygons(pos, labels);
    EXPECT_TRUE(polys.empty());
}

TEST(Cov3Data_Boundary, AutoWithoutPerZReturnsSingleLayer) {
    Eigen::MatrixXd pos(2, 4);
    pos << 0.0, 2.0, 2.0, 0.0,
           0.0, 0.0, 2.0, 2.0;
    std::vector<int> labels = {1, 1, 1, 1};

    auto [joined, stack] = baysor::boundary_polygons_auto(
        pos, labels, /*estimate_per_z=*/true, /*cell_names=*/nullptr, /*verbose=*/true);

    ASSERT_EQ(stack.size(), 1u);
    EXPECT_EQ(stack[0].first, "2d");
    ASSERT_EQ(joined.size(), 1u);
    EXPECT_EQ(joined.count("1"), 1u);
}

TEST(Cov3Data_Boundary, AdmixtureInsideBoundaryTriangleCutsTheBorder) {
    // Cell 1: square with an interior point; cell 2: a single intruder that
    // sits inside the boundary triangle (c1, c2, center).
    Eigen::MatrixXd with_intruder(2, 6);
    with_intruder << 0.0, 4.0, 4.0, 0.0, 2.0, 2.0,
                     0.0, 0.0, 4.0, 4.0, 1.8, 0.9;
    std::vector<int> labels_with = {1, 1, 1, 1, 1, 2};
    auto polys_with = baysor::boundary_polygons(with_intruder, labels_with);

    // Without the intruder the border is the plain 4-gon of the square.
    Eigen::MatrixXd clean(2, 5);
    clean << 0.0, 4.0, 4.0, 0.0, 2.0,
             0.0, 0.0, 4.0, 4.0, 1.8;
    std::vector<int> labels_clean(5, 1);
    auto polys_clean = baysor::boundary_polygons(clean, labels_clean);

    ASSERT_EQ(polys_clean.count("1"), 1u);
    ASSERT_EQ(polys_clean.at("1").cols(), 4);
    EXPECT_NEAR(baysor::polygon_area(polys_clean.at("1")), 16.0, 1e-9);

    // With the admixture point the boundary triangle is excluded and the
    // border reroutes through the interior center: a 5-gon of area 12.4.
    ASSERT_EQ(polys_with.count("1"), 1u);
    ASSERT_EQ(polys_with.at("1").cols(), 5);
    EXPECT_NEAR(baysor::polygon_area(polys_with.at("1")), 12.4, 1e-9);
    // The intruder cell keeps its 4-point star polygon.
    ASSERT_EQ(polys_with.count("2"), 1u);
    EXPECT_EQ(polys_with.at("2").cols(), 4);
}

TEST(Cov3Data_Boundary, AdmixtureSkipGuardKeepsConvexPentagonBorder) {
    // Convex pentagon: its "middle" triangle has exactly one border edge and
    // both ends of every internal edge are hull vertices with two border
    // edges each -> the skip guard fires and the border is unchanged.
    Eigen::MatrixXd with_intruder(2, 6);
    with_intruder << 0.0, 4.0, 4.0, 2.0, 0.0, 2.0,
                     0.0, 0.0, 4.0, 6.0, 4.0, 2.8;
    std::vector<int> labels_with = {1, 1, 1, 1, 1, 2};

    Eigen::MatrixXd clean(2, 5);
    clean << 0.0, 4.0, 4.0, 2.0, 0.0,
             0.0, 0.0, 4.0, 6.0, 4.0;
    std::vector<int> labels_clean(5, 1);

    auto polys_with = baysor::boundary_polygons(with_intruder, labels_with);
    auto polys_clean = baysor::boundary_polygons(clean, labels_clean);

    ASSERT_EQ(polys_clean.count("1"), 1u);
    ASSERT_EQ(polys_with.count("1"), 1u);
    ASSERT_EQ(polys_clean.at("1").cols(), 5);
    ASSERT_EQ(polys_with.at("1").cols(), 5);
    // The skip guard keeps the full pentagon boundary: area 20 either way.
    EXPECT_NEAR(baysor::polygon_area(polys_clean.at("1")), 20.0, 1e-9);
    EXPECT_NEAR(baysor::polygon_area(polys_with.at("1")), 20.0, 1e-9);
}

TEST(Cov3Data_Boundary, GridWithAbsentLabelAndSinglePointLabels) {
    // Label 2 never appears in the grid (borders[1] stays empty -> the
    // <2-points branch), labels 1 and 3 each produce a single border point
    // -> also the <2-points branch.
    Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> grid(2, 2);
    grid << 1u, 3u,
            0u, 0u;

    auto polys = baysor::boundary_polygons_from_grid(grid);
    ASSERT_EQ(polys.size(), 3u);
    for (const auto& p : polys) {
        EXPECT_EQ(p.cols(), 0);
    }

    // Two adjacent pixels of one label: each contributes exactly one border
    // point, so there are 2 border points but only 2 distinct (duplicate)
    // locations -> no Delaunay faces -> empty polygon entry.
    Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic> row(2, 2);
    row << 1u, 1u,
           0u, 0u;
    auto polys_row = baysor::boundary_polygons_from_grid(row);
    ASSERT_EQ(polys_row.size(), 1u);
    EXPECT_EQ(polys_row[0].cols(), 0);
}

TEST(Cov3Data_Boundary, AutoBinnedZStackProducesTenLayers) {
    auto run = [](int n_mols) {
        Eigen::MatrixXd pos(3, n_mols);
        for (int i = 0; i < n_mols; ++i) {
            pos(0, i) = static_cast<double>(i % 7);
            pos(1, i) = static_cast<double>(i / 7);
            pos(2, i) = static_cast<double>(i % 15);   // 15 unique z > 10 slices
        }
        std::vector<int> labels(n_mols, 1);
        return baysor::boundary_polygons_auto(
            pos, labels, /*estimate_per_z=*/true, /*cell_names=*/nullptr, /*verbose=*/true);
    };

    // 41 molecules: the 2.5% quantile index (0.025 * 40 = 1.0) is integral.
    auto [joined41, stack41] = run(41);
    EXPECT_EQ(joined41.count("1"), 1u);
    ASSERT_EQ(stack41.size(), 11u);           // "2d" + 10 binned layers
    EXPECT_EQ(stack41[0].first, "2d");
    std::set<std::string> layer_names;
    for (size_t i = 1; i < stack41.size(); ++i) {
        layer_names.insert(stack41[i].first);
        EXPECT_EQ(stack41[i].first.front(), '[');
        EXPECT_NE(stack41[i].first.find(','), std::string::npos);
        EXPECT_EQ(stack41[i].first.back(), ']');
    }
    EXPECT_EQ(layer_names.size(), 10u);

    // 42 molecules: both quantile indices are fractional (interpolation path).
    auto [joined42, stack42] = run(42);
    EXPECT_EQ(joined42.count("1"), 1u);
    ASSERT_EQ(stack42.size(), 11u);
    EXPECT_EQ(stack42[0].first, "2d");
    EXPECT_EQ(stack42[1].first.front(), '[');
}

TEST(Cov3Data_Boundary, AutoBinnedZStackHonoursCustomSliceLimit) {
    // 20 distinct z values with a 4-point square per layer and a single cell,
    // so every quantile bin is non-empty. More z values than either limit.
    constexpr int n_z = 20;
    constexpr int points_per_z = 4;
    Eigen::MatrixXd pos(3, n_z * points_per_z);
    for (int zi = 0; zi < n_z; ++zi) {
        const double z = 0.5 * static_cast<double>(zi);
        const int o = zi * points_per_z;
        pos.col(o + 0) << 0.0, 0.0, z;
        pos.col(o + 1) << 2.0, 0.0, z;
        pos.col(o + 2) << 2.0, 2.0, z;
        pos.col(o + 3) << 0.0, 2.0, z;
    }
    std::vector<int> labels(n_z * points_per_z, 1);
    std::vector<std::string> cell_names = {"cell_1"};

    auto sink = std::make_shared<baysor_test::CapturingSink>();
    baysor_test::LoggerGuard logger(sink);

    auto count_layers = [&](int max_z_slices) {
        auto [joined, stack] = baysor::boundary_polygons_auto(
            pos, labels, /*estimate_per_z=*/true, &cell_names, /*verbose=*/true,
            max_z_slices);
        EXPECT_EQ(joined.count("cell_1"), 1u);
        EXPECT_EQ(stack[0].first, "2d");
        std::set<std::string> names;
        for (size_t i = 1; i < stack.size(); ++i) {
            names.insert(stack[i].first);
            EXPECT_EQ(stack[i].second.count("cell_1"), 1u);
        }
        EXPECT_EQ(names.size(), stack.size() - 1);  // distinct layer names
        return static_cast<int>(stack.size()) - 1;
    };

    // More z values (20) than either limit, and the two limits give different
    // layer counts.
    EXPECT_EQ(count_layers(5), 5);
    EXPECT_EQ(count_layers(10), 10);

    // The warning names the config option so users can find how to change it.
    const std::string logged = sink->data();
    EXPECT_NE(logged.find("Too many z values"), std::string::npos) << logged;
    EXPECT_NE(logged.find("max_z_slices"), std::string::npos) << logged;
}

TEST(Cov3Data_Boundary, InternalBorderFilterStopsAfterMaxIterations) {
    // The max_iters guard is only visible through its warning; capture it.
    auto sink = std::make_shared<baysor_test::CapturingSink>();
    baysor_test::LoggerGuard logger(sink);

    // Square + interior center. Triangles (vertex indices into pos):
    //   T0 = (0, 1, 4), T1 = (1, 2, 4), T2 = (2, 3, 4), T3 = (3, 0, 4)
    // with the intruder inside T0 only.
    Eigen::MatrixXd pos(2, 5);
    pos << 0.0, 4.0, 4.0, 0.0, 2.0,
           0.0, 0.0, 4.0, 4.0, 1.8;
    Eigen::MatrixXd non_cell(2, 1);
    non_cell << 2.0, 0.9;

    const std::vector<std::array<int, 3>> triangles = {
        {0, 1, 4}, {1, 2, 4}, {2, 3, 4}, {3, 0, 4}};

    // One iteration is enough to exclude exactly the intruder-carrying T0,
    // which is not convergence -> the max_iters warning fires, and the
    // remaining border is {T0's neighbours plus the untouched spokes}.
    auto edges_1 = baysor::internal::find_border_without_admixture(
        triangles, pos, non_cell, /*max_iters=*/1);
    const auto got_1 = cov3_sorted_edges(edges_1);
    const std::vector<std::pair<int, int>> expected_1 = {
        {0, 3}, {0, 4}, {1, 2}, {1, 4}, {2, 3}};
    EXPECT_EQ(got_1, expected_1);
    EXPECT_NE(sink->data().find(
                  "Polygon filtering did not converge within 1 iterations"),
              std::string::npos) << sink->data();
    sink->clear();

    // More iterations converge to the same border (the only exclusion is T0:
    // the remaining triangles either have >1 border edge or no admixture),
    // so the warning must not fire.
    auto edges_full = baysor::internal::find_border_without_admixture(
        triangles, pos, non_cell, /*max_iters=*/100);
    EXPECT_EQ(cov3_sorted_edges(edges_full), expected_1);
    EXPECT_EQ(sink->data().find("did not converge"), std::string::npos)
        << sink->data();

    // Without an admixture point nothing is excluded: the plain hull.
    Eigen::MatrixXd empty_non_cell(2, 0);
    auto edges_clean = baysor::internal::find_border_without_admixture(
        triangles, pos, empty_non_cell, /*max_iters=*/1);
    const std::vector<std::pair<int, int>> expected_hull = {
        {0, 1}, {0, 3}, {1, 2}, {2, 3}};
    EXPECT_EQ(cov3_sorted_edges(edges_clean), expected_hull);
}

TEST(Cov3Data_Boundary, InternalBorderToPolyHonoursLengthCap) {
    // The cap guard is only visible through its warning; capture it.
    auto sink = std::make_shared<baysor_test::CapturingSink>();
    baysor_test::LoggerGuard logger(sink);

    // A simple 4-cycle closes under the default cap ...
    const std::vector<std::pair<int, int>> cycle = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    auto poly = baysor::internal::border_edges_to_poly(cycle);
    EXPECT_EQ(poly.size(), 4u);
    EXPECT_EQ(sink->data().find("Could not build a polygon border"), std::string::npos)
        << sink->data();

    // ... but with a cap of two steps the walk cannot return to the start
    // and the function warns and gives up.
    auto capped = baysor::internal::border_edges_to_poly(cycle, /*max_border_len=*/2);
    EXPECT_TRUE(capped.empty());
    EXPECT_NE(sink->data().find("Could not build a polygon border of size 4"),
              std::string::npos) << sink->data();
    sink->clear();

    // Degenerate input: at most two border edges can never form a polygon.
    EXPECT_TRUE(baysor::internal::border_edges_to_poly({{0, 1}}).empty());
    EXPECT_TRUE(baysor::internal::border_edges_to_poly({}).empty());
}

// ============================================================================
// neighborhood_composition.cpp
// ============================================================================

TEST(Cov3Data_Neighborhood, DistanceFloorRetriesWhenNearestNeighborsCoincide) {
    // Three coincident points force the exact-closest search to widen its
    // k (2 -> 4) before any non-zero distance shows up.
    Eigen::MatrixXd pos(2, 5);
    pos << 0.0, 0.0, 0.0, 3.0, 10.0,
           0.0, 0.0, 0.0, 4.0,  0.0;

    const double floor = baysor::neighborhood_distance_floor(pos);

    // Closest non-zero distances: {5, 5, 5, 5, 10} -> median 5.
    EXPECT_DOUBLE_EQ(floor, 5.0);
}

TEST(Cov3Data_Neighborhood, AutoDistanceFloorMatchesExplicitFloor) {
    Eigen::MatrixXd pos(2, 12);
    for (int i = 0; i < 12; ++i) {
        pos(0, i) = static_cast<double>(i % 4) * 1.5;
        pos(1, i) = static_cast<double>(i / 4) * 2.0 + 0.3 * (i % 3);
    }
    std::vector<int> genes = {1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3};
    const int n_genes = 3;
    const int k = 4;

    Eigen::MatrixXf emb(2, n_genes);
    emb << 1.0f, 0.0f, 0.0f,
           0.0f, 1.0f, 0.0f;

    // distance_floor <= 0 makes stream_projected_neighborhood_vectors compute
    // the median closest distance itself; that must equal the explicit value.
    const double explicit_floor = baysor::neighborhood_distance_floor(pos);
    auto auto_floor = baysor::project_neighborhood_vectors(
        pos, genes, k, emb, n_genes, /*query_ids=*/nullptr, /*confidences=*/nullptr,
        /*normalize_by_dist=*/true, /*normalize=*/true, /*distance_floor=*/-1.0,
        /*log_transform=*/false);
    auto explicit_mat = baysor::project_neighborhood_vectors(
        pos, genes, k, emb, n_genes, /*query_ids=*/nullptr, /*confidences=*/nullptr,
        /*normalize_by_dist=*/true, /*normalize=*/true, explicit_floor,
        /*log_transform=*/false);

    ASSERT_EQ(auto_floor.rows(), explicit_mat.rows());
    ASSERT_EQ(auto_floor.cols(), explicit_mat.cols());
    EXPECT_EQ((auto_floor - explicit_mat).norm(), 0.0f);

    // The floor actually matters: a huge floor changes the weighting.
    auto huge_floor = baysor::project_neighborhood_vectors(
        pos, genes, k, emb, n_genes, /*query_ids=*/nullptr, /*confidences=*/nullptr,
        /*normalize_by_dist=*/true, /*normalize=*/true, 1000.0,
        /*log_transform=*/false);
    EXPECT_GT((auto_floor - huge_floor).norm(), 0.0f);
}

// ============================================================================
// Molecule-graph reuse: the triangulation of the confidence step is reused for
// the segmentation graph.
// ============================================================================

namespace {

// Regular grid, no duplicate coordinates: direct and reused builds must agree
// bit-for-bit (with duplicates the two normalized point sets differ by
// design — see the RNG-parity test below).
baysor::MoleculeData make_reuse_grid() {
    baysor::MoleculeData data;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            data.x.push_back(static_cast<double>(i));
            data.y.push_back(static_cast<double>(j));
            data.gene.push_back(1 + ((i + j) % 2));
        }
    }
    data.gene_names = {"A", "B"};
    return data;
}

// Same grid plus exact coordinate duplicates, so normalize_points consumes a
// non-empty jitter batch of the global RNG.
baysor::MoleculeData make_reuse_duplicates() {
    baysor::MoleculeData data = make_reuse_grid();
    for (int d = 0; d < 4; ++d) {
        data.x.push_back(data.x[d]);
        data.y.push_back(data.y[d]);
        data.gene.push_back(1);
    }
    return data;
}

void expect_same_csr(const baysor::AdjList& a, const baysor::AdjList& b) {
    EXPECT_EQ(a.indptr, b.indptr);
    EXPECT_EQ(a.indices, b.indices);
    ASSERT_EQ(a.weights.size(), b.weights.size());
    for (size_t i = 0; i < a.weights.size(); ++i) {
        ASSERT_EQ(a.weights[i], b.weights[i]) << "weight " << i;
    }
}

} // namespace

TEST(Cov3Data_GraphReuse, PrecomputedEdgesMatchDirectBuild) {
    const auto data = make_reuse_grid();

    const auto direct = baysor::build_molecule_graph(data, /*filter=*/true);
    auto edges = baysor::compute_molecule_adjacency(data);
    // No duplicate coordinates: the reuse trigger holds.
    ASSERT_EQ(edges.normalize_rng_draws, 0);
    const auto reused = baysor::build_molecule_graph(
        data, /*filter=*/true, /*use_local_gene_similarities=*/false,
        baysor::AdjacencyType::Auto, /*composition_neighborhood=*/0, /*n_gene_pcs=*/0,
        std::move(edges));

    ASSERT_GT(direct.nnz(), 0);
    expect_same_csr(reused, direct);

    // The unfiltered variant (the confidence step's MRF view) matches a
    // direct unfiltered build as well.
    const auto direct_unfiltered = baysor::build_molecule_graph(data, /*filter=*/false);
    auto edges2 = baysor::compute_molecule_adjacency(data);
    const auto reused_unfiltered = baysor::build_molecule_graph(
        data, /*filter=*/false, /*use_local_gene_similarities=*/false,
        baysor::AdjacencyType::Auto, /*composition_neighborhood=*/0, /*n_gene_pcs=*/0,
        std::move(edges2));
    ASSERT_GT(direct_unfiltered.nnz(), 0);
    expect_same_csr(reused_unfiltered, direct_unfiltered);
}

TEST(Cov3Data_GraphReuse, PrecomputedEdgesPreserveGlobalRngStream) {
    // With duplicate coordinates the two normalize_points() runs consume
    // different jitter batches, so build_molecule_graph falls back to
    // recomputing: the segmentation graph and the global RNG stream must be
    // exactly what the historical double build produced (confidence build,
    // then segmentation build).
    baysor_test::GlobalRngGuard rng_guard;
    const auto data = make_reuse_duplicates();
    const int n = data.n_molecules();

    // Direct flow: two recomputes (confidence build, then segmentation build).
    baysor::reset_global_xoshiro_rng();
    (void)baysor::build_molecule_graph(data, /*filter=*/false);
    const auto direct_segmentation = baysor::build_molecule_graph(data, /*filter=*/true);
    const double after_direct = baysor::global_xoshiro_rng().rand_float64();

    // Reuse flow: one compute, then the segmentation build from its edges —
    // which detects the duplicates and recomputes instead of reusing.
    baysor::reset_global_xoshiro_rng();
    auto edges = baysor::compute_molecule_adjacency(data);
    ASSERT_GT(edges.normalize_rng_draws, 0);  // the duplicates really jittered
    (void)baysor::build_molecule_graph_from_edges(edges, n);
    const auto reused_segmentation = baysor::build_molecule_graph(
        data, /*filter=*/true, /*use_local_gene_similarities=*/false,
        baysor::AdjacencyType::Auto, /*composition_neighborhood=*/0, /*n_gene_pcs=*/0,
        std::move(edges));
    const double after_reuse = baysor::global_xoshiro_rng().rand_float64();

    ASSERT_EQ(after_direct, after_reuse);
    // The fallback graph is the historical one, bit for bit.
    ASSERT_GT(direct_segmentation.nnz(), 0);
    expect_same_csr(reused_segmentation, direct_segmentation);
}

TEST(Cov3Data_GraphReuse, ConfidenceDetailsAdjacencyBuildsSameGraph) {
    // The confidence step returns the edges it computed; the segmentation
    // graph built from them equals the graph built directly.
    const auto data = make_reuse_grid();

    auto details = baysor::estimate_confidence_details(data, /*nn_id=*/3);
    const auto fresh_edges = baysor::compute_molecule_adjacency(data);
    ASSERT_EQ(details.adjacency.edge_src, fresh_edges.edge_src);
    ASSERT_EQ(details.adjacency.edge_dst, fresh_edges.edge_dst);
    ASSERT_EQ(details.adjacency.edge_dists, fresh_edges.edge_dists);
    EXPECT_EQ(details.adjacency.normalize_rng_draws, fresh_edges.normalize_rng_draws);

    const auto direct = baysor::build_molecule_graph(data, /*filter=*/true);
    const auto reused = baysor::build_molecule_graph(
        data, /*filter=*/true, /*use_local_gene_similarities=*/false,
        baysor::AdjacencyType::Auto, /*composition_neighborhood=*/0, /*n_gene_pcs=*/0,
        std::move(details.adjacency));
    ASSERT_GT(direct.nnz(), 0);
    expect_same_csr(reused, direct);
}
