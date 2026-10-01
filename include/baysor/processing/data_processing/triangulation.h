#pragma once

#include <vector>
#include <Eigen/Dense>

namespace baysor {

/// Result of adjacency computation
struct AdjacencyResult {
    std::vector<int> edge_src;
    std::vector<int> edge_dst;
    std::vector<double> edge_dists;

    /// Number of global-RNG draws normalize_points() spent jittering duplicate
    /// points while producing this result (0 when there are no duplicates).
    /// Two normalize_points() runs over the same points consume different
    /// jitter batches whenever this is > 0, so their triangulations may differ;
    /// a reuser must recompute in that case and only reuse when this is 0
    /// (then the normalization is deterministic and a recomputation would
    /// produce these exact edges).
    int normalize_rng_draws = 0;
};

enum class AdjacencyType { Auto, Triangulation, Knn, Both };

/// Compute adjacency list from point positions using Delaunay triangulation and/or KNN.
/// For 3D data, only KNN is supported.
AdjacencyResult adjacency_list(
    const Eigen::MatrixXd& points,    // dims x n_points
    bool filter = true,
    double n_mads = 2.0,
    int k_adj = 5,
    AdjacencyType type = AdjacencyType::Auto
);

/// Normalize points to [1+eps, 2-eps] range for Delaunay tessellation.
/// If rng_draws is not null, it receives the number of global-RNG draws the
/// duplicate-point jitter consumed (see AdjacencyResult::normalize_rng_draws).
Eigen::MatrixXd normalize_points(const Eigen::MatrixXd& points, int* rng_draws = nullptr);

/// Filter edges longer than median + n_mads * MAD in log-distance
void filter_long_edges(AdjacencyResult& result, double n_mads = 2.0);

} // namespace baysor
