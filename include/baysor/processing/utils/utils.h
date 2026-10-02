#pragma once

#include <Eigen/Sparse>
#include <vector>

namespace baysor {

/// Build a sparse count vector from gene IDs (like Julia's count_array_sparse)
Eigen::SparseVector<float> count_array_sparse(
    const int* values, int n, int total,
    const double* weights = nullptr,
    bool normalize = false
);

/// Parallel KNN query: for each column of query_points, find k nearest neighbors in tree_points.
/// Row i of the result lives at [i*k, (i+1)*k) of `indices` / `distances`.
struct KnnResult {
    std::vector<int> indices;       // n * k
    std::vector<double> distances;  // n * k, Euclidean distances
    int n = 0;                      // query rows (0 for an empty tree or query, or k <= 0)
    int k = 0;                      // neighbors per row, clamped to the number of tree points

    const int* idx_row(int i) const { return indices.data() + static_cast<std::size_t>(i) * k; }
    const double* dist_row(int i) const { return distances.data() + static_cast<std::size_t>(i) * k; }
};

KnnResult knn_parallel(
    const Eigen::MatrixXd& tree_points,
    const Eigen::MatrixXd& query_points,
    int k,
    bool sorted = false
);

/// Distance from each point to its kth (0-based) neighbor among `points`
/// itself (the point is neighbor 0), as in knn_parallel(points, points, k,
/// true), but computed in blocks without the n x k result. k is clamped to
/// the number of points, kth to k - 1; returns {} for no points or k <= 0.
std::vector<double> knn_kth_distances(
    const Eigen::MatrixXd& points,
    int k,
    int kth
);

/// Queries per k-NN block, bounding the block's neighbor indices and
/// distances to about 32 MiB so that memory does not grow with k.
int knn_block_size(int k);

/// Orders a k-NN row that is sorted by distance (as nanoflann returns it) by
/// (distance, index): only the indices inside runs of equal distances move.
void sort_tied_runs_by_index(int* indices, const double* distances, int k);

} // namespace baysor
