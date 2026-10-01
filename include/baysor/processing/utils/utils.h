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
/// Returns (indices, distances), each n_points x k.
///
/// Flat row-major storage (review-bmm.md C7-1): row i of the result lives at
/// [i*k, (i+1)*k) of `indices` / `distances`, so the whole result is two
/// allocations instead of two heap vectors per query. k is clamped to the
/// number of tree points; `n == 0` marks "no results" (empty tree, empty
/// query, or k <= 0), matching the former empty vector<vector> return.
struct KnnResult {
    std::vector<int> indices;       // n * k
    std::vector<double> distances;  // n * k, real (sqrt'ed) distances
    int n = 0;                      // number of query rows (0 = no results)
    int k = 0;                      // neighbors per row (clamped to n_tree)

    const int* idx_row(int i) const {
        return indices.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(k);
    }
    int* idx_row(int i) {
        return indices.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(k);
    }
    const double* dist_row(int i) const {
        return distances.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(k);
    }
    double* dist_row(int i) {
        return distances.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(k);
    }
};

KnnResult knn_parallel(
    const Eigen::MatrixXd& tree_points,
    const Eigen::MatrixXd& query_points,
    int k,
    bool sorted = false
);

/// Distance from each point to its kth (0-based) neighbor among `points`
/// itself (the point counts as neighbor 0), i.e. the row
/// `knn_parallel(points, points, k, true).distances[i][kth]` without keeping
/// the n x k result: the kd-tree is built once and queries run in blocks, so
/// peak memory is O(block * k) instead of O(n * k) (REPORT.md 6.4, the 3 GiB
/// confidence kNN result at 10.6M molecules). k is clamped to the number of
/// points, kth to k - 1; returns {} for no points or k <= 0.
std::vector<double> knn_kth_distances(
    const Eigen::MatrixXd& points,
    int k,
    int kth
);

} // namespace baysor
