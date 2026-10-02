#include "baysor/processing/utils/utils.h"
#include "baysor/utils/thread_pool.h"
#include <third_party/nanoflann.hpp>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace baysor {

// ============================================================================
// count_array_sparse
// ============================================================================

Eigen::SparseVector<float> count_array_sparse(
    const int* values, int n, int total,
    const double* weights,
    bool normalize
) {
    if (n == 0) return Eigen::SparseVector<float>(total);

    // Sort indices by gene ID
    std::vector<int> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(), [&](int a, int b) {
        return values[a] < values[b];
    });

    // Accumulate counts grouped by gene ID
    std::vector<int> indices;
    std::vector<float> counts;
    constexpr double min_val = 1e-5;

    int last_id = values[perm[0]];
    double cnt = 0.0;
    int id = 0;

    for (int pi = 0; pi < n; ++pi) {
        int i = perm[pi];
        id = values[i];
        if (id != last_id) {
            if (cnt > min_val) {
                counts.push_back(static_cast<float>(cnt));
                indices.push_back(last_id - 1); // 1-based gene IDs -> 0-based sparse index
            }
            cnt = 0.0;
            last_id = id;
        }
        cnt += weights ? weights[i] : 1.0;
    }
    if (cnt > min_val) {
        counts.push_back(static_cast<float>(cnt));
        indices.push_back(id - 1); // 1-based -> 0-based
    }

    if (normalize && !counts.empty()) {
        float s = 0.0f;
        for (float c : counts) s += c;
        if (s > 0.0f) {
            for (float& c : counts) c /= s;
        }
    }

    // Build SparseVector
    Eigen::SparseVector<float> sv(total);
    sv.reserve(static_cast<int>(indices.size()));
    for (size_t i = 0; i < indices.size(); ++i) {
        if (indices[i] >= 0 && indices[i] < total) {
            sv.insert(indices[i]) = counts[i];
        }
    }
    return sv;
}

// ============================================================================
// knn_parallel (nanoflann-based)
// ============================================================================

// Adaptor for dims x N Eigen matrix (columns = points)
struct EigenColMajorAdaptor {
    const Eigen::MatrixXd& mat;
    EigenColMajorAdaptor(const Eigen::MatrixXd& m) : mat(m) {}

    inline size_t kdtree_get_point_count() const { return static_cast<size_t>(mat.cols()); }

    inline double kdtree_get_pt(const size_t idx, const size_t dim) const {
        return mat(static_cast<Eigen::Index>(dim), static_cast<Eigen::Index>(idx));
    }

    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, EigenColMajorAdaptor>,
    EigenColMajorAdaptor,
    -1,  // dynamic dimensionality
    int  // index type
>;

int knn_block_size(int k) {
    constexpr std::size_t budget_bytes = std::size_t(32) << 20;
    constexpr std::size_t min_block = 2048;  // keeps the 256-query parallel chunks fed
    constexpr std::size_t max_block = 32768;
    const std::size_t per_query = static_cast<std::size_t>(std::max(k, 1)) * (sizeof(int) + sizeof(double));
    return static_cast<int>(std::clamp(budget_bytes / per_query, min_block, max_block));
}

void sort_tied_runs_by_index(int* indices, const double* distances, int k) {
    for (int begin = 0, end; begin < k; begin = end) {
        for (end = begin + 1; end < k && distances[end] == distances[begin]; ++end) {}
        if (end - begin > 1) std::sort(indices + begin, indices + end);
    }
}

KnnResult knn_parallel(
    const Eigen::MatrixXd& tree_points,
    const Eigen::MatrixXd& query_points,
    int k,
    bool sorted
) {
    const int n_dims = static_cast<int>(tree_points.rows());
    const int n_tree = static_cast<int>(tree_points.cols());
    const int n_query = static_cast<int>(query_points.cols());

    if (n_tree == 0 || n_query == 0 || k <= 0) {
        return {};
    }

    // Clamp k to available points
    k = std::min(k, n_tree);

    KnnResult result;
    result.n = n_query;
    result.k = k;
    result.indices.resize(static_cast<std::size_t>(n_query) * k);
    result.distances.resize(static_cast<std::size_t>(n_query) * k);

    // Build KD-tree
    EigenColMajorAdaptor adaptor(tree_points);
    KDTree tree(n_dims, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(/* max_leaf = */ 10));

    parallel_for(0, n_query, 256, [&](int i) {
        int* row_indices = result.indices.data() + static_cast<std::size_t>(i) * k;
        double* row_distances = result.distances.data() + static_cast<std::size_t>(i) * k;

        nanoflann::KNNResultSet<double, int> resultSet(k);
        resultSet.init(row_indices, row_distances);
        tree.findNeighbors(
            resultSet,
            query_points.col(i).data(),
            nanoflann::SearchParameters(/*eps=*/0.0f, /*sorted=*/sorted)
        );

        // Keep sorted=true deterministic even when the backend does not define
        // a stable tie order for equal-distance neighbors.
        if (sorted) {
            sort_tied_runs_by_index(row_indices, row_distances, k);
        }

        // nanoflann returns squared distances; convert to actual distances.
        for (int j = 0; j < k; ++j) {
            row_distances[j] = std::sqrt(row_distances[j]);
        }
    });

    return result;
}

std::vector<double> knn_kth_distances(
    const Eigen::MatrixXd& points,
    int k,
    int kth
) {
    const int n = static_cast<int>(points.cols());
    if (n == 0 || k <= 0) {
        return {};
    }
    k = std::min(k, n);
    kth = std::min(kth, k - 1);

    EigenColMajorAdaptor adaptor(points);
    KDTree tree(static_cast<int>(points.rows()), adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(/* max_leaf = */ 10));

    const int block = std::min(knn_block_size(k), n);
    std::vector<int> scratch_indices(static_cast<std::size_t>(block) * k);
    std::vector<double> scratch_distances(static_cast<std::size_t>(block) * k);
    std::vector<double> out(n);

    for (int block_begin = 0; block_begin < n; block_begin += block) {
        parallel_for(block_begin, std::min(block_begin + block, n), 256, [&](int i) {
            const std::size_t offset = static_cast<std::size_t>(i - block_begin) * k;
            nanoflann::KNNResultSet<double, int> resultSet(k);
            resultSet.init(scratch_indices.data() + offset, scratch_distances.data() + offset);
            tree.findNeighbors(
                resultSet,
                points.col(i).data(),
                nanoflann::SearchParameters(/*eps=*/0.0f, /*sorted=*/true)
            );
            // The kth distance does not depend on the order of tied neighbors.
            out[i] = std::sqrt(scratch_distances[offset + kth]);
        });
    }

    return out;
}

} // namespace baysor
