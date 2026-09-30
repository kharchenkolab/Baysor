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
    result.indices.resize(n_query);
    result.distances.resize(n_query);

    // Build KD-tree
    EigenColMajorAdaptor adaptor(tree_points);
    KDTree tree(n_dims, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(/* max_leaf = */ 10));

    parallel_for(0, n_query, 256, [&](int i) {
        result.indices[i].resize(k);
        result.distances[i].resize(k);

        nanoflann::KNNResultSet<double, int> resultSet(k);
        resultSet.init(result.indices[i].data(), result.distances[i].data());
        tree.findNeighbors(
            resultSet,
            query_points.col(i).data(),
            nanoflann::SearchParameters(/*eps=*/0.0f, /*sorted=*/sorted)
        );

        // Keep sorted=true deterministic even when the backend does not define
        // a stable tie order for equal-distance neighbors.
        if (sorted) {
            std::vector<int> order(k);
            std::iota(order.begin(), order.end(), 0);
            std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
                if (result.distances[i][a] != result.distances[i][b]) {
                    return result.distances[i][a] < result.distances[i][b];
                }
                return result.indices[i][a] < result.indices[i][b];
            });

            std::vector<int> sorted_indices(k);
            std::vector<double> sorted_distances(k);
            for (int j = 0; j < k; ++j) {
                sorted_indices[j] = result.indices[i][order[j]];
                sorted_distances[j] = result.distances[i][order[j]];
            }
            result.indices[i].swap(sorted_indices);
            result.distances[i].swap(sorted_distances);
        }

        // nanoflann returns squared distances; convert to actual distances.
        for (double& d : result.distances[i]) {
            d = std::sqrt(d);
        }
    });

    return result;
}

} // namespace baysor
