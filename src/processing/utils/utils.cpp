#include "baysor/processing/utils/utils.h"
#include "baysor/utils/thread_pool.h"
#include <third_party/nanoflann.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
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

KDTree make_kdtree(const Eigen::MatrixXd& tree_points, const EigenColMajorAdaptor& adaptor) {
    return KDTree(
        static_cast<int>(tree_points.rows()), adaptor,
        nanoflann::KDTreeSingleIndexAdaptorParams(/* max_leaf = */ 10)
    );
}

// nanoflann returns neighbours in nondecreasing distance order with
// equal-distance neighbours in traversal order (KNNResultSet::addPoint keeps
// the row sorted and inserts ties after existing ones). For sorted results we
// reorder every tied run by index, which is exactly the stable sort by
// (distance, index) done previously over the whole row: distances outside tied
// runs keep their order, and distances within a run are equal, so only the
// indices move. Tied runs are scanned (C7-2) instead of sorting every row, so
// the common no-tie case allocates nothing and touches the row once.
void sort_tied_runs_by_index(int* indices, const double* squared_distances, int k) {
    int run_begin = 0;
    while (run_begin < k) {
        int run_end = run_begin + 1;
        while (run_end < k && squared_distances[run_end] == squared_distances[run_begin]) {
            ++run_end;
        }
        if (run_end - run_begin > 1) {
            // Insertion sort the tied run by index; squared distances are equal
            // inside the run, so the distance row needs no permutation.
            for (int a = run_begin + 1; a < run_end; ++a) {
                const int v = indices[a];
                int b = a;
                while (b > run_begin && indices[b - 1] > v) {
                    indices[b] = indices[b - 1];
                    --b;
                }
                indices[b] = v;
            }
        }
        run_begin = run_end;
    }
}

KnnResult knn_parallel(
    const Eigen::MatrixXd& tree_points,
    const Eigen::MatrixXd& query_points,
    int k,
    bool sorted
) {
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
    KDTree tree = make_kdtree(tree_points, adaptor);

    parallel_for(0, n_query, 256, [&](int i) {
        int* row_indices = result.idx_row(i);
        double* row_distances = result.dist_row(i);

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

    // One kd-tree over all points (same tree knn_parallel builds), queries in
    // blocks so the scratch n_block x k buffers stay bounded: the whole-slide
    // confidence step reads only the kth distance per molecule, so keeping the
    // full n x k result would peak at n * k * 12 bytes (3 GiB at 10.6M).
    EigenColMajorAdaptor adaptor(points);
    KDTree tree = make_kdtree(points, adaptor);

    constexpr std::int64_t scratch_budget_bytes = 32 << 20;
    std::int64_t block = scratch_budget_bytes / (static_cast<std::int64_t>(k) * 12);
    block = std::max<std::int64_t>(block, 4096);
    block = std::min<std::int64_t>(block, 65536);
    block = std::min<std::int64_t>(block, n);

    std::vector<int> scratch_indices(static_cast<std::size_t>(block) * k);
    std::vector<double> scratch_distances(static_cast<std::size_t>(block) * k);
    std::vector<double> out(n);

    for (int block_begin = 0; block_begin < n; block_begin += static_cast<int>(block)) {
        const int block_end = static_cast<int>(
            std::min<std::int64_t>(block_begin + block, n));
        parallel_for(block_begin, block_end, 256, [&](int i) {
            const int local = i - block_begin;
            int* row_indices = scratch_indices.data() + static_cast<std::size_t>(local) * k;
            double* row_distances = scratch_distances.data() + static_cast<std::size_t>(local) * k;

            nanoflann::KNNResultSet<double, int> resultSet(k);
            resultSet.init(row_indices, row_distances);
            tree.findNeighbors(
                resultSet,
                points.col(i).data(),
                nanoflann::SearchParameters(/*eps=*/0.0f, /*sorted=*/true)
            );
            // Only the distance value is read here: tied runs may keep any
            // order, the kth smallest distance is the same either way.
            out[i] = std::sqrt(row_distances[kth]);
        });
    }

    return out;
}

} // namespace baysor
