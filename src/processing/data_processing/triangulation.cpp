#include "baysor/processing/data_processing/triangulation.h"
#include "baysor/processing/utils/utils.h"
#include "baysor/utils/general.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <unordered_set>

namespace baysor {

// ============================================================================
// normalize_points
// ============================================================================

Eigen::MatrixXd normalize_points(const Eigen::MatrixXd& points) {
    const int dims = static_cast<int>(points.rows());
    const int n = static_cast<int>(points.cols());
    if (n == 0) return points;

    Eigen::MatrixXd out = points;

    // Shift: subtract per-row minimum
    for (int d = 0; d < dims; ++d) {
        double mn = out.row(d).minCoeff();
        out.row(d).array() -= mn;
    }

    // Scale: divide by global max * 1.1, then shift to [1.01, ~1.9]
    double mx = out.maxCoeff();
    if (mx > 0) {
        out /= mx * 1.1;
    }
    out.array() += 1.01;

    // Jitter duplicate points (KNN with distance=0 creates self-loops)
    if (n > 1) {
        auto knn = knn_parallel(out, out, 2, true);
        auto& rng = global_xoshiro_rng();
        for (int i = 0; i < n; ++i) {
            if (knn.k >= 2 && knn.dist_row(i)[1] < 1e-6) {
                for (int d = 0; d < dims; ++d) {
                    out(d, i) += (rng.rand_float64() - 0.5) * 2e-5;
                }
            }
        }
    }

    return out;
}

// ============================================================================
// filter_long_edges
// ============================================================================

void filter_long_edges(AdjacencyResult& result, double n_mads) {
    int n = static_cast<int>(result.edge_dists.size());
    if (n == 0) return;

    // Compute log10 distances
    std::vector<double> log_dists(n);
    for (int i = 0; i < n; ++i) {
        log_dists[i] = std::log10(std::max(result.edge_dists[i], 1e-30));
    }

    // Median
    std::vector<double> sorted_log = log_dists;
    std::sort(sorted_log.begin(), sorted_log.end());
    double median;
    if (n % 2 == 0) {
        median = (sorted_log[n / 2 - 1] + sorted_log[n / 2]) / 2.0;
    } else {
        median = sorted_log[n / 2];
    }

    // MAD (median absolute deviation), normalized
    std::vector<double> abs_devs(n);
    for (int i = 0; i < n; ++i) {
        abs_devs[i] = std::abs(log_dists[i] - median);
    }
    std::sort(abs_devs.begin(), abs_devs.end());
    double mad;
    if (n % 2 == 0) {
        mad = (abs_devs[n / 2 - 1] + abs_devs[n / 2]) / 2.0;
    } else {
        mad = abs_devs[n / 2];
    }
    mad *= 1.4826; // normalize=true (consistent estimator of std)

    double threshold = median + n_mads * mad;

    // Filter
    AdjacencyResult filtered;
    for (int i = 0; i < n; ++i) {
        if (log_dists[i] < threshold) {
            filtered.edge_src.push_back(result.edge_src[i]);
            filtered.edge_dst.push_back(result.edge_dst[i]);
            filtered.edge_dists.push_back(result.edge_dists[i]);
        }
    }
    result = std::move(filtered);
}

// ============================================================================
// adjacency_list
// ============================================================================

AdjacencyResult adjacency_list(
    const Eigen::MatrixXd& points,
    bool filter,
    double n_mads,
    int k_adj,
    AdjacencyType type
) {
    const int dims = static_cast<int>(points.rows());
    const int n = static_cast<int>(points.cols());

    if (n <= 1) return {};

    // Auto-select type
    if (type == AdjacencyType::Auto) {
        type = (dims == 3) ? AdjacencyType::Knn : AdjacencyType::Triangulation;
    }
    if (dims == 3 && type == AdjacencyType::Triangulation) {
        type = AdjacencyType::Knn; // 3D only supports KNN
    }

    Eigen::MatrixXd norm_pts = normalize_points(points);

    std::vector<std::pair<int, int>> tri_edges;
    std::vector<std::pair<int, int>> knn_edges;

    // --- Triangulation (2D only, using CGAL with vertex info) ---
    if (dims == 2 && (type == AdjacencyType::Triangulation || type == AdjacencyType::Both)) {
        // We'll use CGAL with vertex info to store our original index
        using K = CGAL::Exact_predicates_inexact_constructions_kernel;

        // Use Triangulation_vertex_base_with_info to store the original index
        using Vb = CGAL::Triangulation_vertex_base_with_info_2<int, K>;
        using Fb = CGAL::Triangulation_face_base_2<K>;
        using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
        using DT = CGAL::Delaunay_triangulation_2<K, Tds>;
        using Point = K::Point_2;
        using PointWithInfo = std::pair<Point, int>;

        std::vector<PointWithInfo> indexed_pts(n);
        for (int i = 0; i < n; ++i) {
            indexed_pts[i] = {Point(norm_pts(0, i), norm_pts(1, i)), i};
        }

        DT dt;
        dt.insert(indexed_pts.begin(), indexed_pts.end());

        // Extract edges
        for (auto eit = dt.finite_edges_begin(); eit != dt.finite_edges_end(); ++eit) {
            auto face = eit->first;
            int idx = eit->second;
            auto v1 = face->vertex((idx + 1) % 3);
            auto v2 = face->vertex((idx + 2) % 3);
            int i1 = v1->info();
            int i2 = v2->info();
            int lo = std::min(i1, i2);
            int hi = std::max(i1, i2);
            tri_edges.push_back({lo, hi});
        }
    }

    // --- KNN ---
    if (type == AdjacencyType::Knn || type == AdjacencyType::Both) {
        auto knn = knn_parallel(norm_pts, norm_pts, k_adj + 1, true);
        for (int i = 0; i < n; ++i) {
            // Skip self (index 0 is the point itself when sorted)
            const int* row = knn.idx_row(i);
            for (int j = 1; j < knn.k; ++j) {
                int nb = row[j];
                int lo = std::min(i, nb);
                int hi = std::max(i, nb);
                knn_edges.push_back({lo, hi});
            }
        }
    }

    std::vector<std::pair<int, int>> ordered_edges;
    if (type == AdjacencyType::Triangulation) {
        ordered_edges = std::move(tri_edges);
    } else if (type == AdjacencyType::Knn) {
        ordered_edges = std::move(knn_edges);
    } else {
        ordered_edges.reserve(knn_edges.size() + tri_edges.size());
        ordered_edges.insert(ordered_edges.end(), knn_edges.begin(), knn_edges.end());
        ordered_edges.insert(ordered_edges.end(), tri_edges.begin(), tri_edges.end());
    }

    // TODO(parity): Julia keeps first-occurrence ordering of the incoming
    // edge stream. Before dedup we canonically sort the edges so the result
    // is a pure function of the edge *set* (see sort rationale below), which
    // for a sorted stream is still well-defined "first occurrence" ordering.
    //
    // Rationale for the sort: CGAL's finite_edges iterator emits an edge from
    // whichever of its two adjacent faces has the LOWER HEAP ADDRESS
    // (Triangulation_ds_iterators_2.h: associated_edge() compares raw
    // Face_handle pointers). The emission order therefore depends on where
    // the triangulation's blocks happen to land in the heap, which shifts
    // with the length of the `-o` output path (early std::string chunk sizes)
    // and with allocator history (e.g. concurrent parquet decoding threads).
    // That order flows into the CSR adjacency lists and hence into the
    // floating-point summation order of the MRF E-step, where last-bit weight
    // differences flip stochastic assignments and make 1-thread runs
    // depend on the output path length. Sorting removes the layout
    // dependence entirely.
    std::sort(ordered_edges.begin(), ordered_edges.end());

    // Julia keeps the first occurrence of each undirected edge.
    std::unordered_set<std::uint64_t> seen;
    seen.reserve(ordered_edges.size() * 2 + 1);

    AdjacencyResult result;
    result.edge_src.reserve(ordered_edges.size());
    result.edge_dst.reserve(ordered_edges.size());
    result.edge_dists.reserve(ordered_edges.size());

    for (const auto& edge : ordered_edges) {
        const int lo = edge.first;
        const int hi = edge.second;
        const std::uint64_t key =
            (static_cast<std::uint64_t>(static_cast<std::uint32_t>(lo)) << 32)
            | static_cast<std::uint32_t>(hi);
        if (!seen.insert(key).second) {
            continue;
        }
        double dist = (norm_pts.col(lo) - norm_pts.col(hi)).norm();
        result.edge_src.push_back(lo);
        result.edge_dst.push_back(hi);
        result.edge_dists.push_back(dist);
    }

    if (filter) {
        filter_long_edges(result, n_mads);
    }

    return result;
}

} // namespace baysor
