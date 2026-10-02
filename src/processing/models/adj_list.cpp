#include "baysor/processing/models/adj_list.h"
#include <algorithm>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace baysor {

AdjList AdjList::from_edge_list(
    const int* edges_src, const int* edges_dst,
    const double* edge_weights, int n_edges, int n_verts
) {
    // indptr has n_verts + 1 entries and there are 2 * n_edges directed
    // edges, all indexed with int
    constexpr int kMax = std::numeric_limits<int>::max();
    if (n_verts < 0)
        throw std::invalid_argument("AdjList::from_edge_list: n_verts must be non-negative");
    if (n_edges < 0)
        throw std::invalid_argument("AdjList::from_edge_list: n_edges must be non-negative");
    if (n_verts == kMax)
        throw std::length_error("AdjList::from_edge_list: n_verts is too large for 32-bit indices");
    if (n_edges > kMax / 2)
        throw std::length_error("AdjList::from_edge_list: n_edges is too large for 32-bit indices");

    if (n_edges == 0 || n_verts == 0) {
        AdjList adj; // GCOVR_EXCL_LINE: exception-cleanup block for this declaration is only reachable if `indptr.assign` throws inside this frame (allocation failure); no portable test input can induce that
        adj.indptr.assign(n_verts + 1, 0);
        return adj;
    }

    // Match Julia's convert_edge_list_to_adj_list ordering:
    // for each vertex, first edges where it is in edge_src (in input order),
    // then edges where it is in edge_dst (also in input order).
    AdjList adj;
    adj.indptr.assign(n_verts + 1, 0);
    for (int i = 0; i < n_edges; ++i) {
        adj.indptr[edges_src[i] + 1]++;
        adj.indptr[edges_dst[i] + 1]++;
    }

    for (int i = 1; i <= n_verts; ++i) {
        adj.indptr[i] += adj.indptr[i - 1];
    }

    const int n_directed = n_edges * 2;
    adj.indices.resize(n_directed);
    adj.weights.resize(n_directed);

    std::vector<int32_t> write_pos = adj.indptr;

    for (int i = 0; i < n_edges; ++i) {
        const int src = edges_src[i];
        const int dst = edges_dst[i];
        const double w = edge_weights[i];
        const int pos = write_pos[src]++;
        adj.indices[pos] = dst;
        adj.weights[pos] = w;
    }

    for (int i = 0; i < n_edges; ++i) {
        const int src = edges_dst[i];
        const int dst = edges_src[i];
        const double w = edge_weights[i];
        const int pos = write_pos[src]++;
        adj.indices[pos] = dst;
        adj.weights[pos] = w;
    }

    return adj;
}

} // namespace baysor
