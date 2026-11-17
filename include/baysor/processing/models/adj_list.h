#pragma once

#include <vector>
#include <cstdint>

namespace baysor {

/// CSR-format adjacency list for the molecule graph (MRF).
///
/// Replaces Julia's Vector{Vector{Int}} with a flat, cache-friendly layout.
/// For molecule i, neighbors are indices[indptr[i]..indptr[i+1]),
/// with corresponding weights[indptr[i]..indptr[i+1]).
///
/// This layout is directly usable as a sparse matrix for GPU SpMM operations.
struct AdjList {
    // Row offsets are 64-bit: a molecule graph has 2 x n_edges directed edges
    // (~6 x n_molecules for a 2D Delaunay graph), which exceeds INT32_MAX at
    // ~358M molecules. Molecule ids (indices) still fit in int32.
    std::vector<int64_t> indptr;    // size n_molecules + 1
    std::vector<int32_t> indices;   // size nnz (neighbor molecule IDs)
    std::vector<double> weights;    // size nnz (edge weights)

    int n_molecules() const { return static_cast<int>(indptr.size()) - 1; }
    int64_t nnz() const { return static_cast<int64_t>(indices.size()); }

    // Access neighbors of molecule i
    int neighbor_count(int i) const { return static_cast<int>(indptr[i + 1] - indptr[i]); }
    const int32_t* neighbor_ids(int i) const { return indices.data() + indptr[i]; }
    const double* neighbor_weights(int i) const { return weights.data() + indptr[i]; }

    /// Build from edge list (2 x n_edges) and edge weights
    static AdjList from_edge_list(
        const int* edges_src, const int* edges_dst,
        const double* edge_weights, int64_t n_edges, int n_verts
    );
};

} // namespace baysor
