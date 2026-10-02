#pragma once

#include "baysor/processing/models/bmm_data.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/processing/distributions/mv_normal.h"
#include "baysor/processing/data_processing/triangulation.h"
#include <Eigen/Dense>
#include <optional>
#include <vector>

namespace baysor {

struct MoleculeData;

/// Initial cell parameters (centers + covariances + assignment)
template<int N>
struct InitialParams {
    using Mat = Eigen::Matrix<double, N, N>;
    Eigen::MatrixXd centers;             // n_cells x N
    std::vector<Mat> covs;               // per-cell covariance
    std::vector<int> assignment;          // per-molecule initial cell ID
};

/// Place initial cell centers uniformly in space
template<int N>
InitialParams<N> cell_centers_uniformly(
    const Eigen::MatrixXd& pos_data,
    int n_clusters,
    const std::vector<double>* confidences = nullptr,
    double scale = -1.0
);

/// Compute the molecule-graph edge list once (normalized points, Delaunay
/// triangulation / kNN, canonical dedup), unfiltered. The result can be passed
/// to build_molecule_graph() to avoid recomputing the triangulation.
AdjacencyResult compute_molecule_adjacency(const MoleculeData& data);

/// Build the weighted molecule graph (CSR AdjList) from edges produced by
/// compute_molecule_adjacency() or adjacency_list().
AdjList build_molecule_graph_from_edges(const AdjacencyResult& edges, int n_molecules);

/// Build the molecule adjacency graph (MRF)
AdjList build_molecule_graph(
    const MoleculeData& data,
    bool filter = true,
    bool use_local_gene_similarities = false,
    AdjacencyType type = AdjacencyType::Auto,
    int composition_neighborhood = 0,
    int n_gene_pcs = 0,
    /// Unfiltered edges from compute_molecule_adjacency() to reuse (then
    /// filtered here); ignored when the data has duplicate coordinates.
    std::optional<AdjacencyResult> precomputed_edges = std::nullopt
);

/// Full initialization of BmmData from molecule data (adj_list is moved in)
template<int N>
BmmData<N> initialize_bmm_data(
    const MoleculeData& mol_data,
    AdjList adj_list,
    int n_cells_init,
    double scale,
    const std::string& scale_std = "25%",
    double prior_seg_confidence = 0.5,
    int min_molecules_per_cell = 3,
    bool verbose = true
);

} // namespace baysor
