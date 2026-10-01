#pragma once

#include "baysor/processing/models/bmm_data.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/processing/distributions/mv_normal.h"
#include "baysor/processing/data_processing/triangulation.h"
#include <Eigen/Dense>
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

/// Build the molecule adjacency graph (MRF)
AdjList build_molecule_graph(
    const MoleculeData& data,
    bool filter = true,
    bool use_local_gene_similarities = false,
    AdjacencyType type = AdjacencyType::Auto,
    int composition_neighborhood = 0,
    int n_gene_pcs = 0
);

/// Full initialization of BmmData from molecule data. The adjacency list is
/// taken by value and moved into BmmData: pass an rvalue when the caller no
/// longer needs it, so that only one copy of the molecule graph exists.
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
