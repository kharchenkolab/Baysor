#pragma once

#include "baysor/processing/models/adj_list.h"
#include "baysor/reporting/color_utils.h"
#include "baysor/utils/options.h"
#include <Eigen/Dense>
#include <memory>
#include <vector>

namespace baysor {

/// Result of MRF-based molecule clustering
struct ClusteringResult {
    Eigen::MatrixXd exprs;              // n_clusters × n_genes expression profiles
    std::vector<int> assignment;         // per-molecule cluster ID (1-based)
    Eigen::MatrixXd assignment_probs;    // n_clusters × n_molecules soft assignments
    std::vector<double> diffs;           // max_diff convergence trace
    std::vector<double> change_fracs;    // fraction-changed convergence trace
    std::shared_ptr<NcvProjectedModel> ncv_projected_model;
};

struct GraphClusteringSummary {
    int micro_clusters = 0;
    int final_clusters = 0;
    double chosen_resolution = 1.0;
    std::vector<double> move_fracs;
};

struct ClusteringOptions {
    ClusterMethod method = ClusterMethod::Mrf;
    int n_clusters = 4;
    double resolution = 1.0;
    int graph_k = 15;
    int spatial_k = 0;
    int n_dims = 20;
    int basis_sample_size = 100000;
    double tol = 0.01;
    double mrf_weight = 1.0;
    int max_iters = -1;
};

ClusteringResult cluster_molecules(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    const AdjList& adj_list,
    const std::vector<double>& confidence,
    const ClusteringOptions& options,
    bool verbose = true
);

std::vector<int> louvain_partition(
    const AdjList& graph,
    double resolution = 1.0,
    int max_passes = 100
);

std::vector<int> leiden_partition(
    const AdjList& graph,
    double resolution = 1.0,
    int max_passes = 100
);

AdjList build_knn_similarity_graph(
    const Eigen::MatrixXf& mol_vecs,
    const std::vector<double>& confidence,
    int k
);

std::vector<int> graph_partition_to_target(
    const AdjList& graph,
    const Eigen::MatrixXf& mol_vecs,
    const std::vector<double>& confidence,
    ClusterMethod method,
    int target_clusters,
    double resolution_seed = 1.0,
    int max_passes = 100,
    GraphClusteringSummary* summary = nullptr
);

/// Core MRF-based categorical EM for molecule clustering.
/// exprs_init: optional pre-computed n_clusters × n_genes expression profiles
///             (e.g. from ICA). Pass nullptr to use deterministic hash fallback.
ClusteringResult cluster_molecules_on_mrf(
    const std::vector<int>& genes,      // 1-based gene IDs
    const AdjList& adj_list,
    const std::vector<double>& confidence,
    int n_clusters,
    double tol       = 0.01,
    double mrf_weight = 1.0,
    int max_iters    = -1,              // -1 = auto (max(10000, n_mols/200))
    bool verbose     = true,
    const Eigen::MatrixXd* exprs_init = nullptr  // nullptr → hash-based fallback
);

/// Gene panels up to this size whiten the ICA input with the dense
/// covariance eigen-decomposition (Julia's algorithm, O(genes^3): ~1.5 min of
/// CPU and ~0.4 GB at 3,000 genes); larger panels use the top-k eigenpairs
/// from a truncated SVD (irlba). Both give the same eigenpairs up to the
/// eigenvector signs, which are arbitrary in the dense solver and fixed in the
/// truncated one, so the threshold also decides which panels keep their
/// previous FastICA start: on a 2,672-gene Xenium 5K crop the truncated path
/// changes the segmentation by about as much as a different thread count.
constexpr int ica_dense_whitening_max_genes = 3000;

/// ICA-initialized molecule clustering (matches Julia's DataFrame wrapper).
/// Computes pairwise gene spatial co-occurrence → FastICA → initializes EM.
/// Falls back to hash initialization if ICA fails to converge.
ClusteringResult cluster_molecules_ica(
    const std::vector<int>& genes,
    const AdjList& adj_list,
    const std::vector<double>& confidence,
    int n_clusters,
    double tol       = 0.01,
    double mrf_weight = 1.0,
    int max_iters    = -1,
    bool verbose     = true,
    int dense_whitening_max_genes = ica_dense_whitening_max_genes
);

namespace detail {

/// Top-k eigenpairs of the row covariance of an ICA input matrix, in
/// decreasing eigenvalue order.
struct IcaWhitening {
    Eigen::VectorXd eigenvalues;   // k
    Eigen::MatrixXd eigenvectors;  // n_features × k
};

/// Dense covariance + full symmetric eigen-decomposition (Julia's whitening).
IcaWhitening ica_whitening_dense(const Eigen::MatrixXd& X, int k);

/// Truncated SVD of the implicitly centred X; each eigenvector's
/// largest-magnitude entry is positive.
IcaWhitening ica_whitening_truncated(const Eigen::MatrixXd& X, int k);

} // namespace detail

ClusteringResult cluster_molecules_louvain(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    const AdjList& adj_list,
    const std::vector<double>& confidence,
    double resolution = 1.0,
    int graph_k = 15,
    int spatial_k = 0,
    int target_clusters = 4,
    int n_dims = 20,
    int basis_sample_size = 100000,
    bool verbose = true
);

ClusteringResult cluster_molecules_leiden(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    const AdjList& adj_list,
    const std::vector<double>& confidence,
    double resolution = 1.0,
    int graph_k = 15,
    int spatial_k = 0,
    int target_clusters = 4,
    int n_dims = 20,
    int basis_sample_size = 100000,
    bool verbose = true
);

} // namespace baysor
