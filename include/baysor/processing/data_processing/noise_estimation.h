#pragma once

#include "baysor/processing/models/adj_list.h"
#include "baysor/processing/data_processing/triangulation.h"
#include <vector>
#include <Eigen/Dense>

namespace baysor {

struct MoleculeData;

/// Result of noise probability fitting
struct NoiseFitResult {
    Eigen::MatrixXd assignment_probs;  // n_molecules x 2 (signal, noise)
    std::vector<int> assignment;        // 1 = signal, 2 = noise
    double signal_mu, signal_sigma;
    double noise_mu, noise_sigma;
    std::vector<double> diffs;
};

/// Full details of confidence estimation, including per-molecule distances
/// used by the noise model and the fitted posterior.
struct ConfidenceEstimationDetails {
    std::vector<double> edge_lengths;
    NoiseFitResult fit_result;
    int nn_id = 0;
    /// Unfiltered molecule-graph edges of the noise-model MRF, for reuse by
    /// build_molecule_graph()
    AdjacencyResult adjacency;
};

/// Fit two-component mixture (signal vs noise) on KNN distances using MRF-regularized EM
NoiseFitResult fit_noise_probabilities(
    const std::vector<double>& edge_lengths,
    const AdjList& adj_list,
    const std::vector<double>* min_confidence = nullptr,
    int max_iters = 10000,
    double tol = 0.005,
    bool verbose = false
);

/// Estimate per-molecule confidence scores and append to MoleculeData
void append_confidence(
    MoleculeData& data,
    int nn_id,
    double prior_confidence = 0.5
);

/// Estimate confidence and return the full noise-model diagnostics.
ConfidenceEstimationDetails estimate_confidence_details(
    const MoleculeData& data,
    int nn_id,
    double prior_confidence = 0.5
);

} // namespace baysor
