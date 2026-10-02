#pragma once

#include "baysor/data_loading/data.h"
#include "baysor/processing/bmm_algorithm/molecule_clustering.h"
#include "baysor/processing/data_processing/noise_estimation.h"
#include "baysor/processing/data_processing/boundary_estimation.h"
#include "baysor/reporting/color_utils.h"
#include "baysor/utils/options.h"

#include <Eigen/Dense>
#include <string>
#include <unordered_map>
#include <vector>

namespace baysor {

/// Generate the post-segmentation diagnostic report used by `run --plot`.
std::string generate_run_diagnostic_html(
    const MoleculeData& data,
    const std::vector<double>& edge_lengths,
    const NoiseFitResult& noise_result,
    int confidence_nn_id,
    const std::vector<int>& assignment,
    const std::vector<std::unordered_map<int, int>>& n_components_trace,
    const std::vector<double>& assignment_confidence,
    const ClusteringResult* clustering_result,
    const NcvReportEmbedding* ncv_report,
    const Eigen::MatrixXd& cell_stats,
    const std::vector<std::string>& cell_stat_col_names,
    const PriorInputOptions& prior_opts,
    double scale,
    const std::string& scale_std
);

/// Generate the post-segmentation molecule/border report used by `run --plot`.
/// The molecule images are at most max_plot_size pixels on their longer side
/// (`[plotting] max_plot_size`).
std::string generate_run_segmentation_html(
    const MoleculeData& data,
    const std::vector<int>& assignment,
    const std::vector<std::string>& ncv_color,
    const std::vector<int>* molecule_clusters,
    const PolygonCollection* polygons,
    int max_plot_size = PlottingOptions{}.max_plot_size
);

} // namespace baysor
