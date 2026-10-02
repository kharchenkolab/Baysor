#pragma once

#include "baysor/data_loading/data.h"
#include "baysor/processing/data_processing/boundary_estimation.h"
#include "baysor/processing/data_processing/noise_estimation.h"
#include "baysor/reporting/color_utils.h"
#include "baysor/utils/options.h"
#include <nlohmann/json.hpp>
#include <cstdint>
#include <string>
#include <vector>

namespace baysor {

/// Width in pixels of a molecule scatter of (x, y) whose longer side is
/// max_size_px (height follows the data aspect ratio as in rasterize_scatter,
/// so tall data gets a narrower image). Values < 1 mean the default
/// `[plotting] max_plot_size`.
int scatter_width_for_max_size(
    const std::vector<double>& x,
    const std::vector<double>& y,
    int max_size_px
);

/// RGB raster of a scatter plot: row-major, 3 bytes per pixel.
struct ScatterRaster {
    std::vector<uint8_t> pixels;
    int width_px = 0;
    int height_px = 0;

    bool empty() const { return pixels.empty(); }
};

/// Rasterize all molecules, coloured by the given hex strings, with optional
/// polygon outlines. width_px controls the output width; height is derived
/// from the data aspect ratio (at most 4x the width). Empty input gives an
/// empty raster.
ScatterRaster rasterize_scatter(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<std::string>& colors,
    const PolygonCollection* polygons = nullptr,
    int width_px = 6000,
    int point_radius_px = 0
);

/// Rasterize molecules coloured by confidence (blue-orange gradient).
ScatterRaster rasterize_confidence(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& confidence,
    int width_px = 6000,
    int point_radius_px = 0
);

/// Encode rasters as base64 PNG data URIs ("" for an empty raster), one
/// thread-pool task per image.
std::vector<std::string> encode_png_data_uris(const std::vector<ScatterRaster>& rasters);

/// Generate Vega-Lite spec: noise estimation histogram + fitted PDFs
nlohmann::json vega_noise_histogram(
    const std::vector<double>& edge_lengths,
    const std::vector<double>& confidence,
    double signal_mu, double signal_sigma,
    double noise_mu, double noise_sigma,
    int nn_id = 0,
    int n_bins = 50
);

/// Generate Vega-Lite spec: gene frequency stacked bar chart
nlohmann::json vega_gene_frequency(
    const std::vector<int>& genes,
    const std::vector<double>& confidence,
    const std::vector<std::string>& gene_names
);

/// Generate Vega-Lite spec: gene structure scatter (UMAP of genes by spatial co-occurrence)
nlohmann::json vega_gene_structure(const GeneStructureEmbedding& emb);

/// Generate complete HTML preview report. The molecule images are at most
/// max_plot_size pixels on their longer side.
std::string generate_preview_html(
    const MoleculeData& data,
    const std::vector<std::string>& gene_colors,
    const std::vector<double>& edge_lengths,
    const NoiseFitResult& noise_result,
    int confidence_nn_id,
    const GeneStructureEmbedding* gene_structure = nullptr,
    int max_plot_size = PlottingOptions{}.max_plot_size
);

} // namespace baysor
