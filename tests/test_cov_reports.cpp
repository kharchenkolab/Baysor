// Coverage tests for src/reporting/{preview_report,run_report,color_utils,gene_structure}.cpp (COV-4).
//
// Verifies the generated HTML reports contain their key sections and plot
// specs, exercises the scatter/confidence PNG renderers and the Vega-Lite
// spec builders, and covers the NCV colour-embedding variants plus the
// gene-structure analysis.

#include <gtest/gtest.h>

#include <nlohmann/json.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "baysor/data_loading/data.h"
#include "baysor/processing/bmm_algorithm/molecule_clustering.h"
#include "baysor/processing/data_processing/noise_estimation.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/reporting/color_utils.h"
#include "baysor/reporting/preview_report.h"
#include "baysor/reporting/run_report.h"
#include "baysor/utils/options.h"

namespace {

// Deliberately wide, flat coordinates: the report renderers default to a
// 6000 px width, so a large x-range keeps the raster tiny and the tests fast.
baysor::MoleculeData wide_data() {
    baysor::MoleculeData d;
    d.x = {0.0, 300.0, 600.0, 900.0};
    d.y = {0.0, 0.5, 1.0, 0.5};
    d.gene = {1, 1, 2, 2};
    d.gene_names = {"A", "B"};
    d.confidence = {0.005, 0.9, 0.9, 0.9};
    return d;
}

baysor::NoiseFitResult noise_fit() {
    baysor::NoiseFitResult nr;
    nr.assignment_probs = Eigen::MatrixXd::Zero(4, 2);
    nr.assignment = {1, 1, 2, 2};
    nr.signal_mu = 0.1;
    nr.signal_sigma = 0.02;
    nr.noise_mu = 1.1;
    nr.noise_sigma = 0.1;
    nr.diffs = {0.1};
    return nr;
}

baysor::AdjList chain_adj_list(int n_molecules) {
    std::vector<int> src, dst;
    std::vector<double> wts;
    for (int i = 0; i < n_molecules - 1; ++i) {
        src.push_back(i);
        dst.push_back(i + 1);
        wts.push_back(1.0);
    }
    return baysor::AdjList::from_edge_list(src.data(), dst.data(), wts.data(),
                                           static_cast<int>(src.size()), n_molecules);
}

Eigen::MatrixXf random_vecs(int rows, int cols, unsigned seed) {
    Eigen::MatrixXf m(rows, cols);
    std::mt19937 rng(seed);
    std::normal_distribution<float> dist(0.0f, 1.0f);
    for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c) m(r, c) = dist(rng);
    return m;
}

// Minimal base64 decoder for the "data:image/png;base64,..." URIs returned
// by the PNG renderers.
std::string base64_decode(const std::string& in) {
    static const char* kAlphabet =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out;
    out.reserve(in.size() * 3 / 4);
    int val = 0;
    int valb = -8;
    for (unsigned char c : in) {
        if (c == '=') break;  // padding: end of payload
        const char* p = std::strchr(kAlphabet, static_cast<char>(c));
        if (p == nullptr || c == '\0') break;
        val = (val << 6) + static_cast<int>(p - kAlphabet);
        valb += 6;
        if (valb >= 0) {
            out.push_back(static_cast<char>((val >> valb) & 0xFF));
            valb -= 8;
        }
    }
    return out;
}

// Verify the PNG signature of a base64 data URI and return the dimensions
// from the IHDR chunk (big-endian); {0, 0} plus a failure on any problem.
std::pair<uint32_t, uint32_t> png_ihdr(const std::string& data_uri) {
    const std::string prefix = "data:image/png;base64,";
    if (data_uri.rfind(prefix, 0) != 0) {
        ADD_FAILURE() << "missing PNG data-URI prefix: " << data_uri.substr(0, 64);
        return {0, 0};
    }
    const std::string bytes = base64_decode(data_uri.substr(prefix.size()));
    static const uint8_t kSig[8] = {0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
    if (bytes.size() < 24 || std::memcmp(bytes.data(), kSig, sizeof(kSig)) != 0) {
        ADD_FAILURE() << "decoded payload is not a PNG (bad signature or truncated, "
                      << bytes.size() << " bytes)";
        return {0, 0};
    }
    if (bytes.compare(12, 4, "IHDR") != 0) {
        ADD_FAILURE() << "first PNG chunk is not IHDR";
        return {0, 0};
    }
    auto be32 = [&](size_t off) {
        return (static_cast<uint32_t>(static_cast<uint8_t>(bytes[off])) << 24) |
               (static_cast<uint32_t>(static_cast<uint8_t>(bytes[off + 1])) << 16) |
               (static_cast<uint32_t>(static_cast<uint8_t>(bytes[off + 2])) << 8) |
               static_cast<uint32_t>(static_cast<uint8_t>(bytes[off + 3]));
    };
    return {be32(16), be32(20)};
}

} // namespace

// ============================================================================
// render_scatter_png / render_confidence_png
// ============================================================================

TEST(Cov4PreviewRender, ScatterPngWithPolygonsAndInvalidHexDigits) {
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> y = {0.0, 1.0, 0.0, 1.0};
    // 'G', 'z' are not hex digits -> parsed as 0 (exercises the fallback digit).
    std::vector<std::string> colors = {"#FF0000", "#GG0000", "#00ffzz", "#123456"};

    baysor::PolygonCollection polygons;
    Eigen::MatrixXd quad_cols(2, 4);
    quad_cols << 0.0, 1.0, 1.0, 0.0,
                 0.0, 0.0, 1.0, 1.0;
    polygons["c1"] = quad_cols;
    Eigen::MatrixXd quad_rows(4, 2);
    quad_rows << 0.0, 1.0,
                 1.0, 0.0,
                 0.0, 0.0,
                 1.0, 1.0;
    polygons["c2"] = quad_rows;
    Eigen::MatrixXd degenerate(3, 3);  // neither 2xN nor Nx2 -> skipped
    degenerate.setZero();
    polygons["c3"] = degenerate;

    auto png = baysor::render_scatter_png(x, y, colors, &polygons,
                                          /*width_px=*/64);
    ASSERT_GE(png.size(), 22u);
    EXPECT_EQ(png.substr(0, 22), "data:image/png;base64,");
    // Decode and check the PNG signature plus the IHDR dimensions:
    // 64 px wide, height = 64 * (yrange=1) / (xrange=3) = 21 px.
    const auto dims = png_ihdr(png);
    EXPECT_EQ(dims.first, 64u);
    EXPECT_EQ(dims.second, 21u);

    // Empty input renders nothing at all.
    EXPECT_EQ(baysor::render_scatter_png({}, {}, {}), "");
}

TEST(Cov4PreviewRender, ConfidencePngCoversBothBranchesAndClamping) {
    std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> y = {0.0, 1.0, 0.0, 1.0, 0.5};
    // 0.0/0.4 -> orange branch, -0.5 -> clamped low, 0.6 -> blue branch,
    // 1.5 -> clamped high.
    std::vector<double> confidence = {0.0, 0.4, -0.5, 0.6, 1.5};

    auto png = baysor::render_confidence_png(x, y, confidence,
                                             /*width_px=*/64);
    ASSERT_GE(png.size(), 22u);
    EXPECT_EQ(png.substr(0, 22), "data:image/png;base64,");
    // Decode and check the PNG signature plus the IHDR dimensions:
    // 64 px wide, height = 64 * (yrange=1) / (xrange=4) = 16 px.
    const auto dims = png_ihdr(png);
    EXPECT_EQ(dims.first, 64u);
    EXPECT_EQ(dims.second, 16u);

    EXPECT_EQ(baysor::render_confidence_png({}, {}, {}), "");
}

// ============================================================================
// Vega-Lite spec builders
// ============================================================================

TEST(Cov4PreviewSpec, GeneFrequencyCountsSortsAndFilters) {
    std::vector<std::string> names = {"B", "A", "Rare"};
    std::vector<int> genes;
    for (int i = 0; i < 60; ++i) genes.push_back(1);   // B
    for (int i = 0; i < 50; ++i) genes.push_back(2);   // A
    genes.push_back(3);   // Rare: 1 molecule < 1% of 113 -> filtered out
    genes.push_back(0);   // invalid gene id 0 -> skipped
    genes.push_back(4);   // out-of-range gene id -> skipped
    std::vector<double> confidence(genes.size(), 0.9);
    for (int i = 40; i < 60; ++i) confidence[i] = 0.2;  // 20 B molecules are noise

    auto spec = baysor::vega_gene_frequency(genes, confidence, names);

    EXPECT_EQ(spec["title"], "Gene frequency");
    EXPECT_EQ(spec["mark"], "bar");
    const auto& values = spec["data"]["values"];
    ASSERT_EQ(values.size(), 4u);
    // Sorted alphabetically: A before B; Rare filtered; ids 0/4 skipped.
    EXPECT_EQ(values[0]["gene"], "A");
    EXPECT_EQ(values[0]["count"], 50);
    EXPECT_EQ(values[0]["type"], "Real");
    EXPECT_EQ(values[1]["gene"], "A");
    EXPECT_EQ(values[1]["count"], 0);
    EXPECT_EQ(values[1]["type"], "Noise");
    EXPECT_EQ(values[2]["gene"], "B");
    EXPECT_EQ(values[2]["count"], 40);
    EXPECT_EQ(values[2]["type"], "Real");
    EXPECT_EQ(values[3]["gene"], "B");
    EXPECT_EQ(values[3]["count"], 20);
    EXPECT_EQ(values[3]["type"], "Noise");
}

TEST(Cov4PreviewSpec, NoiseHistogramEmptyInputIsNull) {
    auto spec = baysor::vega_noise_histogram({}, {}, 0.1, 0.05, 1.0, 0.3, 6, 10);
    EXPECT_TRUE(spec.is_null());
}

TEST(Cov4PreviewSpec, NoiseHistogramLayersAndSeries) {
    std::vector<double> edges;
    for (int i = 0; i < 20; ++i) edges.push_back(0.1 + 0.1 * i);
    std::vector<double> confidence(20, 0.8);

    auto spec = baysor::vega_noise_histogram(edges, confidence,
                                             /*signal_mu=*/0.5, /*signal_sigma=*/0.2,
                                             /*noise_mu=*/1.5, /*noise_sigma=*/0.4,
                                             /*nn_id=*/6, /*n_bins=*/10);

    EXPECT_EQ(spec["title"], "Noise estimation");
    ASSERT_TRUE(spec.contains("layer"));
    EXPECT_EQ(spec["layer"].size(), 3u);
    EXPECT_EQ(spec["data"]["values"].size(), 30u);  // 3 series x 10 bins
    EXPECT_EQ(spec["layer"][0]["encoding"]["x"]["title"],
              "Distance to 6th nearest neighbor");
    std::set<std::string> types;
    for (const auto& v : spec["data"]["values"]) types.insert(v["type"]);
    EXPECT_EQ(types, (std::set<std::string>{"Observed", "Intracellular", "Background"}));
    EXPECT_EQ(spec["resolve"]["legend"]["color"], "shared");
}

TEST(Cov4PreviewSpec, GeneStructureSpecClampsMarkerSizes) {
    baysor::GeneStructureEmbedding emb;
    emb.x = {0.1, 0.25};
    emb.y = {0.4, 0.2};
    emb.gene_names = {"G1", "G2"};
    emb.marker_sizes = {0.5, 2.5};  // 0.5 is clamped up to 1.0

    auto spec = baysor::vega_gene_structure(emb);

    EXPECT_EQ(spec["title"], "Gene structure");
    const auto& values = spec["data"]["values"];
    ASSERT_EQ(values.size(), 2u);
    EXPECT_EQ(values[0]["gene"], "G1");
    EXPECT_DOUBLE_EQ(values[0]["size"], 1.0);
    EXPECT_EQ(values[1]["gene"], "G2");
    EXPECT_DOUBLE_EQ(values[1]["size"], 2.5);
    ASSERT_TRUE(spec.contains("layer"));
    EXPECT_EQ(spec["layer"].size(), 2u);  // points + text labels
}

// ============================================================================
// generate_preview_html
// ============================================================================

TEST(Cov4PreviewHtml, ContainsSectionsPlotsAndNoiseStats) {
    auto data = wide_data();
    auto noise = noise_fit();
    std::vector<double> edges = {0.5, 1.0, 1.5, 2.0};
    std::vector<std::string> colors = {"#111111", "#222222", "#333333", "#444444"};

    baysor::GeneStructureEmbedding emb;
    emb.x = {0.1, 0.25};
    emb.y = {0.4, 0.2};
    emb.gene_names = {"A", "B"};
    emb.marker_sizes = {1.0, 1.0};

    auto html = baysor::generate_preview_html(data, colors, edges, noise,
                                              /*confidence_nn_id=*/10, &emb);

    EXPECT_NE(html.find("Baysor Preview Report"), std::string::npos);
    EXPECT_NE(html.find("id=\"scatter_img\""), std::string::npos);
    EXPECT_NE(html.find("id=\"conf_img\""), std::string::npos);
    EXPECT_NE(html.find("data:image/png;base64,"), std::string::npos);
    EXPECT_NE(html.find("openFullRes"), std::string::npos);

    // Noise statistics: 1 of 4 molecules below 0.01 confidence.
    EXPECT_NE(html.find("Minimal noise level: 25.0%"), std::string::npos);
    EXPECT_NE(html.find("Expected noise level: 32.4%"), std::string::npos);

    EXPECT_NE(html.find("vegaEmbed('#vg_noise_dist'"), std::string::npos);
    EXPECT_NE(html.find("Noise estimation"), std::string::npos);
    EXPECT_NE(html.find("vegaEmbed('#vg_gene_freq'"), std::string::npos);
    EXPECT_NE(html.find("Gene frequency"), std::string::npos);
    EXPECT_NE(html.find("id=\"vg_gene_structure\""), std::string::npos);
    EXPECT_NE(html.find("vegaEmbed('#vg_gene_structure'"), std::string::npos);
    EXPECT_NE(html.find("Gene structure"), std::string::npos);
}

TEST(Cov4PreviewHtml, GeneStructureSectionOmittedWhenAbsent) {
    auto data = wide_data();
    auto noise = noise_fit();
    std::vector<double> edges = {0.5, 1.0};
    std::vector<std::string> colors = {"#111111", "#222222", "#333333", "#444444"};

    auto without_ptr = baysor::generate_preview_html(data, colors, edges, noise, 10,
                                                     nullptr);
    EXPECT_EQ(without_ptr.find("vegaEmbed('#vg_gene_structure'"), std::string::npos);
    EXPECT_EQ(without_ptr.find("id=\"vg_gene_structure\""), std::string::npos);
    EXPECT_NE(without_ptr.find("vegaEmbed('#vg_gene_freq'"), std::string::npos);

    baysor::GeneStructureEmbedding empty_emb;  // x is empty -> section omitted
    auto with_empty = baysor::generate_preview_html(data, colors, edges, noise, 10,
                                                    &empty_emb);
    EXPECT_EQ(with_empty.find("id=\"vg_gene_structure\""), std::string::npos);
}

TEST(Cov4PreviewHtml, EmptyMoleculeSetReportsZeroNoise) {
    baysor::MoleculeData data;  // no molecules at all
    auto noise = noise_fit();
    auto html = baysor::generate_preview_html(data, {}, {}, noise, 10, nullptr);

    EXPECT_NE(html.find("Minimal noise level: 0.0%"), std::string::npos);
    EXPECT_NE(html.find("Expected noise level: 0.0%"), std::string::npos);
    EXPECT_NE(html.find("vegaEmbed('#vg_noise_dist', null"), std::string::npos);
    EXPECT_EQ(html.find("data:image/png;base64,"), std::string::npos);
}

// ============================================================================
// generate_run_diagnostic_html
// ============================================================================

namespace {

std::string diagnostic_html(const baysor::MoleculeData& data,
                            const std::vector<double>& edges,
                            const baysor::NoiseFitResult& noise,
                            const std::vector<int>& assignment,
                            const std::vector<std::unordered_map<int, int>>& trace,
                            const std::vector<double>& assign_conf,
                            const baysor::ClusteringResult* clustering,
                            const baysor::NcvReportEmbedding* ncv,
                            const Eigen::MatrixXd& cell_stats,
                            const std::vector<std::string>& cols,
                            const baysor::PriorInputOptions& prior,
                            const std::string& scale_std) {
    return baysor::generate_run_diagnostic_html(
        data, edges, noise, /*confidence_nn_id=*/10, assignment, trace, assign_conf,
        clustering, ncv, cell_stats, cols, prior, /*scale=*/4.5, scale_std);
}

} // namespace

TEST(Cov4RunHtml, DiagnosticHtmlEscapesSpecialCharacters) {
    auto data = wide_data();
    auto noise = noise_fit();
    std::vector<double> edges = {0.1, 0.12, 1.0, 1.2};
    std::vector<int> assignment = {1, 1, 0, 2};
    std::vector<std::unordered_map<int, int>> trace = {
        {{1, 2}, {2, 1}},
        {{1, 2}, {2, 2}},
    };
    std::vector<double> assign_conf = {1.0, 1.0, 0.0, 1.0};
    Eigen::MatrixXd stats(2, 3);
    stats << 1.0, 2.0, 3.0,
             4.0, 5.0, 6.0;
    std::vector<std::string> cols = {"area", "density", "elongation"};

    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Column;
    prior.column_name = R"(c1&c2<c3>"c4")";

    auto html = diagnostic_html(data, edges, noise, assignment, trace, assign_conf,
                                nullptr, nullptr, stats, cols, prior,
                                R"(s&t<d>"q")");

    EXPECT_NE(html.find("Prior type: column"), std::string::npos);
    EXPECT_NE(html.find("Prior column: c1&amp;c2&lt;c3&gt;&quot;c4&quot;"),
              std::string::npos);
    EXPECT_EQ(html.find("c1&c2<c3>"), std::string::npos);
    EXPECT_NE(html.find("scale_std=s&amp;t&lt;d&gt;&quot;q&quot;"), std::string::npos);
    EXPECT_EQ(html.find("s&t<d>"), std::string::npos);
}

TEST(Cov4RunHtml, DiagnosticHtmlPriorTypesNoneImageBoundary) {
    auto data = wide_data();
    auto noise = noise_fit();
    std::vector<double> edges = {0.1, 1.0};
    std::vector<int> assignment = {1, 1, 0, 2};
    std::vector<std::unordered_map<int, int>> trace = {{{1, 2}, {2, 2}}};
    std::vector<double> assign_conf = {1.0, 1.0, 0.0, 1.0};
    Eigen::MatrixXd stats(2, 3);
    stats << 1.0, 2.0, 3.0,
             4.0, 5.0, 6.0;
    std::vector<std::string> cols = {"area", "density", "elongation"};

    // PriorInputType::None: no prior lines at all, empty prior_segmentation.
    baysor::PriorInputOptions none_prior;
    auto none_html = diagnostic_html(data, edges, noise, assignment, trace, assign_conf,
                                     nullptr, nullptr, stats, cols, none_prior, "1.2");
    EXPECT_NE(none_html.find("Prior type: none"), std::string::npos);
    EXPECT_EQ(none_html.find("Prior source:"), std::string::npos);
    EXPECT_EQ(none_html.find("Prior column:"), std::string::npos);
    EXPECT_EQ(none_html.find("Molecules with prior label"), std::string::npos);

    // Image prior prints (escaped) source path.
    baysor::PriorInputOptions image_prior;
    image_prior.type = baysor::PriorInputType::Image;
    image_prior.path = R"(/masks/m&1.tif)";
    auto image_html = diagnostic_html(data, edges, noise, assignment, trace, assign_conf,
                                      nullptr, nullptr, stats, cols, image_prior, "1.2");
    EXPECT_NE(image_html.find("Prior type: image"), std::string::npos);
    EXPECT_NE(image_html.find("Prior source: /masks/m&amp;1.tif"), std::string::npos);

    // Boundary prior.
    baysor::PriorInputOptions boundary_prior;
    boundary_prior.type = baysor::PriorInputType::Boundary;
    boundary_prior.path = "boundaries.csv";
    auto boundary_html = diagnostic_html(data, edges, noise, assignment, trace,
                                         assign_conf, nullptr, nullptr, stats, cols,
                                         boundary_prior, "1.2");
    EXPECT_NE(boundary_html.find("Prior type: boundary"), std::string::npos);
    EXPECT_NE(boundary_html.find("Prior source: boundaries.csv"), std::string::npos);
}

TEST(Cov4RunHtml, DiagnosticHtmlConvergenceTraceAndClustering) {
    auto data = wide_data();
    auto noise = noise_fit();
    std::vector<double> edges = {0.1, 1.0};
    std::vector<int> assignment = {1, 1, 0, 2};
    std::vector<double> assign_conf = {1.0, 1.0, 0.0, 1.0};
    Eigen::MatrixXd stats(2, 3);
    stats << 1.0, 2.0, 3.0,
             4.0, 5.0, 6.0;
    std::vector<std::string> cols = {"area", "density", "elongation"};
    baysor::PriorInputOptions prior;

    // Empty trace -> vega_convergence_trace returns null json.
    auto empty_trace_html = diagnostic_html(
        data, edges, noise, assignment, /*trace=*/{}, assign_conf,
        nullptr, nullptr, stats, cols, prior, "1.2");
    EXPECT_NE(empty_trace_html.find("vegaEmbed('#vg_seg_conv', null"), std::string::npos);
    EXPECT_EQ(empty_trace_html.find("vg_clust_conv"), std::string::npos);

    // Non-empty trace + clustering with matching convergence vectors.
    std::vector<std::unordered_map<int, int>> trace = {
        {{1, 2}, {2, 1}},
        {{1, 2}, {2, 2}},
    };
    baysor::ClusteringResult clustering;
    clustering.assignment = {1, 1, 2, 2};
    clustering.diffs = {0.5, 0.1};
    clustering.change_fracs = {0.4, 0.05};
    auto html = diagnostic_html(data, edges, noise, assignment, trace, assign_conf,
                                &clustering, nullptr, stats, cols, prior, "1.2");
    EXPECT_NE(html.find("Segmentation convergence"), std::string::npos);
    EXPECT_NE(html.find("vegaEmbed('#vg_clust_conv'"), std::string::npos);
    EXPECT_NE(html.find("Molecule clustering convergence"), std::string::npos);
    EXPECT_NE(html.find("Max prob. difference"), std::string::npos);
    EXPECT_NE(html.find("Molecules changed"), std::string::npos);

    // Mismatched convergence vectors -> clustering plot skipped.
    baysor::ClusteringResult mismatched;
    mismatched.diffs = {0.5};
    mismatched.change_fracs = {};
    auto mismatch_html = diagnostic_html(data, edges, noise, assignment, trace,
                                         assign_conf, &mismatched, nullptr, stats, cols,
                                         prior, "1.2");
    EXPECT_EQ(mismatch_html.find("vg_clust_conv"), std::string::npos);
}

TEST(Cov4RunHtml, DiagnosticHtmlHistogramEdgeCases) {
    auto data = wide_data();
    data.confidence = {0.9, 0.9, 0.9, 0.9};  // all-equal -> single-bin histograms
    auto noise = noise_fit();
    std::vector<double> edges = {0.1, 1.0};
    std::vector<int> assignment = {1, 2, 3, 4};  // every cell has 1 molecule
    std::vector<std::unordered_map<int, int>> trace = {{{1, 2}, {2, 4}}};
    std::vector<double> assign_conf = {0.9, 0.9, 0.9, 0.9};

    // Only "area" exists (with a NaN entry); density/elongation are missing
    // -> column_index returns -1 and extract_col bails out.
    Eigen::MatrixXd stats(2, 1);
    stats << std::numeric_limits<double>::quiet_NaN(),
             1.0;
    std::vector<std::string> cols = {"area"};
    baysor::PriorInputOptions prior;

    auto html = diagnostic_html(data, edges, noise, assignment, trace, assign_conf,
                                nullptr, nullptr, stats, cols, prior, "1.2");

    // Single-bin histograms: cell sizes all 1.0 and confidence all 0.9.
    EXPECT_NE(html.find(R"({"count":4,"x":0.5,"x2":1.5})"), std::string::npos);
    EXPECT_NE(html.find(R"({"count":4,"x":0.4,"x2":1.4})"), std::string::npos);
    // NaN is dropped from the area column -> single observation.
    EXPECT_NE(html.find(R"({"count":1,"x":0.5,"x2":1.5})"), std::string::npos);
    EXPECT_NE(html.find("Cell area"), std::string::npos);
    EXPECT_NE(html.find("Cell density"), std::string::npos);
    EXPECT_NE(html.find("Cell elongation"), std::string::npos);
    EXPECT_NE(html.find("Assignment confidence"), std::string::npos);
}

TEST(Cov4RunHtml, DiagnosticHtmlSkipsManifoldWhenSampleSizesDisagree) {
    auto data = wide_data();
    auto noise = noise_fit();
    std::vector<double> edges = {0.1, 1.0};
    std::vector<int> assignment = {1, 1, 0, 2};
    std::vector<std::unordered_map<int, int>> trace = {{{1, 2}, {2, 2}}};
    std::vector<double> assign_conf = {1.0, 1.0, 0.0, 1.0};
    Eigen::MatrixXd stats(2, 3);
    stats << 1.0, 2.0, 3.0,
             4.0, 5.0, 6.0;
    std::vector<std::string> cols = {"area", "density", "elongation"};
    baysor::PriorInputOptions prior;

    baysor::NcvReportEmbedding ncv;
    ncv.colors = {"#ff0000", "#00ff00", "#0000ff", "#ff00ff"};
    ncv.sample_ids = {0, 2, 3};
    ncv.sample_umap_x = {0.0};  // deliberately inconsistent with sample_ids
    ncv.sample_umap_y = {0.0, 1.0, 0.5};

    auto html = diagnostic_html(data, edges, noise, assignment, trace, assign_conf,
                                nullptr, &ncv, stats, cols, prior, "1.2");
    EXPECT_EQ(html.find("NCV / clustering manifold"), std::string::npos);
    // The confidence section is always present.
    EXPECT_NE(html.find("id=\"vg_noise\""), std::string::npos);
}

// ============================================================================
// generate_run_segmentation_html
// ============================================================================

TEST(Cov4RunHtml, SegmentationHtmlMinimalWithoutOptionalLayers) {
    auto data = wide_data();
    std::vector<int> assignment = {1, 1, 0, 2};
    std::vector<int> empty_clusters;

    auto html = baysor::generate_run_segmentation_html(data, assignment,
                                                       /*ncv_color=*/{},
                                                       &empty_clusters,
                                                       /*polygons=*/nullptr);
    EXPECT_NE(html.find("Baysor Segmentation Plot"), std::string::npos);
    EXPECT_NE(html.find("Final cell assignment"), std::string::npos);
    EXPECT_NE(html.find("data:image/png;base64,"), std::string::npos);
    EXPECT_EQ(html.find("Local expression similarity (NCV)"), std::string::npos);
    EXPECT_EQ(html.find("Molecule clustering"), std::string::npos);
}

TEST(Cov4RunHtml, SegmentationHtmlRendersClusterLayerWithZeroIds) {
    auto data = wide_data();
    std::vector<int> assignment = {1, 1, 0, 2};
    std::vector<int> clusters = {0, 1, 0, 2};  // 0 -> black, others -> palette

    auto html = baysor::generate_run_segmentation_html(data, assignment,
                                                       /*ncv_color=*/{},
                                                       &clusters,
                                                       /*polygons=*/nullptr);
    EXPECT_NE(html.find("Molecule clustering"), std::string::npos);
    EXPECT_NE(html.find("#clusters"), std::string::npos);
    EXPECT_EQ(html.find("Local expression similarity (NCV)"), std::string::npos);
}

// ============================================================================
// color_utils: normalize_embedding_to_lab_range
// ============================================================================

TEST(Cov4ColorUtils, NormalizeEmbeddingWithLogColors) {
    Eigen::MatrixXd emb(3, 60);
    std::mt19937 rng(42);
    std::normal_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 60; ++j) emb(i, j) = dist(rng);

    baysor::normalize_embedding_to_lab_range(emb, 10.0, 90.0, 0.0125, /*log_colors=*/true);

    for (int j = 0; j < 60; ++j) {
        EXPECT_GE(emb(0, j), 10.0 - 1e-6);
        EXPECT_LE(emb(0, j), 90.0 + 1e-6);
        EXPECT_GE(emb(1, j), -100.0 - 1e-6);
        EXPECT_LE(emb(1, j), 100.0 + 1e-6);
        EXPECT_GE(emb(2, j), -100.0 - 1e-6);
        EXPECT_LE(emb(2, j), 100.0 + 1e-6);
    }
    EXPECT_TRUE(emb.allFinite());
}

TEST(Cov4ColorUtils, NormalizeEmbeddingDegenerateAndEmptyInputs) {
    // All-zero embedding: max quantile is 0 and per-row log range is 0.
    Eigen::MatrixXd zeros = Eigen::MatrixXd::Zero(3, 10);
    baysor::normalize_embedding_to_lab_range(zeros, 10.0, 90.0, 0.0125, true);
    EXPECT_TRUE(zeros.allFinite());
    for (int j = 0; j < 10; ++j) {
        EXPECT_DOUBLE_EQ(zeros(0, j), 10.0);
        EXPECT_DOUBLE_EQ(zeros(1, j), -100.0);
        EXPECT_DOUBLE_EQ(zeros(2, j), -100.0);
    }

    Eigen::MatrixXd empty(3, 0);
    baysor::normalize_embedding_to_lab_range(empty, 10.0, 90.0, 0.0125, true);
    EXPECT_EQ(empty.cols(), 0);
}

// ============================================================================
// color_utils: NCV colour embedding variants
// ============================================================================

TEST(Cov4ColorUtils, ReportEmbeddingReturnsAnchorDiagnostics) {
    Eigen::MatrixXf mol_vecs = random_vecs(8, 60, /*seed=*/42);
    std::vector<double> confidence(60, 0.99);

    auto res = baysor::gene_composition_report_embedding(
        mol_vecs, confidence, /*sample_size=*/30, /*seed=*/5,
        /*n_pca_dims=*/4, /*graph_k=*/10);

    ASSERT_EQ(res.colors.size(), 60u);
    for (const auto& c : res.colors) {
        ASSERT_EQ(c.size(), 7u);
        EXPECT_EQ(c[0], '#');
    }
    ASSERT_EQ(res.sample_ids.size(), 30u);
    ASSERT_EQ(res.sample_umap_x.size(), res.sample_ids.size());
    ASSERT_EQ(res.sample_umap_y.size(), res.sample_ids.size());
    EXPECT_EQ(res.anchor_count, 60);
    EXPECT_DOUBLE_EQ(res.chosen_threshold, 0.95);
    for (int id : res.sample_ids) {
        EXPECT_GE(id, 0);
        EXPECT_LT(id, 60);
    }
    for (double v : res.sample_umap_x) EXPECT_TRUE(std::isfinite(v));
    for (double v : res.sample_umap_y) EXPECT_TRUE(std::isfinite(v));
}

TEST(Cov4ColorUtils, ColorEmbeddingFallsBackWithTooFewComponents) {
    // n_components < 3 -> UMAP cannot run, everything falls back to grey.
    Eigen::MatrixXf mol_vecs = random_vecs(2, 10, /*seed=*/7);
    std::vector<double> confidence(10, 0.99);

    auto colors = baysor::gene_composition_color_embedding(mol_vecs, confidence,
                                                           /*sample_size=*/5);
    ASSERT_EQ(colors.size(), 10u);
    for (const auto& c : colors) EXPECT_EQ(c, "#808080");
}

TEST(Cov4ColorUtils, ReportEmbeddingStreamsWithFreshBasisFit) {
    constexpr int groups = 3;
    constexpr int per_group = 24;
    constexpr int n = groups * per_group;

    Eigen::MatrixXd pos(2, n);
    std::vector<int> genes(n, 1);
    std::vector<double> confidence(n, 0.99);
    for (int g = 0; g < groups; ++g) {
        for (int i = 0; i < per_group; ++i) {
            const int idx = g * per_group + i;
            pos(0, idx) = 30.0 * g + 0.5 * i;
            pos(1, idx) = 5.0 * (i % 6);
            genes[idx] = 1 + g;
        }
    }

    auto res = baysor::gene_composition_report_embedding_streaming(
        pos, genes, /*n_genes=*/3, confidence,
        /*k_neighbors=*/8,
        /*basis_sample_size=*/48,
        /*sample_size=*/24,
        /*seed=*/7,
        /*n_pca_dims=*/3,
        /*graph_k=*/15);

    ASSERT_EQ(res.colors.size(), static_cast<size_t>(n));
    std::set<std::string> unique(res.colors.begin(), res.colors.end());
    EXPECT_GT(unique.size(), 1u);
    ASSERT_EQ(res.sample_ids.size(), res.sample_umap_x.size());
    ASSERT_EQ(res.sample_ids.size(), res.sample_umap_y.size());
    EXPECT_GT(res.sample_ids.size(), 0u);
    EXPECT_GT(res.anchor_count, 0);
    for (double v : res.sample_umap_x) EXPECT_TRUE(std::isfinite(v));
    for (double v : res.sample_umap_y) EXPECT_TRUE(std::isfinite(v));
}

TEST(Cov4ColorUtils, ReportEmbeddingStreamsFallsBackWhenAnchorsBelowThreshold) {
    constexpr int n = 4;
    Eigen::MatrixXd pos(2, n);
    pos << 0.0, 1.0, 0.0, 1.0,
           0.0, 0.0, 1.0, 1.0;
    std::vector<int> genes = {1, 2, 1, 2};
    std::vector<double> confidence(n, 0.1);  // nothing passes the 0.5 floor

    baysor::NcvProjectedModel model;
    model.basis.spatial_k = 6;
    model.basis.basis_ids = {0, 1};
    model.basis.basis_vecs = Eigen::MatrixXf::Ones(3, 2);
    model.basis.gene_emb_t = Eigen::MatrixXf::Ones(3, 3);
    model.basis.distance_floor = 1.0;
    model.mol_vecs = Eigen::MatrixXf::Ones(3, n);

    auto res = baysor::gene_composition_report_embedding_streaming(
        pos, genes, /*n_genes=*/2, confidence,
        /*k_neighbors=*/6,
        /*basis_sample_size=*/100,
        /*sample_size=*/20,
        /*seed=*/1,
        /*n_pca_dims=*/3,
        /*graph_k=*/10,
        &model);

    ASSERT_EQ(res.colors.size(), static_cast<size_t>(n));
    for (const auto& c : res.colors) EXPECT_EQ(c, "#808080");
    EXPECT_TRUE(res.sample_ids.empty());
    EXPECT_EQ(res.anchor_count, 0);
    EXPECT_DOUBLE_EQ(res.chosen_threshold, 0.5);
}

TEST(Cov4ColorUtils, ReportEmbeddingStreamsRefitsWhenPrecomputedKDiffers) {
    constexpr int groups = 3;
    constexpr int per_group = 24;
    constexpr int n = groups * per_group;

    Eigen::MatrixXd pos(2, n);
    std::vector<int> genes(n, 1);
    std::vector<double> confidence(n, 0.99);
    for (int g = 0; g < groups; ++g) {
        for (int i = 0; i < per_group; ++i) {
            const int idx = g * per_group + i;
            pos(0, idx) = 30.0 * g + 0.5 * i;
            pos(1, idx) = 5.0 * (i % 6);
            genes[idx] = 1 + g;
        }
    }

    auto model = baysor::fit_ncv_projected_model(
        pos, genes, /*n_genes=*/3, confidence,
        /*k_neighbors=*/8,
        /*basis_sample_size=*/48,
        /*n_components=*/20,
        /*include_full_projection=*/true);
    ASSERT_GT(model.basis.basis_ids.size(), 1u);
    ASSERT_EQ(model.basis.spatial_k, 8);

    // Requested k differs from the precomputed model -> it is ignored and a
    // fresh basis is fitted.
    auto res = baysor::gene_composition_report_embedding_streaming(
        pos, genes, /*n_genes=*/3, confidence,
        /*k_neighbors=*/10,
        /*basis_sample_size=*/48,
        /*sample_size=*/24,
        /*seed=*/7,
        /*n_pca_dims=*/3,
        /*graph_k=*/15,
        &model);

    ASSERT_EQ(res.colors.size(), static_cast<size_t>(n));
    std::set<std::string> unique(res.colors.begin(), res.colors.end());
    EXPECT_GT(unique.size(), 1u);
    EXPECT_GT(res.anchor_count, 0);
}

TEST(Cov4ColorUtils, FitNcvBasisModelDownsamplesSpatialGridSelection) {
    // 60 candidates spread over only 3 of 4 grid bins (top-right is empty):
    // the grid pass overshoots basis_sample_size=4 (6 ids), so the explicit
    // downsampling pass runs.
    constexpr int n = 60;
    Eigen::MatrixXd pos(2, n);
    std::vector<int> genes(n);
    std::vector<double> confidence(n, 0.99);
    for (int i = 0; i < n; ++i) {
        const int q = i % 3;
        if (q == 0) {          // top-left: x < 5, y >= 5
            pos(0, i) = (i % 5) * 0.9;
            pos(1, i) = 5.0 + (i % 4) * 1.1;
        } else if (q == 1) {   // bottom-left: x < 5, y < 5
            pos(0, i) = (i % 5) * 0.9;
            pos(1, i) = (i % 4) * 1.1;
        } else {               // bottom-right: x >= 5, y < 5
            pos(0, i) = 5.0 + (i % 5) * 0.9;
            pos(1, i) = (i % 4) * 1.1;
        }
        genes[i] = 1 + (i % 3);
    }

    auto model = baysor::fit_ncv_basis_model(pos, genes, /*n_genes=*/3, confidence,
                                             /*k_neighbors=*/2,
                                             /*basis_sample_size=*/4,
                                             /*n_components=*/6);

    EXPECT_EQ(model.basis_ids.size(), 4u);  // exactly the requested budget
    EXPECT_EQ(model.spatial_k, 2);
    EXPECT_GT(model.distance_floor, 0.0);
    EXPECT_GT(model.gene_emb_t.size(), 0);
    std::set<int> unique(model.basis_ids.begin(), model.basis_ids.end());
    EXPECT_EQ(unique.size(), model.basis_ids.size());  // no duplicates
}

// ============================================================================
// gene_structure.cpp
// ============================================================================

TEST(Cov4GeneStructure, PairwiseCorrelationRespectsConfidenceThreshold) {
    constexpr int n = 8;
    std::vector<int> genes = {1, 2, 1, 2, 1, 2, 1, 2};
    auto adj = chain_adj_list(n);

    std::vector<double> high(n, 0.99);
    auto cor = baysor::pairwise_gene_spatial_cor(genes, high, adj);
    ASSERT_EQ(cor.rows(), 2);
    ASSERT_EQ(cor.cols(), 2);
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 2; ++c) {
            EXPECT_GE(cor(r, c), 0.0);
            EXPECT_LE(cor(r, c), 1.0 + 1e-9);
        }
    }
    // On a strictly alternating chain each gene only neighbours the other one.
    EXPECT_GT(cor(0, 1), 0.0);
    EXPECT_GT(cor(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(cor(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(cor(1, 1), 0.0);

    // Everything below the 0.95 default threshold is skipped.
    std::vector<double> low(n, 0.9);
    auto zero = baysor::pairwise_gene_spatial_cor(genes, low, adj);
    ASSERT_EQ(zero.rows(), 2);
    EXPECT_DOUBLE_EQ(zero.sum(), 0.0);
}

TEST(Cov4GeneStructure, EstimateEmbeddingWithPositiveCorrelations) {
    constexpr int n = 12;
    std::vector<int> genes = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4};
    std::vector<std::string> names = {"A", "B", "C", "D"};
    std::vector<double> confidence(n, 1.0);
    auto adj = chain_adj_list(n);

    auto emb = baysor::estimate_gene_structure_embedding(genes, names, confidence, adj,
                                                         /*seed=*/7);

    EXPECT_EQ(emb.gene_names, names);
    ASSERT_EQ(emb.x.size(), 4u);
    ASSERT_EQ(emb.y.size(), 4u);
    ASSERT_EQ(emb.marker_sizes.size(), 4u);
    const double log3 = std::log(3.0);  // each gene has exactly 3 molecules
    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(std::isfinite(emb.x[i]));
        EXPECT_TRUE(std::isfinite(emb.y[i]));
        EXPECT_NEAR(emb.marker_sizes[i], log3, 1e-12);
    }
}

TEST(Cov4GeneStructure, EstimateEmbeddingWithoutPositiveCorrelations) {
    constexpr int n = 12;
    // Last molecule carries gene id 0 (ignored by the counting loop) and all
    // confidences are below the 0.95 threshold, so no positive correlations
    // exist and min_cor/max_cor keep their defaults.
    std::vector<int> genes = {1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 0};
    std::vector<std::string> names = {"A", "B", "C", "D"};
    std::vector<double> confidence(n, 0.5);
    auto adj = chain_adj_list(n);

    auto emb = baysor::estimate_gene_structure_embedding(genes, names, confidence, adj,
                                                         /*seed=*/7);

    EXPECT_EQ(emb.gene_names, names);
    ASSERT_EQ(emb.x.size(), 4u);
    ASSERT_EQ(emb.y.size(), 4u);
    ASSERT_EQ(emb.marker_sizes.size(), 4u);
    const double log3 = std::log(3.0);
    const double log2 = std::log(2.0);  // gene D only has 2 valid molecules
    EXPECT_NEAR(emb.marker_sizes[0], log3, 1e-12);
    EXPECT_NEAR(emb.marker_sizes[1], log3, 1e-12);
    EXPECT_NEAR(emb.marker_sizes[2], log3, 1e-12);
    EXPECT_NEAR(emb.marker_sizes[3], log2, 1e-12);
    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(std::isfinite(emb.x[i]));
        EXPECT_TRUE(std::isfinite(emb.y[i]));
    }
}
