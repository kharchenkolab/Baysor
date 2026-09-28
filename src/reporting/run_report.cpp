#include "baysor/reporting/run_report.h"

#include "baysor/reporting/preview_report.h"
#include "baysor/utils/general.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <unordered_map>

namespace baysor {

namespace {

std::string html_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        switch (c) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            default: out += c; break;
        }
    }
    return out;
} // GCOVR_EXCL_LINE: gcov: only an exception-unwind block is attributed to this closing brace; normal exit is reported on the preceding line

std::string prior_type_name(PriorInputType t) {
    switch (t) {
        case PriorInputType::None: return "none";
        case PriorInputType::Column: return "column";
        case PriorInputType::Image: return "image";
        case PriorInputType::Boundary: return "boundary";
    }
    return "unknown"; // GCOVR_EXCL_LINE: unreachable defensive return after exhaustive switch over the enum
}

std::vector<int> count_molecules_per_cell(const std::vector<int>& assignment) {
    int max_cell = 0;
    for (int a : assignment) max_cell = std::max(max_cell, a);
    std::vector<int> counts(max_cell, 0);
    for (int a : assignment) {
        if (a > 0) counts[a - 1]++;
    }
    return counts;
}

nlohmann::json vega_generic_histogram(
    const std::vector<double>& values,
    const std::string& title,
    const std::string& x_title,
    int n_bins = 40
) {
    std::vector<double> finite_vals;
    finite_vals.reserve(values.size());
    for (double v : values) {
        if (std::isfinite(v)) finite_vals.push_back(v);
    }

    nlohmann::json vals = nlohmann::json::array();
    if (finite_vals.empty()) {
        return {
            {"$schema", "https://vega.github.io/schema/vega-lite/v5.json"},
            {"title", title},
            {"width", 500}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
            {"height", 250}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
            {"data", {{"values", vals}}}
        };
    }

    double min_v = *std::min_element(finite_vals.begin(), finite_vals.end());
    double max_v = *std::max_element(finite_vals.begin(), finite_vals.end());
    if (max_v <= min_v) {
        vals.push_back({{"x", min_v - 0.5}, {"x2", min_v + 0.5}, {"count", static_cast<int>(finite_vals.size())}});
    } else {
        double width = (max_v - min_v) / n_bins;
        std::vector<int> counts(n_bins, 0);
        for (double v : finite_vals) {
            int bin = std::min(n_bins - 1, static_cast<int>((v - min_v) / width));
            counts[bin]++;
        }
        for (int i = 0; i < n_bins; ++i) {
            double x0 = min_v + i * width;
            double x1 = (i == n_bins - 1) ? max_v : (x0 + width);
            vals.push_back({{"x", x0}, {"x2", x1}, {"count", counts[i]}});
        }
    }

    return {
        {"$schema", "https://vega.github.io/schema/vega-lite/v5.json"},
        {"title", title},
        {"width", 500}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
        {"height", 250}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
        {"data", {{"values", vals}}},
        {"mark", "bar"},
        {"encoding", {
            {"x", {{"field", "x"}, {"type", "quantitative"}, {"title", x_title}}},
            {"x2", {{"field", "x2"}}},
            {"y", {{"field", "count"}, {"type", "quantitative"}, {"title", "Count"}}}
        }}
    };
}

nlohmann::json vega_convergence_trace(
    const std::vector<std::unordered_map<int, int>>& trace
) {
    nlohmann::json vals = nlohmann::json::array();
    if (trace.empty()) {
        return {};
    }

    std::vector<int> thresholds;
    for (const auto& [thr, _] : trace.front()) thresholds.push_back(thr);
    std::sort(thresholds.begin(), thresholds.end());

    for (int iter = 0; iter < static_cast<int>(trace.size()); ++iter) {
        for (int thr : thresholds) {
            auto it = trace[iter].find(thr);
            if (it == trace[iter].end()) continue;
            vals.push_back({
                {"iteration", iter},
                {"threshold", std::to_string(thr)},
                {"n_cells", it->second}
            });
        }
    }

    return {
        {"$schema", "https://vega.github.io/schema/vega-lite/v5.json"},
        {"title", "Segmentation convergence"},
        {"width", 520}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
        {"height", 280}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
        {"data", {{"values", vals}}},
        {"mark", {{"type", "line"}, {"point", false}}},
        {"encoding", {
            {"x", {{"field", "iteration"}, {"type", "quantitative"}, {"title", "Iteration"}}},
            {"y", {{"field", "n_cells"}, {"type", "quantitative"}, {"title", "Num. cells"}}},
            {"color", {{"field", "threshold"}, {"type", "nominal"}, {"title", "Min #molecules"}}}
        }}
    };
}

nlohmann::json vega_clustering_convergence(
    const std::vector<double>& diffs,
    const std::vector<double>& change_fracs,
    const std::string& title
) {
    nlohmann::json vals = nlohmann::json::array();
    int n = std::min(diffs.size(), change_fracs.size());
    for (int i = 0; i < n; ++i) {
        vals.push_back({{"iteration", i + 1}, {"value", 100.0 * diffs[i]}, {"metric", "Max prob. difference"}});
        vals.push_back({{"iteration", i + 1}, {"value", 100.0 * change_fracs[i]}, {"metric", "Molecules changed"}});
    }

    return {
        {"$schema", "https://vega.github.io/schema/vega-lite/v5.json"},
        {"title", title},
        {"width", 420}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
        {"height", 250}, // GCOVR_EXCL_LINE: gcov: only exception-unwind blocks are attributed to this line; the executed code is reported on the neighbouring lines
        {"data", {{"values", vals}}},
        {"mark", {{"type", "line"}, {"point", false}}},
        {"encoding", {
            {"x", {{"field", "iteration"}, {"type", "quantitative"}, {"title", "Iteration"}}},
            {"y", {{"field", "value"}, {"type", "quantitative"}, {"title", "Change, %"}}},
            {"color", {{"field", "metric"}, {"type", "nominal"}, {"title", "Metric"}}}
        }}
    };
}

std::vector<std::string> assignment_colors(const std::vector<int>& assignment) {
    static const std::vector<std::string> palette = {
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
    };

    std::vector<std::string> colors(assignment.size(), "#000000");
    for (int i = 0; i < static_cast<int>(assignment.size()); ++i) {
        int a = assignment[i];
        if (a <= 0) {
            colors[i] = "#000000";
        } else {
            colors[i] = palette[(a - 1) % palette.size()];
        }
    }
    return colors;
} // GCOVR_EXCL_LINE: gcov: only an exception-unwind block is attributed to this closing brace; normal exit is reported on the preceding line

std::vector<std::string> cluster_colors(const std::vector<int>& clusters) {
    static const std::vector<std::string> palette = {
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"
    };

    std::vector<std::string> colors(clusters.size(), "#000000");
    for (int i = 0; i < static_cast<int>(clusters.size()); ++i) {
        int c = clusters[i];
        colors[i] = (c > 0) ? palette[(c - 1) % palette.size()] : "#000000";
    }
    return colors;
} // GCOVR_EXCL_LINE: gcov: only an exception-unwind block is attributed to this closing brace; normal exit is reported on the preceding line

std::vector<std::string> subset_colors(
    const std::vector<std::string>& colors,
    const std::vector<int>& ids
) {
    std::vector<std::string> out;
    out.reserve(ids.size());
    for (int idx : ids) out.push_back(colors[idx]);
    return out;
} // GCOVR_EXCL_LINE: gcov: only an exception-unwind block is attributed to this closing brace; normal exit is reported on the preceding line

std::vector<int> subset_ints(
    const std::vector<int>& vals,
    const std::vector<int>& ids
) {
    std::vector<int> out;
    out.reserve(ids.size());
    for (int idx : ids) out.push_back(vals[idx]);
    return out;
} // GCOVR_EXCL_LINE: gcov: only an exception-unwind block is attributed to this closing brace; normal exit is reported on the preceding line

void normalize_to_unit_square(
    std::vector<double>& x,
    std::vector<double>& y
) {
    if (x.empty() || y.empty()) return;

    auto [xmin_it, xmax_it] = std::minmax_element(x.begin(), x.end());
    auto [ymin_it, ymax_it] = std::minmax_element(y.begin(), y.end());
    double xmin = *xmin_it, xmax = *xmax_it;
    double ymin = *ymin_it, ymax = *ymax_it;
    double xrange = std::max(xmax - xmin, 1e-12);
    double yrange = std::max(ymax - ymin, 1e-12);

    for (double& v : x) v = (v - xmin) / xrange;
    for (double& v : y) v = (v - ymin) / yrange;
}

int column_index(const std::vector<std::string>& names, const std::string& target) {
    for (int i = 0; i < static_cast<int>(names.size()); ++i) {
        if (names[i] == target) return i;
    }
    return -1;
}

std::vector<double> extract_col(const Eigen::MatrixXd& mat, int idx) {
    std::vector<double> out;
    if (idx < 0 || idx >= mat.cols()) return out;
    out.reserve(mat.rows());
    for (int i = 0; i < mat.rows(); ++i) {
        double v = mat(i, idx);
        if (std::isfinite(v)) out.push_back(v);
    }
    return out;
} // GCOVR_EXCL_LINE: gcov: only an exception-unwind block is attributed to this closing brace; normal exit is reported on the preceding line

} // namespace

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
) {
    auto noise_spec = vega_noise_histogram(
        edge_lengths, data.confidence,
        noise_result.signal_mu, noise_result.signal_sigma,
        noise_result.noise_mu, noise_result.noise_sigma,
        confidence_nn_id
    );
    auto conv_spec = vega_convergence_trace(n_components_trace);

    std::vector<int> cell_sizes_i = count_molecules_per_cell(assignment);
    std::vector<double> cell_sizes(cell_sizes_i.begin(), cell_sizes_i.end());
    auto cell_size_spec = vega_generic_histogram(cell_sizes, "Number of molecules per cell", "Num. molecules");
    auto mol_conf_spec = vega_generic_histogram(data.confidence, "Molecule confidence", "Confidence");
    auto assign_conf_spec = vega_generic_histogram(assignment_confidence, "Assignment confidence", "Assignment confidence");

    int area_idx = column_index(cell_stat_col_names, "area");
    int dens_idx = column_index(cell_stat_col_names, "density");
    int elong_idx = column_index(cell_stat_col_names, "elongation");
    auto area_spec = vega_generic_histogram(extract_col(cell_stats, area_idx), "Cell area", "Area");
    auto dens_spec = vega_generic_histogram(extract_col(cell_stats, dens_idx), "Cell density", "Density");
    auto elong_spec = vega_generic_histogram(extract_col(cell_stats, elong_idx), "Cell elongation", "Elongation");

    std::string ncv_umap_png;
    std::string cluster_umap_png;
    if (ncv_report && !ncv_report->sample_ids.empty() &&
        ncv_report->sample_ids.size() == ncv_report->sample_umap_x.size() &&
        ncv_report->sample_ids.size() == ncv_report->sample_umap_y.size()) {
        auto umap_x = ncv_report->sample_umap_x;
        auto umap_y = ncv_report->sample_umap_y;
        normalize_to_unit_square(umap_x, umap_y);
        ncv_umap_png = render_scatter_png(
            umap_x,
            umap_y,
            subset_colors(ncv_report->colors, ncv_report->sample_ids),
            nullptr,
            1540,
            2
        );
        if (clustering_result && clustering_result->assignment.size() == data.n_molecules()) {
            auto sampled_clusters = subset_ints(clustering_result->assignment, ncv_report->sample_ids);
            cluster_umap_png = render_scatter_png(
                umap_x,
                umap_y,
                cluster_colors(sampled_clusters),
                nullptr,
                1540,
                2
            );
        }
    }

    int n_noise = 0;
    for (int a : assignment) n_noise += (a == 0);
    int n_with_prior = 0;
    for (int s : data.prior_segmentation) n_with_prior += (s > 0);
    int max_prior_seg = 0;
    for (int s : data.prior_segmentation) max_prior_seg = std::max(max_prior_seg, s);

    std::ostringstream html;
    html << R"(<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Baysor Run Report</title>
<script src="https://cdn.jsdelivr.net/npm/vega@5"></script>
<script src="https://cdn.jsdelivr.net/npm/vega-lite@5"></script>
<script src="https://cdn.jsdelivr.net/npm/vega-embed@6"></script>
<style>
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 20px; background: #fff; }
h1 { border-bottom: 2px solid #333; padding-bottom: 5px; }
.grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(520px, 1fr)); gap: 24px; align-items: start; }
.plot-container { margin: 0; }
.png-plot { display: block; max-width: 100%; height: auto; border: 1px solid #ddd; }
.stats { background: #f5f5f5; padding: 10px 15px; border-radius: 5px; margin: 10px 0 20px 0; font-size: 14px; line-height: 1.7; }
ul { line-height: 2; }
a { color: #1f77b4; }
</style>
</head>
<body>
<h2>Contents</h2>
<ul>
<li><a href="#summary">Summary</a></li>
<li><a href="#convergence">Convergence</a></li>
)";
    if (!ncv_umap_png.empty()) {
        html << "<li><a href=\"#ncv_manifold\">NCV / clustering manifold</a></li>\n";
    }
    html << R"(
<li><a href="#confidence">Confidence</a></li>
<li><a href="#cells">Cell statistics</a></li>
</ul>

<hr>
<h1 id="summary">Summary</h1>
<div class="stats">
)";
    html << "Molecules: " << data.n_molecules() << "<br>\n";
    html << "Genes: " << data.n_genes() << "<br>\n";
    html << "Final cells: " << *std::max_element(assignment.begin(), assignment.end()) << "<br>\n";
    html << "Noise molecules: " << n_noise << " (" << std::fixed << std::setprecision(1)
         << (100.0 * n_noise / std::max<int>(1, data.n_molecules())) << "%)<br>\n";
    html << "Prior type: " << html_escape(prior_type_name(prior_opts.type)) << "<br>\n";
    if (prior_opts.type == PriorInputType::Column) {
        html << "Prior column: " << html_escape(prior_opts.column_name) << "<br>\n";
    } else if (prior_opts.type != PriorInputType::None) {
        html << "Prior source: " << html_escape(prior_opts.path) << "<br>\n";
    }
    if (!data.prior_segmentation.empty()) {
        html << "Molecules with prior label: " << n_with_prior << " / " << data.n_molecules()
             << " (" << std::fixed << std::setprecision(1)
             << (100.0 * n_with_prior / std::max<int>(1, data.n_molecules())) << "%)<br>\n";
        html << "Prior segments represented: " << max_prior_seg << "<br>\n";
    }
    html << "Scale: " << std::fixed << std::setprecision(2) << scale
         << " (scale_std=" << html_escape(scale_std) << ")\n";
    html << R"(
</div>

<hr>
<h1 id="convergence">Convergence</h1>
<div class="grid">
<div id="vg_seg_conv"></div>
)";
    if (clustering_result && !clustering_result->diffs.empty() &&
        clustering_result->diffs.size() == clustering_result->change_fracs.size()) {
        html << R"(<div id="vg_clust_conv"></div>
)";
    }
    html << R"(</div>
)";
    if (!ncv_umap_png.empty()) {
        html << R"(

<hr>
<h1 id="ncv_manifold">NCV / clustering manifold</h1>
<div class="grid">
<div class="plot-container">
<h3>Colored by NCV</h3>
<img class="png-plot" src=")";
        html << ncv_umap_png;
        html << R"(" alt="NCV manifold colored by NCV">
</div>
)";
        if (!cluster_umap_png.empty()) {
            html << R"(
<div class="plot-container">
<h3>Colored by cluster assignment</h3>
<img class="png-plot" src=")";
            html << cluster_umap_png;
            html << R"(" alt="NCV manifold colored by cluster assignment">
</div>
)";
        }
        html << R"(</div>
)";
    }
    html << R"(

<hr>
<h1 id="confidence">Confidence</h1>
<div class="grid">
<div id="vg_noise"></div>
<div id="vg_mol_conf"></div>
)";
    if (!assignment_confidence.empty()) {
        html << R"(<div id="vg_assign_conf"></div>
)";
    }
    html << R"(</div>

<hr>
<h1 id="cells">Cell statistics</h1>
<div class="grid">
<div id="vg_cell_size"></div>
<div id="vg_area"></div>
<div id="vg_density"></div>
<div id="vg_elongation"></div>
</div>

<script>
)";
    html << "vegaEmbed('#vg_seg_conv', " << conv_spec.dump() << ", {mode: 'vega-lite'});\n";
    if (clustering_result && !clustering_result->diffs.empty() &&
        clustering_result->diffs.size() == clustering_result->change_fracs.size()) {
        auto clust_spec = vega_clustering_convergence(
            clustering_result->diffs, clustering_result->change_fracs, "Molecule clustering convergence");
        html << "vegaEmbed('#vg_clust_conv', " << clust_spec.dump() << ", {mode: 'vega-lite'});\n";
    }
    html << "vegaEmbed('#vg_noise', " << noise_spec.dump() << ", {mode: 'vega-lite'});\n";
    html << "vegaEmbed('#vg_mol_conf', " << mol_conf_spec.dump() << ", {mode: 'vega-lite'});\n";
    if (!assignment_confidence.empty()) {
        html << "vegaEmbed('#vg_assign_conf', " << assign_conf_spec.dump() << ", {mode: 'vega-lite'});\n";
    }
    html << "vegaEmbed('#vg_cell_size', " << cell_size_spec.dump() << ", {mode: 'vega-lite'});\n";
    html << "vegaEmbed('#vg_area', " << area_spec.dump() << ", {mode: 'vega-lite'});\n";
    html << "vegaEmbed('#vg_density', " << dens_spec.dump() << ", {mode: 'vega-lite'});\n";
    html << "vegaEmbed('#vg_elongation', " << elong_spec.dump() << ", {mode: 'vega-lite'});\n";
    html << R"(
</script>
</body>
</html>
)";
    return html.str();
}

std::string generate_run_segmentation_html(
    const MoleculeData& data,
    const std::vector<int>& assignment,
    const std::vector<std::string>& ncv_color,
    const std::vector<int>* molecule_clusters,
    const PolygonCollection* polygons
) {
    std::string assign_png = render_scatter_png(data.x, data.y, assignment_colors(assignment), polygons);
    std::string ncv_png;
    if (!ncv_color.empty()) {
        ncv_png = render_scatter_png(data.x, data.y, ncv_color, polygons);
    }

    std::string cluster_png;
    if (molecule_clusters && !molecule_clusters->empty()) {
        cluster_png = render_scatter_png(data.x, data.y, cluster_colors(*molecule_clusters), polygons);
    }

    std::ostringstream html;
    html << R"(<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Baysor Segmentation Plot</title>
<style>
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; margin: 20px; background: #fff; }
h1 { border-bottom: 2px solid #333; padding-bottom: 5px; }
.plot-container { margin: 20px 0 40px 0; }
.png-plot { display: block; max-width: 100%; height: auto; cursor: zoom-in; border: 1px solid #ddd; }
ul { line-height: 2; }
a { color: #1f77b4; }
</style>
<script>
function openFullRes(src) {
    var parts = src.split(',');
    var mime = parts[0].match(/:(.*?);/)[1];
    var raw = atob(parts[1]);
    var arr = new Uint8Array(raw.length);
    for (var i = 0; i < raw.length; i++) arr[i] = raw.charCodeAt(i);
    var blob = new Blob([arr], {type: mime});
    window.open(URL.createObjectURL(blob), '_blank');
}
</script>
</head>
<body>
<h2>Contents</h2>
<ul>
<li><a href="#assignment">Final cell assignment</a></li>
)";
    if (!ncv_png.empty()) {
        html << "<li><a href=\"#ncv\">Local expression similarity (NCV)</a></li>\n";
    }
    if (!cluster_png.empty()) {
        html << "<li><a href=\"#clusters\">Molecule clustering</a></li>\n";
    }
    html << R"(</ul>
<p>Click image to open full resolution in a new tab.</p>

<hr>
<h1 id="assignment">Final cell assignment</h1>
<div class="plot-container">
<img class="png-plot" src=")";
    html << assign_png;
    html << R"(" alt="Final cell assignment">
</div>
)";

    if (!ncv_png.empty()) {
        html << R"(
<hr>
<h1 id="ncv">Local expression similarity (NCV)</h1>
<div class="plot-container">
<img class="png-plot" src=")";
        html << ncv_png;
        html << "\" alt=\"Local expression similarity (NCV)\">\n"
                "</div>\n";
    }

    if (!cluster_png.empty()) {
        html << R"(
<hr>
<h1 id="clusters">Molecule clustering</h1>
<div class="plot-container">
<img class="png-plot" src=")";
        html << cluster_png;
        html << R"(" alt="Molecule clustering">
</div>
)";
    }

    html << R"(
<script>
document.querySelectorAll('.png-plot').forEach(function(img) {
    img.addEventListener('click', function() { openFullRes(this.src); });
});
</script>
</body>
</html>
)";
    return html.str();
}

} // namespace baysor
