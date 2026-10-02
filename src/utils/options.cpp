#include "baysor/utils/options.h"
#include "baysor/utils/general.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <climits>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <functional>

namespace baysor {

ClusterMethod parse_cluster_method(const std::string& method_raw) {
    std::string method = strip(method_raw);
    std::transform(method.begin(), method.end(), method.begin(), ::tolower);
    if (method.empty() || method == "mrf" || method == "ica_mrf" || method == "ica-mrf" || method == "ica") {
        return ClusterMethod::Mrf;
    }
    if (method == "none") {
        return ClusterMethod::None;
    }
    if (method == "louvain") {
        return ClusterMethod::Louvain;
    }
    if (method == "leiden") {
        return ClusterMethod::Leiden;
    }
    throw std::runtime_error(
        "cluster_method must be one of 'mrf', 'louvain', 'leiden', or 'none' "
        "(legacy alias 'ica_mrf' is also accepted)"
    );
}

std::string cluster_method_to_string(ClusterMethod method) {
    switch (method) {
        case ClusterMethod::None: return "none";
        case ClusterMethod::Mrf: return "mrf";
        case ClusterMethod::Louvain: return "louvain";
        case ClusterMethod::Leiden: return "leiden";
    }
    return "mrf"; // GCOVR_EXCL_LINE: unreachable
}

int default_cluster_count(ClusterMethod method) {
    switch (method) {
        case ClusterMethod::Louvain:
        case ClusterMethod::Leiden:
            return 10;
        case ClusterMethod::None:
            return 0;
        case ClusterMethod::Mrf:
        default:
            return 4;
    }
}

double parse_scale_std(const std::string& scale_std_raw, double scale) {
    std::string s = strip(scale_std_raw);
    if (s.empty()) return 0.25 * scale;

    if (s.back() == '%') {
        double pct = std::stod(s.substr(0, s.size() - 1));
        return pct * scale / 100.0;
    }
    return std::stod(s);
}

int default_param_value(const std::string& param, int min_molecules_per_cell,
                        int n_molecules, int n_genes) {
    if (min_molecules_per_cell <= 0) {
        throw std::runtime_error("Either '" + param + "' or 'min_molecules_per_cell' must be provided");
    }
    min_molecules_per_cell = imax(min_molecules_per_cell, 3);

    if (param == "min_molecules_per_segment") {
        // Use integer division (floor) to match Julia's round-half-to-even behaviour:
        // round(10/4) = round(2.5) → 2 in Julia (banker's rounding), 3 in C++ std::round.
        return imax(min_molecules_per_cell / 4, 2);
    }
    if (param == "confidence_nn_id") {
        return imax(min_molecules_per_cell / 2 + 1, 5);
    }
    if (param == "composition_neighborhood") {
        if (n_genes <= 0) return imax(min_molecules_per_cell, 3);
        return imax(imax(n_genes / 10, min_molecules_per_cell), 3);
    }
    if (param == "n_gene_pcs") {
        if (n_genes <= 0) throw std::runtime_error("Either 'n_gene_pcs' or 'n_genes' must be provided");
        return imin(imax(n_genes / 3, 30), imin(100, n_genes));
    }
    if (param == "n_cells_init") {
        if (n_molecules <= 0) throw std::runtime_error("Either 'n_cells_init' or 'n_molecules' must be provided");
        return (n_molecules / min_molecules_per_cell) * 2;
    }
    throw std::runtime_error("Unknown parameter: " + param);
}

void fill_and_check_molecule_input_options(MoleculeInputOptions& opts) {
    if (opts.min_molecules_per_cell <= 0) {
        throw std::runtime_error("'min_molecules_per_cell' must be positive");
    }
    if (opts.confidence_nn_id == 0) {
        opts.confidence_nn_id = default_param_value("confidence_nn_id",
                                                     opts.min_molecules_per_cell);
    }
}

void fill_and_check_prior_input_options(PriorInputOptions& opts, int min_molecules_per_cell) {
    if (opts.type == PriorInputType::Column && opts.column_name.empty()) {
        throw std::runtime_error("Prior column input requires a non-empty column_name");
    }
    if ((opts.type == PriorInputType::Image || opts.type == PriorInputType::Boundary) &&
        opts.path.empty()) {
        throw std::runtime_error("Prior image/boundary input requires a non-empty path");
    }
    if (opts.type != PriorInputType::None && opts.min_molecules_per_segment == 0) {
        opts.min_molecules_per_segment = default_param_value(
            "min_molecules_per_segment", min_molecules_per_cell);
    }
}

void fill_and_check_plotting_options(PlottingOptions& opts, int min_molecules_per_cell,
                                     int n_genes) {
    if (opts.gene_composition_neighborhood <= 0) {
        opts.gene_composition_neighborhood = default_param_value(
            "composition_neighborhood", min_molecules_per_cell, -1, n_genes
        );
    }
    if (opts.ncv_method != "ri" && opts.ncv_method != "dense" && opts.ncv_method != "sparse") {
        throw std::runtime_error("ncv_method must be one of 'ri', 'dense', or 'sparse'");
    }
    if (opts.max_z_slices < 1) {
        throw std::runtime_error("max_z_slices must be at least 1");
    }
}

// ============================================================================
// Minimal TOML parser — handles the flat [section] key = value format Baysor uses.
// No dependency on toml++; the config format is simple enough.
// ============================================================================

namespace {

struct TomlValue {
    std::string raw; // unquoted string
};

using TomlSection = std::unordered_map<std::string, TomlValue>;
using TomlDoc = std::unordered_map<std::string, TomlSection>;

TomlDoc parse_toml_simple(const std::string& path) {
    TomlDoc doc;
    std::ifstream f(path);
    if (!f.is_open()) {
        throw std::runtime_error("Cannot open config file: " + path);
    }

    std::string current_section;
    std::string line;
    while (std::getline(f, line)) {
        // Strip comment
        auto comment_pos = line.find('#');
        if (comment_pos != std::string::npos) line = line.substr(0, comment_pos);
        line = strip(line);
        if (line.empty()) continue;

        // Section header
        if (line.front() == '[' && line.back() == ']') {
            current_section = strip(line.substr(1, line.size() - 2));
            // Normalize to lowercase
            std::transform(current_section.begin(), current_section.end(),
                           current_section.begin(), ::tolower);
            continue;
        }

        // Key = value
        auto eq_pos = line.find('=');
        if (eq_pos == std::string::npos) continue;

        std::string key = strip(line.substr(0, eq_pos));
        std::string val = strip(line.substr(eq_pos + 1));

        // Remove quotes from string values
        if (val.size() >= 2 &&
            ((val.front() == '"' && val.back() == '"') ||
             (val.front() == '\'' && val.back() == '\''))) {
            val = val.substr(1, val.size() - 2);
        }

        doc[current_section][key] = {val};
    }
    return doc;
}

std::string toml_get(const TomlSection& sec, const std::string& key, const std::string& def) {
    auto it = sec.find(key);
    return it != sec.end() ? it->second.raw : def;
}

// Remove TOML digit separators (`1_000` -> `1000`). Other underscores are
// kept, so such values fail to parse as in TOML.jl.
std::string strip_digit_separators(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '_' && i > 0 && i + 1 < s.size() &&
            std::isdigit(static_cast<unsigned char>(s[i - 1])) &&
            std::isdigit(static_cast<unsigned char>(s[i + 1]))) {
            continue;
        }
        out.push_back(s[i]);
    }
    return out;
}

// Like Julia (TOML.jl, Configurations.jl), reject malformed or wrong-typed
// values instead of silently keeping the default.
[[noreturn]] void throw_invalid_value(const std::string& raw, const std::string& key,
                                      const std::string& expected) {
    throw std::runtime_error(
        "Invalid value '" + raw + "' for config key '" + key + "': expected " + expected);
}

// The whole value as a number, or nullopt if it does not parse.
std::optional<double> parse_number(const std::string& raw) {
    const std::string s = strip_digit_separators(raw);
    try {
        size_t pos = 0;
        const double value = std::stod(s, &pos);
        if (pos == s.size()) return value;
    } catch (const std::exception&) {}
    return std::nullopt;
}

int toml_get_int(const TomlSection& sec, const std::string& key, int def) {
    auto it = sec.find(key);
    if (it == sec.end()) return def;
    // Integral floats (`50.0`, `1e2`) are accepted, as Configurations.jl
    // converts them with `Base.convert(Int, ...)`.
    const auto value = parse_number(it->second.raw);
    if (!value || *value != std::trunc(*value) || *value < INT_MIN || *value > INT_MAX) {
        throw_invalid_value(it->second.raw, key, "an integer");
    }
    return static_cast<int>(*value);
}

double toml_get_double(const TomlSection& sec, const std::string& key, double def) {
    auto it = sec.find(key);
    if (it == sec.end()) return def;
    const auto value = parse_number(it->second.raw);
    if (!value) throw_invalid_value(it->second.raw, key, "a number");
    return *value;
}

bool toml_get_bool(const TomlSection& sec, const std::string& key, bool def) {
    auto it = sec.find(key);
    if (it == sec.end()) return def;
    std::string v = it->second.raw;
    std::transform(v.begin(), v.end(), v.begin(), ::tolower);
    if (v == "true" || v == "1") return true;
    if (v == "false" || v == "0") return false;
    throw_invalid_value(it->second.raw, key, "a boolean (true or false)");
}

} // anonymous namespace

RunOptions load_config(const std::string& path) {
    RunOptions opts;

    if (path.empty()) return opts;

    auto doc = parse_toml_simple(path);

    // Top-level keys (before any [section])
    if (doc.count("")) {
        opts.threads = toml_get_int(doc[""], "threads", opts.threads);
    }

    auto apply_molecule_section = [&](const TomlSection& sec) {
        opts.molecules.x_col = toml_get(sec, "x", opts.molecules.x_col);
        opts.molecules.y_col = toml_get(sec, "y", opts.molecules.y_col);
        opts.molecules.z_col = toml_get(sec, "z", opts.molecules.z_col);
        opts.molecules.gene_col = toml_get(sec, "gene", opts.molecules.gene_col);
        opts.molecules.qv_col = toml_get(sec, "qv", opts.molecules.qv_col);
        opts.molecules.force_2d = toml_get_bool(sec, "force_2d", opts.molecules.force_2d);
        opts.molecules.min_molecules_per_gene = toml_get_int(
            sec, "min_molecules_per_gene", opts.molecules.min_molecules_per_gene);
        opts.molecules.exclude_genes = toml_get(sec, "exclude_genes", opts.molecules.exclude_genes);
        opts.molecules.min_molecules_per_cell = toml_get_int(
            sec, "min_molecules_per_cell", opts.molecules.min_molecules_per_cell);
        opts.molecules.confidence_nn_id = toml_get_int(
            sec, "confidence_nn_id", opts.molecules.confidence_nn_id);
        opts.molecules.min_qv = toml_get_double(sec, "min_qv", opts.molecules.min_qv);
        opts.molecules.x_min = toml_get_double(sec, "x_min", opts.molecules.x_min);
        opts.molecules.x_max = toml_get_double(sec, "x_max", opts.molecules.x_max);
        opts.molecules.y_min = toml_get_double(sec, "y_min", opts.molecules.y_min);
        opts.molecules.y_max = toml_get_double(sec, "y_max", opts.molecules.y_max);
        opts.molecules.z_min = toml_get_double(sec, "z_min", opts.molecules.z_min);
        opts.molecules.z_max = toml_get_double(sec, "z_max", opts.molecules.z_max);
        opts.prior.min_molecules_per_segment = toml_get_int(
            sec, "min_molecules_per_segment", opts.prior.min_molecules_per_segment);
    };

    // [molecules] (preferred), [data] (backward compatible)
    if (doc.count("molecules")) {
        apply_molecule_section(doc["molecules"]);
    }
    if (doc.count("data")) {
        apply_molecule_section(doc["data"]);
    }

    // [segmentation]
    if (doc.count("segmentation")) {
        auto& sec = doc["segmentation"];
        opts.segmentation.scale = toml_get_double(sec, "scale", opts.segmentation.scale);
        opts.segmentation.scale_std = toml_get(sec, "scale_std", opts.segmentation.scale_std);
        opts.segmentation.cluster_method =
            parse_cluster_method(toml_get(sec, "cluster_method",
                                         cluster_method_to_string(opts.segmentation.cluster_method)));
        opts.segmentation.n_clusters = toml_get_int(sec, "n_clusters", opts.segmentation.n_clusters);
        opts.segmentation.cluster_resolution = toml_get_double(
            sec, "cluster_resolution", opts.segmentation.cluster_resolution);
        opts.segmentation.cluster_graph_k = toml_get_int(
            sec, "cluster_graph_k", opts.segmentation.cluster_graph_k);
        opts.segmentation.cluster_n_dims = toml_get_int(
            sec, "cluster_n_dims", opts.segmentation.cluster_n_dims);
        opts.segmentation.cluster_basis_sample_size = toml_get_int(
            sec, "cluster_basis_sample_size", opts.segmentation.cluster_basis_sample_size);
        opts.segmentation.prior_segmentation_confidence = toml_get_double(
            sec, "prior_segmentation_confidence", opts.segmentation.prior_segmentation_confidence);
        opts.segmentation.iters = toml_get_int(sec, "iters", opts.segmentation.iters);
        opts.segmentation.tol   = toml_get_double(sec, "tol", opts.segmentation.tol);
        opts.segmentation.n_cells_init = toml_get_int(sec, "n_cells_init", opts.segmentation.n_cells_init);
        opts.segmentation.nuclei_genes = toml_get(sec, "nuclei_genes", opts.segmentation.nuclei_genes);
        opts.segmentation.cyto_genes = toml_get(sec, "cyto_genes", opts.segmentation.cyto_genes);
        // Backward-compatible prior keys formerly stored under [segmentation]
        opts.prior.estimate_scale_from_prior = toml_get_bool(
            sec, "estimate_scale_from_centers", opts.prior.estimate_scale_from_prior);
        opts.prior.unassigned_label = toml_get(
            sec, "unassigned_prior_label", opts.prior.unassigned_label);
    }
    if (opts.segmentation.n_clusters <= 0) {
        opts.segmentation.n_clusters = default_cluster_count(opts.segmentation.cluster_method);
    }

    // [prior]
    if (doc.count("prior")) {
        auto& sec = doc["prior"];
        std::string type_raw = toml_get(sec, "type", "");
        std::string type = type_raw;
        std::transform(type.begin(), type.end(), type.begin(), ::tolower);
        if (type == "none" || type.empty()) opts.prior.type = PriorInputType::None;
        else if (type == "column") opts.prior.type = PriorInputType::Column;
        else if (type == "image") opts.prior.type = PriorInputType::Image;
        else if (type == "boundary") opts.prior.type = PriorInputType::Boundary;
        else throw_invalid_value(type_raw, "type", "one of 'none', 'column', 'image', 'boundary'");
        opts.prior.path = toml_get(sec, "path", opts.prior.path);
        opts.prior.column_name = toml_get(sec, "column_name", opts.prior.column_name);
        opts.prior.unassigned_label = toml_get(
            sec, "unassigned_label", opts.prior.unassigned_label);
        opts.prior.min_molecules_per_segment = toml_get_int(
            sec, "min_molecules_per_segment", opts.prior.min_molecules_per_segment);
        opts.prior.estimate_scale_from_prior = toml_get_bool(
            sec, "estimate_scale_from_prior", opts.prior.estimate_scale_from_prior);
    }

    // [plotting]
    if (doc.count("plotting")) {
        auto& sec = doc["plotting"];
        opts.plotting.gene_composition_neighborhood = toml_get_int(
            sec, "gene_composition_neigborhood",  // note: Julia has the typo "neigborhood"
            opts.plotting.gene_composition_neighborhood);
        // Also accept the correct spelling
        opts.plotting.gene_composition_neighborhood = toml_get_int(
            sec, "gene_composition_neighborhood",
            opts.plotting.gene_composition_neighborhood);
        opts.plotting.min_pixels_per_cell = toml_get_int(sec, "min_pixels_per_cell",
                                                          opts.plotting.min_pixels_per_cell);
        opts.plotting.max_plot_size = toml_get_int(sec, "max_plot_size", opts.plotting.max_plot_size);
        opts.plotting.max_z_slices = toml_get_int(sec, "max_z_slices", opts.plotting.max_z_slices);
        opts.plotting.ncv_method = toml_get(sec, "ncv_method", opts.plotting.ncv_method);
    }

    return opts;
}

// ============================================================================
// save_params_toml — serialise RunOptions to a TOML dump (matches Julia format)
// ============================================================================

void save_params_toml(const RunOptions& opts, const std::string& cli_cmd,
                      const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("save_params_toml: cannot open " + path);

    f << "# CLI params: `" << cli_cmd << "`\n";
    f << "threads = " << opts.threads << "\n";

    f << "[molecules]\n";
    f << "x = \"" << opts.molecules.x_col << "\"\n";
    f << "y = \"" << opts.molecules.y_col << "\"\n";
    f << "z = \"" << opts.molecules.z_col << "\"\n";
    f << "gene = \"" << opts.molecules.gene_col << "\"\n";
    f << "qv = \"" << opts.molecules.qv_col << "\"\n";
    f << "force_2d = " << (opts.molecules.force_2d ? "true" : "false") << "\n";
    f << "min_molecules_per_gene = " << opts.molecules.min_molecules_per_gene << "\n";
    f << "exclude_genes = \"" << opts.molecules.exclude_genes << "\"\n";
    f << "min_molecules_per_cell = " << opts.molecules.min_molecules_per_cell << "\n";
    f << "confidence_nn_id = " << opts.molecules.confidence_nn_id << "\n";
    f << "min_qv = " << opts.molecules.min_qv << "\n";
    f << "x_min = " << opts.molecules.x_min << "\n";
    f << "x_max = " << opts.molecules.x_max << "\n";
    f << "y_min = " << opts.molecules.y_min << "\n";
    f << "y_max = " << opts.molecules.y_max << "\n";
    f << "z_min = " << opts.molecules.z_min << "\n";
    f << "z_max = " << opts.molecules.z_max << "\n";

    f << "\n[prior]\n";
    const char* prior_type = "none";
    switch (opts.prior.type) {
        case PriorInputType::None: prior_type = "none"; break;
        case PriorInputType::Column: prior_type = "column"; break;
        case PriorInputType::Image: prior_type = "image"; break;
        case PriorInputType::Boundary: prior_type = "boundary"; break;
    }
    f << "type = \"" << prior_type << "\"\n";
    f << "path = \"" << opts.prior.path << "\"\n";
    f << "column_name = \"" << opts.prior.column_name << "\"\n";
    f << "unassigned_label = \"" << opts.prior.unassigned_label << "\"\n";
    f << "min_molecules_per_segment = " << opts.prior.min_molecules_per_segment << "\n";
    f << "estimate_scale_from_prior = "
      << (opts.prior.estimate_scale_from_prior ? "true" : "false") << "\n";

    f << "\n[segmentation]\n";
    f << "scale = " << opts.segmentation.scale << "\n";
    f << "scale_std = \"" << opts.segmentation.scale_std << "\"\n";
    f << "cluster_method = \"" << cluster_method_to_string(opts.segmentation.cluster_method) << "\"\n";
    f << "n_clusters = " << opts.segmentation.n_clusters << "\n";
    f << "cluster_resolution = " << opts.segmentation.cluster_resolution << "\n";
    f << "cluster_graph_k = " << opts.segmentation.cluster_graph_k << "\n";
    f << "cluster_n_dims = " << opts.segmentation.cluster_n_dims << "\n";
    f << "cluster_basis_sample_size = " << opts.segmentation.cluster_basis_sample_size << "\n";
    f << "prior_segmentation_confidence = " << opts.segmentation.prior_segmentation_confidence << "\n";
    f << "iters = " << opts.segmentation.iters << "\n";
    f << "tol = " << opts.segmentation.tol << "\n";
    f << "n_cells_init = " << opts.segmentation.n_cells_init << "\n";
    f << "nuclei_genes = \"" << opts.segmentation.nuclei_genes << "\"\n";
    f << "cyto_genes = \"" << opts.segmentation.cyto_genes << "\"\n";

    f << "\n[plotting]\n";
    f << "gene_composition_neigborhood = " << opts.plotting.gene_composition_neighborhood << "\n";
    f << "min_pixels_per_cell = " << opts.plotting.min_pixels_per_cell << "\n";
    f << "max_plot_size = " << opts.plotting.max_plot_size << "\n";
    f << "max_z_slices = " << opts.plotting.max_z_slices << "\n";
    f << "ncv_method = \"" << opts.plotting.ncv_method << "\"\n";
}

} // namespace baysor
