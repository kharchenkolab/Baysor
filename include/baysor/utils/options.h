#pragma once

#include <string>
#include <stdexcept>
#include <limits>

namespace baysor {

enum class PriorInputType {
    None,
    Column,
    Image,
    Boundary
};

enum class ClusterMethod {
    None,
    Mrf,
    Louvain,
    Leiden
};

struct MoleculeInputOptions {
    std::string x_col = "x";
    std::string y_col = "y";
    std::string z_col = "z";
    std::string gene_col = "gene";
    std::string qv_col = "qv";
    bool force_2d = false;
    int min_molecules_per_gene = 1;
    std::string exclude_genes;       // comma-separated, supports regex (e.g. "Blank*,MALAT1")
    int min_molecules_per_cell = 0;
    int confidence_nn_id = 0;
    double min_qv = -1.0;
    double x_min = -std::numeric_limits<double>::infinity();
    double x_max =  std::numeric_limits<double>::infinity();
    double y_min = -std::numeric_limits<double>::infinity();
    double y_max =  std::numeric_limits<double>::infinity();
    double z_min = -std::numeric_limits<double>::infinity();
    double z_max =  std::numeric_limits<double>::infinity();
};

struct PriorInputOptions {
    PriorInputType type = PriorInputType::None;
    std::string path;
    std::string column_name;
    std::string unassigned_label = "0";
    int min_molecules_per_segment = 0;
    bool estimate_scale_from_prior = true;
};

struct SegmentationOptions {
    double scale = -1.0;
    std::string scale_std = "25%";
    ClusterMethod cluster_method = ClusterMethod::Mrf;
    int n_clusters = 0;
    double cluster_resolution = 1.0;
    int cluster_graph_k = 15;
    int cluster_n_dims = 20;
    int cluster_basis_sample_size = 100000;
    double prior_segmentation_confidence = 0.2;
    int iters = 500;
    double tol = 0.0;     ///< convergence tolerance for main BMM loop (0 = run all iters, matching Julia)
    int n_cells_init = 0;
    std::string nuclei_genes;
    std::string cyto_genes;
};

ClusterMethod parse_cluster_method(const std::string& method);
std::string cluster_method_to_string(ClusterMethod method);
int default_cluster_count(ClusterMethod method);

struct PlottingOptions {
    int gene_composition_neighborhood = 0;
    int min_pixels_per_cell = 15;
    int max_plot_size = 3000;
    int max_z_slices = 10;  ///< z-layers for 3D polygons; larger z-stacks are binned
    std::string ncv_method = "ri";
};

struct RunOptions {
    MoleculeInputOptions molecules;
    PriorInputOptions prior;
    SegmentationOptions segmentation;
    PlottingOptions plotting;
    int threads = 0;  ///< thread-pool size; 0 = OMP_NUM_THREADS or physical cores
};

/// Parse scale_std: "25%" relative to scale, or absolute number like "5.0"
double parse_scale_std(const std::string& scale_std, double scale);

/// Fill in default parameter values derived from min_molecules_per_cell.
/// Throws if required info not available.
int default_param_value(const std::string& param, int min_molecules_per_cell,
                        int n_molecules = -1, int n_genes = -1);

/// Fill in zero-valued derived options from min_molecules_per_cell
void fill_and_check_molecule_input_options(MoleculeInputOptions& opts);
void fill_and_check_prior_input_options(PriorInputOptions& opts, int min_molecules_per_cell);

/// Fill plotting defaults
void fill_and_check_plotting_options(PlottingOptions& opts, int min_molecules_per_cell,
                                     int n_genes = -1);

/// Load RunOptions from a TOML file. Missing keys retain defaults.
RunOptions load_config(const std::string& path);

/// Serialize RunOptions to a TOML dump file (like Julia's segmentation_params.dump.toml).
/// cli_cmd is the full CLI invocation string recorded as a comment.
void save_params_toml(const RunOptions& opts, const std::string& cli_cmd,
                      const std::string& path);

} // namespace baysor
