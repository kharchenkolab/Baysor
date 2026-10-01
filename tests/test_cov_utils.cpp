// Coverage tests for src/utils/{general,options,xenium}.cpp and
// include/baysor/utils/{general,julia_int_dict,xoshiro}.h (task COV-1).
//
// File-local helpers live in an anonymous namespace; every suite name is
// prefixed with Cov1 so it cannot clash with other coverage test files.

#include <gtest/gtest.h>

#include "baysor/utils/general.h"
#include "baysor/utils/julia_int_dict.h"
#include "baysor/utils/options.h"
#include "baysor/utils/xenium.h"
#include "baysor/utils/xoshiro.h"
#include "baysor/processing/bmm_algorithm/bmm_algorithm.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/processing/models/bmm_data.h"
#include "baysor/processing/models/component.h"
#include "baysor/processing/distributions/mv_normal.h"
#include "baysor/processing/distributions/categorical_smoothed.h"

#include <atomic>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "test_cov_helpers.h"

namespace {

// Portable RAII temp directory (see tests/test_cov_helpers.h): unique via a
// counter plus a random suffix, no getpid()/POSIX.
using TempDir = baysor_test::TempDir;

std::string write_file(const TempDir& dir, const std::string& name,
                       const std::string& content) {
    const std::string p = dir.file(name);
    std::ofstream f(p);
    f << content;
    return p;
}

} // namespace

// ============================================================================
// src/utils/general.cpp
// ============================================================================

TEST(Cov1Utils_General, CountArrayWeightedBranches) {
    // Empty input -> empty result.
    EXPECT_TRUE(baysor::count_array_weighted({}, {}).empty());

    // All-zero values with inferred max_value <= 0 -> empty result.
    EXPECT_TRUE(baysor::count_array_weighted({0, 0, 0}, {1.0, 1.0, 1.0}).empty());

    // Inferred max_value = max(values) = 5; zeros skipped; both explicit and
    // inferred paths are exercised.
    std::vector<int> values{1, 4, 0, 4, 5};
    std::vector<double> weights{1.0, 2.0, 100.0, 3.0, 0.5};
    auto counts = baysor::count_array_weighted(values, weights);
    ASSERT_EQ(counts.size(), 5u);
    EXPECT_DOUBLE_EQ(counts[0], 1.0);   // value 1
    EXPECT_DOUBLE_EQ(counts[1], 0.0);   // value 2 absent
    EXPECT_DOUBLE_EQ(counts[2], 0.0);   // value 3 absent
    EXPECT_DOUBLE_EQ(counts[3], 5.0);   // value 4: 2.0 + 3.0
    EXPECT_DOUBLE_EQ(counts[4], 0.5);   // value 5

    // Explicit max smaller than the data truncates (and ignores out-of-range).
    auto truncated = baysor::count_array_weighted(values, weights, 2);
    ASSERT_EQ(truncated.size(), 2u);
    EXPECT_DOUBLE_EQ(truncated[0], 1.0);
    EXPECT_DOUBLE_EQ(truncated[1], 0.0);

    // All-negative values with inferred max <= 0 -> empty result.
    EXPECT_TRUE(baysor::count_array_weighted({-5, -1}, {1.0, 1.0}).empty());
}

TEST(Cov1Utils_General, FsampleMt19937) {
    double w0[] = {0.0, 1.0};
    std::mt19937 rng(42);

    // n == 0 -> -1
    EXPECT_EQ(baysor::fsample(w0, 0, rng), -1);

    // weights[0] == 0 -> the cumulative sum starts at 0, so t > 0 forces the
    // loop body to run and index 1 to be selected.
    EXPECT_EQ(baysor::fsample(w0, 2, rng), 1);

    // Deterministic draws stay in range for generic weights.
    double w1[] = {10.0, 1.0, 1.0};
    for (int i = 0; i < 20; ++i) {
        int idx = baysor::fsample(w1, 3, rng);
        ASSERT_GE(idx, 0);
        ASSERT_LT(idx, 3);
    }

    // arr overloads from the header (both rng flavours).
    const int arr[] = {7, 8, 9};
    int picked = baysor::fsample(arr, w1, 3, rng);
    EXPECT_GE(picked, 7);
    EXPECT_LE(picked, 9);

    baysor::Xoshiro256pp xrng(1);
    picked = baysor::fsample(arr, w1, 3, xrng);
    EXPECT_GE(picked, 7);
    EXPECT_LE(picked, 9);
}

TEST(Cov1Utils_General, FsampleXoshiro) {
    double w0[] = {0.0, 1.0};
    baysor::Xoshiro256pp rng(1);

    EXPECT_EQ(baysor::fsample(w0, 0, rng), -1);

    // total = 1 and weights[0] = 0, so the while loop body is entered and the
    // second bin (the only one with weight) is returned.
    EXPECT_EQ(baysor::fsample(w0, 2, rng), 1);
}

TEST(Cov1Utils_General, GlobalXoshiroResetIsDeterministic) {
    // Reseeds the global RNG: restore the default stream afterwards so later
    // tests cannot observe this state (order independence).
    baysor_test::GlobalRngGuard rng_guard;
    baysor::reset_global_xoshiro_rng(123);
    double first = baysor::global_xoshiro_rng().rand_float64();
    baysor::reset_global_xoshiro_rng(123);
    double second = baysor::global_xoshiro_rng().rand_float64();
    EXPECT_EQ(first, second);
    EXPECT_GE(first, 0.0);
    EXPECT_LT(first, 1.0);
}

TEST(Cov1Utils_General, WeightedMean) {
    double v[] = {1.0, 2.0, 3.0};
    double w[] = {1.0, 1.0, 2.0};
    // (1 + 2 + 6) / 4 = 2.25
    EXPECT_DOUBLE_EQ(baysor::wmean(v, w, 3), 2.25);

    // Zero total weight -> 0.0; zero count -> 0.0.
    double zw[] = {0.0, 0.0, 0.0};
    EXPECT_DOUBLE_EQ(baysor::wmean(v, zw, 3), 0.0);
    EXPECT_DOUBLE_EQ(baysor::wmean(v, w, 0), 0.0);
}

TEST(Cov1Utils_General, EstimateDifferenceL0Branches) {
    // Two 2x2 column-major matrices.
    double m1[] = {1.0, 2.0, 3.0, 4.0};
    double m2[] = {1.0, 2.5, 3.0, 9.0};

    // No weights: max |diff| = 5.0, both columns change above 1e-7.
    auto r1 = baysor::estimate_difference_l0(m1, m2, 2, 2);
    EXPECT_DOUBLE_EQ(r1.max_diff, 5.0);
    EXPECT_DOUBLE_EQ(r1.change_frac, 1.0);

    // Weight the first column down below the threshold: only 1 of 2 changes.
    double cw[] = {0.1, 1.0};
    auto r2 = baysor::estimate_difference_l0(m1, m2, 2, 2, cw, 1.0);
    EXPECT_DOUBLE_EQ(r2.max_diff, 5.0);
    EXPECT_DOUBLE_EQ(r2.change_frac, 0.5);

    // Zero columns -> zero change fraction.
    auto r3 = baysor::estimate_difference_l0(m1, m2, 2, 0);
    EXPECT_DOUBLE_EQ(r3.change_frac, 0.0);
}

// ============================================================================
// src/utils/options.cpp
// ============================================================================

TEST(Cov1Utils_Options, ParseClusterMethodAllValuesAndError) {
    EXPECT_EQ(baysor::parse_cluster_method("none"), baysor::ClusterMethod::None);
    EXPECT_EQ(baysor::parse_cluster_method(" MRf "), baysor::ClusterMethod::Mrf);
    EXPECT_EQ(baysor::parse_cluster_method("ica-mrf"), baysor::ClusterMethod::Mrf);
    EXPECT_EQ(baysor::parse_cluster_method("ICA"), baysor::ClusterMethod::Mrf);
    EXPECT_EQ(baysor::parse_cluster_method("louvain"), baysor::ClusterMethod::Louvain);
    EXPECT_EQ(baysor::parse_cluster_method("LOUVAIN"), baysor::ClusterMethod::Louvain);
    EXPECT_EQ(baysor::parse_cluster_method("leiden"), baysor::ClusterMethod::Leiden);

    EXPECT_THROW(baysor::parse_cluster_method("kmeans"), std::runtime_error);
    try {
        baysor::parse_cluster_method("kmeans");
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("cluster_method must be one of"),
                  std::string::npos);
    }
}

TEST(Cov1Utils_Options, ClusterMethodToStringAllValues) {
    EXPECT_EQ(baysor::cluster_method_to_string(baysor::ClusterMethod::None), "none");
    EXPECT_EQ(baysor::cluster_method_to_string(baysor::ClusterMethod::Mrf), "mrf");
    EXPECT_EQ(baysor::cluster_method_to_string(baysor::ClusterMethod::Louvain), "louvain");
    EXPECT_EQ(baysor::cluster_method_to_string(baysor::ClusterMethod::Leiden), "leiden");
}

TEST(Cov1Utils_Options, DefaultClusterCountAllMethods) {
    EXPECT_EQ(baysor::default_cluster_count(baysor::ClusterMethod::None), 0);
    EXPECT_EQ(baysor::default_cluster_count(baysor::ClusterMethod::Mrf), 4);
    EXPECT_EQ(baysor::default_cluster_count(baysor::ClusterMethod::Louvain), 10);
    EXPECT_EQ(baysor::default_cluster_count(baysor::ClusterMethod::Leiden), 10);
}

TEST(Cov1Utils_Options, DefaultParamValueErrors) {
    // min_molecules_per_cell <= 0 is a hard error.
    try {
        baysor::default_param_value("confidence_nn_id", 0);
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("min_molecules_per_cell"), std::string::npos);
    }

    // n_gene_pcs requires n_genes.
    EXPECT_THROW(baysor::default_param_value("n_gene_pcs", 10, -1, 0), std::runtime_error);
    // n_cells_init requires n_molecules.
    EXPECT_THROW(baysor::default_param_value("n_cells_init", 10, 0, -1), std::runtime_error);
    // Unknown parameter name.
    try {
        baysor::default_param_value("no_such_param", 10, 100, 10);
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("Unknown parameter: no_such_param"),
                  std::string::npos);
    }
}

TEST(Cov1Utils_Options, DefaultParamValueGeneDrivenBranches) {
    // composition_neighborhood without genes falls back to min_molecules_per_cell.
    EXPECT_EQ(baysor::default_param_value("composition_neighborhood", 10, -1, 0), 10);
    // ... and with genes takes max(n_genes/10, min_molecules_per_cell, 3).
    EXPECT_EQ(baysor::default_param_value("composition_neighborhood", 10, -1, 250), 25);
    EXPECT_EQ(baysor::default_param_value("composition_neighborhood", 30, -1, 50), 30);

    // n_gene_pcs value path: clamp(n_genes/3, 30..100, n_genes).
    EXPECT_EQ(baysor::default_param_value("n_gene_pcs", 10, -1, 90), 30);
    EXPECT_EQ(baysor::default_param_value("n_gene_pcs", 10, -1, 400), 100);
    EXPECT_EQ(baysor::default_param_value("n_gene_pcs", 10, -1, 12), 12);

    // n_cells_init value path.
    EXPECT_EQ(baysor::default_param_value("n_cells_init", 10, 100, -1), 20);
}

TEST(Cov1Utils_Options, FillAndCheckMoleculeInputRejectsNonPositive) {
    baysor::MoleculeInputOptions opts;
    opts.min_molecules_per_cell = 0;
    try {
        baysor::fill_and_check_molecule_input_options(opts);
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("'min_molecules_per_cell' must be positive"),
                  std::string::npos);
    }
}

TEST(Cov1Utils_Options, FillAndCheckPriorInputValidation) {
    // Column type without a column name.
    baysor::PriorInputOptions col;
    col.type = baysor::PriorInputType::Column;
    try {
        baysor::fill_and_check_prior_input_options(col, 10);
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("non-empty column_name"), std::string::npos);
    }

    // Image type without a path.
    baysor::PriorInputOptions img;
    img.type = baysor::PriorInputType::Image;
    EXPECT_THROW(baysor::fill_and_check_prior_input_options(img, 10), std::runtime_error);

    // Boundary type without a path.
    baysor::PriorInputOptions bnd;
    bnd.type = baysor::PriorInputType::Boundary;
    try {
        baysor::fill_and_check_prior_input_options(bnd, 10);
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("non-empty path"), std::string::npos);
    }

    // Valid column fills the derived min_molecules_per_segment
    // (default_param_value("min_molecules_per_segment", 10) == max(10/4, 2) == 2).
    baysor::PriorInputOptions ok;
    ok.type = baysor::PriorInputType::Column;
    ok.column_name = "cell_id";
    baysor::fill_and_check_prior_input_options(ok, 10);
    EXPECT_EQ(ok.min_molecules_per_segment, 2);

    // None never fills min_molecules_per_segment.
    baysor::PriorInputOptions none;
    baysor::fill_and_check_prior_input_options(none, 10);
    EXPECT_EQ(none.min_molecules_per_segment, 0);
}

TEST(Cov1Utils_Options, FillAndCheckPlottingOptions) {
    baysor::PlottingOptions opts;
    opts.ncv_method = "dense";
    baysor::fill_and_check_plotting_options(opts, 10, 250);
    EXPECT_EQ(opts.gene_composition_neighborhood, 25);

    // Without genes the neighbourhood falls back to min_molecules_per_cell.
    baysor::PlottingOptions opts2;
    opts2.ncv_method = "sparse";
    baysor::fill_and_check_plotting_options(opts2, 12, -1);
    EXPECT_EQ(opts2.gene_composition_neighborhood, 12);

    // Invalid ncv_method.
    baysor::PlottingOptions bad;
    bad.ncv_method = "nope";
    try {
        baysor::fill_and_check_plotting_options(bad, 10, 10);
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("ncv_method"), std::string::npos);
    }

    // A non-positive min_molecules_per_cell makes the default lookup throw.
    baysor::PlottingOptions zero;
    EXPECT_THROW(baysor::fill_and_check_plotting_options(zero, 0, 10), std::runtime_error);

    // max_z_slices must be at least 1.
    baysor::PlottingOptions bad_z;
    bad_z.max_z_slices = 0;
    try {
        baysor::fill_and_check_plotting_options(bad_z, 10, 10);
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("max_z_slices"), std::string::npos);
    }
}

// ----------------------------------------------------------------------------
// load_config / save_params_toml
// ----------------------------------------------------------------------------

TEST(Cov1Utils_Options, LoadConfigMissingFileThrows) {
    try {
        baysor::load_config("/nonexistent/dir/config.toml");
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("Cannot open config file"),
                  std::string::npos);
    }
    // Empty path returns defaults without touching the filesystem.
    auto defaults = baysor::load_config("");
    EXPECT_EQ(defaults.molecules.x_col, "x");
}

TEST(Cov1Utils_Options, LoadConfigDataSectionAndTypedValues) {
    TempDir dir("cov1_utils");
    const auto path = write_file(dir, "cfg.toml",
        "# comment line\n"
        "\n"
        "line_without_equals\n"
        "[DATA]\n"
        "x = 'pos_x'\n"
        "y = \"pos_y\"\n"
        "gene = \"feature\"\n"
        "min_molecules_per_gene = 2\n"
        "min_molecules_per_cell = 5\n"
        "min_qv = 10.5\n"
        "force_2d = 1\n"
        "\n"
        "[segmentation]\n"
        "scale = 25.0\n"
        "scale_std = 10%\n"
        "cluster_method = \"louvain\"\n"
        "estimate_scale_from_centers = true\n"
        "unassigned_prior_label = \"EMPTY\"\n"
        "\n"
        "[prior]\n"
        "type = \"column\"\n"
        "column_name = \"cell\"\n"
        "estimate_scale_from_prior = false\n"
        "\n"
        "[plotting]\n"
        "gene_composition_neigborhood = 15\n"
        "gene_composition_neighborhood = 25\n"
        "min_pixels_per_cell = 20\n"
        "max_plot_size = 500\n"
        "max_z_slices = 20\n"
        "ncv_method = \"dense\"\n");

    auto opts = baysor::load_config(path);

    // [data] section (header is lowercased) applies molecule keys.
    EXPECT_EQ(opts.molecules.x_col, "pos_x");
    EXPECT_EQ(opts.molecules.y_col, "pos_y");
    EXPECT_EQ(opts.molecules.gene_col, "feature");
    EXPECT_EQ(opts.molecules.min_molecules_per_gene, 2);
    // Well-formed int/float values parse (unparsable ones now raise an error,
    // see Bug3_ConfigErrors in tests/test_bugfix_correctness.cpp).
    EXPECT_EQ(opts.molecules.min_molecules_per_cell, 5);
    EXPECT_DOUBLE_EQ(opts.molecules.min_qv, 10.5);
    // force_2d = "1" parses as true.
    EXPECT_TRUE(opts.molecules.force_2d);

    EXPECT_DOUBLE_EQ(opts.segmentation.scale, 25.0);
    EXPECT_EQ(opts.segmentation.scale_std, "10%");
    EXPECT_EQ(opts.segmentation.cluster_method, baysor::ClusterMethod::Louvain);
    // n_clusters unset -> method-specific default.
    EXPECT_EQ(opts.segmentation.n_clusters, 10);
    // Backward-compatible keys from [segmentation]: the value maps onto
    // estimate_scale_from_prior, which the [prior] section below then
    // overwrites to false.
    EXPECT_EQ(opts.prior.unassigned_label, "EMPTY");
    // [prior] section wins afterwards.
    EXPECT_EQ(opts.prior.type, baysor::PriorInputType::Column);
    EXPECT_EQ(opts.prior.column_name, "cell");
    EXPECT_FALSE(opts.prior.estimate_scale_from_prior);

    // [plotting]: both spellings parsed, correct one wins.
    EXPECT_EQ(opts.plotting.gene_composition_neighborhood, 25);
    EXPECT_EQ(opts.plotting.min_pixels_per_cell, 20);
    EXPECT_EQ(opts.plotting.max_plot_size, 500);
    EXPECT_EQ(opts.plotting.max_z_slices, 20);
    EXPECT_EQ(opts.plotting.ncv_method, "dense");
}

TEST(Cov1Utils_Options, LoadConfigPriorTypes) {
    TempDir dir("cov1_utils");

    const auto img_path = write_file(dir, "img.toml",
        "[prior]\n"
        "type = \"IMAGE\"\n"
        "path = \"mask.tiff\"\n");
    auto img = baysor::load_config(img_path);
    EXPECT_EQ(img.prior.type, baysor::PriorInputType::Image);
    EXPECT_EQ(img.prior.path, "mask.tiff");

    const auto bnd_path = write_file(dir, "bnd.toml",
        "[prior]\n"
        "type = \"boundary\"\n"
        "path = \"b.csv\"\n"
        "min_molecules_per_segment = 7\n");
    auto bnd = baysor::load_config(bnd_path);
    EXPECT_EQ(bnd.prior.type, baysor::PriorInputType::Boundary);
    EXPECT_EQ(bnd.prior.min_molecules_per_segment, 7);

    const auto none_path = write_file(dir, "none.toml",
        "[prior]\n"
        "type = \"none\"\n");
    auto none = baysor::load_config(none_path);
    EXPECT_EQ(none.prior.type, baysor::PriorInputType::None);
}

TEST(Cov1Utils_Options, LoadConfigBooleanParsing) {
    TempDir dir("cov1_utils");
    // "true"/"false" cover both accepting branches of the TOML bool parser;
    // unrecognised bool text now raises an error instead of silently keeping
    // the default (see Bug3_ConfigErrors in tests/test_bugfix_correctness.cpp).
    const auto path = write_file(dir, "bools.toml",
        "[molecules]\n"
        "force_2d = true\n"
        "[segmentation]\n"
        "estimate_scale_from_centers = false\n");
    auto opts = baysor::load_config(path);
    EXPECT_TRUE(opts.molecules.force_2d); // "true" branch of the bool parser
    // "false" branch: [segmentation]'s backward-compatible
    // estimate_scale_from_centers key maps onto estimate_scale_from_prior.
    EXPECT_FALSE(opts.prior.estimate_scale_from_prior);
}

TEST(Cov1Utils_Options, SaveParamsTomlRoundtripAndError) {
    baysor::RunOptions opts;
    opts.molecules.x_col = "x_location";
    opts.molecules.y_col = "y_location";
    opts.molecules.z_col = "z_location";
    opts.molecules.gene_col = "feature_name";
    opts.molecules.qv_col = "qv";
    opts.molecules.force_2d = true;
    opts.molecules.min_molecules_per_gene = 3;
    opts.molecules.exclude_genes = "Blank*";
    opts.molecules.min_molecules_per_cell = 12;
    opts.molecules.confidence_nn_id = 7;
    opts.molecules.min_qv = 20.0;
    opts.molecules.x_min = 0.0;
    opts.molecules.x_max = 100.0;
    opts.molecules.y_min = -5.0;
    opts.molecules.y_max = 50.0;
    opts.molecules.z_min = -1.0;
    opts.molecules.z_max = 1.0;
    opts.prior.type = baysor::PriorInputType::Column;
    opts.prior.path = "some_path";
    opts.prior.column_name = "cell_id";
    opts.prior.unassigned_label = "UNASSIGNED";
    opts.prior.min_molecules_per_segment = 4;
    opts.prior.estimate_scale_from_prior = false;
    opts.segmentation.scale = 30.0;
    opts.segmentation.scale_std = "50%";
    opts.segmentation.cluster_method = baysor::ClusterMethod::Leiden;
    opts.segmentation.n_clusters = 12;
    opts.segmentation.cluster_resolution = 0.5;
    opts.segmentation.cluster_graph_k = 11;
    opts.segmentation.cluster_n_dims = 7;
    opts.segmentation.cluster_basis_sample_size = 1234;
    opts.segmentation.prior_segmentation_confidence = 0.4;
    opts.segmentation.iters = 250;
    opts.segmentation.n_cells_init = 60;
    opts.segmentation.nuclei_genes = "G1,G2";
    opts.segmentation.cyto_genes = "G3";
    opts.plotting.gene_composition_neighborhood = 33;
    opts.plotting.min_pixels_per_cell = 17;
    opts.plotting.max_plot_size = 444;
    opts.plotting.max_z_slices = 42;
    opts.plotting.ncv_method = "sparse";

    TempDir dir("cov1_utils");
    const auto out = dir.file("params.toml");
    baysor::save_params_toml(opts, "baysor run -d data.csv", out);

    std::ifstream in(out);
    ASSERT_TRUE(in.good());
    std::stringstream ss;
    ss << in.rdbuf();
    const std::string content = ss.str();

    EXPECT_NE(content.find("# CLI params: `baysor run -d data.csv`"), std::string::npos);
    EXPECT_NE(content.find("[molecules]"), std::string::npos);
    EXPECT_NE(content.find("x = \"x_location\""), std::string::npos);
    EXPECT_NE(content.find("force_2d = true"), std::string::npos);
    EXPECT_NE(content.find("exclude_genes = \"Blank*\""), std::string::npos);
    EXPECT_NE(content.find("min_qv = 20"), std::string::npos);
    EXPECT_NE(content.find("type = \"column\""), std::string::npos);
    EXPECT_NE(content.find("unassigned_label = \"UNASSIGNED\""), std::string::npos);
    EXPECT_NE(content.find("estimate_scale_from_prior = false"), std::string::npos);
    EXPECT_NE(content.find("cluster_method = \"leiden\""), std::string::npos);
    EXPECT_NE(content.find("scale_std = \"50%\""), std::string::npos);
    EXPECT_NE(content.find("ncv_method = \"sparse\""), std::string::npos);

    // Roundtrip: everything we wrote comes back through load_config.
    auto back = baysor::load_config(out);
    EXPECT_EQ(back.molecules.x_col, "x_location");
    EXPECT_EQ(back.molecules.gene_col, "feature_name");
    EXPECT_TRUE(back.molecules.force_2d);
    EXPECT_EQ(back.molecules.min_molecules_per_gene, 3);
    EXPECT_EQ(back.molecules.min_molecules_per_cell, 12);
    EXPECT_DOUBLE_EQ(back.molecules.min_qv, 20.0);
    EXPECT_DOUBLE_EQ(back.molecules.x_max, 100.0);
    EXPECT_DOUBLE_EQ(back.molecules.y_min, -5.0);
    EXPECT_EQ(back.prior.type, baysor::PriorInputType::Column);
    EXPECT_EQ(back.prior.column_name, "cell_id");
    EXPECT_EQ(back.prior.unassigned_label, "UNASSIGNED");
    EXPECT_EQ(back.prior.min_molecules_per_segment, 4);
    EXPECT_FALSE(back.prior.estimate_scale_from_prior);
    EXPECT_DOUBLE_EQ(back.segmentation.scale, 30.0);
    EXPECT_EQ(back.segmentation.cluster_method, baysor::ClusterMethod::Leiden);
    EXPECT_EQ(back.segmentation.n_clusters, 12);
    EXPECT_EQ(back.segmentation.iters, 250);
    EXPECT_EQ(back.plotting.max_z_slices, 42);
    EXPECT_EQ(back.plotting.ncv_method, "sparse");

    // The prior-type switch in save_params_toml covers every enum value.
    const struct {
        baysor::PriorInputType type;
        const char* label;
    } kTypes[] = {
        {baysor::PriorInputType::None, "none"},
        {baysor::PriorInputType::Image, "image"},
        {baysor::PriorInputType::Boundary, "boundary"},
    };
    int i = 0;
    for (const auto& t : kTypes) {
        baysor::RunOptions o;
        o.prior.type = t.type;
        const auto p = dir.file("type_" + std::to_string(i++) + ".toml");
        baysor::save_params_toml(o, "cmd", p);
        std::ifstream f(p);
        std::stringstream s2;
        s2 << f.rdbuf();
        EXPECT_NE(s2.str().find(std::string("type = \"") + t.label + "\""),
                  std::string::npos) << t.label;
    }

    // Error path: unwritable destination.
    try {
        baysor::save_params_toml(opts, "cmd", dir.file("no_such_dir/params.toml"));
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("save_params_toml: cannot open"),
                  std::string::npos);
    }
}

// ============================================================================
// src/utils/xenium.cpp
// ============================================================================

TEST(Cov1Utils_Xenium, ManifestResolution) {
    TempDir dir("cov1_utils");

    // Missing manifest -> open error.
    try {
        baysor::load_xenium_manifest_context(dir.file("experiment.xenium"));
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("Could not open Xenium manifest"),
                  std::string::npos);
    }

    // Manifest exists but no transcripts file next to it.
    write_file(dir, "experiment.xenium", "{}\n");
    try {
        baysor::load_xenium_manifest_context(dir.file("experiment.xenium"));
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("Could not locate transcripts.parquet"),
                  std::string::npos);
    }

    // Only transcripts.csv.gz present -> fallback is used.
    write_file(dir, "transcripts.csv.gz", "stub");
    auto ctx = baysor::load_xenium_manifest_context(dir.file("experiment.xenium"));
    EXPECT_EQ(ctx.manifest_path, dir.file("experiment.xenium"));
    EXPECT_EQ(ctx.dataset_dir, dir.path.string());
    EXPECT_EQ(ctx.transcripts_path, dir.file("transcripts.csv.gz"));

    // transcripts.parquet wins when both exist.
    write_file(dir, "transcripts.parquet", "stub");
    auto ctx2 = baysor::load_xenium_manifest_context(dir.file("experiment.xenium"));
    EXPECT_EQ(ctx2.transcripts_path, dir.file("transcripts.parquet"));
}

// ============================================================================
// include/baysor/utils/julia_int_dict.h
// ============================================================================

TEST(Cov1Utils_JuliaIntDict, GrowAndCollisionProbing) {
    baysor::JuliaIntDoubleDict dict;
    EXPECT_EQ(dict.size(), 0);

    // 100 distinct keys must force multiple grow() cycles and, by the
    // pigeonhole principle over the 16-slot initial table, probe past
    // occupied slots (linear probing) along the way.
    for (int i = 0; i < 100; ++i) dict.add(i, 1.0);
    EXPECT_EQ(dict.size(), 100);

    int seen = 0;
    double total = 0.0;
    dict.for_each([&](int key, double value) {
        ++seen;
        total += value;
        EXPECT_GE(key, 0);
        EXPECT_LT(key, 100);
    });
    EXPECT_EQ(seen, 100);
    EXPECT_DOUBLE_EQ(total, 100.0);

    // Adding an existing key accumulates instead of inserting.
    dict.add(50, 2.5);
    EXPECT_EQ(dict.size(), 100);
    double v50 = -1.0;
    dict.for_each([&](int key, double value) {
        if (key == 50) v50 = value;
    });
    EXPECT_DOUBLE_EQ(v50, 3.5);

    // Negative keys work too (julia hash on sign-extended values).
    dict.add(-7, 4.0);
    EXPECT_EQ(dict.size(), 101);

    dict.clear();
    EXPECT_EQ(dict.size(), 0);
    seen = 0;
    dict.for_each([&](int, double) { ++seen; });
    EXPECT_EQ(seen, 0);
}

TEST(Cov1Utils_JuliaIntDict, EStep3DCoversBmmInstance) {
    // The 3-D E-step instantiates JuliaIntDoubleDict::for_each inside
    // expect_dirichlet_spatial<3>; exercise it through the public template.
    using baysor::AdjList;
    using baysor::BmmData;
    using baysor::CategoricalSmoothed;
    using baysor::Component;
    using baysor::MvNormal;
    using baysor::ShapePrior;

    BmmData<3> data;
    data.position_data.resize(3, 6);
    data.position_data <<
        0.0, 0.1, 10.0, 10.1, 20.0, 20.1,
        0.0, 0.0,  0.0,  0.0,  0.0,  0.0,
        0.0, 0.1,  0.0,  0.1,  0.0,  0.1;
    data.composition_data = {0, 0, 1, 1, 0, 1};
    data.confidence.assign(6, 1.0);

    const int edge_src[] = {0, 1, 3, 4};
    const int edge_dst[] = {1, 2, 4, 5};
    const double edge_wt[] = {1.0, 1.0, 1.0, 1.0};
    data.adj_list = AdjList::from_edge_list(edge_src, edge_dst, edge_wt, 4, 6);

    ShapePrior<3> prior;
    prior.std_values << 0.25, 0.25, 0.25;
    prior.std_value_stds << 0.05, 0.05, 0.05;
    prior.n_samples = 3;

    const Eigen::Matrix3d sigma = Eigen::Matrix3d::Identity() * 0.05;
    const Eigen::Vector3d centers[] = {
        (Eigen::Vector3d() << 0.05, 0.0, 0.05).finished(),
        (Eigen::Vector3d() << 10.05, 0.0, 0.05).finished(),
        (Eigen::Vector3d() << 20.05, 0.0, 0.05).finished()
    };
    for (int ci = 0; ci < 3; ++ci) {
        MvNormal<3> pos_params(centers[ci], sigma);
        CategoricalSmoothed comp_params(2, 1.0);
        comp_params.set_uniform_counts(1.0f);
        data.components.emplace_back(pos_params, comp_params, prior, ci + 1);
    }
    data.assignment = {1, 1, 2, 2, 3, 3};
    data.max_component_guid = 3;
    data.noise_position_density = 1e-6;
    data.noise_density = 1e-6;
    data.prior_seg_confidence = 0.2;
    data.cluster_penalty_mult = 0.25;
    data.use_gene_smoothing = true;
    data.min_nuclei_frac = 0.1;
    data.mrf_strength = 0.1;
    data.real_edge_weight = 1.0;

    const auto before = data.assignment;
    auto stats = baysor::expect_dirichlet_spatial<3>(data, /*stochastic=*/false);

    // Candidate components come only from direct graph neighbours. Molecules
    // 3 and 4 (0-based 2, 3) sit at the end of their chains, so their only
    // candidate is the far-away component across the gap; its position
    // density underflows to zero and they fall to noise. Everything else
    // stays with its own spatial group.
    EXPECT_EQ(data.assignment, (std::vector<int>{1, 1, 0, 0, 3, 3}));
    EXPECT_EQ(stats.n_changed, 2);
    // n_changed must equal the number of assignments that actually changed,
    // and every assignment must stay in the valid component-id range.
    std::int64_t expected_changed = 0;
    for (size_t i = 0; i < before.size(); ++i) {
        if (before[i] != data.assignment[i]) ++expected_changed;
    }
    EXPECT_EQ(stats.n_changed, expected_changed);
    for (int a : data.assignment) {
        EXPECT_GE(a, 0);
        EXPECT_LE(a, 3);
    }
}
