#include <third_party/CLI11.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "baysor/utils/options.h"
#include "baysor/data_loading/data.h"
#include "baysor/data_loading/prior_segmentation.h"
#include "baysor/processing/data_processing/noise_estimation.h"
#include "baysor/processing/data_processing/neighborhood_composition.h"
#include "baysor/processing/data_processing/initialization.h"
#include "baysor/processing/bmm_algorithm/bmm_algorithm.h"
#include "baysor/processing/bmm_algorithm/molecule_clustering.h"
#include "baysor/processing/utils/utils.h"
#include "baysor/processing/utils/convex_hull.h"
#include "baysor/processing/bmm_algorithm/tracing.h"
#include "baysor/reporting/color_utils.h"
#include "baysor/reporting/output.h"
#include "baysor/reporting/preview_report.h"
#include "baysor/reporting/run_report.h"

#include "baysor/utils/general.h"
#include "baysor/utils/thread_pool.h"
#include "baysor/utils/xenium.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <thread>
#include <type_traits>
#include <unordered_map>

using namespace baysor;

// Resolve the worker-thread count for the global thread pool:
// explicit --threads (or config `threads`) > BAYSOR_NUM_THREADS >
// OMP_NUM_THREADS (kept for backward compatibility with scripts and the
// benchmark harness) > std::thread::hardware_concurrency().
static int resolve_thread_count(int requested) {
    if (requested > 0) return requested;
    for (const char* var : {"BAYSOR_NUM_THREADS", "OMP_NUM_THREADS"}) {
        if (const char* env = std::getenv(var)) {
            // OMP_NUM_THREADS may be a comma-separated list; take the first.
            try {
                int n = std::stoi(env);
                if (n > 0) return n;
            } catch (...) {
                // ignore malformed values and fall through
            }
        }
    }
    return default_thread_count();
}

// ============================================================================
// Subcommand: run
// ============================================================================

int cmd_run(
    const std::string& coordinates,
    RunOptions& opts,
    const std::string& output,
    OutputStyle output_style,
    bool plot,
    bool skip_ncv_color,
    const std::string& polygon_format,
    const std::string& count_matrix_format,
    const std::string& cli_cmd
) {
    spdlog::info("Loading data from '{}'...", coordinates);

    fill_and_check_molecule_input_options(opts.molecules);
    fill_and_check_prior_input_options(opts.prior, opts.molecules.min_molecules_per_cell);

    auto data = load_molecules(
        coordinates, opts.molecules, opts.prior);
    spdlog::info("Loaded {} transcripts, {} genes.", data.n_molecules(), data.n_genes());

    fill_and_check_plotting_options(opts.plotting, opts.molecules.min_molecules_per_cell, data.n_genes());
    const bool infer_n_cells_init = opts.segmentation.n_cells_init <= 0;

    // Load prior segmentation if provided
    if (opts.prior.type != PriorInputType::None) {
        auto [scale, scale_std] = load_prior_segmentation(
            data, opts.prior, opts.molecules.min_molecules_per_cell
        );

        if (opts.prior.estimate_scale_from_prior && scale > 0) {
            opts.segmentation.scale = scale;
            opts.segmentation.scale_std = std::to_string(scale_std);
        }
    }

    if (infer_n_cells_init) {
        int inferred_n_cells_init = default_param_value(
            "n_cells_init", opts.molecules.min_molecules_per_cell, data.n_molecules());

        if (!data.prior_segmentation.empty()) {
            int max_label = *std::max_element(data.prior_segmentation.begin(), data.prior_segmentation.end());
            if (max_label > 0) {
                std::vector<int> seg_counts(max_label + 1, 0);
                int n_unassigned = 0;
                for (int lab : data.prior_segmentation) {
                    if (lab > 0) {
                        seg_counts[lab]++;
                    } else {
                        ++n_unassigned;
                    }
                }

                int n_active_prior_segments = 0;
                for (int lab = 1; lab <= max_label; ++lab) {
                    if (seg_counts[lab] > 0) ++n_active_prior_segments;
                }

                if (n_active_prior_segments > 0) {
                    constexpr double prior_segment_multiplier = 2.25;
                    constexpr double unassigned_multiplier = 2.0;
                    int prior_based_n_cells_init = static_cast<int>(std::ceil(
                        prior_segment_multiplier * static_cast<double>(n_active_prior_segments) +
                        unassigned_multiplier * static_cast<double>(n_unassigned) /
                            std::max(opts.molecules.min_molecules_per_cell, 1)
                    ));
                    prior_based_n_cells_init = std::max(prior_based_n_cells_init, n_active_prior_segments);
                    inferred_n_cells_init = std::min(inferred_n_cells_init, prior_based_n_cells_init);

                    spdlog::info( // GCOVR_EXCL_LINE: gcov attributes only this multi-line call's exception-cleanup block to its first line; the call itself is counted on the following lines
                        "Using prior-aware n_cells_init={} (active prior segments={}, unassigned molecules={}, "
                        "default without prior would be {}).",
                        inferred_n_cells_init, n_active_prior_segments, n_unassigned,
                        default_param_value("n_cells_init", opts.molecules.min_molecules_per_cell, data.n_molecules())
                    );
                }
            }
        }

        opts.segmentation.n_cells_init = inferred_n_cells_init;
    }

    if (opts.segmentation.scale <= 0) {
        spdlog::error("Scale could not be determined. Either provide prior_segmentation or set --scale.");
        return 1;
    }

    spdlog::info("Using scale={:.2f}, scale_std={}", // GCOVR_EXCL_LINE: gcov attributes only this multi-line call's exception-cleanup block to its first line; the call itself is counted on the following lines
                 opts.segmentation.scale, opts.segmentation.scale_std);

    double psc = opts.segmentation.prior_segmentation_confidence;

    std::vector<double> noise_edge_lengths;
    NoiseFitResult noise_fit;
    int confidence_nn_id = opts.molecules.confidence_nn_id;
    spdlog::info("Estimating confidence...");
    if (plot) {
        auto conf_details = estimate_confidence_details(data, opts.molecules.confidence_nn_id, psc);
        confidence_nn_id = conf_details.nn_id;
        noise_edge_lengths = std::move(conf_details.edge_lengths);
        noise_fit = std::move(conf_details.fit_result);
        data.confidence.resize(data.n_molecules());
        for (int i = 0; i < data.n_molecules(); ++i) {
            data.confidence[i] = noise_fit.assignment_probs(i, 0);
        }
    } else {
        append_confidence(data, opts.molecules.confidence_nn_id, psc);
    }

    // Build molecule adjacency graph (MRF)
    spdlog::info("Building molecule graph...");
    auto adj_list = build_molecule_graph(data);

    // Create output directory
    {
        std::string mkdir_cmd = "mkdir -p \"" + output + "\"";
        if (std::system(mkdir_cmd.c_str()) != 0) {
            spdlog::warn("Could not create output directory '{}'", output);
        }
    }
    auto out_paths = get_output_paths(output, output_style, count_matrix_format);

    // Set up dual logger: console + log file (matches Julia's setup_logger)
    {
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        auto file_sink    = std::make_shared<spdlog::sinks::basic_file_sink_mt>(
                                out_paths.log_file, /*truncate=*/true);
        auto logger = std::make_shared<spdlog::logger>(
                          "baysor", spdlog::sinks_init_list{console_sink, file_sink});
        logger->set_level(spdlog::level::info);
        logger->flush_on(spdlog::level::info);
        spdlog::set_default_logger(logger);
    }

    // Guard unimplemented features
    if (!opts.segmentation.nuclei_genes.empty() || !opts.segmentation.cyto_genes.empty()) {
        spdlog::error("--nuclei-genes / --cyto-genes compartment segmentation is not yet implemented.");
        return 1;
    }

    int min_mols     = opts.molecules.min_molecules_per_cell;
    int n_iters   = opts.segmentation.iters;
    double scale  = opts.segmentation.scale;
    int n_cells   = opts.segmentation.n_cells_init;
    const std::string& scale_std = opts.segmentation.scale_std;

    // Optional molecule clustering (pre-segmentation cell type assignment).
    // Produces a coarse compatibility prior for segmentation; the default
    // ICA/MRF path matches Julia, while alternative methods plug into the
    // same dispatcher and downstream contract.
    std::vector<int> mol_clusters;
    std::optional<ClusteringResult> clustering_result;
    {
        ClusteringOptions clustering_opts;
        clustering_opts.method = opts.segmentation.cluster_method;
        clustering_opts.n_clusters = opts.segmentation.n_clusters;
        clustering_opts.resolution = opts.segmentation.cluster_resolution;
        clustering_opts.graph_k = opts.segmentation.cluster_graph_k;
        clustering_opts.spatial_k = opts.plotting.gene_composition_neighborhood;
        clustering_opts.n_dims = opts.segmentation.cluster_n_dims;
        clustering_opts.basis_sample_size = opts.segmentation.cluster_basis_sample_size;

        const bool run_clustering =
            clustering_opts.method == ClusterMethod::Louvain ||
            clustering_opts.method == ClusterMethod::Leiden ||
            (clustering_opts.method == ClusterMethod::Mrf && clustering_opts.n_clusters > 1);

        if (run_clustering) {
            if (clustering_opts.method == ClusterMethod::Mrf) {
                spdlog::info("Clustering molecules into {} types (ICA init)...",
                             clustering_opts.n_clusters);
            } else if (clustering_opts.method == ClusterMethod::Louvain) {
                spdlog::info("Clustering molecules with Louvain on NCV kNN graph (k={})...",
                             clustering_opts.graph_k);
            } else {
                spdlog::info("Clustering molecules with Leiden on NCV kNN graph (k={})...",
                             clustering_opts.graph_k);
            }

            clustering_result = cluster_molecules(
                (clustering_opts.method == ClusterMethod::Louvain ||
                 clustering_opts.method == ClusterMethod::Leiden)
                    ? data.position_matrix() : Eigen::MatrixXd(),
                data.gene, adj_list, data.confidence,
                clustering_opts, /*verbose=*/true
            );
            if (clustering_result && !clustering_result->assignment.empty()) {
                mol_clusters = clustering_result->assignment;  // 1-based cluster IDs
                spdlog::info("Molecule clustering complete.");
            }
        }
    }

    // Dispatch on dimensionality
    auto run_segmentation = [&](auto tag) {
        constexpr int N = decltype(tag)::value;

        spdlog::info("Initializing BmmData ({}D)...", N);
        auto bm_data = initialize_bmm_data<N>(
            data, adj_list, n_cells, scale, scale_std, psc, min_mols, /*verbose=*/true);

        // Wire molecule clusters into BmmData
        if (!mol_clusters.empty()) {
            bm_data.cluster_per_molecule = mol_clusters;
        }

        // History depth: match Julia's round(iters * 0.1)
        int history_depth = std::max(1, n_iters / 10);

        spdlog::info("Running segmentation ({} iters, history_depth={}, tol={})...", // GCOVR_EXCL_LINE: gcov attributes only this multi-line call's exception-cleanup block to its first line; the call itself is counted on the following lines
                     n_iters, history_depth, opts.segmentation.tol);
        // Julia hardcodes min_n_samples=2 in drop_unused_components! — match that exactly.
        // min_mols = min_molecules_per_cell = display threshold only.
        bmm(bm_data, /*min_molecules_drop=*/2, n_iters,
            history_depth,
            /*verbose=*/true,
            /*component_split_step=*/3,
            /*refine=*/true,
            /*freeze_composition=*/false,
            /*freeze_position=*/false,
            /*freeze_components=*/false,
            opts.segmentation.tol,
            /*min_molecules_display=*/min_mols);

        int n_cells_final = bm_data.n_components();
        spdlog::info("Segmentation complete: {} cells.", n_cells_final);

        std::vector<std::string> ncv_color;
        std::optional<NcvReportEmbedding> ncv_report;
        if (!skip_ncv_color) {
            spdlog::info("Computing neighborhood composition colors...");
            int ncv_spatial_k = opts.plotting.gene_composition_neighborhood;
            int ncv_graph_k = opts.segmentation.cluster_graph_k;
            auto pos = data.position_matrix();
            const NcvProjectedModel* shared_ncv_model =
                (clustering_result && clustering_result->ncv_projected_model)
                    ? clustering_result->ncv_projected_model.get()
                    : nullptr;
            if (plot) {
                ncv_report = gene_composition_report_embedding_streaming(
                    pos, data.gene, data.n_genes(), data.confidence, ncv_spatial_k,
                    100000, 20000, 42, 10, ncv_graph_k, shared_ncv_model
                );
                ncv_color = ncv_report->colors;
            } else {
                ncv_color = gene_composition_color_embedding_streaming(
                    pos, data.gene, data.n_genes(), data.confidence, ncv_spatial_k,
                    100000, 20000, 42, 10, ncv_graph_k, shared_ncv_model
                );
            }
        }

        // Save per-molecule segmentation table
        spdlog::info("Saving segmented molecule table...");
        const std::vector<double>* ac_ptr = bm_data.assignment_confidence.empty()
                                            ? nullptr : &bm_data.assignment_confidence;
        if (output_style == OutputStyle::Parquet) {
            save_segmented_df_parquet(
                data, bm_data.assignment, data.gene_names, out_paths.segmented_df,
                &ncv_color, ac_ptr, mol_clusters.empty() ? nullptr : &mol_clusters
            );
        } else {
            save_segmented_df(
                data, bm_data.assignment, data.gene_names, out_paths.segmented_df,
                &ncv_color, ac_ptr, mol_clusters.empty() ? nullptr : &mol_clusters
            );
        }

        // Save per-cell stats CSV
        spdlog::info("Saving cell stats...");
        Eigen::MatrixXd cell_stats_mat;
        std::vector<std::string> cell_stat_col_names;
        std::vector<std::string> cell_names(n_cells_final);
        for (int i = 0; i < n_cells_final; ++i) {
            cell_names[i] = "cell_" + std::to_string(i + 1);
        }
        {
            auto ids_by_cell = split_ids(bm_data.assignment, n_cells_final, true);

            // Precompute lifespan map (guid → lifespan)
            std::unordered_map<int,int> lifespan_map;
            if (!bm_data.assignment_history.empty()) {
                lifespan_map = estimate_component_lifespan(bm_data.assignment_history);
            }

            // Column names — order matches Julia:
            // cell, x, y, [z,] [cluster,] n_transcripts, density, elongation, area,
            // avg_confidence, [avg_assignment_confidence,] [max_cluster_frac,] [lifespan]
            bool has_cluster = !bm_data.cluster_per_cell.empty();
            bool has_ac      = !bm_data.assignment_confidence.empty();
            bool has_clmol   = !bm_data.cluster_per_molecule.empty();
            bool has_lspan   = !lifespan_map.empty();

            cell_stat_col_names = {"x", "y"};
            if (N == 3) cell_stat_col_names.push_back("z");
            if (has_cluster) cell_stat_col_names.push_back("cluster");    // matches Julia position
            cell_stat_col_names.push_back("n_transcripts");
            cell_stat_col_names.push_back("density");
            cell_stat_col_names.push_back("elongation");
            cell_stat_col_names.push_back("area");
            cell_stat_col_names.push_back("avg_confidence");
            if (has_ac)    cell_stat_col_names.push_back("avg_assignment_confidence");
            if (has_clmol) cell_stat_col_names.push_back("max_cluster_frac");
            if (has_lspan) cell_stat_col_names.push_back("lifespan");

            // Build name→index map so writing order is independent of names order
            std::unordered_map<std::string, int> ci_map;
            for (int i = 0; i < static_cast<int>(cell_stat_col_names.size()); ++i)
                ci_map[cell_stat_col_names[i]] = i;

            int n_cols = static_cast<int>(cell_stat_col_names.size());
            cell_stats_mat.resize(n_cells_final, n_cols);
            cell_stats_mat.fill(std::numeric_limits<double>::quiet_NaN());

            for (int ci = 0; ci < n_cells_final; ++ci) {
                const auto& ids = ids_by_cell[ci];
                int np = static_cast<int>(ids.size());

                // --- Position means ---
                double sx = 0, sy = 0, sz = 0;
                double sum_conf = 0.0, sum_ac = 0.0;
                for (int mol : ids) {
                    sx += data.x[mol]; sy += data.y[mol];
                    if (N == 3 && !data.z.empty()) sz += data.z[mol];
                    sum_conf += data.confidence[mol];
                    if (has_ac) sum_ac += bm_data.assignment_confidence[mol];
                }
                double denom = np > 0 ? np : 1.0;

                cell_stats_mat(ci, ci_map["x"]) = sx / denom;
                cell_stats_mat(ci, ci_map["y"]) = sy / denom;
                if (N == 3) cell_stats_mat(ci, ci_map["z"]) = sz / denom;

                cell_stats_mat(ci, ci_map["n_transcripts"]) = np;
                cell_stats_mat(ci, ci_map["avg_confidence"]) = sum_conf / denom;
                if (has_ac)    cell_stats_mat(ci, ci_map["avg_assignment_confidence"]) = sum_ac / denom;
                if (has_cluster) cell_stats_mat(ci, ci_map["cluster"]) = bm_data.cluster_per_cell[ci];

                // --- Convex hull metrics (only for cells with > 2 molecules) ---
                if (np > 2) {
                    Eigen::MatrixXd pos2d(2, np);
                    for (int j = 0; j < np; ++j) {
                        pos2d(0, j) = data.x[ids[j]];
                        pos2d(1, j) = data.y[ids[j]];
                    }
                    auto hull = convex_hull(pos2d);
                    double area = polygon_area(hull);
                    cell_stats_mat(ci, ci_map["area"])    = area;
                    cell_stats_mat(ci, ci_map["density"]) = (area > 0) ? (np / area)
                                                               : std::numeric_limits<double>::quiet_NaN();

                    // Elongation: eigenvalue ratio of 2D sample covariance
                    Eigen::Matrix2d cov = Eigen::Matrix2d::Zero();
                    double mx = sx / denom, my = sy / denom;
                    for (int j = 0; j < np; ++j) {
                        double dx = data.x[ids[j]] - mx;
                        double dy = data.y[ids[j]] - my;
                        cov(0,0) += dx*dx; cov(0,1) += dx*dy;
                        cov(1,0) += dx*dy; cov(1,1) += dy*dy;
                    }
                    cov /= np;
                    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(cov);
                    auto ev = eig.eigenvalues();  // ascending order
                    cell_stats_mat(ci, ci_map["elongation"]) = (ev(0) > 1e-10)
                        ? (ev(1) / ev(0)) : std::numeric_limits<double>::quiet_NaN();
                }

                // max_cluster_frac
                if (has_clmol && np > 0) {
                    std::unordered_map<int,int> cl_cnt;
                    for (int mol : ids) cl_cnt[bm_data.cluster_per_molecule[mol]]++;
                    int mc = 0;
                    for (auto& [k, v] : cl_cnt) if (v > mc) mc = v;
                    cell_stats_mat(ci, ci_map["max_cluster_frac"]) = static_cast<double>(mc) / np;
                }

                // lifespan
                if (has_lspan) {
                    int guid = bm_data.components[ci].guid;
                    auto it2 = lifespan_map.find(guid);
                    cell_stats_mat(ci, ci_map["lifespan"]) = (it2 != lifespan_map.end()) ? it2->second : -1;
                }
            }

            if (output_style == OutputStyle::Parquet) {
                save_cell_stat_df_parquet(cell_stats_mat, cell_names, cell_stat_col_names, out_paths.cell_stats);
            } else {
                save_cell_stat_df(cell_stats_mat, cell_names, cell_stat_col_names, out_paths.cell_stats);
            }
        }

        // Save cell polygons as GeoJSON
        PolygonCollection poly_joined;
        PolygonStack poly_stack;
        bool have_polygons = false;
        if (output_style == OutputStyle::Parquet || polygon_format != "none" || plot) {
            auto pos = data.position_matrix();
            auto polys = boundary_polygons_auto(
                pos, bm_data.assignment, /*estimate_per_z=*/(N == 3), &cell_names, /*verbose=*/true);
            poly_joined = std::move(polys.first);
            poly_stack = std::move(polys.second);
            have_polygons = true;
        }

        if (output_style == OutputStyle::Parquet || polygon_format != "none") {
            spdlog::info("Saving cell polygons...");
            if (output_style == OutputStyle::Parquet) {
                if (N == 3) {
                    save_polygon_stack_geoparquet(poly_stack, out_paths);
                } else {
                    save_polygons_geoparquet(poly_joined, out_paths.polygons_2d);
                }
            } else {
                if (N == 3) {
                    save_polygon_stack_geojson(poly_stack, out_paths, polygon_format);
                } else {
                    save_polygons_geojson(poly_joined, out_paths.polygons_2d, polygon_format);
                }
            }
        }

        // Save count matrix
        spdlog::info("Saving count matrix...");
        {
            int ng = data.n_genes();
            std::vector<std::string> cell_names_cm(n_cells_final);
            for (int i = 0; i < n_cells_final; ++i)
                cell_names_cm[i] = "cell_" + std::to_string(i + 1);

            // Build sparse count matrix: n_cells x n_genes
            // Using composition_data which is 0-based internally
            std::vector<Eigen::Triplet<float>> trips;
            auto ids_by_cell_cm = split_ids(bm_data.assignment, n_cells_final, true);
            for (int ci = 0; ci < n_cells_final; ++ci) {
                std::unordered_map<int,float> cnt;
                for (int mol : ids_by_cell_cm[ci]) {
                    int g = bm_data.composition_data[mol];  // 0-based
                    if (g >= 0 && g < ng) cnt[g] += 1.0f;
                }
                for (auto& [g, v] : cnt) {
                    trips.emplace_back(ci, g, v);
                }
            }
            if (output_style == OutputStyle::Parquet) {
                Eigen::SparseMatrix<float> cm(n_cells_final, ng);
                cm.setFromTriplets(trips.begin(), trips.end());
                save_matrix_to_10x_h5(cm, data.gene_names, cell_names_cm, out_paths.counts);
            } else if (count_matrix_format == "tsv") {
                Eigen::SparseMatrix<float> cm(n_cells_final, ng);
                cm.setFromTriplets(trips.begin(), trips.end());
                Eigen::SparseMatrix<double> cm_d = cm.cast<double>();
                save_matrix_to_tsv(cm_d, data.gene_names, cell_names_cm, out_paths.counts);
            } else {
                Eigen::SparseMatrix<float, Eigen::RowMajor> cm(n_cells_final, ng);
                cm.setFromTriplets(trips.begin(), trips.end());
                save_matrix_to_loom(cm, data.gene_names, cell_names_cm, out_paths.counts);
            }
        }

        // Save parameters dump
        save_params_toml(opts, cli_cmd, out_paths.params_dump);

        if (plot) {
            spdlog::info("Generating HTML run report...");
            auto diagnostic_html = generate_run_diagnostic_html(
                data,
                noise_edge_lengths,
                noise_fit,
                confidence_nn_id,
                bm_data.assignment,
                bm_data.n_components_trace,
                bm_data.assignment_confidence,
                clustering_result ? &*clustering_result : nullptr,
                ncv_report ? &*ncv_report : nullptr,
                cell_stats_mat,
                cell_stat_col_names,
                opts.prior,
                opts.segmentation.scale,
                opts.segmentation.scale_std
            );
            std::ofstream diag_file(out_paths.diagnostic_report);
            if (!diag_file) {
                spdlog::warn("Could not write diagnostic report '{}'", out_paths.diagnostic_report);
            } else {
                diag_file << diagnostic_html;
            }

            auto segmentation_html = generate_run_segmentation_html(
                data,
                bm_data.assignment,
                ncv_color,
                mol_clusters.empty() ? nullptr : &mol_clusters,
                have_polygons ? &poly_joined : nullptr
            );
            std::ofstream seg_plot_file(out_paths.molecule_plot);
            if (!seg_plot_file) {
                spdlog::warn("Could not write segmentation plot '{}'", out_paths.molecule_plot);
            } else {
                seg_plot_file << segmentation_html;
            }
        }

        spdlog::info("Results saved to '{}' in {} style", output, to_string(output_style));
    };

    if (data.is_3d()) {
        run_segmentation(std::integral_constant<int, 3>{});
    } else {
        run_segmentation(std::integral_constant<int, 2>{});
    }

    return 0;
}

// ============================================================================
// Subcommand: preview
// ============================================================================

int cmd_preview(
    const std::string& coordinates,
    RunOptions& opts,
    const std::string& output
) {
    spdlog::info("Loading data from '{}'...", coordinates);

    fill_and_check_molecule_input_options(opts.molecules);

    auto data = load_molecules(coordinates, opts.molecules);
    spdlog::info("Loaded {} transcripts, {} genes.", data.n_molecules(), data.n_genes());

    fill_and_check_plotting_options(opts.plotting, opts.molecules.min_molecules_per_cell, data.n_genes());

    // Confidence estimation — done once; reuse knn, adj_list, and noise_result
    // everywhere below instead of recomputing them.
    spdlog::info("Estimating noise level...");
    int nn_id = opts.molecules.confidence_nn_id;
    if (nn_id <= 0) nn_id = std::max(data.n_genes() / 10, 10);

    auto pos = data.position_matrix();
    // Block-wise kth-neighbour distances (same values the former full kNN
    // result held; only these distances are read).
    std::vector<double> edge_lengths = knn_kth_distances(pos, nn_id + 1, nn_id);

    auto adj_list    = build_molecule_graph(data, false);
    auto noise_result = fit_noise_probabilities(edge_lengths, adj_list, nullptr, 100, 0.005, true);

    data.confidence.resize(data.n_molecules());
    for (int i = 0; i < data.n_molecules(); ++i)
        data.confidence[i] = noise_result.assignment_probs(i, 0);

    spdlog::info("Done. Noise estimation complete.");

    // Gene composition colors
    spdlog::info("Estimating local colors...");
    int ncv_spatial_k = opts.plotting.gene_composition_neighborhood;
    int ncv_graph_k = opts.segmentation.cluster_graph_k;
    auto gene_colors = gene_composition_color_embedding_streaming(
        pos, data.gene, data.n_genes(), data.confidence, ncv_spatial_k,
        100000, 20000, 42, 10, ncv_graph_k
    );
    spdlog::info("Done.");

    // Gene structure embedding (reuses adj_list already built above)
    spdlog::info("Estimating gene structure...");
    auto gene_structure = estimate_gene_structure_embedding(
        data.gene, data.gene_names, data.confidence, adj_list);
    spdlog::info("Done.");

    // Generate HTML report
    spdlog::info("Generating HTML report...");
    auto html = generate_preview_html(data, gene_colors, edge_lengths, noise_result, nn_id, &gene_structure);

    std::ofstream out_file(output);
    if (!out_file) {
        spdlog::error("Could not write to '{}'", output);
        return 1;
    }
    out_file << html;
    out_file.close();

    spdlog::info("Preview saved to '{}'", output);
    return 0;
}

// ============================================================================
// Subcommand: segfree
// ============================================================================

int cmd_segfree(
    const std::string& coordinates,
    RunOptions& opts,
    int k_neighbors,
    const std::string& output
) {
    spdlog::info("Loading data from '{}'...", coordinates);

    fill_and_check_molecule_input_options(opts.molecules);

    auto data = load_molecules(coordinates, opts.molecules);
    spdlog::info("Loaded {} transcripts, {} genes.", data.n_molecules(), data.n_genes());

    if (k_neighbors <= 0) {
        k_neighbors = default_param_value(
            "composition_neighborhood", opts.molecules.min_molecules_per_cell, -1, data.n_genes());
    }
    spdlog::info("Using k={} neighbors for NCV composition.", k_neighbors);

    // Neighborhood count matrix: n_genes × n_mols sparse
    spdlog::info("Estimating neighborhoods...");
    auto pos = data.position_matrix();
    auto neighb_cm = neighborhood_count_matrix(pos, data.gene, k_neighbors, data.n_genes());

    // Log-transform (matching Julia and preview pipeline)
    for (int k = 0; k < neighb_cm.outerSize(); ++k) {
        for (Eigen::SparseMatrix<float>::InnerIterator it(neighb_cm, k); it; ++it) {
            it.valueRef() = static_cast<float>(std::log(it.value() * 10000.0f + 1e-5f));
        }
    }

    // Per-molecule confidence (noise model)
    spdlog::info("Estimating molecule confidences...");
    append_confidence(data, opts.molecules.confidence_nn_id);

    // Gene vectors via randomized indexing, then UMAP color embedding
    spdlog::info("Estimating gene colors...");
    auto mol_vecs = estimate_gene_vectors(neighb_cm, data.gene, 20, "ri", true);
    auto gene_colors = gene_composition_color_embedding(mol_vecs, data.confidence);

    // Cell/NCV names: "V{i}" (1-based), matching Julia's get_cell_name(:ncv)
    int n = data.n_molecules();
    std::vector<std::string> ncv_names(n);
    for (int i = 0; i < n; ++i) ncv_names[i] = "V" + std::to_string(i + 1);

    // Save Loom file.
    // Julia stores neighb_cm transposed: n_mols × n_genes (rows = cells, cols = genes).
    spdlog::info("Saving results to '{}'...", output);
    Eigen::SparseMatrix<float> ncv_mat = neighb_cm.transpose();

    LoomColAttrs col_attrs;
    col_attrs["ncv_color"]  = gene_colors;  // vector<string>
    col_attrs["confidence"] = std::vector<double>(data.confidence.begin(),
                                                   data.confidence.end());

    save_matrix_to_loom(ncv_mat, data.gene_names, ncv_names, output, col_attrs);

    spdlog::info("Done.");
    return 0;
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char* argv[]) {
    CLI::App app{"Baysor — Bayesian cell segmentation of spatial transcriptomics data"};
    app.set_version_flag("--version", std::string("baysor ") + BAYSOR_VERSION,
                         "Print the Baysor version and exit");
    app.require_subcommand(1);

    // Pre-scan argv for -c/--config so we can load config before CLI11 registers
    // option defaults. This way config values serve as defaults and explicit CLI
    // flags override them (proper config-then-CLI precedence).
    std::string config_path;
    for (int i = 1; i < argc - 1; ++i) {
        std::string arg = argv[i];
        if (arg == "-c" || arg == "--config") {
            config_path = argv[i + 1];
            break;
        }
    }

    RunOptions opts;
    if (!config_path.empty()) {
        try {
            opts = load_config(config_path);
            spdlog::info("Loaded config from '{}'", config_path);
        } catch (const std::exception& e) {
            spdlog::error("Failed to load config '{}': {}", config_path, e.what());
            return 1;
        }
    }

    // ---- run ----
    auto* run = app.add_subcommand("run", "Run cell segmentation");

    std::string run_coordinates, run_prior_seg;
    std::string run_output = "segmentation";
    std::string run_output_style = "legacy";
    std::string run_polygon_format = "FeatureCollection";
    std::string run_count_format = "loom";
    std::string run_cluster_method = cluster_method_to_string(opts.segmentation.cluster_method);
    bool run_plot = false;
    bool run_skip_ncv_color = false;

    run->add_option("coordinates", run_coordinates,
        "CSV/Parquet molecule table, or a Xenium experiment.xenium manifest")
        ->required();
    run->add_option("prior_segmentation", run_prior_seg,
        "Prior segmentation as image mask, boundary CSV/Parquet, or ':column_name' in the coordinates file");

    run->add_option("-c,--config", config_path,
        "TOML file with configuration");
    run->add_option("-x,--x-column", opts.molecules.x_col,
        "Name of x column (default: x)");
    run->add_option("-y,--y-column", opts.molecules.y_col,
        "Name of y column (default: y)");
    run->add_option("-z,--z-column", opts.molecules.z_col,
        "Name of z column (default: z)");
    run->add_option("-g,--gene-column", opts.molecules.gene_col,
        "Name of gene column (default: gene)");
    run->add_option("--qv-column", opts.molecules.qv_col,
        "Name of quality-value column used by --min-qv (default: qv)");
    run->add_option("-m,--min-molecules-per-cell", opts.molecules.min_molecules_per_cell,
        "Minimal number of molecules for a cell to be considered as real");
    run->add_option("-s,--scale", opts.segmentation.scale,
        "Approximate cell radius. Sets estimate-scale-from-centers to false");
    run->add_option("--scale-std", opts.segmentation.scale_std,
        "Std of scale across cells. Number or 'N%' relative to scale (default: 25%)");
    run->add_option("--cluster-method", run_cluster_method,
        "Molecule clustering prior: mrf, louvain, leiden, or none (default: mrf; legacy alias: ica_mrf)");
    run->add_option("--n-clusters", opts.segmentation.n_clusters,
        "Target number of molecule clusters / major cell types (exact for mrf; merged target for louvain/leiden; default: 4 for mrf, 10 for louvain/leiden)");
    run->add_option("--cluster-resolution", opts.segmentation.cluster_resolution,
        "Advanced overclustering resolution for cluster-method=louvain or leiden (default: 1.0)");
    run->add_option("--cluster-graph-k", opts.segmentation.cluster_graph_k,
        "Number of NCV nearest neighbors used for graph clustering and NCV UMAPs (default: 15)");
    run->add_option("--cluster-n-dims", opts.segmentation.cluster_n_dims,
        "Number of NCV dimensions used by cluster-method=louvain or leiden (default: 20)");
    run->add_option("--cluster-basis-sample-size", opts.segmentation.cluster_basis_sample_size,
        "Maximum number of basis anchors used by cluster-method=louvain or leiden (default: 100000)");
    run->add_option("--prior-segmentation-confidence", opts.segmentation.prior_segmentation_confidence,
        "Confidence of prior segmentation results, in [0,1] (default: 0.2)");
    run->add_option("--min-molecules-per-gene", opts.molecules.min_molecules_per_gene,
        "Minimal number of molecules per gene (default: 1)");
    run->add_option("--exclude-genes", opts.molecules.exclude_genes,
        "Comma-separated list of genes or patterns to exclude (e.g. 'Blank*,MALAT1')");
    run->add_option("--min-qv", opts.molecules.min_qv,
        "Drop molecules with qv below this threshold during input loading");
    run->add_option("--x-min", opts.molecules.x_min,
        "Minimum x coordinate to keep during input loading");
    run->add_option("--x-max", opts.molecules.x_max,
        "Maximum x coordinate to keep during input loading");
    run->add_option("--y-min", opts.molecules.y_min,
        "Minimum y coordinate to keep during input loading");
    run->add_option("--y-max", opts.molecules.y_max,
        "Maximum y coordinate to keep during input loading");
    run->add_option("--z-min", opts.molecules.z_min,
        "Minimum z coordinate to keep during input loading");
    run->add_option("--z-max", opts.molecules.z_max,
        "Maximum z coordinate to keep during input loading");
    run->add_option("--nuclei-genes", opts.segmentation.nuclei_genes,
        "Comma-separated list of nuclei-specific genes");
    run->add_option("--cyto-genes", opts.segmentation.cyto_genes,
        "Comma-separated list of cytoplasm-specific genes");
    run->add_option("-o,--output", run_output,
        "Output directory (default: segmentation)");
    run->add_option("--output-style", run_output_style,
        "Output bundle style: legacy or parquet (default: legacy)");
    run->add_option("--polygon-format", run_polygon_format,
        "Polygon output format: FeatureCollection, GeometryCollection, or none (default: FeatureCollection)");
    run->add_option("--count-matrix-format", run_count_format,
        "Count matrix format: loom or tsv (default: loom)");
    run->add_flag("-p,--plot", run_plot,
        "Save an HTML diagnostic plot");
    run->add_flag("--skip-ncv-color", run_skip_ncv_color,
        "Skip neighborhood composition color embedding to speed up development runs");
    run->add_flag("--force-2d", opts.molecules.force_2d,
        "Ignore z-column in the data");
    run->add_option("--iters", opts.segmentation.iters,
        "Maximum number of algorithm iterations (default: 500)");
    run->add_option("--tol", opts.segmentation.tol,
        "Convergence tolerance: stop when <tol fraction of molecules change assignment "
        "over 20 consecutive iterations. 0 = always run all --iters (default: 0)");
    run->add_option("--n-cells-init", opts.segmentation.n_cells_init,
        "Initial number of cells (default: auto)");
    run->add_option("--unassigned-prior-label", opts.prior.unassigned_label,
        "Label for unassigned cells in prior segmentation (default: 0)");
    run->add_option("-t,--threads", opts.threads,
        "Number of worker threads (default: BAYSOR_NUM_THREADS, then OMP_NUM_THREADS, "
        "then the number of physical CPU cores)");

    // ---- preview ----
    auto* preview = app.add_subcommand("preview", "Plot a dataset preview");

    std::string prev_coordinates;
    std::string prev_output = "preview.html";

    preview->add_option("coordinates", prev_coordinates,
        "CSV/Parquet molecule table, or a Xenium experiment.xenium manifest")
        ->required();
    preview->add_option("-c,--config", config_path,
        "TOML file with configuration");
    preview->add_option("-x,--x-column", opts.molecules.x_col,
        "Name of x column (default: x)");
    preview->add_option("-y,--y-column", opts.molecules.y_col,
        "Name of y column (default: y)");
    preview->add_option("-z,--z-column", opts.molecules.z_col,
        "Name of z column (default: z)");
    preview->add_option("-g,--gene-column", opts.molecules.gene_col,
        "Name of gene column (default: gene)");
    preview->add_option("--qv-column", opts.molecules.qv_col,
        "Name of quality-value column used by --min-qv (default: qv)");
    preview->add_option("-m,--min-molecules-per-cell", opts.molecules.min_molecules_per_cell,
        "Minimal number of molecules for a cell to be considered as real");
    preview->add_option("--min-qv", opts.molecules.min_qv,
        "Drop molecules with qv below this threshold during input loading");
    preview->add_option("--x-min", opts.molecules.x_min,
        "Minimum x coordinate to keep during input loading");
    preview->add_option("--x-max", opts.molecules.x_max,
        "Maximum x coordinate to keep during input loading");
    preview->add_option("--y-min", opts.molecules.y_min,
        "Minimum y coordinate to keep during input loading");
    preview->add_option("--y-max", opts.molecules.y_max,
        "Maximum y coordinate to keep during input loading");
    preview->add_option("--z-min", opts.molecules.z_min,
        "Minimum z coordinate to keep during input loading");
    preview->add_option("--z-max", opts.molecules.z_max,
        "Maximum z coordinate to keep during input loading");
    preview->add_option("-o,--output", prev_output,
        "Output HTML file (default: preview.html)");
    preview->add_flag("--force-2d", opts.molecules.force_2d,
        "Ignore z-column in the data");
    preview->add_option("-t,--threads", opts.threads,
        "Number of worker threads (default: BAYSOR_NUM_THREADS, then OMP_NUM_THREADS, "
        "then the number of physical CPU cores)");

    // ---- segfree ----
    auto* segfree = app.add_subcommand("segfree", "Extract Neighborhood Composition Vectors (NCVs)");

    std::string sf_coordinates;
    std::string sf_output = "ncvs.loom";
    int sf_k_neighbors = 0;

    segfree->add_option("coordinates", sf_coordinates,
        "CSV/Parquet molecule table, or a Xenium experiment.xenium manifest")
        ->required();
    segfree->add_option("-c,--config", config_path,
        "TOML file with configuration");
    segfree->add_option("-x,--x-column", opts.molecules.x_col,
        "Name of x column (default: x)");
    segfree->add_option("-y,--y-column", opts.molecules.y_col,
        "Name of y column (default: y)");
    segfree->add_option("-z,--z-column", opts.molecules.z_col,
        "Name of z column (default: z)");
    segfree->add_option("-g,--gene-column", opts.molecules.gene_col,
        "Name of gene column (default: gene)");
    segfree->add_option("--qv-column", opts.molecules.qv_col,
        "Name of quality-value column used by --min-qv (default: qv)");
    segfree->add_option("-m,--min-molecules-per-cell", opts.molecules.min_molecules_per_cell,
        "Minimal number of molecules for a cell to be considered as real");
    segfree->add_option("--min-qv", opts.molecules.min_qv,
        "Drop molecules with qv below this threshold during input loading");
    segfree->add_option("--x-min", opts.molecules.x_min,
        "Minimum x coordinate to keep during input loading");
    segfree->add_option("--x-max", opts.molecules.x_max,
        "Maximum x coordinate to keep during input loading");
    segfree->add_option("--y-min", opts.molecules.y_min,
        "Minimum y coordinate to keep during input loading");
    segfree->add_option("--y-max", opts.molecules.y_max,
        "Maximum y coordinate to keep during input loading");
    segfree->add_option("--z-min", opts.molecules.z_min,
        "Minimum z coordinate to keep during input loading");
    segfree->add_option("--z-max", opts.molecules.z_max,
        "Maximum z coordinate to keep during input loading");
    segfree->add_option("-k,--k-neighbors", sf_k_neighbors,
        "Number of neighbors for segmentation-free pseudo-cells (default: inferred)");
    segfree->add_option("-o,--output", sf_output,
        "Output .loom file (default: ncvs.loom)");
    segfree->add_flag("--force-2d", opts.molecules.force_2d,
        "Ignore z-column in the data");
    segfree->add_option("-t,--threads", opts.threads,
        "Number of worker threads (default: BAYSOR_NUM_THREADS, then OMP_NUM_THREADS, "
        "then the number of physical CPU cores)");

    // ---- Parse ----
    CLI11_PARSE(app, argc, argv);

    // Configure the global thread pool once, before any parallel work.
    int n_threads = resolve_thread_count(opts.threads);
    set_thread_pool_size(n_threads);
    spdlog::info("Using {} threads", n_threads);

    // Reconstruct CLI command string for params dump
    std::string cli_cmd;
    for (int i = 0; i < argc; ++i) {
        if (i > 0) cli_cmd += " ";
        cli_cmd += argv[i];
    }

    auto resolve_xenium_input = [](const std::string& coordinates) -> std::string {
        if (!is_xenium_manifest_path(coordinates)) {
            return coordinates;
        }
        auto xenium_ctx = load_xenium_manifest_context(coordinates);
        return xenium_ctx.transcripts_path;
    };

    // Dispatch
    try {
        if (run->parsed()) {
            try {
                opts.segmentation.cluster_method = parse_cluster_method(run_cluster_method);
            } catch (const std::exception& e) {
                spdlog::error("{}", e.what());
                return 1;
            }
            if (opts.segmentation.n_clusters <= 0) {
                opts.segmentation.n_clusters = default_cluster_count(opts.segmentation.cluster_method);
            }
            OutputStyle output_style;
            try {
                output_style = parse_output_style(run_output_style);
            } catch (const std::exception& e) {
                spdlog::error("{}", e.what());
                return 1;
            }
            std::string resolved_run_input = run_coordinates;
            try {
                resolved_run_input = resolve_xenium_input(run_coordinates);
            } catch (const std::exception& e) {
                spdlog::error("{}", e.what());
                return 1;
            }
            if (output_style != OutputStyle::Legacy) {
                if (run_polygon_format != "FeatureCollection") {
                    spdlog::warn("--polygon-format is ignored for output style '{}'", run_output_style);
                }
                if (run_count_format != "loom") {
                    spdlog::warn("--count-matrix-format is ignored for output style '{}'", run_output_style);
                }
            }
            if (!run_prior_seg.empty()) {
                apply_prior_input_spec(opts.prior, run_prior_seg);
            }
            if (opts.segmentation.scale > 0) {
                opts.prior.estimate_scale_from_prior = false;
            }
            if (opts.prior.type == PriorInputType::None && opts.segmentation.scale <= 0) {
                spdlog::error("Either prior_segmentation or --scale must be provided.");
                return 1;
            }
            return cmd_run(resolved_run_input, opts, run_output, output_style,
                           run_plot, run_skip_ncv_color, run_polygon_format, run_count_format, cli_cmd);
        }

        if (preview->parsed()) {
            return cmd_preview(resolve_xenium_input(prev_coordinates), opts, prev_output);
        }

        if (segfree->parsed()) {
            return cmd_segfree(resolve_xenium_input(sf_coordinates), opts, sf_k_neighbors, sf_output);
        }
    } catch (const std::exception& e) {
        spdlog::error("{}", e.what());
        return 1;
    } catch (...) {
        spdlog::error("Unknown error"); // GCOVR_EXCL_LINE: unreachable, every exception thrown by Baysor or its dependencies derives from std::exception
        return 1; // GCOVR_EXCL_LINE: unreachable, only reachable via the never-entered catch-all above
    }

    return 0; // GCOVR_EXCL_LINE: unreachable, require_subcommand(1) guarantees exactly one subcommand and run/preview/segfree are all dispatched above
}
