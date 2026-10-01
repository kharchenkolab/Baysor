#include "baysor/reporting/color_utils.h"
#include "baysor/processing/data_processing/neighborhood_composition.h"
#include "baysor/processing/data_processing/umap_wrappers.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/processing/utils/utils.h"
#include "baysor/utils/thread_pool.h"

#include <Eigen/SVD>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <cstdio>

namespace baysor {

static constexpr const char* NCV_FALLBACK_COLOR = "#808080";

// ============================================================================
// Quantile helper (for internal use)
// ============================================================================

static double quantile_vec(const std::vector<double>& v, double p) {
    if (v.empty()) return 0.0;
    std::vector<double> s = v;
    std::sort(s.begin(), s.end());
    double idx = p * (s.size() - 1);
    int lo = static_cast<int>(std::floor(idx));
    int hi = std::min(lo + 1, static_cast<int>(s.size()) - 1);
    double frac = idx - lo;
    return s[lo] * (1.0 - frac) + s[hi] * frac;
}

// ============================================================================
// normalize_embedding_to_lab_range
// ============================================================================

void normalize_embedding_to_lab_range(
    Eigen::MatrixXd& embedding,
    double l_min, double l_max,
    double trim_frac,
    bool log_colors
) {
    int n = static_cast<int>(embedding.cols());
    if (n == 0) return;

    // Per-row: subtract quantile(row, trim_frac), clamp to 0
    for (int r = 0; r < 3; ++r) {
        std::vector<double> row_vals(n);
        for (int i = 0; i < n; ++i) row_vals[i] = embedding(r, i);
        double q_lo = quantile_vec(row_vals, trim_frac);
        for (int i = 0; i < n; ++i) {
            embedding(r, i) = std::max(embedding(r, i) - q_lo, 0.0);
        }
    }

    // Global max quantile
    std::vector<double> all_vals;
    all_vals.reserve(3 * n);
    for (int r = 0; r < 3; ++r) {
        for (int i = 0; i < n; ++i) {
            all_vals.push_back(embedding(r, i));
        }
    }
    double max_val = quantile_vec(all_vals, 1.0 - trim_frac);
    if (max_val > 0) {
        embedding /= max_val;
    }
    embedding = embedding.cwiseMin(1.0);

    if (log_colors) {
        // Compute global 5th percentile
        all_vals.clear();
        for (int r = 0; r < 3; ++r) {
            for (int i = 0; i < n; ++i) {
                all_vals.push_back(embedding(r, i));
            }
        }
        double q05 = std::max(quantile_vec(all_vals, 0.05), 1e-3);
        for (int r = 0; r < 3; ++r) {
            for (int i = 0; i < n; ++i) {
                embedding(r, i) = std::log10(embedding(r, i) + q05);
            }
        }
        // Per-row normalize to [0, 1]
        for (int r = 0; r < 3; ++r) {
            double mn = embedding.row(r).minCoeff();
            embedding.row(r).array() -= mn;
            double mx = embedding.row(r).maxCoeff();
            if (mx > 0) embedding.row(r) /= mx;
        }
    }

    // L channel: scale to [l_min, l_max]
    embedding.row(0) *= (l_max - l_min);
    embedding.row(0).array() += l_min;

    // a, b channels: scale to [-100, 100]
    embedding.row(1).array() -= 0.5;
    embedding.row(1) *= 200.0;
    embedding.row(2).array() -= 0.5;
    embedding.row(2) *= 200.0;
}

// ============================================================================
// LAB → sRGB → hex conversion
// ============================================================================

// Standard LAB → XYZ → linear RGB → sRGB pipeline
// Reference white: D65 (X=0.95047, Y=1.0, Z=1.08883)

static inline double lab_f_inv(double t) {
    constexpr double delta = 6.0 / 29.0;
    if (t > delta) return t * t * t;
    return 3.0 * delta * delta * (t - 4.0 / 29.0);
}

static inline double srgb_gamma(double c) {
    if (c <= 0.0031308) return 12.92 * c;
    return 1.055 * std::pow(c, 1.0 / 2.4) - 0.055;
}

static inline int clamp_byte(double v) {
    int b = static_cast<int>(std::round(v * 255.0));
    return std::max(0, std::min(255, b));
}

static std::string lab_to_hex(double L, double a, double b) {
    // LAB → XYZ
    double fy = (L + 16.0) / 116.0;
    double fx = a / 500.0 + fy;
    double fz = fy - b / 200.0;

    double X = 0.95047 * lab_f_inv(fx);
    double Y = 1.00000 * lab_f_inv(fy);
    double Z = 1.08883 * lab_f_inv(fz);

    // XYZ → linear RGB (sRGB D65)
    double rlin =  3.2404542 * X - 1.5371385 * Y - 0.4985314 * Z;
    double glin = -0.9692660 * X + 1.8760108 * Y + 0.0415560 * Z;
    double blin =  0.0556434 * X - 0.2040259 * Y + 1.0572252 * Z;

    // Linear RGB → sRGB with gamma
    int R = clamp_byte(srgb_gamma(rlin));
    int G = clamp_byte(srgb_gamma(glin));
    int B = clamp_byte(srgb_gamma(blin));

    char buf[8];
    std::snprintf(buf, sizeof(buf), "#%02X%02X%02X", R, G, B);
    return std::string(buf);
}

std::vector<std::string> embedding_to_hex(const Eigen::MatrixXd& lab_embedding) {
    int n = static_cast<int>(lab_embedding.cols());
    std::vector<std::string> colors(n);
    for (int i = 0; i < n; ++i) {
        colors[i] = lab_to_hex(lab_embedding(0, i), lab_embedding(1, i), lab_embedding(2, i));
    }
    return colors;
}

// ============================================================================
// gene_composition_color_embedding (UMAP-based)
// ============================================================================

namespace {

struct NcvSampleSelection {
    std::vector<int> sample_ids;
    double chosen_threshold = 0.95;
    int anchor_count = 0;
};

struct LabNormalizationParams {
    double l_min = 10.0;
    double l_max = 90.0;
    double trim_frac = 0.0125;
    bool log_colors = false;
    std::array<double, 3> row_q_lo{0.0, 0.0, 0.0};
    double max_val = 1.0;
    double q05 = 1e-3;
    std::array<double, 3> log_row_min{0.0, 0.0, 0.0};
    std::array<double, 3> log_row_scale{1.0, 1.0, 1.0};
};

struct NcvInterpolationModel {
    Eigen::VectorXf sample_mean;
    Eigen::MatrixXf pca_basis;
    Eigen::MatrixXd anchor_pca;
    Eigen::MatrixXd anchor_emb;
    int interp_k = 5;
};

NcvSampleSelection select_ncv_sample_ids(
    const Eigen::MatrixXf& mol_vecs,
    const std::vector<double>& confidence,
    int sample_size
);

LabNormalizationParams fit_lab_normalization_params(
    const Eigen::MatrixXd& embedding,
    double l_min,
    double l_max,
    double trim_frac,
    bool log_colors
);

void apply_lab_normalization_params(
    Eigen::MatrixXd& embedding,
    const LabNormalizationParams& params
);

NcvInterpolationModel fit_ncv_interpolation_model(
    const Eigen::MatrixXf& anchor_vecs,
    const Eigen::MatrixXd& anchor_emb,
    int n_pca_dims,
    int graph_k
);

Eigen::MatrixXd interpolate_ncv_embedding(
    const NcvInterpolationModel& model,
    const Eigen::MatrixXf& query_vecs
);

static void fill_colors_from_projected_vectors(
    NcvReportEmbedding& result,
    const Eigen::MatrixXf& basis_vecs,
    const std::vector<int>& basis_ids,
    const std::vector<double>& basis_conf,
    const Eigen::MatrixXf* all_mol_vecs,
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    const Eigen::MatrixXf& gene_emb_t,
    int n_genes,
    int k_neighbors,
    double dist_floor,
    int sample_size,
    int seed,
    int n_pca_dims,
    int graph_k,
    bool include_report_umap
) {
    NcvSampleSelection selection = select_ncv_sample_ids(basis_vecs, basis_conf, sample_size);
    result.chosen_threshold = selection.chosen_threshold;
    result.anchor_count = selection.anchor_count;
    if (static_cast<int>(selection.sample_ids.size()) <= 1) {
        spdlog::warn("NCV color embedding fallback: insufficient sampled anchors after basis selection.");
        return;
    }

    const int n_sample = static_cast<int>(selection.sample_ids.size());
    spdlog::info(
        "NCV color embedding: selected {} anchors with confidence >= {:.2f}; sampling {} for UMAP fit.",
        selection.anchor_count, selection.chosen_threshold, n_sample
    );

    Eigen::MatrixXf sample_vecs(basis_vecs.rows(), n_sample);
    for (int i = 0; i < n_sample; ++i) sample_vecs.col(i) = basis_vecs.col(selection.sample_ids[i]);
    // The 3-D colour UMAP and the 2-D report UMAP share one kNN graph.
    const UmapGraph umap_graph = umap_fuzzy_graph(sample_vecs.cast<double>(), graph_k);
    Eigen::MatrixXd sample_emb = umap_embed_graph(umap_graph, 3, 200, seed, 2.0);

    if (include_report_umap) {
        Eigen::MatrixXd sample_umap2d = umap_embed_graph(umap_graph, 2, 200, seed, 2.0);
        result.sample_ids.resize(n_sample);
        result.sample_umap_x.resize(n_sample);
        result.sample_umap_y.resize(n_sample);
        for (int i = 0; i < n_sample; ++i) {
            result.sample_ids[i] = basis_ids[selection.sample_ids[i]];
            result.sample_umap_x[i] = sample_umap2d(0, i);
            result.sample_umap_y[i] = sample_umap2d(1, i);
        }
    }

    NcvInterpolationModel sample_model = fit_ncv_interpolation_model(
        sample_vecs, sample_emb, n_pca_dims, graph_k
    );
    Eigen::MatrixXd basis_emb = interpolate_ncv_embedding(sample_model, basis_vecs);
    LabNormalizationParams lab_params = fit_lab_normalization_params(basis_emb, 10.0, 90.0, 0.0125, false);

    if (all_mol_vecs) {
        constexpr int block_size = 32768;
        for (int block_start = 0; block_start < all_mol_vecs->cols(); block_start += block_size) {
            int block_n = std::min(block_size, static_cast<int>(all_mol_vecs->cols()) - block_start);
            Eigen::MatrixXf block_vecs = all_mol_vecs->middleCols(block_start, block_n);
            Eigen::MatrixXd block_emb = interpolate_ncv_embedding(sample_model, block_vecs);
            apply_lab_normalization_params(block_emb, lab_params);
            auto block_colors = embedding_to_hex(block_emb);
            for (int i = 0; i < block_n; ++i) {
                result.colors[static_cast<size_t>(block_start + i)] = block_colors[static_cast<size_t>(i)];
            }
        }
        return;
    }

    stream_projected_neighborhood_vectors(
        pos_data, genes, k_neighbors, gene_emb_t, n_genes, nullptr, nullptr,
        true, true, dist_floor, true, 32768,
        [&](int, const std::vector<int>& block_query_ids, const Eigen::MatrixXf& block_vecs) {
            Eigen::MatrixXd block_emb = interpolate_ncv_embedding(sample_model, block_vecs);
            apply_lab_normalization_params(block_emb, lab_params);
            auto block_colors = embedding_to_hex(block_emb);
            for (size_t i = 0; i < block_query_ids.size(); ++i) {
                result.colors[block_query_ids[i]] = block_colors[i];
            }
        }
    );
}

NcvSampleSelection select_ncv_sample_ids(
    const Eigen::MatrixXf& mol_vecs,
    const std::vector<double>& confidence,
    int sample_size
) {
    const int n_mols = static_cast<int>(mol_vecs.cols());
    sample_size = std::min(sample_size, n_mols);

    const int min_anchor_count = std::min(sample_size, n_mols);
    const std::array<double, 8> threshold_ladder = {0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.60, 0.50};

    std::vector<int> high_conf_ids;
    double chosen_threshold = threshold_ladder.front();
    for (double thr : threshold_ladder) {
        high_conf_ids.clear();
        high_conf_ids.reserve(n_mols);
        for (int i = 0; i < n_mols; ++i) {
            if (confidence[i] >= thr) high_conf_ids.push_back(i);
        }
        chosen_threshold = thr;
        if (static_cast<int>(high_conf_ids.size()) >= min_anchor_count) {
            break;
        }
    }

    if (static_cast<int>(high_conf_ids.size()) < sample_size) {
        sample_size = static_cast<int>(high_conf_ids.size());
    }

    NcvSampleSelection result;
    result.chosen_threshold = chosen_threshold;
    result.anchor_count = static_cast<int>(high_conf_ids.size());
    if (sample_size <= 1) return result;

    int hc = static_cast<int>(high_conf_ids.size());
    std::vector<std::pair<float, int>> sum_ids(hc);
    for (int i = 0; i < hc; ++i) {
        float s = mol_vecs.col(high_conf_ids[i]).sum();
        sum_ids[i] = {s, high_conf_ids[i]};
    }
    std::sort(sum_ids.begin(), sum_ids.end());
    for (int i = 0; i < hc; ++i) high_conf_ids[i] = sum_ids[i].second;

    result.sample_ids.reserve(sample_size);
    for (int i = 0; i < sample_size; ++i) {
        int idx = static_cast<int>(std::round(static_cast<double>(i) * (hc - 1) / (sample_size - 1)));
        result.sample_ids.push_back(high_conf_ids[idx]);
    }
    return result;
}

NcvReportEmbedding compute_ncv_embedding(
    const Eigen::MatrixXf& mol_vecs,
    const std::vector<double>& confidence,
    int sample_size,
    int seed,
    int n_pca_dims,
    int graph_k,
    bool include_report_umap
) {
    const int n_components = static_cast<int>(mol_vecs.rows());
    const int n_mols = static_cast<int>(mol_vecs.cols());

    NcvReportEmbedding result;
    result.colors.assign(n_mols, NCV_FALLBACK_COLOR);

    if (n_mols == 0 || n_components < 3) {
        return result;
    }

    NcvSampleSelection selection = select_ncv_sample_ids(mol_vecs, confidence, sample_size);
    result.chosen_threshold = selection.chosen_threshold;
    result.anchor_count = selection.anchor_count;
    result.sample_ids = selection.sample_ids;
    sample_size = static_cast<int>(selection.sample_ids.size());

    // Large runs can legitimately have no molecules above the hard 0.95
    // confidence cutoff historically used for NCV-color anchor selection. In that case
    // segmentation is still valid, but the UMAP fit/interpolation path below
    // would hit divisions by (sample_size - 1) and invalid KNN sizes.
    if (sample_size <= 1) {
        spdlog::warn( // GCOVR_EXCL_LINE: gcov: only an exception-unwind block is attributed here; the executed warn() call is reported on the following lines
            "NCV color embedding fallback: insufficient anchor molecules after adaptive thresholding "
            "(max_conf={:.4f}, threshold={:.2f}, anchors={}, sample_size={}).",
            confidence.empty() ? 0.0 : *std::max_element(confidence.begin(), confidence.end()),
            selection.chosen_threshold,
            selection.anchor_count,
            sample_size
        );
        return result;
    }

    spdlog::info(
        "NCV color embedding: selected {} anchors with confidence >= {:.2f}; sampling {} for UMAP fit.",
        selection.anchor_count, selection.chosen_threshold, sample_size
    );

    // Build float sample matrix (n_components x sample_size) — no cast needed.
    Eigen::MatrixXf sample_mat(n_components, sample_size);
    for (int i = 0; i < sample_size; ++i)
        sample_mat.col(i) = mol_vecs.col(selection.sample_ids[i]);

    // UMAP fit: use spread=2.0 to match Julia's UmapFit defaults, which produce
    // better colour separation than the umappp default of spread=1.0.
    // The 3-D colour UMAP and the 2-D report UMAP share one kNN graph.
    const UmapGraph umap_graph = umap_fuzzy_graph(sample_mat.cast<double>(), /*n_neighbors=*/graph_k);
    Eigen::MatrixXd sample_emb = umap_embed_graph(umap_graph, 3, /*n_epochs=*/200, seed, /*spread=*/2.0);

    if (include_report_umap) {
        Eigen::MatrixXd sample_umap2d = umap_embed_graph(umap_graph, 2, /*n_epochs=*/200, seed, /*spread=*/2.0);
        result.sample_umap_x.resize(sample_size);
        result.sample_umap_y.resize(sample_size);
        for (int i = 0; i < sample_size; ++i) {
            result.sample_umap_x[i] = sample_umap2d(0, i);
            result.sample_umap_y[i] = sample_umap2d(1, i);
        }
    }

    // PCA-reduce sample from n_components → n_pca dimensions before building the
    // interpolation KNN index.  In 20D the VPtree degenerates to exhaustive search
    // (curse of dimensionality), making interpolation O(N_mols × N_sample).
    // In 3D the VPtree search is O(log N_sample), giving a ~1000× speedup.
    // The first few PCs capture the dominant variance in gene-expression space,
    // so neighbour quality changes are minimal for colour interpolation.
    Eigen::VectorXf sample_mean = sample_mat.rowwise().mean();
    Eigen::MatrixXf sample_centered = sample_mat.colwise() - sample_mean;

    // Thin U of (n_components × sample_size): columns are principal components.
    Eigen::BDCSVD<Eigen::MatrixXf> svd(sample_centered, Eigen::ComputeThinU);
    // Thin U only has min(rows, cols) columns, so clamp to min(n_pca_dims,
    // n_components, sample_size): with fewer sampled anchors than n_pca_dims the
    // old clamp asked for more columns than exist (Eigen assertion in Debug,
    // out-of-bounds read in Release). Julia's gene_composition_color_embedding
    // performs no PCA truncation at all (it interpolates on the raw sampled
    // columns), so using every available component matches it as closely as
    // possible.
    const int n_pca = std::min({n_pca_dims, n_components, sample_size});
    Eigen::MatrixXf pca_basis = svd.matrixU().leftCols(n_pca);  // n_components × n_pca

    // Project sample to PCA space: n_pca × sample_size
    Eigen::MatrixXf sample_pca =
        (pca_basis.transpose() * sample_centered).eval();

    // Project all molecules to PCA space: n_pca × n_mols
    Eigen::MatrixXd all_pca =
        (pca_basis.transpose() * (mol_vecs.colwise() - sample_mean))
        .cast<double>();

    // KNN in PCA-3D space: nanoflann KD-tree via knn_parallel (already OMP-parallel).
    // tree = sample_pca (n_pca × sample_size), query = all_pca (n_pca × n_mols).
    // graph_k <= 0 must not reach knn_parallel, which then returns no results
    // and the interpolation loop below reads past the empty index lists.
    int k_interp = std::max(1, std::min(graph_k, sample_size - 1));
    auto knn = knn_parallel(sample_pca.cast<double>(), all_pca, k_interp);

    // Weighted interpolation of UMAP coordinates.
    Eigen::MatrixXd emb(3, n_mols);
    constexpr double dist_offset = 1e-10;

    parallel_for_static(0, n_mols, [&](int i) {
        double w_sum = 0.0;
        Eigen::Vector3d weighted = Eigen::Vector3d::Zero();
        const int* nn_indices = knn.idx_row(i);
        const double* nn_distances = knn.dist_row(i);
        for (int j = 0; j < knn.k; ++j) {
            double w = 1.0 / (nn_distances[j] + dist_offset);
            weighted += w * sample_emb.col(nn_indices[j]);
            w_sum    += w;
        }
        emb.col(i) = weighted / w_sum;
    });

    normalize_embedding_to_lab_range(emb);
    result.colors = embedding_to_hex(emb);
    return result;
}

LabNormalizationParams fit_lab_normalization_params(
    const Eigen::MatrixXd& embedding,
    double l_min = 10.0, double l_max = 90.0,
    double trim_frac = 0.0125,
    bool log_colors = false
) {
    LabNormalizationParams params;
    params.l_min = l_min;
    params.l_max = l_max;
    params.trim_frac = trim_frac;
    params.log_colors = log_colors;

    const int n = static_cast<int>(embedding.cols());
    if (n == 0) return params;

    Eigen::MatrixXd work = embedding;
    for (int r = 0; r < 3; ++r) {
        std::vector<double> row_vals(n);
        for (int i = 0; i < n; ++i) row_vals[i] = work(r, i);
        params.row_q_lo[r] = quantile_vec(row_vals, trim_frac);
        for (int i = 0; i < n; ++i) {
            work(r, i) = std::max(work(r, i) - params.row_q_lo[r], 0.0);
        }
    }

    std::vector<double> all_vals;
    all_vals.reserve(3 * n);
    for (int r = 0; r < 3; ++r)
        for (int i = 0; i < n; ++i)
            all_vals.push_back(work(r, i));
    params.max_val = quantile_vec(all_vals, 1.0 - trim_frac);
    if (params.max_val <= 0.0) params.max_val = 1.0;

    work /= params.max_val;
    work = work.cwiseMin(1.0);

    if (log_colors) {
        // GCOVR_EXCL_START: unreachable: the sole caller of fit_lab_normalization_params hardcodes log_colors=false
        all_vals.clear();
        for (int r = 0; r < 3; ++r)
            for (int i = 0; i < n; ++i)
                all_vals.push_back(work(r, i));
        params.q05 = std::max(quantile_vec(all_vals, 0.05), 1e-3);
        for (int r = 0; r < 3; ++r) {
            for (int i = 0; i < n; ++i) work(r, i) = std::log10(work(r, i) + params.q05);
            params.log_row_min[r] = work.row(r).minCoeff();
            work.row(r).array() -= params.log_row_min[r];
            params.log_row_scale[r] = work.row(r).maxCoeff();
            if (params.log_row_scale[r] <= 0.0) params.log_row_scale[r] = 1.0;
        // GCOVR_EXCL_STOP: unreachable: the sole caller of fit_lab_normalization_params hardcodes log_colors=false
        }
    }

    return params;
}

void apply_lab_normalization_params(Eigen::MatrixXd& embedding, const LabNormalizationParams& params) {
    const int n = static_cast<int>(embedding.cols());
    if (n == 0) return;

    for (int r = 0; r < 3; ++r) {
        for (int i = 0; i < n; ++i) {
            embedding(r, i) = std::max(embedding(r, i) - params.row_q_lo[r], 0.0);
        }
    }

    embedding /= params.max_val;
    embedding = embedding.cwiseMin(1.0);

    if (params.log_colors) {
        // GCOVR_EXCL_START: unreachable: params.log_colors can only be false, see fit_lab_normalization_params
        for (int r = 0; r < 3; ++r) {
            for (int i = 0; i < n; ++i) {
                embedding(r, i) = std::log10(embedding(r, i) + params.q05);
            }
            embedding.row(r).array() -= params.log_row_min[r];
            embedding.row(r) /= params.log_row_scale[r];
        // GCOVR_EXCL_STOP: unreachable: params.log_colors can only be false, see fit_lab_normalization_params
        }
    }

    embedding.row(0) *= (params.l_max - params.l_min);
    embedding.row(0).array() += params.l_min;
    embedding.row(1).array() -= 0.5;
    embedding.row(1) *= 200.0;
    embedding.row(2).array() -= 0.5;
    embedding.row(2) *= 200.0;
}

std::vector<int> select_basis_anchor_ids(
    const Eigen::MatrixXd& pos_data,
    const std::vector<double>& confidence,
    int basis_sample_size
) {
    const int n = static_cast<int>(confidence.size());
    basis_sample_size = std::min(basis_sample_size, n);
    const std::array<double, 8> threshold_ladder = {0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.60, 0.50};
    std::vector<int> candidates;
    double chosen_threshold = threshold_ladder.front();
    for (double thr : threshold_ladder) {
        candidates.clear();
        candidates.reserve(n);
        for (int i = 0; i < n; ++i) if (confidence[i] >= thr) candidates.push_back(i);
        chosen_threshold = thr;
        if (static_cast<int>(candidates.size()) >= basis_sample_size) break;
    }
    (void)chosen_threshold;
    if (static_cast<int>(candidates.size()) <= basis_sample_size) return candidates;

    double xmin = std::numeric_limits<double>::infinity();
    double xmax = -std::numeric_limits<double>::infinity();
    double ymin = std::numeric_limits<double>::infinity();
    double ymax = -std::numeric_limits<double>::infinity();
    for (int idx : candidates) {
        xmin = std::min(xmin, pos_data(0, idx));
        xmax = std::max(xmax, pos_data(0, idx));
        ymin = std::min(ymin, pos_data(1, idx));
        ymax = std::max(ymax, pos_data(1, idx));
    }
    const double width = std::max(xmax - xmin, 1.0);
    const double height = std::max(ymax - ymin, 1.0);
    int nx = std::max(1, static_cast<int>(std::round(std::sqrt(basis_sample_size * width / height))));
    int ny = std::max(1, static_cast<int>(std::ceil(static_cast<double>(basis_sample_size) / nx)));
    std::vector<std::vector<int>> bins(static_cast<size_t>(nx * ny));
    for (int idx : candidates) {
        int bx = std::min(nx - 1, std::max(0, static_cast<int>(((pos_data(0, idx) - xmin) / width) * nx)));
        int by = std::min(ny - 1, std::max(0, static_cast<int>(((pos_data(1, idx) - ymin) / height) * ny)));
        bins[static_cast<size_t>(by) * static_cast<size_t>(nx) + static_cast<size_t>(bx)].push_back(idx);
    }
    int occupied = 0;
    for (const auto& bin : bins) if (!bin.empty()) occupied++;
    int per_bin = std::max(1, static_cast<int>(std::ceil(static_cast<double>(basis_sample_size) / std::max(1, occupied))));
    std::vector<int> selected;
    selected.reserve(basis_sample_size + per_bin);
    for (const auto& bin : bins) {
        if (bin.empty()) continue;
        if (static_cast<int>(bin.size()) <= per_bin) {
            selected.insert(selected.end(), bin.begin(), bin.end());
        } else {
            for (int i = 0; i < per_bin; ++i) {
                int idx = static_cast<int>(std::round(static_cast<double>(i) * (bin.size() - 1) / std::max(1, per_bin - 1)));
                selected.push_back(bin[idx]);
            }
        }
    }
    if (static_cast<int>(selected.size()) <= basis_sample_size) return selected;
    std::vector<int> downsampled;
    downsampled.reserve(basis_sample_size);
    for (int i = 0; i < basis_sample_size; ++i) {
        int idx = static_cast<int>(std::round(static_cast<double>(i) * (selected.size() - 1) / std::max(1, basis_sample_size - 1)));
        downsampled.push_back(selected[idx]);
    }
    return downsampled;
}

NcvInterpolationModel fit_ncv_interpolation_model(
    const Eigen::MatrixXf& anchor_vecs,
    const Eigen::MatrixXd& anchor_emb,
    int n_pca_dims,
    int graph_k
) {
    const int n_components = static_cast<int>(anchor_vecs.rows());
    const int n_anchors = static_cast<int>(anchor_vecs.cols());
    // Thin U of (n_components × n_anchors) has min(rows, cols) columns, so the
    // truncation must be clamped to the anchor count as well: tiny datasets can
    // have fewer anchors than n_pca_dims, and leftCols() would then read past
    // the end of matrixU() (assertion in Debug, out-of-bounds read in Release).
    const int n_pca = std::min({n_pca_dims, n_components, n_anchors});
    NcvInterpolationModel model;
    model.sample_mean = anchor_vecs.rowwise().mean();
    Eigen::MatrixXf centered = anchor_vecs.colwise() - model.sample_mean;
    Eigen::BDCSVD<Eigen::MatrixXf> svd(centered, Eigen::ComputeThinU);
    model.pca_basis = svd.matrixU().leftCols(n_pca);
    model.anchor_pca = (model.pca_basis.transpose() * centered).cast<double>();
    model.anchor_emb = anchor_emb;
    model.interp_k = std::max(1, graph_k);
    return model;
}

Eigen::MatrixXd interpolate_ncv_embedding(
    const NcvInterpolationModel& model,
    const Eigen::MatrixXf& query_vecs
) {
    const int n_query = static_cast<int>(query_vecs.cols());
    if (n_query == 0) return Eigen::MatrixXd(3, 0);
    Eigen::MatrixXd query_pca =
        (model.pca_basis.transpose() * (query_vecs.colwise() - model.sample_mean)).cast<double>();
    int k_interp = std::min(model.interp_k, static_cast<int>(model.anchor_pca.cols()) - 1);
    if (k_interp < 1) k_interp = 1;
    auto knn = knn_parallel(model.anchor_pca, query_pca, k_interp);
    Eigen::MatrixXd emb(3, n_query);
    constexpr double dist_offset = 1e-10;
    parallel_for_static(0, n_query, [&](int i) {
        double w_sum = 0.0;
        Eigen::Vector3d weighted = Eigen::Vector3d::Zero();
        const int* nn_indices = knn.idx_row(i);
        const double* nn_distances = knn.dist_row(i);
        for (int j = 0; j < knn.k; ++j) {
            double w = 1.0 / (nn_distances[j] + dist_offset);
            weighted += w * model.anchor_emb.col(nn_indices[j]);
            w_sum += w;
        }
        emb.col(i) = weighted / w_sum;
    });
    return emb;
}

NcvReportEmbedding compute_ncv_embedding_streaming(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int n_genes,
    const std::vector<double>& confidence,
    int k_neighbors,
    int basis_sample_size,
    int sample_size,
    int seed,
    int n_pca_dims,
    int graph_k,
    bool include_report_umap,
    const NcvProjectedModel* precomputed_model
) {
    const int n_mols = static_cast<int>(genes.size());
    NcvReportEmbedding result;
    result.colors.assign(n_mols, NCV_FALLBACK_COLOR);
    if (n_mols == 0 || n_genes <= 0) return result;

    if (precomputed_model &&
        precomputed_model->basis.spatial_k == k_neighbors &&
        static_cast<int>(precomputed_model->basis.basis_ids.size()) > 1) {
        std::vector<double> basis_conf(precomputed_model->basis.basis_ids.size());
        for (size_t i = 0; i < precomputed_model->basis.basis_ids.size(); ++i) {
            basis_conf[i] = confidence[precomputed_model->basis.basis_ids[i]];
        }
        fill_colors_from_projected_vectors(
            result,
            precomputed_model->basis.basis_vecs,
            precomputed_model->basis.basis_ids,
            basis_conf,
            precomputed_model->mol_vecs.cols() == n_mols ? &precomputed_model->mol_vecs : nullptr,
            pos_data,
            genes,
            precomputed_model->basis.gene_emb_t,
            n_genes,
            k_neighbors,
            precomputed_model->basis.distance_floor,
            sample_size,
            seed,
            n_pca_dims,
            graph_k,
            include_report_umap
        );
        return result;
    }
    if (precomputed_model &&
        precomputed_model->basis.spatial_k != k_neighbors &&
        static_cast<int>(precomputed_model->basis.basis_ids.size()) > 1) {
        spdlog::info(
            "Ignoring precomputed NCV model because spatial neighborhood k differs "
            "(model={}, requested={}).",
            precomputed_model->basis.spatial_k, k_neighbors
        );
    }

    auto basis_model = fit_ncv_basis_model(
        pos_data, genes, n_genes, confidence, k_neighbors, basis_sample_size, 20
    );
    if (static_cast<int>(basis_model.basis_ids.size()) <= 1) return result;

    std::vector<double> basis_conf(basis_model.basis_ids.size());
    for (size_t i = 0; i < basis_model.basis_ids.size(); ++i) basis_conf[i] = confidence[basis_model.basis_ids[i]];
    fill_colors_from_projected_vectors(
        result,
        basis_model.basis_vecs,
        basis_model.basis_ids,
        basis_conf,
        nullptr,
        pos_data,
        genes,
        basis_model.gene_emb_t,
        n_genes,
        k_neighbors,
        basis_model.distance_floor,
        sample_size,
        seed,
        n_pca_dims,
        graph_k,
        include_report_umap
    );

    return result;
}

} // namespace

std::vector<std::string> gene_composition_color_embedding(
    const Eigen::MatrixXf& mol_vecs,
    const std::vector<double>& confidence,
    int sample_size,
    int seed,
    int n_pca_dims,
    int graph_k
) {
    return compute_ncv_embedding(mol_vecs, confidence, sample_size, seed, n_pca_dims, graph_k, false).colors;
}

NcvReportEmbedding gene_composition_report_embedding(
    const Eigen::MatrixXf& mol_vecs,
    const std::vector<double>& confidence,
    int sample_size,
    int seed,
    int n_pca_dims,
    int graph_k
) {
    return compute_ncv_embedding(mol_vecs, confidence, sample_size, seed, n_pca_dims, graph_k, true);
}

NcvBasisModel fit_ncv_basis_model(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int n_genes,
    const std::vector<double>& confidence,
    int k_neighbors,
    int basis_sample_size,
    int n_components
) {
    NcvBasisModel model;
    model.spatial_k = k_neighbors;
    model.basis_ids = select_basis_anchor_ids(pos_data, confidence, basis_sample_size);
    if (static_cast<int>(model.basis_ids.size()) <= 1) return model;

    spdlog::info("NCV basis anchors: selected {} molecules for basis learning.", model.basis_ids.size());
    model.distance_floor = neighborhood_distance_floor(pos_data);
    auto anchor_cm = neighborhood_count_matrix_subset(
        pos_data, genes, model.basis_ids, k_neighbors, n_genes, nullptr, true, true, model.distance_floor
    );
    for (int k = 0; k < anchor_cm.outerSize(); ++k) {
        for (Eigen::SparseMatrix<float>::InnerIterator it(anchor_cm, k); it; ++it) {
            it.valueRef() = static_cast<float>(std::log(it.value() * 10000.0f + 1e-5f));
        }
    }

    model.gene_emb_t = estimate_gene_vectors(anchor_cm, genes, n_components, "ri", false);
    model.basis_vecs = project_gene_vectors(model.gene_emb_t, anchor_cm);
    return model;
}

NcvProjectedModel fit_ncv_projected_model(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int n_genes,
    const std::vector<double>& confidence,
    int k_neighbors,
    int basis_sample_size,
    int n_components,
    bool include_full_projection
) {
    NcvProjectedModel model;
    model.basis = fit_ncv_basis_model(
        pos_data, genes, n_genes, confidence, k_neighbors, basis_sample_size, n_components
    );
    if (!include_full_projection || model.basis.gene_emb_t.cols() == 0) return model;

    model.mol_vecs = project_neighborhood_vectors(
        pos_data, genes, k_neighbors, model.basis.gene_emb_t, n_genes,
        nullptr, nullptr, true, true, model.basis.distance_floor, true
    );
    return model;
}

std::vector<std::string> gene_composition_color_embedding_streaming(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int n_genes,
    const std::vector<double>& confidence,
    int k_neighbors,
    int basis_sample_size,
    int sample_size,
    int seed,
    int n_pca_dims,
    int graph_k,
    const NcvProjectedModel* precomputed_model
) {
    return compute_ncv_embedding_streaming(
        pos_data, genes, n_genes, confidence, k_neighbors,
        basis_sample_size, sample_size, seed, n_pca_dims, graph_k, false, precomputed_model
    ).colors;
}

NcvReportEmbedding gene_composition_report_embedding_streaming(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& genes,
    int n_genes,
    const std::vector<double>& confidence,
    int k_neighbors,
    int basis_sample_size,
    int sample_size,
    int seed,
    int n_pca_dims,
    int graph_k,
    const NcvProjectedModel* precomputed_model
) {
    return compute_ncv_embedding_streaming(
        pos_data, genes, n_genes, confidence, k_neighbors,
        basis_sample_size, sample_size, seed, n_pca_dims, graph_k, true, precomputed_model
    );
}

} // namespace baysor
