#include "baysor/processing/bmm_algorithm/molecule_clustering.h"
#include "baysor/reporting/color_utils.h"
#include "baysor/utils/general.h"
#include "baysor/utils/thread_pool.h"

#include <Eigen/Dense>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>
#include <stdexcept>

namespace baysor {

// ============================================================================
// fast_ica — symmetric deflation FastICA (tanh nonlinearity)
//
// Input:  X  — data matrix (n_features × n_samples)
//         n_components — number of independent components
// Returns: unmixing matrix W (n_features × n_components)
//          whose columns correspond to the components used by Julia's
//          MultivariateStats.ICA fit via ica_fit.W.
//
// Mirrors Julia's MultivariateStats.ICA usage in cluster_molecules_on_mrf,
// including its argument check `k <= min(m, n) || error(...)`: throws
// std::invalid_argument instead of silently clamping, so the caller's
// try/catch falls back exactly like Julia's wrapper does.
// ============================================================================
static Eigen::MatrixXd fast_ica(
    const Eigen::MatrixXd& X,
    int n_components,
    int max_iter = 1000,
    double tol   = 1e-5,
    unsigned int seed = 42
) {
    int n_features = static_cast<int>(X.rows());
    int n_samples  = static_cast<int>(X.cols());
    // Julia's fit(ICA, X, k) rejects k > min(m, n) with an error that
    // cluster_molecules_on_mrf's wrapper catches to fall back to hash/random
    // init. Clamping here instead would return fewer columns than the caller
    // expects and make it index W out of bounds.
    if (n_components > std::min(n_features, n_samples))
        throw std::invalid_argument("k must not exceed min(m, n).");

    // 1. Center rows to match Julia's preprocess_mean/centralize path.
    Eigen::VectorXd mean_vec = X.rowwise().mean();
    Eigen::MatrixXd Xc = X.colwise() - mean_vec;

    // 2. Whiten with the top-k covariance eigenvectors, matching:
    //    C = Z * Z' / (n - 1); W0 = P * Diagonal(1 ./ sqrt.(v)); Z = W0' * Z
    Eigen::MatrixXd cov =
        (Xc * Xc.transpose()) / static_cast<double>(std::max(1, n_samples - 1));
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(cov);
    Eigen::VectorXd lambdas(n_components);
    Eigen::MatrixXd P(n_features, n_components);
    for (int i = 0; i < n_components; ++i) {
        int src_idx = n_features - 1 - i;
        lambdas(i) = eig.eigenvalues()(src_idx);
        P.col(i) = eig.eigenvectors().col(src_idx);
    }

    // Clamp tiny eigenvalues to avoid division by ~0
    for (int i = 0; i < n_components; ++i)
        if (lambdas(i) < 1e-10) lambdas(i) = 1e-10;

    Eigen::MatrixXd W0 =
        P * lambdas.cwiseSqrt().cwiseInverse().asDiagonal(); // m × k
    Eigen::MatrixXd Z = W0.transpose() * Xc;                // k × n

    std::mt19937 rng(seed);
    std::normal_distribution<double> ndist(0.0, 1.0);
    Eigen::MatrixXd W(n_components, n_components);
    for (int c = 0; c < n_components; ++c) {
        for (int r = 0; r < n_components; ++r) {
            W(r, c) = ndist(rng);
        }
        double norm = W.col(c).norm();
        if (norm > 0.0) W.col(c) /= norm;
    }

    Eigen::MatrixXd W_prev(n_components, n_components);
    Eigen::MatrixXd U(n_samples, n_components);
    Eigen::MatrixXd Y(n_components, n_components);
    Eigen::VectorXd E1(n_components);

    for (int iter = 0; iter < max_iter; ++iter) {
        W_prev = W;

        // U <- W' * X, stored as (n_samples × n_components)
        U.noalias() = Z.transpose() * W;

        // Tanh nonlinearity with a=1 and mean derivative per component.
        for (int c = 0; c < n_components; ++c) {
            double deriv_sum = 0.0;
            for (int i = 0; i < n_samples; ++i) {
                double t = std::tanh(U(i, c));
                U(i, c) = t;
                deriv_sum += 1.0 - t * t;
            }
            E1(c) = deriv_sum / static_cast<double>(n_samples);
        }

        // Y <- E{x g(w'x)}
        Y.noalias() = (Z * U) / static_cast<double>(n_samples);

        // W <- Y - E{g'(w'x)} * W
        for (int c = 0; c < n_components; ++c) {
            W.col(c) = Y.col(c) - E1(c) * W.col(c);
        }

        // Symmetric decorrelation: W <- W * (W'W)^(-1/2)
        Eigen::MatrixXd gram = W.transpose() * W;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> gram_eig(gram);
        Eigen::VectorXd gram_vals = gram_eig.eigenvalues();
        for (int i = 0; i < gram_vals.size(); ++i) {
            if (gram_vals(i) < 1e-10) gram_vals(i) = 1e-10;
        }
        Eigen::MatrixXd invsqrt =
            gram_eig.eigenvectors() *
            gram_vals.cwiseSqrt().cwiseInverse().asDiagonal() *
            gram_eig.eigenvectors().transpose();
        W = W * invsqrt;

        Eigen::MatrixXd prod = W * W_prev.transpose();
        double chg = 0.0;
        for (int i = 0; i < n_components; ++i) {
            chg = std::max(chg, std::abs(std::abs(prod(i, i)) - 1.0));
        }
        if (chg < tol) break;
    }

    // 4. Map the whitened solution back to the original feature space.
    return W0 * W;
}

// ============================================================================
// cluster_molecules_on_mrf — categorical variant
// ============================================================================
//
// Port of Julia's cluster_molecules_on_mrf (molecule_clustering.jl)
// using the CatMixture (categorical gene-expression) mixture model.
//
// Layout conventions (matching Julia):
//   exprs            : n_clusters × n_genes    (row = cluster, col = gene)
//   assignment_probs : n_clusters × n_molecules (row = cluster, col = mol)
//   genes            : 1-based (0 = missing/unknown gene, skip)
//
// Optional exprs_init: pre-computed n_clusters × n_genes expression profiles
// (e.g. from ICA). When nullptr, falls back to deterministic hash perturbation.

ClusteringResult cluster_molecules_on_mrf(
    const std::vector<int>& genes,   // 1-based gene IDs
    const AdjList& adj_list,
    const std::vector<double>& confidence,
    int n_clusters,
    double tol,
    double mrf_weight,
    int max_iters,
    bool verbose,
    const Eigen::MatrixXd* exprs_init
) {
    int n_mols  = static_cast<int>(genes.size());
    if (n_mols == 0 || n_clusters <= 1) return {};

    // Infer n_genes from max gene ID (genes are 1-based)
    int n_genes = 0;
    for (int g : genes) if (g > n_genes) n_genes = g;
    if (n_genes == 0) return {};

    if (max_iters <= 0)
        max_iters = std::max(10000, n_mols / 200);

    constexpr int n_iters_without_update = 20;

    // ------------------------------------------------------------------
    // Step 0: Pre-multiply adj weights by neighbor confidence
    //   adj_weights_conf[i][j] = original_weight[i][j] * confidence[neighbor[j]]
    // Stored as a flat mirror of adj_list (same indptr).
    // ------------------------------------------------------------------
    std::vector<double> adj_w_conf(adj_list.nnz());
    for (int i = 0; i < n_mols; ++i) {
        int start = adj_list.indptr[i];
        int end   = adj_list.indptr[i + 1];
        const int32_t* nb  = adj_list.indices.data() + start;
        const double*  wt  = adj_list.weights.data()  + start;
        double*        out = adj_w_conf.data()         + start;
        for (int j = 0; j < end - start; ++j) {
            out[j] = wt[j] * confidence[nb[j]];
        }
    }

    // ------------------------------------------------------------------
    // Step 1: Initialize expression profiles
    //
    // If exprs_init is provided (e.g. from ICA), use it directly (after
    // normalizing rows to sum to 1).  Otherwise, fall back to Julia's
    // deterministic hash-perturbed gene-frequency initialization
    // (init_cell_type_exprs with init_mod=10000).
    // ------------------------------------------------------------------

    // Global gene frequency (0-based genes) — needed for fallback init
    std::vector<double> gene_freq(n_genes, 0.0);
    int n_valid = 0;
    for (int g1b : genes) {
        int g0 = g1b - 1;
        if (g0 >= 0 && g0 < n_genes) { gene_freq[g0] += 1.0; ++n_valid; }
    }
    if (n_valid > 0) for (double& f : gene_freq) f /= n_valid;

    Eigen::MatrixXd exprs(n_clusters, n_genes);

    if (exprs_init && exprs_init->rows() == n_clusters
                   && exprs_init->cols() == n_genes) {
        // Use provided initialization directly, matching Julia's
        // init_cell_type_exprs(..., cell_type_exprs=...) fast path.
        exprs = *exprs_init;
    } else {
        // Fallback: deterministic hash perturbation of global gene frequency
        // Mirrors Julia's init_cell_type_exprs (init_mod=10000 path)
        constexpr int init_mod = 10000;
        for (int k = 0; k < n_clusters; ++k) {
            double row_sum = 0.0;
            for (int g = 0; g < n_genes; ++g) {
                std::size_t h = std::hash<long long>{}(
                    static_cast<long long>(g + 1) *
                    static_cast<long long>((k + 1) * (k + 1)));
                double noise = static_cast<double>(h % init_mod) / 100000.0;
                exprs(k, g) = gene_freq[g] * (0.95 + noise);
                row_sum += exprs(k, g);
            }
            // Pseudocount smoothing matching Julia: (x + 1) / (sum(x) + 1)
            double norm = row_sum + 1.0;
            for (int g = 0; g < n_genes; ++g)
                exprs(k, g) = (exprs(k, g) + 1.0) / norm;
        }
    }

    // assignment_probs[k][i] = exprs[k][gene[i]], then column-normalize
    // (matches Julia's init_assignment_probs_inner)
    Eigen::MatrixXd probs(n_clusters, n_mols);
    for (int i = 0; i < n_mols; ++i) {
        int g0 = genes[i] - 1;
        double col_sum = 0.0;
        for (int k = 0; k < n_clusters; ++k) {
            double p = (g0 >= 0) ? exprs(k, g0) : (1.0 / n_clusters);
            probs(k, i) = p;
            col_sum += p;
        }
        if (col_sum > 1e-100)
            for (int k = 0; k < n_clusters; ++k) probs(k, i) /= col_sum;
        else
            for (int k = 0; k < n_clusters; ++k) probs(k, i) = 1.0 / n_clusters;
    }

    Eigen::MatrixXd prev_probs(n_clusters, n_mols);

    std::vector<double> max_diffs;
    std::vector<double> change_fracs;
    max_diffs.reserve(max_iters);
    change_fracs.reserve(max_iters);

    // The convergence statistics are computed inside the parallel E-step, per
    // chunk of molecules. Chunk boundaries are fixed and the reductions (a
    // maximum and a count) are exact, so the result does not depend on the
    // thread count and equals the former serial pass.
    constexpr std::int64_t mol_chunk = 512;
    const std::int64_t n_chunks = (n_mols + mol_chunk - 1) / mol_chunk;
    std::vector<double> chunk_max_diff(static_cast<size_t>(n_chunks));
    std::vector<int> chunk_n_changed(static_cast<size_t>(n_chunks));

    // ------------------------------------------------------------------
    // Step 2: EM loop
    // ------------------------------------------------------------------
    int n_iters_done = 0;
    for (int iter = 0; iter < max_iters; ++iter) {
        n_iters_done = iter + 1;
        // prev_probs takes the last iteration's probabilities; the E-step
        // overwrites every entry of probs.
        probs.swap(prev_probs);

        // ---- E-step and convergence statistics (parallel over molecules) ----
        run_parallel_chunks(0, n_mols, mol_chunk, Scheduling::Dynamic,
                            [&](std::int64_t chunk_begin, std::int64_t chunk_end, int) {
            double chunk_max = 0.0;
            int chunk_changed = 0;
            for (int i = static_cast<int>(chunk_begin); i < static_cast<int>(chunk_end); ++i) {
                int  g0    = genes[i] - 1;  // 0-based gene (< 0 if missing)
                int  start = adj_list.indptr[i];
                int  end   = adj_list.indptr[i + 1];
                const int32_t* nb_ids = adj_list.indices.data() + start;
                const double*  nb_wt  = adj_w_conf.data()        + start;
                int  n_nb  = end - start;

                double col_sum = 0.0;
                for (int k = 0; k < n_clusters; ++k) {
                    // MRF term: weighted sum of neighbor probabilities for cluster k
                    double c_d = 0.0;
                    for (int j = 0; j < n_nb; ++j) {
                        double a_p = prev_probs(k, nb_ids[j]);
                        if (a_p > 1e-5) c_d += nb_wt[j] * a_p;
                    }
                    double mrf_prior = std::exp(mrf_weight * c_d);

                    // Expression likelihood (skip if gene unknown)
                    double expr_ll = (g0 >= 0) ? exprs(k, g0) : 1.0;
                    probs(k, i) = expr_ll * mrf_prior;
                    col_sum    += probs(k, i);
                }

                // Normalize column
                if (col_sum > 1e-100) {
                    for (int k = 0; k < n_clusters; ++k) probs(k, i) /= col_sum;
                } else {
                    for (int k = 0; k < n_clusters; ++k) probs(k, i) = 1.0 / n_clusters;
                }

                // Largest confidence-weighted probability change of this molecule
                double conf = confidence[i];
                double mol_max = 0.0;
                for (int k = 0; k < n_clusters; ++k) {
                    double d = std::abs(probs(k, i) - prev_probs(k, i)) * conf;
                    if (d > mol_max) mol_max = d;
                }
                if (mol_max > chunk_max) chunk_max = mol_max;
                if (mol_max > 1e-7) ++chunk_changed;
            }
            const size_t chunk = static_cast<size_t>(chunk_begin / mol_chunk);
            chunk_max_diff[chunk] = chunk_max;
            chunk_n_changed[chunk] = chunk_changed;
        });

        // ---- M-step with pseudocount ----
        exprs.setZero();
        for (int i = 0; i < n_mols; ++i) {
            int g0 = genes[i] - 1;
            if (g0 < 0) continue;
            double conf = confidence[i];
            for (int k = 0; k < n_clusters; ++k) {
                exprs(k, g0) += conf * probs(k, i);
            }
        }
        for (int k = 0; k < n_clusters; ++k) {
            double row_sum = exprs.row(k).sum();
            double norm = row_sum + 1.0;
            for (int g = 0; g < n_genes; ++g) {
                exprs(k, g) = (exprs(k, g) + 1.0) / norm;
            }
        }

        // ---- Convergence check ----
        double max_diff = 0.0;
        int n_changed = 0;
        for (std::int64_t c = 0; c < n_chunks; ++c) {
            if (chunk_max_diff[c] > max_diff) max_diff = chunk_max_diff[c];
            n_changed += chunk_n_changed[c];
        }
        max_diffs.push_back(max_diff);
        change_fracs.push_back(static_cast<double>(n_changed) / n_mols);

        if (verbose && (iter % 100 == 0 || iter < 5)) {
            spdlog::info("  Clustering iter {:4d}: max_diff={:.4f}, change_frac={:.4f}", // GCOVR_EXCL_LINE: gcov exception-cleanup artifact: the call is counted on the following line; this line only runs when an exception unwinds through the statement
                         iter + 1, max_diff, change_fracs.back());
        }

        // Stop if last n_iters_without_update all below tol
        if (iter + 1 > n_iters_without_update) {
            int look_back = std::min(n_iters_without_update + 1,
                                     static_cast<int>(max_diffs.size()));
            double worst = 0.0;
            for (int t = static_cast<int>(max_diffs.size()) - look_back;
                 t < static_cast<int>(max_diffs.size()); ++t) {
                if (max_diffs[t] > worst) worst = max_diffs[t];
            }
            if (worst < tol) {
                if (verbose)
                    spdlog::info("Clustering converged after {} iterations. Max diff: {:.4f}", // GCOVR_EXCL_LINE: gcov exception-cleanup artifact: the call is counted on the following line; this line only runs when an exception unwinds through the statement
                                 iter + 1, max_diff);
                break;
            }
        }
    }

    if (verbose && !max_diffs.empty()) {
        spdlog::info("Clustering stopped after {} iterations. Max diff: {:.4f}. Converged: {}",
                     n_iters_done, max_diffs.back(), max_diffs.back() < tol ? "true" : "false");
    }

    // ---- Final M-step without pseudocount ----
    exprs.setZero();
    for (int i = 0; i < n_mols; ++i) {
        int g0 = genes[i] - 1;
        if (g0 < 0) continue;
        double conf = confidence[i];
        for (int k = 0; k < n_clusters; ++k) {
            exprs(k, g0) += conf * probs(k, i);
        }
    }
    for (int k = 0; k < n_clusters; ++k) {
        double row_sum = exprs.row(k).sum();
        if (row_sum > 0) exprs.row(k) /= row_sum;
    }

    // ---- Hard assignment: 1-based cluster IDs ----
    std::vector<int> assignment(n_mols);
    for (int i = 0; i < n_mols; ++i) {
        int best_k = 0;
        double best_p = probs(0, i);
        for (int k = 1; k < n_clusters; ++k) {
            if (probs(k, i) > best_p) { best_p = probs(k, i); best_k = k; }
        }
        assignment[i] = best_k + 1;  // 1-based
    }

    return ClusteringResult{
        std::move(exprs),
        std::move(assignment),
        std::move(probs),
        std::move(max_diffs),
        std::move(change_fracs)
    };
}

// ============================================================================
// cluster_molecules_ica
//
// Port of Julia's DataFrame-based wrapper for cluster_molecules_on_mrf:
//   1. Compute pairwise gene spatial co-occurrence (already in color_utils)
//   2. Run FastICA on the correlation matrix to get n_clusters gene profiles
//   3. Call core EM with ICA init; fall back to hash init if ICA throws
// ============================================================================
ClusteringResult cluster_molecules_ica(
    const std::vector<int>& genes,
    const AdjList& adj_list,
    const std::vector<double>& confidence,
    int n_clusters,
    double tol,
    double mrf_weight,
    int max_iters,
    bool verbose
) {
    if (n_clusters <= 1) return {};

    // Infer n_genes
    int n_genes = 0;
    for (int g : genes) if (g > n_genes) n_genes = g;

    // 1. Gene spatial correlation matrix (n_genes × n_genes)
    Eigen::MatrixXd cor_mat = pairwise_gene_spatial_cor(genes, confidence, adj_list);
    // cor_mat is indexed 0-based (genes 0..n_genes-1 from 1-based input)

    // 2. FastICA on the correlation matrix → unmixing matrix (n_genes × n_clusters)
    std::unique_ptr<Eigen::MatrixXd> exprs_init_ptr;
    try {
        Eigen::MatrixXd W = fast_ica(cor_mat, n_clusters);
        // Convert mixing matrix to expression profiles:
        //   ct_exprs_init[k][g] = abs(W[g][k]) / sum_g'(abs(W[g'][k]))
        // Matches Julia: (abs.(ica_fit.W) ./ sum(abs.(ica_fit.W), dims=1))'
        Eigen::MatrixXd exprs(n_clusters, n_genes);
        for (int k = 0; k < n_clusters; ++k) {
            double col_sum = W.col(k).cwiseAbs().sum();
            if (col_sum < 1e-10) col_sum = 1.0;
            for (int g = 0; g < n_genes; ++g)
                exprs(k, g) = std::abs(W(g, k)) / col_sum;
        }
        exprs_init_ptr = std::make_unique<Eigen::MatrixXd>(std::move(exprs));
        if (verbose) spdlog::info("ICA initialization succeeded ({} components).", n_clusters);
    } catch (const std::exception& e) {
        // Discard any partially built init so *every* failure path falls back
        // to hash initialization, consistent with the log message (and with
        // Julia, which leaves ct_exprs_init = nothing on any exception).
        exprs_init_ptr.reset();
        spdlog::warn("ICA did not converge ({}), falling back to hash initialization.", e.what());
    } catch (...) {
        exprs_init_ptr.reset(); // GCOVR_EXCL_LINE: this handler is entered only by a non-std exception, and no Baysor or dependency code throws one
        spdlog::warn("ICA failed, falling back to hash initialization."); // GCOVR_EXCL_LINE: reachable only if a non-std exception escapes the ICA try block; same justification as main.cpp's catch-all (no Baysor or dependency exception lacks std::exception)
    }

    // 3. Core EM with ICA init (or nullptr → hash fallback)
    return cluster_molecules_on_mrf(
        genes, adj_list, confidence,
        n_clusters, tol, mrf_weight, max_iters, verbose,
        exprs_init_ptr.get()
    );
}

} // namespace baysor
