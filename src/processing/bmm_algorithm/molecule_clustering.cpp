#include "baysor/processing/bmm_algorithm/molecule_clustering.h"
#include "baysor/utils/general.h"
#include "baysor/utils/thread_pool.h"

#include <Eigen/Dense>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

namespace baysor {

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

    // Molecules grouped by gene in increasing molecule order (CSR), so that the
    // M-step runs in parallel over genes with the sums of a sequential pass.
    std::vector<int> gene_mol_offsets(static_cast<size_t>(n_genes) + 1, 0);
    for (int g1b : genes) {
        if (g1b >= 1) ++gene_mol_offsets[static_cast<size_t>(g1b)];
    }
    std::partial_sum(gene_mol_offsets.begin(), gene_mol_offsets.end(), gene_mol_offsets.begin());
    std::vector<int> gene_mols(static_cast<size_t>(gene_mol_offsets.back()));
    {
        std::vector<int> next(gene_mol_offsets.begin(), gene_mol_offsets.end() - 1);
        for (int i = 0; i < n_mols; ++i) {
            const int g0 = genes[i] - 1;
            if (g0 >= 0) gene_mols[static_cast<size_t>(next[g0]++)] = i;
        }
    }
    // exprs(k, g) = sum over the molecules of gene g of confidence * probs(k, i)
    auto accumulate_exprs = [&]() {
        parallel_for(0, n_genes, 16, [&](int g) {
            for (int k = 0; k < n_clusters; ++k) exprs(k, g) = 0.0;
            for (int j = gene_mol_offsets[g]; j < gene_mol_offsets[g + 1]; ++j) {
                const int i = gene_mols[static_cast<size_t>(j)];
                const double conf = confidence[i];
                for (int k = 0; k < n_clusters; ++k) {
                    exprs(k, g) += conf * probs(k, i);
                }
            }
        });
    };

    Eigen::MatrixXd prev_probs(n_clusters, n_mols);

    std::vector<double> max_diffs;
    std::vector<double> change_fracs;
    max_diffs.reserve(max_iters);
    change_fracs.reserve(max_iters);

    // Convergence statistics per E-step chunk: the boundaries are fixed and
    // the reductions (a maximum and a count) exact, so they do not depend on
    // the thread count.
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
        probs.swap(prev_probs);  // the E-step overwrites every entry of probs

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
        accumulate_exprs();
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
            spdlog::info("  Clustering iter {:4d}: max_diff={:.4f}, change_frac={:.4f}", iter + 1, max_diff, change_fracs.back());
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
                if (verbose) spdlog::info("Clustering converged after {} iterations. Max diff: {:.4f}", iter + 1, max_diff);
                break;
            }
        }
    }

    if (verbose && !max_diffs.empty()) {
        spdlog::info("Clustering stopped after {} iterations. Max diff: {:.4f}. Converged: {}",
                     n_iters_done, max_diffs.back(), max_diffs.back() < tol ? "true" : "false");
    }

    // ---- Final M-step without pseudocount ----
    accumulate_exprs();
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

} // namespace baysor
