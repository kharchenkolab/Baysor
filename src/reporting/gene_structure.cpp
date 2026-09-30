#include "baysor/reporting/color_utils.h"
#include "baysor/processing/data_processing/umap_wrappers.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/utils/thread_pool.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace baysor {

Eigen::MatrixXd pairwise_gene_spatial_cor(
    const std::vector<int>& genes,
    const std::vector<double>& confidence,
    const AdjList& adj_list,
    double confidence_threshold
) {
    int n = static_cast<int>(genes.size());
    int n_genes = *std::max_element(genes.begin(), genes.end());  // 1-based max

    Eigen::MatrixXd cor_mat = Eigen::MatrixXd::Zero(n_genes, n_genes);
    std::vector<int> counts_per_gene(n_genes, 0);
    for (int gi = 0; gi < n; ++gi) {
        if (confidence[gi] < confidence_threshold) continue;
        const int g = genes[gi] - 1;
        if (g >= 0 && g < n_genes) counts_per_gene[g]++;
    }

    std::vector<std::vector<int>> mols_by_gene(n_genes);
    for (int g = 0; g < n_genes; ++g) mols_by_gene[g].reserve(counts_per_gene[g]);
    for (int gi = 0; gi < n; ++gi) {
        if (confidence[gi] < confidence_threshold) continue;
        const int g = genes[gi] - 1;
        if (g >= 0 && g < n_genes) mols_by_gene[g].push_back(gi);
    }

    // Build the correlation matrix row-by-row. This avoids keeping one dense
    // n_genes x n_genes accumulation matrix per worker.
    std::vector<Eigen::VectorXd> row_buf(thread_pool_size(), Eigen::VectorXd::Zero(n_genes));
    parallel_for(0, n_genes, 1, [&](int g2, int ti) {
        auto& row = row_buf[ti];
        row.setZero();
        for (int gi : mols_by_gene[g2]) {
            const int nc = adj_list.neighbor_count(gi);
            const int32_t* nb_ids = adj_list.neighbor_ids(gi);
            const double* nb_wts = adj_list.neighbor_weights(gi);

            for (int ai = 0; ai < nc; ++ai) {
                const int nb = nb_ids[ai];
                if (confidence[nb] < confidence_threshold) continue;
                const int g1 = genes[nb] - 1;
                if (g1 < 0 || g1 >= n_genes) continue;
                row[g1] += nb_wts[ai];
            }
        }
        cor_mat.row(g2) = row.transpose();
    });

    std::vector<double> sum_weight(n_genes, 0.0);
    parallel_for_static(0, n_genes, [&](int g) {
        sum_weight[g] = cor_mat.row(g).sum() + cor_mat.col(g).sum();
    });

    for (int ci = 0; ci < n_genes; ++ci) {
        for (int ri = 0; ri < n_genes; ++ri) {
            double denom = std::sqrt(sum_weight[ri] * sum_weight[ci]);
            cor_mat(ri, ci) /= std::max(denom, 0.1);
        }
    }

    return cor_mat;
}

GeneStructureEmbedding estimate_gene_structure_embedding(
    const std::vector<int>& genes,
    const std::vector<std::string>& gene_names,
    const std::vector<double>& confidence,
    const AdjList& adj_list,
    int seed
) {
    int n_genes = static_cast<int>(gene_names.size());

    Eigen::MatrixXd cor_mat = pairwise_gene_spatial_cor(genes, confidence, adj_list);

    // Convert correlation to distance, matching Julia:
    //   p_dists = 1 - clamp(cor, min_cor, max_cor) / max_cor
    std::vector<double> flat(cor_mat.data(), cor_mat.data() + cor_mat.size());
    std::vector<double> pos_vals;
    for (double v : flat) if (v > 0) pos_vals.push_back(v);
    std::sort(pos_vals.begin(), pos_vals.end());

    double min_cor = 0.0, max_cor = 1.0;
    if (!pos_vals.empty()) {
        int idx001 = static_cast<int>(0.001 * (pos_vals.size() - 1));
        int idx99  = static_cast<int>(0.99  * (flat.size()    - 1));
        std::nth_element(flat.begin(), flat.begin() + idx99, flat.end());
        min_cor = pos_vals[idx001] / 5.0;
        max_cor = std::max(flat[idx99], 1e-10);
    }

    Eigen::MatrixXd dist_mat(n_genes, n_genes);
    for (int ci = 0; ci < n_genes; ++ci) {
        for (int ri = 0; ri < n_genes; ++ri) {
            double c = std::max(std::min(cor_mat(ri, ci), max_cor), min_cor);
            dist_mat(ri, ci) = 1.0 - c / max_cor;
        }
    }
    dist_mat.diagonal().setZero();

    int k = std::max(std::min(15, n_genes / 2), 2);
    Eigen::MatrixXd emb = umap_embed_precomputed(dist_mat, 2, k, 5000, seed);

    // Molecule counts per gene (for marker size, matching Julia's log(count_array)).
    std::vector<int> counts(n_genes, 0);
    for (int g : genes) {
        int gi = g - 1;
        if (gi >= 0 && gi < n_genes) counts[gi]++;
    }

    GeneStructureEmbedding result;
    result.x.resize(n_genes);
    result.y.resize(n_genes);
    result.gene_names = gene_names;
    result.marker_sizes.resize(n_genes);
    for (int i = 0; i < n_genes; ++i) {
        result.x[i]            = emb(0, i);
        result.y[i]            = emb(1, i);
        result.marker_sizes[i] = std::log(static_cast<double>(std::max(counts[i], 1)));
    }
    return result;
}

} // namespace baysor
