#include "baysor/processing/data_processing/initialization.h"
#include "baysor/data_loading/data.h"
#include "baysor/processing/data_processing/triangulation.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/processing/utils/utils.h"
#include "baysor/utils/options.h"
#include "baysor/utils/general.h"

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace baysor {

// ============================================================================
// Internal helpers
// ============================================================================

static double quantile_vec(std::vector<double> v, double p) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    double idx = p * (v.size() - 1);
    int lo = static_cast<int>(std::floor(idx));
    int hi = static_cast<int>(std::ceil(idx));
    if (lo == hi) return v[lo];
    double frac = idx - lo;
    return v[lo] * (1.0 - frac) + v[hi] * frac;
}

// select_ids_uniformly: evenly-spaced sample from molecules sorted by coordinate-sum.
// confidence_threshold=0.25 for initialization (much lower than for colors).
static std::vector<int> select_ids_uniformly(
    const Eigen::MatrixXd& pos_data,
    int n,
    const std::vector<double>* confidences,
    double confidence_threshold = 0.25
) {
    // Matches Julia: `if n <= 1 error("n must be > 1")` in select_ids_uniformly
    // (initialization.jl). Without this guard the evenly-spaced index step below
    // divides by (n - 1) == 0 and the NaN result indexes out of bounds.
    if (n <= 1) {
        throw std::runtime_error("n must be > 1");
    }

    int total = static_cast<int>(pos_data.cols());

    // Collect high-confidence molecule indices
    std::vector<int> high_conf_ids;
    high_conf_ids.reserve(total);
    if (confidences != nullptr && !confidences->empty()) {
        for (int i = 0; i < total; ++i) {
            if ((*confidences)[i] >= confidence_threshold) high_conf_ids.push_back(i);
        }
    } else {
        high_conf_ids.resize(total);
        std::iota(high_conf_ids.begin(), high_conf_ids.end(), 0);
    }

    if (static_cast<int>(high_conf_ids.size()) < n) {
        spdlog::warn("n={} > high_conf molecules ({}). Using all high-conf molecules.", // GCOVR_EXCL_LINE: dead GCC block; the call executes and is counted on the following line(s)
                     n, high_conf_ids.size());
        n = static_cast<int>(high_conf_ids.size());
    }
    if (n <= 0) return {};

    // Only one high-confidence molecule survived the clamp above (a smaller
    // request already errored above). Julia returns that single id the same way
    // when length(high_conf_ids) < n; take it directly so the evenly-spaced
    // step below never divides by (n - 1) == 0.
    if (n == 1) return {high_conf_ids[0]};

    // Sort high_conf_ids by sum of coordinates
    std::vector<std::pair<double,int>> sum_ids(high_conf_ids.size());
    for (int k = 0; k < static_cast<int>(high_conf_ids.size()); ++k) {
        sum_ids[k] = {pos_data.col(high_conf_ids[k]).sum(), high_conf_ids[k]};
    }
    std::sort(sum_ids.begin(), sum_ids.end());

    int hc = static_cast<int>(sum_ids.size());
    // Pick n evenly-spaced indices (matches Julia's unique(round.(range(1, len, length=n))))
    std::vector<int> result;
    result.reserve(n);
    for (int i = 0; i < n; ++i) {
        int idx = static_cast<int>(std::round(static_cast<double>(i) * (hc - 1) / (n - 1)));
        result.push_back(sum_ids[idx].second);
    }
    return result;
}

// ============================================================================
// build_molecule_graph
// ============================================================================

AdjList build_molecule_graph(
    const MoleculeData& data,
    bool filter,
    bool use_local_gene_similarities,
    AdjacencyType type,
    int composition_neighborhood,
    int n_gene_pcs
) {
    Eigen::MatrixXd pos = data.position_matrix();
    int n = static_cast<int>(pos.cols());

    auto adj_result = adjacency_list(pos, filter, /*n_mads=*/2.0, /*k_adj=*/5, type);

    int n_edges = static_cast<int>(adj_result.edge_src.size());
    if (n_edges == 0) {
        AdjList adj; // GCOVR_EXCL_LINE: dead GCC block; construction counted on the following line
        adj.indptr.assign(n + 1, 0);
        return adj;
    }

    double min_edge_length = quantile_vec(adj_result.edge_dists, 0.3);

    std::vector<double> edge_weights(n_edges);
    for (int i = 0; i < n_edges; ++i) {
        edge_weights[i] = min_edge_length / std::max(adj_result.edge_dists[i], min_edge_length);
    }

    return AdjList::from_edge_list(
        adj_result.edge_src.data(),
        adj_result.edge_dst.data(),
        edge_weights.data(),
        n_edges,
        n
    );
}

// ============================================================================
// cell_centers_uniformly<N>
// ============================================================================

template<int N>
InitialParams<N> cell_centers_uniformly(
    const Eigen::MatrixXd& pos_data,
    int n_clusters,
    const std::vector<double>* confidences,
    double scale
) {
    using Mat = Eigen::Matrix<double, N, N>;

    int n_mols = static_cast<int>(pos_data.cols());
    n_clusters = std::min(n_clusters, n_mols);

    // Select n_clusters initial centers evenly-spaced in coordinate-sum order.
    // Julia (cell_centers_uniformly) passes the clamped value straight into
    // select_ids_uniformly, which errors for n <= 1; there is no lower clamp.
    auto center_ids = select_ids_uniformly(pos_data, n_clusters, confidences);
    n_clusters = static_cast<int>(center_ids.size());
    if (n_clusters == 0) {
        // Degenerate: return empty
        InitialParams<N> result; // GCOVR_EXCL_LINE: dead GCC block; statement counted on the following line
        result.centers = Eigen::MatrixXd::Zero(0, N);
        result.assignment.assign(n_mols, 1);
        return result;
    }

    // Build center matrix: N x n_clusters
    Eigen::MatrixXd center_mat(N, n_clusters);
    for (int k = 0; k < n_clusters; ++k) {
        center_mat.col(k) = pos_data.col(center_ids[k]);
    }

    // Assign each molecule to nearest center via KNN (k=1)
    auto knn = knn_parallel(center_mat, pos_data, 1, /*sorted=*/false);
    std::vector<int> cluster_labels(n_mols);
    for (int i = 0; i < n_mols; ++i) {
        cluster_labels[i] = knn.idx_row(i)[0] + 1;  // 1-based
    }

    // Build covariances
    std::vector<Mat> covs(n_clusters);

    if (scale > 0.0) {
        // Isotropic: scale^2 * Identity
        Mat cov = Mat::Identity() * (scale * scale);
        for (auto& c : covs) c = cov;
    } else {
        // Compute sample covariance per cluster
        auto ids_by_cluster = split_ids(cluster_labels, n_clusters, /*drop_zero=*/true);

        // First pass: compute per-cluster sample covariances
        std::vector<bool> valid(n_clusters, false);
        for (int ci = 0; ci < n_clusters; ++ci) {
            const auto& ids = ids_by_cluster[ci];
            int np = static_cast<int>(ids.size());
            if (np <= N) continue;

            // Mean
            Eigen::Matrix<double, N, 1> mu = Eigen::Matrix<double, N, 1>::Zero();
            for (int i : ids) mu += pos_data.col(i).template head<N>();
            mu /= np;

            // Covariance
            covs[ci].setZero();
            for (int i : ids) {
                auto dx = (pos_data.col(i).template head<N>() - mu).eval();
                covs[ci] += dx * dx.transpose();
            }
            covs[ci] /= np;
            valid[ci] = true;
        }

        // Compute median eigenvalues across valid clusters for fallback
        std::vector<double> all_eigs;
        for (int ci = 0; ci < n_clusters; ++ci) {
            if (!valid[ci]) continue;
            Eigen::SelfAdjointEigenSolver<Mat> es(covs[ci]);
            for (int j = 0; j < N; ++j) {
                all_eigs.push_back(es.eigenvalues()[j]);
            }
        }

        Eigen::Matrix<double, N, 1> median_eigvals = Eigen::Matrix<double, N, 1>::Ones();
        if (!all_eigs.empty()) {
            // Use median of all eigenvalues as fallback variance
            std::sort(all_eigs.begin(), all_eigs.end());
            double med_eig = all_eigs[all_eigs.size() / 2];
            for (int j = 0; j < N; ++j) median_eigvals[j] = std::max(med_eig, 1.0);
        }
        Mat fallback_cov = median_eigvals.asDiagonal();

        // Assign fallback to small clusters
        for (int ci = 0; ci < n_clusters; ++ci) {
            if (!valid[ci]) covs[ci] = fallback_cov;
        }
    }

    InitialParams<N> result;
    result.centers = center_mat.transpose();  // n_clusters x N
    result.covs    = std::move(covs);
    result.assignment = std::move(cluster_labels);
    return result;
}

// ============================================================================
// initialize_bmm_data<N>
// ============================================================================

template<int N>
BmmData<N> initialize_bmm_data(
    const MoleculeData& mol_data,
    const AdjList& adj_list,
    int n_cells_init,
    double scale,
    const std::string& scale_std_str,
    double prior_seg_confidence,
    int min_molecules_per_cell,
    bool verbose
) {
    int n_mols  = mol_data.n_molecules();
    int n_genes = mol_data.n_genes();

    // Parse scale_std
    double scale_std = parse_scale_std(scale_std_str, scale);

    // Shape prior: same std for all spatial dimensions
    Eigen::Matrix<double, N, 1> std_vals, std_stds;
    for (int d = 0; d < N; ++d) {
        std_vals[d] = scale;
        std_stds[d] = scale_std;
    }
    ShapePrior<N> shape_prior_obj{std_vals, std_stds, min_molecules_per_cell};

    if (verbose) {
        spdlog::info("Initializing algorithm. Scale: {:.2f}, scale_std: {:.2f}, "
                     "initial #components: {}, #molecules: {}.",
                     scale, scale_std, n_cells_init, n_mols);
    }

    // Get confidence vector (may be empty before append_confidence is called)
    const std::vector<double>* conf_ptr =
        mol_data.confidence.empty() ? nullptr : &mol_data.confidence;

    // Initial cell placement
    Eigen::MatrixXd pos = mol_data.position_matrix();  // N x n_mols
    auto init = cell_centers_uniformly<N>(pos, n_cells_init, conf_ptr, scale);

    int actual_n_cells = static_cast<int>(init.covs.size());

    // noise_position_density = pdf of N(0, scale^2 * I) at (3*scale, ..., 3*scale)
    // = 1 / ((2π*scale^2)^(N/2)) * exp(-0.5 * N * 9)
    double noise_position_density =
        std::exp(-0.5 * N * 9.0)
        / std::pow(2.0 * kPi * scale * scale, N / 2.0);

    // Build Component objects
    // Each component: MvNormal(center_i, cov_i) + CategoricalSmoothed with all-ones counts
    // shared ShapePrior (non-owning pointer, so we keep it alive in BmmData via a side structure)
    // Each Component owns its ShapePrior by value (std::optional<ShapePrior<N>>),
    // so there are no pointer lifetime concerns when BmmData is moved.
    std::vector<Component<N>> components;
    components.reserve(actual_n_cells);

    auto n_samples_per_cell = count_array(init.assignment, actual_n_cells, /*drop_zero=*/true);

    for (int ci = 0; ci < actual_n_cells; ++ci) {
        Eigen::Matrix<double, N, 1> center = init.centers.row(ci).transpose();
        MvNormal<N> pos_params(center, init.covs[ci]);

        // CategoricalSmoothed with all-ones initial counts (uniform gene prior)
        CategoricalSmoothed component_params(n_genes, 1.0);
        // Set initial counts to 1.0 per gene (matches Julia's dense all-ones gene prior)
        component_params.set_uniform_counts(1.0f);

        Component<N> comp(pos_params, component_params, shape_prior_obj, /*guid=*/ci + 1);
        comp.n_samples = n_samples_per_cell.empty() ? 0 : n_samples_per_cell[ci];
        components.push_back(std::move(comp));
    }

    // Build BmmData
    BmmData<N> bm_data;
    bm_data.position_data = mol_data.position_matrix();
    bm_data.composition_data = mol_data.gene;  // 1-based gene IDs → pass as 0-based below

    // Convert gene IDs from 1-based to 0-based for internal use
    // CategoricalSmoothed uses 0-based gene IDs
    for (int& g : bm_data.composition_data) {
        g = (g > 0) ? g - 1 : -1;
    }

    bm_data.confidence           = mol_data.confidence.empty()
                                    ? std::vector<double>(n_mols, 0.95)
                                    : mol_data.confidence;
    bm_data.cluster_per_molecule = mol_data.cluster;
    bm_data.nuclei_prob_per_molecule = mol_data.nuclei_probs;

    bm_data.adj_list   = adj_list;
    bm_data.components = std::move(components);
    bm_data.assignment = init.assignment;
    bm_data.max_component_guid = actual_n_cells;

    bm_data.noise_position_density = noise_position_density;
    bm_data.noise_density = 0.0;

    bm_data.prior_seg_confidence = prior_seg_confidence;
    bm_data.mrf_strength         = 0.1;
    bm_data.use_gene_smoothing   = true;
    bm_data.min_nuclei_frac      = 0.1;
    bm_data.cluster_penalty_mult = 0.25;
    bm_data.real_edge_weight     = 1.0;

    // Prior segmentation
    if (!mol_data.prior_segmentation.empty()) {
        bm_data.segment_per_molecule = mol_data.prior_segmentation;
        // n_molecules_per_segment: count molecules per segment ID
        int max_seg = *std::max_element(
            mol_data.prior_segmentation.begin(), mol_data.prior_segmentation.end());
        bm_data.n_molecules_per_segment = count_array(
            mol_data.prior_segmentation, max_seg, /*drop_zero=*/true);
        bm_data.update_n_mols_per_segment();
    }

    return bm_data;
}

// Explicit instantiations
template InitialParams<2> cell_centers_uniformly<2>(
    const Eigen::MatrixXd&, int, const std::vector<double>*, double);
template InitialParams<3> cell_centers_uniformly<3>(
    const Eigen::MatrixXd&, int, const std::vector<double>*, double);
template BmmData<2> initialize_bmm_data<2>(
    const MoleculeData&, const AdjList&, int, double,
    const std::string&, double, int, bool);
template BmmData<3> initialize_bmm_data<3>(
    const MoleculeData&, const AdjList&, int, double,
    const std::string&, double, int, bool);

} // namespace baysor
