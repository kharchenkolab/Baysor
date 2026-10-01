#include "baysor/processing/data_processing/noise_estimation.h"
#include "baysor/processing/data_processing/initialization.h"
#include "baysor/processing/utils/utils.h"
#include "baysor/data_loading/data.h"
#include "baysor/utils/general.h"
#include "baysor/utils/thread_pool.h"

#include <spdlog/spdlog.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>

namespace baysor {

// ============================================================================
// Helpers
// ============================================================================

static constexpr double PI = 3.14159265358979323846;

static inline double normal_pdf(double x, double mu, double sigma) {
    if (sigma <= 0.0) return 0.0;
    double z = (x - mu) / sigma;
    return std::exp(-0.5 * z * z) / (sigma * std::sqrt(2.0 * PI));
}

static double quantile_sorted(const std::vector<double>& sorted_v, double p) {
    int n = static_cast<int>(sorted_v.size());
    if (n == 0) return 0.0;
    if (n == 1) return sorted_v[0];
    double idx = p * (n - 1);
    int lo = static_cast<int>(std::floor(idx));
    int hi = std::min(lo + 1, n - 1);
    double frac = idx - lo;
    return sorted_v[lo] * (1.0 - frac) + sorted_v[hi] * frac;
}

static double quantile_vec(const std::vector<double>& v, double p) {
    std::vector<double> s = v;
    std::sort(s.begin(), s.end());
    return quantile_sorted(s, p);
}

// ============================================================================
// E-step: MRF-regularized expectation
// ============================================================================

static void expect_noise_probabilities(
    Eigen::MatrixXd& assignment_probs,       // n x 2
    double mu1, double sigma1,
    double mu2, double sigma2,
    const std::vector<double>& edge_lengths,
    const AdjList& adj_list,
    const std::vector<int>& updating_ids,
    const std::vector<double>* min_confidence
) {
    int n = static_cast<int>(edge_lengths.size());

    // Precompute component densities in parallel.
    std::vector<double> pdf1(n), pdf2(n);
    parallel_for_static(0, n, [&](int i) {
        pdf1[i] = normal_pdf(edge_lengths[i], mu1, sigma1);
        pdf2[i] = normal_pdf(edge_lengths[i], mu2, sigma2);
    });

    // Component sizes. Deterministic reduction: sequential accumulation with
    // 1 thread, fixed buckets merged in index order otherwise.
    double n1 = parallel_reduce<double>(0, n, /*bucket_size=*/1024, 0.0,
        [&](std::int64_t b, std::int64_t e, double& acc) {
            for (std::int64_t i = b; i < e; ++i) acc += assignment_probs(i, 0);
        },
        std::plus<double>());
    double n2 = n - n1;

    int m = static_cast<int>(updating_ids.size());
    for (int ui = 0; ui < m; ++ui) {
        int i = updating_ids[ui];
        double c_d1 = 0.0, c_d2 = 0.0;
        int nc = adj_list.neighbor_count(i);
        const int32_t* nb_ids = adj_list.neighbor_ids(i);
        const double* nb_wts = adj_list.neighbor_weights(i);

        for (int j = 0; j < nc; ++j) {
            // TODO(parity): This mirrors Julia's in-place sweep semantics exactly.
            // Later molecules in the same pass see already-updated neighbors.
            // Revisit after parity work: a Jacobi-style or otherwise explicit
            // EM update would be easier to reason about for stability.
            double ap = assignment_probs(nb_ids[j], 0);
            c_d1 += nb_wts[j] * ap;
            c_d2 += nb_wts[j] * (1.0 - ap);
        }

        double d1_val = n1 * std::exp(c_d1) * pdf1[i];
        double d2_val = n2 * std::exp(c_d2) * pdf2[i];
        double denom = std::max(d1_val + d2_val, 1e-20);
        double p1 = d1_val / denom;

        if (min_confidence) {
            p1 = (*min_confidence)[i] + p1 * (1.0 - (*min_confidence)[i]);
        }

        assignment_probs(i, 0) = p1;
        assignment_probs(i, 1) = 1.0 - p1;
    }
}

// ============================================================================
// fit_noise_probabilities
// ============================================================================

NoiseFitResult fit_noise_probabilities(
    const std::vector<double>& edge_lengths,
    const AdjList& adj_list,
    const std::vector<double>* min_confidence,
    int max_iters,
    double tol,
    bool verbose
) {
    int n = static_cast<int>(edge_lengths.size());

    // Initialization: two-component Gaussian from quantiles
    double q10 = quantile_vec(edge_lengths, 0.1);
    double q90 = quantile_vec(edge_lengths, 0.9);
    double init_std = (q90 - q10) / 4.0;

    double mu1 = q10, sigma1 = init_std;
    double mu2 = q90, sigma2 = init_std;

    // Initial assignment probs from PDFs
    Eigen::MatrixXd assignment_probs(n, 2);
    for (int i = 0; i < n; ++i) {
        assignment_probs(i, 0) = normal_pdf(edge_lengths[i], mu1, sigma1);
        assignment_probs(i, 1) = normal_pdf(edge_lengths[i], mu2, sigma2);
    }

    // Hard-assign outliers to noise
    double outlier_threshold = q90 + 3.0 * init_std;
    std::vector<bool> outlier(n, false);
    std::vector<int> updating_ids;
    updating_ids.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (edge_lengths[i] > outlier_threshold) {
            outlier[i] = true;
            assignment_probs(i, 0) = 0.0;
            assignment_probs(i, 1) = 1.0;
        } else {
            updating_ids.push_back(i);
        }
    }

    // Normalize rows
    for (int i = 0; i < n; ++i) {
        double row_sum = assignment_probs(i, 0) + assignment_probs(i, 1);
        if (row_sum > 0.0) {
            assignment_probs(i, 0) /= row_sum;
            assignment_probs(i, 1) /= row_sum;
        }
    }

    std::vector<double> diffs;

    // EM iterations
    int n_iters = max_iters;
    for (int iter = 0; iter < max_iters; ++iter) {
        // E-step
        expect_noise_probabilities(
            assignment_probs, mu1, sigma1, mu2, sigma2,
            edge_lengths, adj_list, updating_ids, min_confidence
        );

        // M-step: refit Normals from weighted data
        std::vector<double> w1(n), w2(n);
        for (int i = 0; i < n; ++i) {
            w1[i] = assignment_probs(i, 0);
            w2[i] = assignment_probs(i, 1);
        }

        auto [new_mu1, new_sigma1] = wmean_std(edge_lengths.data(), w1.data(), n);
        auto [new_mu2, new_sigma2] = wmean_std(edge_lengths.data(), w2.data(), n);

        // Prevent degenerate sigma
        new_sigma1 = std::max(new_sigma1, 1e-10);
        new_sigma2 = std::max(new_sigma2, 1e-10);

        // Convergence: max relative parameter change
        double param_diff = std::max({ // GCOVR_EXCL_LINE: dead GCC block; statement counted on the following lines
            std::abs(new_mu1 - mu1) / std::max(std::abs(mu1), 1e-20),
            std::abs(new_mu2 - mu2) / std::max(std::abs(mu2), 1e-20),
            std::abs(new_sigma1 - sigma1) / std::max(std::abs(sigma1), 1e-20),
            std::abs(new_sigma2 - sigma2) / std::max(std::abs(sigma2), 1e-20)
        });
        diffs.push_back(param_diff);

        mu1 = new_mu1; sigma1 = new_sigma1;
        mu2 = new_mu2; sigma2 = new_sigma2;

        if (param_diff < tol) {
            n_iters = iter + 1;
            break;
        }
    }

    if (verbose) {
        spdlog::info("Noise estimation stopped after {} iterations. Error: {:.4g}. Converged: {}.",
                     n_iters, diffs.back(), diffs.back() <= tol);
    }

    // Ensure d1 is the signal (lower mean) component
    if (mu1 > mu2) {
        std::swap(mu1, mu2);
        std::swap(sigma1, sigma2);
    }

    // Final posterior computation
    double n1 = 0.0;
    for (int i = 0; i < n; ++i) n1 += assignment_probs(i, 0);
    double n2 = n - n1;

    for (int i = 0; i < n; ++i) {
        assignment_probs(i, 0) = n1 * normal_pdf(edge_lengths[i], mu1, sigma1);
        assignment_probs(i, 1) = n2 * normal_pdf(edge_lengths[i], mu2, sigma2);
    }
    // Re-apply outlier hard assignment
    for (int i = 0; i < n; ++i) {
        if (outlier[i]) {
            assignment_probs(i, 1) = 1.0;
        }
    }
    // Normalize rows
    for (int i = 0; i < n; ++i) {
        double row_sum = assignment_probs(i, 0) + assignment_probs(i, 1);
        if (row_sum > 0.0) {
            assignment_probs(i, 0) /= row_sum;
            assignment_probs(i, 1) /= row_sum;
        }
    }

    // Apply min_confidence floor
    if (min_confidence) {
        for (int i = 0; i < n; ++i) {
            assignment_probs(i, 0) = (*min_confidence)[i] +
                assignment_probs(i, 0) * (1.0 - (*min_confidence)[i]);
            assignment_probs(i, 1) = 1.0 - assignment_probs(i, 0);
        }
    }

    // Assignment: 1 = signal, 2 = noise
    std::vector<int> assignment(n);
    for (int i = 0; i < n; ++i) {
        assignment[i] = (assignment_probs(i, 1) > 0.5) ? 2 : 1;
    }

    NoiseFitResult result;
    result.assignment_probs = std::move(assignment_probs);
    result.assignment = std::move(assignment);
    result.signal_mu = mu1;
    result.signal_sigma = sigma1;
    result.noise_mu = mu2;
    result.noise_sigma = sigma2;
    result.diffs = std::move(diffs);
    return result;
}

// ============================================================================
// append_confidence
// ============================================================================

ConfidenceEstimationDetails estimate_confidence_details(
    const MoleculeData& data, int nn_id, double prior_confidence
) {
    int n = data.n_molecules();
    if (n == 0) return {};

    // Default nn_id
    if (nn_id <= 0) {
        nn_id = std::max(data.n_genes() / 10, 10);
    }

    Eigen::MatrixXd pos = data.position_matrix();

    // KNN distance to the (nn_id+1)-th neighbor (first is self)
    std::vector<double> mean_dists = knn_kth_distances(pos, nn_id + 1, nn_id);

    // Build molecule graph (unfiltered, matching Julia); the edges are
    // returned for reuse by the segmentation graph
    auto adj_edges = compute_molecule_adjacency(data);
    auto adj_list = build_molecule_graph_from_edges(adj_edges, n);

    // Compute min_confidence from prior segmentation if available
    std::vector<double> min_conf;
    const std::vector<double>* min_conf_ptr = nullptr;
    if (!data.prior_segmentation.empty()) {
        min_conf.resize(n);
        double pc2 = prior_confidence * prior_confidence;
        for (int i = 0; i < n; ++i) {
            min_conf[i] = (data.prior_segmentation[i] > 0) ? pc2 : 0.0;
        }
        min_conf_ptr = &min_conf;
    }

    // Fit noise model
    auto result = fit_noise_probabilities(mean_dists, adj_list, min_conf_ptr,
                                          /*max_iters=*/10000, /*tol=*/0.005, /*verbose=*/true);

    spdlog::info("Noise estimation: signal mu={:.4g}, sigma={:.4g}; noise mu={:.4g}, sigma={:.4g}",
                 result.signal_mu, result.signal_sigma, result.noise_mu, result.noise_sigma);

    ConfidenceEstimationDetails details;
    details.edge_lengths = std::move(mean_dists);
    details.fit_result = std::move(result);
    details.nn_id = nn_id;
    details.adjacency = std::move(adj_edges);
    return details;
}

void append_confidence(MoleculeData& data, int nn_id, double prior_confidence) {
    auto details = estimate_confidence_details(data, nn_id, prior_confidence);
    int n = data.n_molecules();
    data.confidence.resize(n);
    for (int i = 0; i < n; ++i) {
        data.confidence[i] = details.fit_result.assignment_probs(i, 0);
    }
}

} // namespace baysor
