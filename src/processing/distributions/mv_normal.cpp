#include "baysor/processing/distributions/mv_normal.h"

#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>

namespace baysor {

template<int N>
MvNormal<N>::MvNormal(const Vec& mu, const Mat& sigma)
    : mu(mu), sigma(sigma)
{
    update_cache();
}

template<int N>
void MvNormal<N>::update_cache() {
    sigma_inv = sigma.inverse();
    // Match Julia's current MvNormalF implementation exactly. In 2D it uses
    // a hardcoded (2π)^3 normalization term rather than (2π)^N.
    pdf_divider = 0.5 * std::log(std::pow(2.0 * kPi, 3) * sigma.determinant());
}

template<int N>
double MvNormal<N>::log_pdf(const double* x) const {
    Vec dx;
    for (int i = 0; i < N; ++i) dx[i] = x[i] - mu[i];
    return -0.5 * dx.dot(sigma_inv * dx) - pdf_divider;
}

template<int N>
void MvNormal<N>::maximize(const double* pos_data, int n_points, int stride,
                            const double* center_probs,
                            const ShapePrior<N>* shape_prior,
                            int n_samples) {
    if (n_points == 0) return;

    // ---- Weighted mean ----
    mu.setZero();
    double sum_w = 0.0;

    for (int i = 0; i < n_points; ++i) {
        double w = (center_probs != nullptr) ? std::max(center_probs[i], 0.01) : 1.0;
        for (int j = 0; j < N; ++j) {
            mu[j] += pos_data[i * stride + j] * w;
        }
        sum_w += w;
    }
    mu /= sum_w;

    // ---- Covariance (requires > N points) ----
    if (n_points <= N) {
        // Not enough points to estimate covariance; keep existing sigma, update cache
        update_cache();
        return;
    }

    // Unweighted sample covariance (matches Julia's estimate_sample_cov!)
    sigma.setZero();
    for (int i = 0; i < n_points; ++i) {
        for (int r = 0; r < N; ++r) {
            double dr = pos_data[i * stride + r] - mu[r];
            for (int c = 0; c < N; ++c) {
                double dc = pos_data[i * stride + c] - mu[c];
                sigma(r, c) += dr * dc;
            }
        }
    }
    sigma /= static_cast<double>(n_points);

    if (shape_prior != nullptr) {
        int ns = (n_samples >= 0) ? n_samples : n_points;
        adjust_cov_by_prior<N>(sigma, *shape_prior, ns);
    } else {
        // Even without shape prior, ensure positive definiteness
        adjust_cov_matrix<N>(sigma);
    }

    update_cache();
}

template<int N>
void MvNormal<N>::maximize_indexed(const Eigen::MatrixXd& pos_data, const int* ids, int n_points,
                                   const std::vector<double>* center_probs,
                                   const ShapePrior<N>* shape_prior,
                                   int n_samples) {
    if (n_points == 0) return;

    mu.setZero();
    double sum_w = 0.0;

    for (int i = 0; i < n_points; ++i) {
        int mol_id = ids[i];
        double w = (center_probs != nullptr) ? std::max((*center_probs)[mol_id], 0.01) : 1.0;
        for (int j = 0; j < N; ++j) {
            mu[j] += pos_data(j, mol_id) * w;
        }
        sum_w += w;
    }
    mu /= sum_w;

    if (n_points <= N) {
        update_cache();
        return;
    }

    sigma.setZero();
    for (int i = 0; i < n_points; ++i) {
        int mol_id = ids[i];
        for (int r = 0; r < N; ++r) {
            double dr = pos_data(r, mol_id) - mu[r];
            for (int c = 0; c < N; ++c) {
                double dc = pos_data(c, mol_id) - mu[c];
                sigma(r, c) += dr * dc;
            }
        }
    }
    sigma /= static_cast<double>(n_points);

    if (shape_prior != nullptr) {
        int ns = (n_samples >= 0) ? n_samples : n_points;
        adjust_cov_by_prior<N>(sigma, *shape_prior, ns);
    } else {
        adjust_cov_matrix<N>(sigma);
    }

    update_cache();
}

// ============================================================================
// adjust_cov_matrix: nudge diagonal until positive-definite
// ============================================================================

template<int N>
void adjust_cov_matrix(Eigen::Matrix<double, N, N>& sigma, double cov_modifier, double tol) {
    using Mat = Eigen::Matrix<double, N, N>;

    // Symmetrize
    sigma = (sigma + sigma.transpose()) * 0.5;

    for (int iter = 0; iter < 100; ++iter) {
        Eigen::LLT<Mat> llt(sigma);
        if (llt.info() == Eigen::Success && sigma.determinant() > tol) return;
        for (int i = 0; i < N; ++i) sigma(i, i) += cov_modifier;
        cov_modifier *= 1.5;
    }
}

// ============================================================================
// adjust_cov_by_prior: Bayesian shrinkage toward expected cell size
// ============================================================================

// Posterior variance formula from Julia MvNormal.jl var_posterior
static double var_posterior(double prior_std, double prior_std_std, double eigen_value,
                             double n_samples, double prior_n_samples) {
    double d_std = std::sqrt(eigen_value) - prior_std;
    double std_adj = prior_std
        + std::copysign(
              (std::sqrt(std::abs(d_std) / prior_std_std + 1.0) - 1.0) * prior_std_std,
              d_std);
    double blended = (prior_n_samples * prior_std + n_samples * std_adj)
                     / (prior_n_samples + n_samples);
    return blended * blended;
}

template<int N>
void adjust_cov_by_prior(Eigen::Matrix<double, N, N>& sigma, const ShapePrior<N>& prior, int n_samples) {
    using Mat = Eigen::Matrix<double, N, N>;

    // Match Julia's current temporary fix exactly. It checks the signed ratio,
    // not the absolute ratio, so any sufficiently negative covariance entry is
    // also zeroed before the eigendecomposition.
    if (N >= 2) {
        double ratio = sigma(1, 0) / std::max(sigma(0, 0), sigma(1, 1));
        if (ratio < 1e-5) {
            sigma(0, 1) = sigma(1, 0) = 0.0;
        }
    }

    Eigen::SelfAdjointEigenSolver<Mat> es(sigma);
    auto eigen_values = es.eigenvalues();
    const Mat& V = es.eigenvectors();

    // Compute posterior variance for each eigenvalue
    Eigen::Matrix<double, N, 1> post_vars;
    for (int i = 0; i < N; ++i) {
        double ev = std::max(eigen_values[i], 0.0);
        post_vars[i] = var_posterior(prior.std_values[i], prior.std_value_stds[i],
                                     ev, static_cast<double>(n_samples),
                                     static_cast<double>(prior.n_samples));
    }

    // Reconstruct: sigma = V * diag(post_vars) * V^{-1}
    // For symmetric PD matrix V is orthogonal so V^{-1} = V^T
    sigma = V * post_vars.asDiagonal() * V.transpose();
    // Symmetrize and ensure pos-def
    sigma = (sigma + sigma.transpose()) * 0.5;
    adjust_cov_matrix<N>(sigma);
}

// Explicit instantiations
template struct MvNormal<2>;
template struct MvNormal<3>;
template void adjust_cov_matrix<2>(Eigen::Matrix2d&, double, double);
template void adjust_cov_matrix<3>(Eigen::Matrix3d&, double, double);
template void adjust_cov_by_prior<2>(Eigen::Matrix2d&, const ShapePrior<2>&, int);
template void adjust_cov_by_prior<3>(Eigen::Matrix3d&, const ShapePrior<3>&, int);

} // namespace baysor
