#pragma once

#include "baysor/utils/general.h"

#include <Eigen/Dense>
#include <cmath>
#include <vector>

namespace baysor {

/// Shape prior for cell size regularization.
/// Template parameter N = spatial dimensionality (2 or 3).
template<int N>
struct ShapePrior {
    Eigen::Matrix<double, N, 1> std_values;
    Eigen::Matrix<double, N, 1> std_value_stds;
    int n_samples;
};

/// Multivariate normal with cached inverse and pdf divider.
/// Fixed-size for 2D or 3D — no heap allocation per component.
template<int N>
struct MvNormal {
    using Vec = Eigen::Matrix<double, N, 1>;
    using Mat = Eigen::Matrix<double, N, N>;

    Vec mu;
    Mat sigma;
    Mat sigma_inv;
    double pdf_divider;  // 0.5 * log((2*pi)^3 * det(sigma))

    MvNormal() : MvNormal(Vec::Zero(), Mat::Identity()) {}

    MvNormal(const Vec& mu, const Mat& sigma);

    /// Recompute sigma_inv and pdf_divider from current sigma
    void update_cache();

    double log_pdf(const double* x) const;
    double pdf(const double* x) const { return std::exp(log_pdf(x)); }

    /// M-step: re-estimate mean and covariance from data points
    void maximize(const double* pos_data, int n_points, int stride,
                  const double* center_probs = nullptr,
                  const ShapePrior<N>* shape_prior = nullptr,
                  int n_samples = -1);
    void maximize_indexed(const Eigen::MatrixXd& pos_data, const int* ids, int n_points,
                          const std::vector<double>* center_probs = nullptr,
                          const ShapePrior<N>* shape_prior = nullptr,
                          int n_samples = -1);
};

// Explicit instantiations declared (defined in .cpp)
extern template struct MvNormal<2>;
extern template struct MvNormal<3>;

/// Adjust covariance matrix to ensure positive-definiteness
template<int N>
void adjust_cov_matrix(Eigen::Matrix<double, N, N>& sigma, double cov_modifier = 1e-4, double tol = 1e-10);

/// Adjust covariance by shape prior (shrinkage toward expected cell size)
template<int N>
void adjust_cov_by_prior(Eigen::Matrix<double, N, N>& sigma, const ShapePrior<N>& prior, int n_samples);

} // namespace baysor
