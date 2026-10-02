// BUG-4: inconsistent normalising constant in MvNormal.
//
// Julia (src/processing/distributions/MvNormal.jl) uses
// norm_pdf_divider(Σ) = 0.5 * log((2π)^3 * det(Σ)) -- power 3 even for 2D --
// in the constructor and in update_cache!. The C++ default constructor used
// (2π)^(N/2) instead, a constant 0.5*log(2π) log-density offset in 2D. Every
// construction path must use (2π)^3.

#include <gtest/gtest.h>

#include "baysor/processing/distributions/mv_normal.h"

#include <cmath>

namespace {

using baysor::MvNormal;

// Independent transcription of Julia's norm_pdf_divider + log_pdf:
//   log_pdf(d, x) = -0.5 * (x-μ)' * inv(Σ) * (x-μ) - 0.5 * log((2π)^3 * det(Σ))
template<int N>
double julia_log_pdf(const Eigen::Matrix<double, N, 1>& mu,
                     const Eigen::Matrix<double, N, N>& sigma,
                     const Eigen::Matrix<double, N, 1>& x) {
    const Eigen::Matrix<double, N, 1> dx = x - mu;
    const double mahalanobis = dx.dot(sigma.inverse() * dx);
    const double divider = 0.5 * std::log(std::pow(2.0 * baysor::kPi, 3) * sigma.determinant());
    return -0.5 * mahalanobis - divider;
}

template<int N>
void expect_julia_log_pdf(const MvNormal<N>& mv, const Eigen::Matrix<double, N, 1>& x) {
    const double expected = julia_log_pdf<N>(mv.mu, mv.sigma, x);
    EXPECT_NEAR(mv.log_pdf(x.data()), expected, 1e-12);
    EXPECT_NEAR(mv.pdf(x.data()), std::exp(expected), 1e-12);
}

template<int N>
void expect_default_is_standard_normal(const Eigen::Matrix<double, N, 1>& x) {
    const MvNormal<N> def;
    EXPECT_TRUE(def.mu.isZero());
    EXPECT_TRUE(def.sigma.isIdentity());
    EXPECT_TRUE(def.sigma_inv.isIdentity());
    // Julia: MvNormalF(μ) with default Σ = I -> divider = 0.5*log((2π)^3)
    EXPECT_NEAR(def.pdf_divider, 1.5 * std::log(2.0 * baysor::kPi), 1e-12);
    expect_julia_log_pdf<N>(def, x);

    // Same as the explicit constructor and as update_cache()
    MvNormal<N> updated;
    updated.update_cache();
    const MvNormal<N> explicit_ctor(def.mu, def.sigma);
    EXPECT_DOUBLE_EQ(def.log_pdf(x.data()), explicit_ctor.log_pdf(x.data()));
    EXPECT_DOUBLE_EQ(def.log_pdf(x.data()), updated.log_pdf(x.data()));
}

} // namespace

// Fails for N=2 without the fix, where the default constructor used (2π)^(N/2).
TEST(Bug4_MvNormalNormaliser, DefaultConstructedAgreesWithExplicitAndUpdated) {
    expect_default_is_standard_normal<2>(Eigen::Vector2d(0.4, -0.7));
    expect_default_is_standard_normal<3>(Eigen::Vector3d(0.4, -0.7, 1.1));
}

// Expected literals computed from Julia's MvNormalF log_pdf formula.
TEST(Bug4_MvNormalNormaliser, JuliaFormulaAgreement) {
    Eigen::Matrix2d sigma2;
    sigma2 << 2.0, 0.5,
              0.5, 1.5;
    const MvNormal<2> mv2(Eigen::Vector2d(1.0, -0.5), sigma2);
    const Eigen::Vector2d x2(1.5, 0.5);
    EXPECT_NEAR(mv2.log_pdf(x2.data()), -3.6035251463623488, 1e-12);
    expect_julia_log_pdf<2>(mv2, x2);

    Eigen::Matrix3d sigma3;
    sigma3 << 2.0, 0.3, 0.1,
              0.3, 1.2, 0.2,
              0.1, 0.2, 0.9;
    const MvNormal<3> mv3(Eigen::Vector3d(0.5, -1.0, 2.0), sigma3);
    const Eigen::Vector3d x3(0.0, 0.0, 1.0);
    EXPECT_NEAR(mv3.log_pdf(x3.data()), -4.426300708163545, 1e-12);
    expect_julia_log_pdf<3>(mv3, x3);
}
