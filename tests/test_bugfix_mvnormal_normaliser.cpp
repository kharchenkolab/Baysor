// BUG-4: inconsistent normalising constant in MvNormal.
//
// Julia reference: src/processing/distributions/MvNormal.jl
//   - norm_pdf_divider(Σ) = 0.5 * log((2π)^3 * det(Σ))   (hardcoded power 3,
//     even for 2D distributions)
//   - the constructor MvNormalF(μ, Σ) uses norm_pdf_divider(Σ)
//   - update_cache!(dist) recomputes norm_pdf_divider(dist.Σ)
//
// The C++ default constructor used to initialise pdf_divider with the
// mathematically conventional (2π)^(N/2) form instead, so a default-constructed
// MvNormal<2> disagreed with the explicit constructor / update_cache() by a
// constant 0.5*log(2π) log-density offset. These tests pin the Julia-parity
// behaviour: every construction path must use (2π)^3.

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

} // namespace

// ---------------------------------------------------------------------------
// Regression: default-constructed vs explicit/updated must agree (fails for
// N=2 without the fix, where the default constructor used (2π)^(N/2)).
// ---------------------------------------------------------------------------
TEST(Bug4_MvNormalNormaliser, DefaultConstructedAgreesWithExplicitAndUpdated) {
    // --- 2D ---
    {
        MvNormal<2> def;
        const Eigen::Vector2d mu = Eigen::Vector2d::Zero();
        MvNormal<2> explicit_ctor(mu, Eigen::Matrix2d::Identity());
        MvNormal<2> updated;
        updated.mu = mu;
        updated.sigma = Eigen::Matrix2d::Identity();
        updated.update_cache();

        const double x[2] = {0.4, -0.7};
        EXPECT_DOUBLE_EQ(def.log_pdf(x), explicit_ctor.log_pdf(x))
            << "default vs explicit constructor (2D)";
        EXPECT_DOUBLE_EQ(def.log_pdf(x), updated.log_pdf(x))
            << "default vs update_cache (2D)";
        // Julia: MvNormalF(μ) with default Σ = I -> divider = 0.5*log((2π)^3)
        EXPECT_NEAR(def.pdf_divider, 1.5 * std::log(2.0 * baysor::kPi), 1e-12);
        EXPECT_NEAR(def.log_pdf(x),
                    -0.5 * (x[0] * x[0] + x[1] * x[1]) - 1.5 * std::log(2.0 * baysor::kPi),
                    1e-12);
    }

    // --- 3D ---
    {
        MvNormal<3> def;
        const Eigen::Vector3d mu = Eigen::Vector3d::Zero();
        MvNormal<3> explicit_ctor(mu, Eigen::Matrix3d::Identity());
        MvNormal<3> updated;
        updated.mu = mu;
        updated.sigma = Eigen::Matrix3d::Identity();
        updated.update_cache();

        const double x[3] = {0.4, -0.7, 1.1};
        EXPECT_DOUBLE_EQ(def.log_pdf(x), explicit_ctor.log_pdf(x))
            << "default vs explicit constructor (3D)";
        EXPECT_DOUBLE_EQ(def.log_pdf(x), updated.log_pdf(x))
            << "default vs update_cache (3D)";
        EXPECT_NEAR(def.pdf_divider, 1.5 * std::log(2.0 * baysor::kPi), 1e-12);
        EXPECT_NEAR(def.log_pdf(x),
                    -0.5 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2])
                        - 1.5 * std::log(2.0 * baysor::kPi),
                    1e-12);
    }
}

// ---------------------------------------------------------------------------
// Numeric agreement with Julia's formula, 2D.
// Expected literals computed from Julia's MvNormalF log_pdf formula with
//   μ = (1.0, -0.5), Σ = [[2.0, 0.5], [0.5, 1.5]], x = (1.5, 0.5).
// ---------------------------------------------------------------------------
TEST(Bug4_MvNormalNormaliser, JuliaFormulaAgreement2D) {
    Eigen::Vector2d mu;
    mu << 1.0, -0.5;
    Eigen::Matrix2d sigma;
    sigma << 2.0, 0.5,
             0.5, 1.5;
    MvNormal<2> mv(mu, sigma);

    Eigen::Vector2d x;
    x << 1.5, 0.5;
    double buf[2] = {x(0), x(1)};

    EXPECT_NEAR(mv.log_pdf(buf), -3.6035251463623488, 1e-12);
    EXPECT_NEAR(mv.log_pdf(buf), julia_log_pdf<2>(mu, sigma, x), 1e-12);
    EXPECT_NEAR(mv.pdf(buf), std::exp(-3.6035251463623488), 1e-12);
}

// ---------------------------------------------------------------------------
// Numeric agreement with Julia's formula, 3D.
// Expected literals computed from Julia's MvNormalF log_pdf formula with
//   μ = (0.5, -1.0, 2.0),
//   Σ = [[2.0, 0.3, 0.1], [0.3, 1.2, 0.2], [0.1, 0.2, 0.9]], x = (0, 0, 1).
// ---------------------------------------------------------------------------
TEST(Bug4_MvNormalNormaliser, JuliaFormulaAgreement3D) {
    Eigen::Vector3d mu;
    mu << 0.5, -1.0, 2.0;
    Eigen::Matrix3d sigma;
    sigma << 2.0, 0.3, 0.1,
             0.3, 1.2, 0.2,
             0.1, 0.2, 0.9;
    MvNormal<3> mv(mu, sigma);

    Eigen::Vector3d x;
    x << 0.0, 0.0, 1.0;
    double buf[3] = {x(0), x(1), x(2)};

    EXPECT_NEAR(mv.log_pdf(buf), -4.426300708163545, 1e-12);
    EXPECT_NEAR(mv.log_pdf(buf), julia_log_pdf<3>(mu, sigma, x), 1e-12);
    EXPECT_NEAR(mv.pdf(buf), std::exp(-4.426300708163545), 1e-12);
}
