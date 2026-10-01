#pragma once

// Baysor's serial UMAP layout optimiser.
//
// This is umappp v2.0.1's internal::optimize_layout (optimize_layout.hpp,
// serial branch) kept inside Baysor instead of patching the FetchContent
// dependency, with two changes:
//  - the embedding dimension is a template parameter for the 2-D and 3-D
//    embeddings Baysor computes, so the distance, clamp and update loops are
//    fully unrolled. The operations and their order are unchanged, so the
//    result is bitwise identical to umappp's optimiser (NDim_ = 0 keeps the
//    run-time dimension for any other size);
//  - optionally, std::pow(d2, b) in the gradient is replaced by fast_pow(), a
//    table-driven approximation without libm calls (relative error < 1e-13).
//    This only changes the last bits of the gradient coefficients.
//
// The header includes umappp's optimize_layout.hpp for EpochData and the
// distance/clamp helpers. Translation units that route umappp's parallel
// loops through the Baysor pool must define UMAPPP_CUSTOM_PARALLEL before
// including any umappp header (see umap_wrappers.cpp).

#include "umappp/optimize_layout.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>

namespace baysor {
namespace umap_detail {

/// Tables of fast_pow, computed at compile time (no libm, so the same on every
/// platform): for the 256 mantissa intervals [1 + i/256, 1 + (i+1)/256),
/// 1/c_i and ln(c_i) at the interval centres c_i, and 2^(j/256), j < 256.
struct FastPowTables {
    double inv_c[256];
    double ln_c[256];
    double exp2_j256[256];
};

constexpr FastPowTables make_fast_pow_tables() {
    FastPowTables t{};
    for (int i = 0; i < 256; ++i) {
        const double c = 1.0 + (i + 0.5) / 256.0;
        t.inv_c[i] = 1.0 / c;
        // ln(c) = 2 atanh(u), u = (c - 1) / (c + 1) < 1/5: odd series to u^41,
        // evaluated by Horner from the smallest term up.
        const double u = (c - 1.0) / (c + 1.0), u2 = u * u;
        double sum = 1.0 / 41;
        for (int k = 19; k >= 0; --k) sum = sum * u2 + 1.0 / (2 * k + 1);
        t.ln_c[i] = 2.0 * u * sum;
    }
    for (int j = 0; j < 256; ++j) {
        // 2^(j/256) = exp(j ln2 / 256), Taylor series of exp on [0, ln2).
        const double z = j * 0.69314718055994531 / 256.0;
        double sum = 1.0;
        for (int k = 25; k >= 1; --k) sum = 1.0 + sum * z / k;
        t.exp2_j256[j] = sum;
    }
    return t;
}

inline constexpr FastPowTables fast_pow_tables = make_fast_pow_tables();

/// x^b for finite x > 0, without libm calls. ln(x) = e ln2 + ln(c_i) +
/// log1p(m / c_i - 1) with the mantissa m in its 1/256-wide interval
/// (|m / c_i - 1| <= 1/512, degree-5 series), and exp(y) = 2^(k/256) e^z with
/// |z| <= ln2 / 512 (degree-4 Taylor polynomial). Both polynomials are
/// evaluated in Estrin form: the UMAP optimiser is one long dependency chain
/// (each update feeds the next distance), so the latency of pow matters more
/// than its instruction count. Maximum relative error ~2e-14 against
/// std::pow; ~20 % lower latency and ~3x fewer instructions than glibc's pow
/// on x86-64 without FMA. (The review's prototype, a degree-9/degree-8 Horner
/// polynomial with ~1e-9 error, had a higher latency than glibc's pow and
/// moved the NCV colours twice as much.) Inputs where the approximation does
/// not apply (x <= 0, subnormal, inf/NaN, |b ln x| >= 708) fall back to
/// std::pow.
inline double fast_pow(double x, double b) {
    const FastPowTables& tab = fast_pow_tables;
    std::uint64_t bits;
    std::memcpy(&bits, &x, sizeof bits);
    // Positive normal finite x: biased exponent in [1, 2046], sign bit clear.
    if (bits - 0x0010000000000000ULL >= 0x7fe0000000000000ULL) {
        return std::pow(x, b);
    }

    const std::int64_t e = static_cast<std::int64_t>(bits >> 52) - 1023;
    const std::uint64_t mbits = (bits & 0x000fffffffffffffULL) | 0x3ff0000000000000ULL; // m in [1, 2)
    double m;
    std::memcpy(&m, &mbits, sizeof m);
    const int i = static_cast<int>((mbits >> 44) & 255); // top 8 mantissa bits
    const double r = m * tab.inv_c[i] - 1.0;
    const double r2 = r * r;
    const double log1p_r = r + r2 * ((-0.5 + r * (1.0 / 3)) + r2 * (-0.25 + r * 0.2));
    const double y = b * ((static_cast<double>(e) * 0.69314718055994531 + tab.ln_c[i]) + log1p_r); // ln(x^b)
    if (!(std::fabs(y) < 708.0)) {
        return std::pow(x, b);
    }

    // y = (k / 256) ln2 + z, k = round(256 y / ln2).
    constexpr double shifter = 6755399441055744.0; // 1.5 * 2^52: round to nearest
    const double kd = (y * (256 * 1.4426950408889634) + shifter) - shifter;
    const std::int64_t k = static_cast<std::int64_t>(kd);
    const double z = y - kd * (0.69314718055994531 / 256);
    const double z2 = z * z;
    const double ez = (1.0 + z) + z2 * ((0.5 + z * (1.0 / 6)) + z2 * (1.0 / 24));
    // 2^(k >> 8), arithmetic shift: floor division for negative k.
    const std::uint64_t scale = static_cast<std::uint64_t>((k >> 8) + 1023) << 52;
    double s;
    std::memcpy(&s, &scale, sizeof s);
    return ez * (tab.exp2_j256[k & 255] * s);
}

template <bool FastPow_, typename Float_>
inline Float_ gradient_pow(Float_ x, Float_ b) {
    if constexpr (FastPow_) {
        static_assert(std::is_same_v<Float_, double>, "fast_pow is implemented for double only");
        return fast_pow(x, b);
    } else {
        return std::pow(x, b);
    }
}

/// umappp::internal::optimize_layout (serial) with a compile-time embedding
/// dimension (NDim_ > 0; NDim_ = 0 uses ndim_rt) and a selectable pow.
template <int NDim_, bool FastPow_, typename Index_, typename Float_, class Rng_>
void optimize_layout(
    std::size_t ndim_rt,
    Float_* embedding,
    umappp::internal::EpochData<Index_, Float_>& setup,
    Float_ a,
    Float_ b,
    Float_ gamma,
    Float_ initial_alpha,
    Rng_& rng,
    int epoch_limit
) {
    using umappp::internal::clamp;
    using umappp::internal::quick_squared_distance;

    const std::size_t ndim = NDim_ > 0 ? static_cast<std::size_t>(NDim_) : ndim_rt;
    auto& n = setup.current_epoch;
    auto num_epochs = setup.total_epochs;
    auto limit_epochs = num_epochs;
    if (epoch_limit > 0) {
        limit_epochs = std::min(epoch_limit, num_epochs);
    }

    const std::size_t num_obs = setup.head.size();
    for (; n < limit_epochs; ++n) {
        const Float_ epoch = n;
        const Float_ alpha = initial_alpha * (1.0 - epoch / num_epochs);

        for (std::size_t i = 0; i < num_obs; ++i) {
            std::size_t start = (i == 0 ? 0 : setup.head[i - 1]), end = setup.head[i];
            Float_* left = embedding + i * ndim;

            for (std::size_t j = start; j < end; ++j) {
                if (setup.epoch_of_next_sample[j] > epoch) {
                    continue;
                }

                {
                    Float_* right = embedding + static_cast<std::size_t>(setup.tail[j]) * ndim;
                    Float_ dist2 = quick_squared_distance(left, right, ndim);
                    const Float_ pd2b = gradient_pow<FastPow_>(dist2, b);
                    const Float_ grad_coef = (-2 * a * b * pd2b) / (dist2 * (a * pd2b + 1.0));

                    for (std::size_t d = 0; d < ndim; ++d) {
                        auto& l = left[d];
                        auto& r = right[d];
                        Float_ gradient = alpha * clamp(grad_coef * (l - r));
                        l += gradient;
                        r -= gradient;
                    }
                }

                // 'epochs_per_negative_sample' is epochs_per_sample[j] / negative_sample_rate,
                // used inline as upstream does to keep the same round-off.
                const std::size_t num_neg_samples = (epoch - setup.epoch_of_next_negative_sample[j]) *
                    setup.negative_sample_rate / setup.epochs_per_sample[j];

                for (std::size_t p = 0; p < num_neg_samples; ++p) {
                    std::size_t sampled = aarand::discrete_uniform(rng, num_obs);
                    if (sampled == i) {
                        continue;
                    }

                    const Float_* right = embedding + sampled * ndim;
                    Float_ dist2 = quick_squared_distance(left, right, ndim);
                    const Float_ grad_coef = 2 * gamma * b / ((0.001 + dist2) * (a * gradient_pow<FastPow_>(dist2, b) + 1.0));

                    for (std::size_t d = 0; d < ndim; ++d) {
                        left[d] += alpha * clamp(grad_coef * (left[d] - right[d]));
                    }
                }

                setup.epoch_of_next_sample[j] += setup.epochs_per_sample[j];
                setup.epoch_of_next_negative_sample[j] = epoch;
            }
        }
    }
}

/// Run optimize_layout with the fixed-dimension instance for 2-D and 3-D
/// embeddings and the run-time dimension otherwise.
template <bool FastPow_, typename Index_, typename Float_, class Rng_>
void optimize_layout_dispatch(
    std::size_t ndim,
    Float_* embedding,
    umappp::internal::EpochData<Index_, Float_>& setup,
    Float_ a,
    Float_ b,
    Float_ gamma,
    Float_ initial_alpha,
    Rng_& rng,
    int epoch_limit
) {
    if (ndim == 3) {
        optimize_layout<3, FastPow_>(ndim, embedding, setup, a, b, gamma, initial_alpha, rng, epoch_limit);
    } else if (ndim == 2) {
        optimize_layout<2, FastPow_>(ndim, embedding, setup, a, b, gamma, initial_alpha, rng, epoch_limit);
    } else {
        optimize_layout<0, FastPow_>(ndim, embedding, setup, a, b, gamma, initial_alpha, rng, epoch_limit);
    }
}

} // namespace umap_detail
} // namespace baysor
