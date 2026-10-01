#pragma once

// Baysor's serial UMAP layout optimiser.
//
// This is umappp v2.0.1's internal::optimize_layout (optimize_layout.hpp,
// serial branch) kept inside Baysor instead of patching the FetchContent
// dependency, with one change:
//  - the embedding dimension is a template parameter for the 2-D and 3-D
//    embeddings Baysor computes, so the distance, clamp and update loops are
//    fully unrolled. The operations and their order are unchanged, so the
//    result is bitwise identical to umappp's optimiser (NDim_ = 0 keeps the
//    run-time dimension for any other size).
//
// The header includes umappp's optimize_layout.hpp for EpochData and the
// distance/clamp helpers. Translation units that route umappp's parallel
// loops through the Baysor pool must define UMAPPP_CUSTOM_PARALLEL before
// including any umappp header (see umap_wrappers.cpp).

#include "umappp/optimize_layout.hpp"

#include <cmath>
#include <cstddef>

namespace baysor {
namespace umap_detail {

/// umappp::internal::optimize_layout (serial) with a compile-time embedding
/// dimension (NDim_ > 0; NDim_ = 0 uses ndim_rt).
template <int NDim_, typename Index_, typename Float_, class Rng_>
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
                    const Float_ pd2b = std::pow(dist2, b);
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
                    const Float_ grad_coef = 2 * gamma * b / ((0.001 + dist2) * (a * std::pow(dist2, b) + 1.0));

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
template <typename Index_, typename Float_, class Rng_>
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
        optimize_layout<3>(ndim, embedding, setup, a, b, gamma, initial_alpha, rng, epoch_limit);
    } else if (ndim == 2) {
        optimize_layout<2>(ndim, embedding, setup, a, b, gamma, initial_alpha, rng, epoch_limit);
    } else {
        optimize_layout<0>(ndim, embedding, setup, a, b, gamma, initial_alpha, rng, epoch_limit);
    }
}

} // namespace umap_detail
} // namespace baysor
