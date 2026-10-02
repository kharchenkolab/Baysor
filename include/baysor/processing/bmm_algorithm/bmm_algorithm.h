#pragma once

#include "baysor/processing/models/bmm_data.h"
#include <cstdint>

namespace baysor {

struct EstepStats {
    std::int64_t n_changed = 0;
};

/// Run the main BMM segmentation algorithm (EM loop).
/// This is the core of Baysor: iterates E-step (molecule assignment) and
/// M-step (component parameter estimation) for n_iters iterations.
///
/// min_molecules_drop:    per-iteration drop threshold (Julia: min_n_samples=2 default)
/// min_molecules_display: threshold used in verbose progress display (Julia: min_molecules_per_cell)
///                        Set to 0 to match min_molecules_drop.
template<int N>
void bmm(BmmData<N>& data,
         int min_molecules_drop = 2,
         int n_iters = 500,
         int assignment_history_depth = 0,
         bool verbose = true,
         int component_split_step = 3,
         bool refine = true,
         bool freeze_composition = false,
         bool freeze_position = false,
         bool freeze_components = false,
         double tol = 0.0,
         int min_molecules_display = 0);  ///< display threshold (0 = same as min_molecules_drop)

/// E-step: reassign molecules to components based on spatial + expression density.
/// `rng_salt` keys the multi-threaded per-chunk RNG streams (pass the iteration
/// index so draws differ between iterations); it has no effect with 1 thread,
/// where the global stream is used in index order.
template<int N>
EstepStats expect_dirichlet_spatial(BmmData<N>& data, bool stochastic = true,
                                    std::uint64_t rng_salt = 0);

/// M-step: re-estimate all component parameters from current assignments
template<int N>
void maximize(BmmData<N>& data, bool freeze_composition = false, bool freeze_position = false);

/// Drop components with fewer than min_n_samples assigned molecules
template<int N>
void drop_unused_components(BmmData<N>& data, int min_n_samples = 2);

/// Split cells that have disconnected spatial components
template<int N>
void split_cells_by_connected_components(BmmData<N>& data);

/// Estimate final assignment by majority vote over history
template<int N>
std::pair<std::vector<int>, std::vector<double>> estimate_assignment_by_history(const BmmData<N>& data);

} // namespace baysor
