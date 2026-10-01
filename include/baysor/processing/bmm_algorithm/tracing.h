#pragma once

#include "baysor/processing/models/bmm_data.h"

namespace baysor {

/// Record number of components above various molecule-count thresholds
template<int N>
void trace_n_components(BmmData<N>& data, int min_molecules_per_cell);

/// Same as trace_n_components, with the per-cell molecule counts taken from
/// `ids_by_comp`, the grouping of the current assignment.
template<int N>
void trace_n_components(BmmData<N>& data, int min_molecules_per_cell,
                        const IdsByComponent& ids_by_comp);

/// Record current assignment (global GUIDs) into history ring buffer
template<int N>
void trace_assignment_history(BmmData<N>& data, int assignment_history_depth);

/// Estimate how long each component has existed in the assignment history
std::unordered_map<int, int> estimate_component_lifespan(
    const std::vector<std::vector<int>>& assignment_history
);

} // namespace baysor
