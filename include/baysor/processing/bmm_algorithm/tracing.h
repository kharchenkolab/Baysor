#pragma once

#include "baysor/processing/models/bmm_data.h"

#include <type_traits>
#include <unordered_map>
#include <vector>

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

/// Same for the delta-encoded history of BmmData (no full rows are built).
std::unordered_map<int, int> estimate_component_lifespan_history(
    const AssignmentHistory& assignment_history
);

/// Overload for AssignmentHistory. A template so that a braced `{}` argument
/// still selects the vector-of-rows overload unambiguously.
template <class History,
          std::enable_if_t<std::is_same_v<History, AssignmentHistory>, int> = 0>
std::unordered_map<int, int> estimate_component_lifespan(const History& assignment_history) {
    return estimate_component_lifespan_history(assignment_history);
}

} // namespace baysor
