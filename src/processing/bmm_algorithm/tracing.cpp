#include "baysor/processing/bmm_algorithm/tracing.h"
#include "baysor/utils/general.h"

#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_set>

namespace baysor {

namespace {

template<int N>
void push_n_components_entry(BmmData<N>& data, int min_molecules_per_cell,
                             const std::vector<int>& n_mols) {
    // Compute unique thresholds: {max(round(0.5*min),1), max(min,1), max(2*min,1), max(5*min,1)}
    std::set<int> thresh_set;
    for (double mult : {0.5, 1.0, 2.0, 5.0}) {
        thresh_set.insert(std::max(static_cast<int>(std::round(mult * min_molecules_per_cell)), 1));
    }

    std::unordered_map<int, int> entry;
    for (int t : thresh_set) {
        int cnt = 0;
        for (int m : n_mols) {
            if (m >= t) cnt++;
        }
        entry[t] = cnt;
    }

    data.n_components_trace.push_back(std::move(entry));
}

} // namespace

template<int N>
void trace_n_components(BmmData<N>& data, int min_molecules_per_cell) {
    push_n_components_entry(data, min_molecules_per_cell, data.num_molecules_per_cell());
}

template<int N>
void trace_n_components(BmmData<N>& data, int min_molecules_per_cell,
                        const IdsByComponent& ids_by_comp) {
    std::vector<int> n_mols(data.n_components());
    for (int c = 0; c < data.n_components(); ++c) n_mols[c] = ids_by_comp.size(c);
    push_n_components_entry(data, min_molecules_per_cell, n_mols);
}

template<int N>
void trace_assignment_history(BmmData<N>& data, int assignment_history_depth) {
    if (assignment_history_depth <= 0) return;

    // Trim to depth - 1 first, so the evicted delta buffer is recycled
    auto& history = data.assignment_history;
    while (!history.empty() && static_cast<int>(history.size()) >= assignment_history_depth) {
        history.pop_front();
    }

    // Global assignment: local 1-based IDs replaced with component GUIDs
    const auto& assignment = data.assignment;
    const auto& components = data.components;
    history.push_back_generated(data.n_molecules(), [&](int i) {
        const int a = assignment[i];
        return (a > 0) ? components[a - 1].guid : 0;
    });
}

std::unordered_map<int, int> estimate_component_lifespan(
    const std::vector<std::vector<int>>& assignment_history
) {
    std::unordered_map<int, int> lifespans;
    int total = static_cast<int>(assignment_history.size());
    if (total == 0) return lifespans;

    // Initialize from the most recent history entry
    std::unordered_set<int> still_tracking;
    for (int guid : assignment_history[total - 1]) {
        if (guid > 0) { still_tracking.insert(guid); lifespans[guid] = 1; }
    }

    // Walk backward: extend lifespan only while the consecutive streak is unbroken
    for (int it = total - 2; it >= 0 && !still_tracking.empty(); --it) {
        // Collect which tracked guids appear in this iteration
        std::unordered_set<int> present;
        for (int guid : assignment_history[it]) {
            if (guid > 0 && still_tracking.count(guid)) present.insert(guid);
        }
        // Keep only guids whose streak continues; drop the rest
        std::unordered_set<int> new_tracking;
        for (int guid : still_tracking) {
            if (present.count(guid)) { lifespans[guid]++; new_tracking.insert(guid); }
        }
        still_tracking = std::move(new_tracking);
    }
    return lifespans;
}

std::unordered_map<int, int> estimate_component_lifespan_history(
    const AssignmentHistory& history
) {
    std::unordered_map<int, int> lifespans;
    const int total = static_cast<int>(history.size());
    if (total == 0) return lifespans;

    // Walk the entries backward from the newest one, keeping the current row
    // and the number of molecules per guid up to date with the deltas.
    std::vector<int> cur = history.newest();
    int max_guid = 0;
    for (int g : cur) max_guid = std::max(max_guid, g);
    for (int t = 0; t + 1 < total; ++t) {
        for (const auto& c : history.changes_after(t)) max_guid = std::max(max_guid, c.value);
    }
    std::vector<int> count(static_cast<size_t>(max_guid) + 1, 0);
    std::vector<int> tracking;
    for (int guid : cur) {
        if (guid < 0) continue;
        if (count[guid]++ == 0 && guid > 0) {
            // Same keys and insertion order as the row-based version
            lifespans[guid] = 1;
            tracking.push_back(guid);
        }
    }

    // Extend the lifespan only while the consecutive streak is unbroken
    for (int it = total - 2; it >= 0 && !tracking.empty(); --it) {
        for (const auto& c : history.changes_after(it)) {
            const int old_g = cur[c.mol];
            if (old_g >= 0) --count[old_g];
            if (c.value >= 0) ++count[c.value];
            cur[c.mol] = c.value;
        }
        size_t kept = 0;
        for (int guid : tracking) {
            if (count[guid] > 0) {
                lifespans[guid]++;
                tracking[kept++] = guid;
            }
        }
        tracking.resize(kept);
    }
    return lifespans;
}

template void trace_n_components<2>(BmmData<2>&, int);
template void trace_n_components<3>(BmmData<3>&, int);
template void trace_n_components<2>(BmmData<2>&, int, const IdsByComponent&);
template void trace_n_components<3>(BmmData<3>&, int, const IdsByComponent&);
template void trace_assignment_history<2>(BmmData<2>&, int);
template void trace_assignment_history<3>(BmmData<3>&, int);

} // namespace baysor
