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

    // Build global assignment: replace local 1-based IDs with component GUIDs
    int n = data.n_molecules();
    // Recycle the buffer of the entry that is about to be evicted instead of
    // allocating a new n-int vector per iteration.
    std::vector<int> global;
    if (static_cast<int>(data.assignment_history.size()) >= assignment_history_depth) {
        global = std::move(data.assignment_history.front());
        data.assignment_history.erase(data.assignment_history.begin());
    }
    global.resize(n);
    for (int i = 0; i < n; ++i) {
        int a = data.assignment[i];
        global[i] = (a > 0) ? data.components[a - 1].guid : 0;
    }

    data.assignment_history.push_back(std::move(global));

    // Trim to depth
    while (static_cast<int>(data.assignment_history.size()) > assignment_history_depth) {
        data.assignment_history.erase(data.assignment_history.begin());
    }
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

template void trace_n_components<2>(BmmData<2>&, int);
template void trace_n_components<3>(BmmData<3>&, int);
template void trace_n_components<2>(BmmData<2>&, int, const IdsByComponent&);
template void trace_n_components<3>(BmmData<3>&, int, const IdsByComponent&);
template void trace_assignment_history<2>(BmmData<2>&, int);
template void trace_assignment_history<3>(BmmData<3>&, int);

} // namespace baysor
