#include "baysor/processing/bmm_algorithm/tracing.h"

#include <algorithm>
#include <cmath>
#include <set>

namespace baysor {

template<int N>
void trace_n_components(BmmData<N>& data, int min_molecules_per_cell,
                        const IdsByComponent& ids_by_comp) {
    // Compute unique thresholds: {max(round(0.5*min),1), max(min,1), max(2*min,1), max(5*min,1)}
    std::set<int> thresh_set;
    for (double mult : {0.5, 1.0, 2.0, 5.0}) {
        thresh_set.insert(std::max(static_cast<int>(std::round(mult * min_molecules_per_cell)), 1));
    }

    std::unordered_map<int, int> entry;
    for (int t : thresh_set) {
        int cnt = 0;
        for (int c = 0; c < data.n_components(); ++c) {
            if (ids_by_comp.size(c) >= t) cnt++;
        }
        entry[t] = cnt;
    }

    data.n_components_trace.push_back(std::move(entry));
}

template<int N>
void trace_assignment_history(BmmData<N>& data, int assignment_history_depth) {
    if (assignment_history_depth <= 0) return;

    // Global assignment: local 1-based IDs replaced with component GUIDs,
    // pushed after trimming the history to depth - 1 entries
    const int n = data.n_molecules();
    auto& history = data.assignment_history;
    history.begin_push(n, static_cast<size_t>(assignment_history_depth));
    parallel_region([&](ParallelRegion& region) {
        history.push_rows(region, n, [&](int i) {
            const int a = data.assignment[i];
            return (a > 0) ? data.components[a - 1].guid : 0;
        });
    });
    history.end_push();
}

std::unordered_map<int, int> estimate_component_lifespan(
    const AssignmentHistory& history
) {
    std::unordered_map<int, int> lifespans;
    if (history.empty()) return lifespans;

    // Walk the entries backward from the newest one, keeping the current row
    // and the number of molecules per guid up to date with the deltas.
    std::vector<int> cur = history.back();
    int max_guid = 0;
    for (int g : cur) max_guid = std::max(max_guid, g);
    for (std::size_t t = 0; t + 1 < history.size(); ++t) {
        for (const auto& c : history.changes_after(t)) max_guid = std::max(max_guid, c.value);
    }
    std::vector<int> count(static_cast<size_t>(max_guid) + 1, 0);
    std::vector<int> tracking;  // guids whose streak is still unbroken
    for (int guid : cur) {
        if (count[guid]++ == 0 && guid > 0) {
            lifespans[guid] = 1;
            tracking.push_back(guid);
        }
    }

    // Extend the lifespan only while the consecutive streak is unbroken
    for (int t = static_cast<int>(history.size()) - 2; t >= 0 && !tracking.empty(); --t) {
        for (const auto& c : history.changes_after(t)) {
            --count[cur[c.mol]];
            ++count[c.value];
            cur[c.mol] = c.value;
        }
        tracking.erase(std::remove_if(tracking.begin(), tracking.end(),
                                      [&](int guid) { return count[guid] == 0; }),
                       tracking.end());
        for (int guid : tracking) lifespans[guid]++;
    }
    return lifespans;
}

template void trace_n_components<2>(BmmData<2>&, int, const IdsByComponent&);
template void trace_n_components<3>(BmmData<3>&, int, const IdsByComponent&);
template void trace_assignment_history<2>(BmmData<2>&, int);
template void trace_assignment_history<3>(BmmData<3>&, int);

} // namespace baysor
