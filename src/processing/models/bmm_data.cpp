#include "baysor/processing/models/bmm_data.h"

namespace baysor {

template<int N>
void BmmData<N>::assign(int mol_id, int component_id) {
    int old_id = assignment[mol_id];
    if (old_id == component_id) return;

    // Update per-segment bookkeeping (only when prior segmentation is active)
    if (!segment_per_molecule.empty()) {
        int seg_id = segment_per_molecule[mol_id];
        if (seg_id > 0) {
            if (component_id > 0) {
                components[component_id - 1].n_molecules_per_segment[seg_id]++;
            }
            if (old_id > 0) {
                auto& seg_map = components[old_id - 1].n_molecules_per_segment;
                auto it = seg_map.find(seg_id);
                if (it != seg_map.end()) {
                    it->second--;
                    if (it->second <= 0) seg_map.erase(it);
                }
            }
        }
    }

    assignment[mol_id] = component_id;
}

template<int N>
void BmmData<N>::update_n_mols_per_segment() {
    if (segment_per_molecule.empty()) return;

    // Clear per-component segment counts
    for (auto& comp : components) {
        comp.n_molecules_per_segment.clear();
    }

    // Accumulate
    int n = n_molecules();
    for (int i = 0; i < n; ++i) {
        int c_cell = assignment[i];
        int c_seg  = segment_per_molecule[i];
        if (c_cell <= 0 || c_seg <= 0) continue;
        components[c_cell - 1].n_molecules_per_segment[c_seg]++;
    }

    // Compute main_segment_per_cell: segment with highest fraction of its molecules in this cell
    int nc = n_components();
    main_segment_per_cell.assign(nc, 0);
    for (int ci = 0; ci < nc; ++ci) {
        update_main_segment(ci);
    }
}

template<int N>
void BmmData<N>::update_n_mols_per_segment_of(int ci, const IdsByComponent& ids_by_comp) {
    auto& seg_map = components[ci].n_molecules_per_segment;
    seg_map.clear();
    const int* ids = ids_by_comp.begin(ci);
    const int np = ids_by_comp.size(ci);
    for (int k = 0; k < np; ++k) {
        int c_seg = segment_per_molecule[ids[k]];
        if (c_seg <= 0) continue;
        seg_map[c_seg]++;
    }
    main_segment_per_cell[ci] = 0;
    update_main_segment(ci);
}

template<int N>
void BmmData<N>::update_main_segment(int ci) {
    const auto& seg_map = components[ci].n_molecules_per_segment;
    if (seg_map.empty()) return;

    int   best_seg   = 0;
    double best_frac = 0.0;
    int   best_size  = 0;

    for (const auto& [si, nms] : seg_map) {
        if (si <= 0 || si > static_cast<int>(n_molecules_per_segment.size())) continue;
        int seg_size = n_molecules_per_segment[si - 1];
        if (seg_size <= 0) continue;
        double frac = static_cast<double>(nms) / seg_size;

        // Matches Julia: (frac > best_frac + 1e-10) or (nms == seg_size && seg_size > best_size)
        if (frac > best_frac + 1e-10
            || (nms == seg_size && seg_size > best_size)) {
            best_frac = frac;
            best_seg  = si;
            best_size = seg_size;
        }
    }
    main_segment_per_cell[ci] = best_seg;
}

template struct BmmData<2>;
template struct BmmData<3>;

} // namespace baysor
