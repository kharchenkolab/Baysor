#pragma once

#include "baysor/utils/julia_int_dict.h"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace baysor {

/// Molecule ids grouped by component in CSR form: the ids of component c
/// (0-based) are ids[offsets[c] .. offsets[c + 1]), in ascending order (the
/// same content and order as split_ids(assignment, nc, drop_zero=true)).
/// Molecules with assignment 0 (noise) are not stored.
struct IdsByComponent {
    std::vector<int> offsets;  // size n_groups + 1
    std::vector<int> ids;

    int n_groups() const { return offsets.empty() ? 0 : static_cast<int>(offsets.size()) - 1; }
    const int* begin(int c) const { return ids.data() + offsets[c]; }
    int size(int c) const { return offsets[c + 1] - offsets[c]; }
};

/// Scratch state of the BMM loop, kept across iterations so that the E-step,
/// the M-step and the connected-component split do not allocate per call.
/// Holds no algorithm state: every buffer is (re)initialized before use.
struct BmmWorkspace {
    // E-step: Jacobi target and per-worker candidate buffers.
    std::vector<int> new_assignment;
    std::vector<JuliaIntDoubleDict>  component_weights;
    std::vector<std::vector<int>>    adj_classes;
    std::vector<std::vector<double>> adj_weights;
    std::vector<std::vector<double>> denses;

    // Grouping of molecules by component. After maximize() it is the M-step
    // grouping of the current assignment, which stays valid until the next
    // E-step changes the assignment.
    IdsByComponent ids_by_comp;
    std::vector<int> group_hist;      // per-worker histograms of the counting sort

    // Applying the E-step result
    std::vector<std::int64_t> worker_count;
    std::int64_t n_changed = 0;

    // M-step: per-worker arenas for the cluster-mode count maps
    std::vector<std::vector<std::byte>> cluster_mode_arena;

    // drop_unused_components
    std::vector<int> id_map;
    int n_kept = 0;

    // Connected-component split: position of each molecule inside its
    // cell's id list (molecule-indexed, written per cell), and per-worker
    // BFS scratch plus the molecules to reset to noise.
    std::vector<int> mol_pos;
    struct SplitScratch {
        std::vector<int> label;
        std::vector<int> queue;
        std::vector<int> cc_size;
        std::vector<int> dropped;
    };
    std::vector<SplitScratch> split_scratch;
};

} // namespace baysor
