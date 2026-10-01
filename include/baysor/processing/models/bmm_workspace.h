#pragma once

#include "baysor/utils/julia_int_dict.h"

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
    std::vector<int> group_fill;

    // Connected-component split: position of each molecule inside its
    // cell's id list (molecule-indexed, written per cell).
    std::vector<int> mol_pos;
};

} // namespace baysor
