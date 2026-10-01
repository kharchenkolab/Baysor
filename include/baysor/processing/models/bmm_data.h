#pragma once

#include "baysor/processing/models/component.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/processing/models/assignment_history.h"
#include "baysor/processing/models/bmm_workspace.h"
#include <vector>
#include <string>
#include <any>
#include <unordered_map>
#include <Eigen/Dense>

namespace baysor {

/// Core data structure for the Bayesian Mixture Model.
/// Holds all molecule data, components (cells), and algorithm state.
///
/// Template parameter N = spatial dimensionality (2 or 3).
template<int N>
struct BmmData {
    // --- Static molecule data ---
    Eigen::MatrixXd position_data;                  // N x n_molecules
    std::vector<int> composition_data;              // gene IDs per molecule (or -1 for missing)
    std::vector<double> confidence;                 // per-molecule confidence
    std::vector<int> cluster_per_molecule;           // MRF cluster assignment (empty if unused)
    std::vector<int> segment_per_molecule;           // prior segmentation labels (empty if unused)
    std::vector<double> nuclei_prob_per_molecule;    // nuclei probability (empty if unused)

    // --- MRF graph ---
    AdjList adj_list;

    // --- Components (cells) ---
    std::vector<Component<N>> components;
    std::vector<int> assignment;                     // molecule -> component index (0 = noise)
    int max_component_guid = 0;

    // --- Noise model ---
    double noise_position_density = 0.0;
    double noise_density = 0.0;

    // --- Per-cell metadata ---
    std::vector<int> cluster_per_cell;

    // --- Prior segmentation bookkeeping ---
    std::vector<int> n_molecules_per_segment;
    std::vector<int> main_segment_per_cell;

    // --- Algorithm parameters ---
    double prior_seg_confidence = 0.5;
    double cluster_penalty_mult = 0.25;
    bool use_gene_smoothing = true;
    double min_nuclei_frac = 0.1;
    double mrf_strength = 0.1;
    double real_edge_weight = 1.0;

    // --- Tracing ---
    AssignmentHistory assignment_history;   // global GUIDs per molecule, delta-encoded
    std::vector<std::unordered_map<int, int>> n_components_trace;

    // --- Output: per-molecule assignment confidence (fraction of history agreeing with final) ---
    std::vector<double> assignment_confidence;

    // --- Scratch buffers of the BMM loop (no algorithm state) ---
    BmmWorkspace workspace;

    // --- Accessors ---
    int n_molecules() const { return static_cast<int>(position_data.cols()); }
    int n_components() const { return static_cast<int>(components.size()); }
    int n_genes() const { return components.empty() ? 0 : components[0].composition_params.size(); }

    /// Assign molecule to component, updating segment bookkeeping
    void assign(int mol_id, int component_id);

    /// Count molecules per cell
    std::vector<int> num_molecules_per_cell() const;

    /// Update n_molecules_per_segment and main_segment_per_cell
    void update_n_mols_per_segment();

    /// The result of update_n_mols_per_segment() for component ci (0-based),
    /// computed from `ids_by_comp`, the grouping of the current assignment
    /// (e.g. workspace.ids_by_comp right after an M-step). The map receives
    /// the operation sequence of the serial molecule loop (clear, then ++ in
    /// ascending molecule order), so its state, iteration order and the
    /// main-segment tie-break are the same. main_segment_per_cell must
    /// already have n_components() entries. Different components may be
    /// updated concurrently.
    void update_n_mols_per_segment_of(int ci, const IdsByComponent& ids_by_comp);

private:
    void update_main_segment(int ci);
public:
};

extern template struct BmmData<2>;
extern template struct BmmData<3>;

} // namespace baysor
