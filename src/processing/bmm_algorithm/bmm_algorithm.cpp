#include "baysor/processing/bmm_algorithm/bmm_algorithm.h"
#include "baysor/processing/bmm_algorithm/tracing.h"
#include "baysor/utils/general.h"
#include "baysor/utils/julia_int_dict.h"
#include "baysor/utils/thread_pool.h"

#include <spdlog/spdlog.h>
#include <spdlog/fmt/fmt.h>

#include <algorithm>
#include <optional>
#include <sstream>
#include <cmath>
#include <numeric>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace baysor {

// ============================================================================
// Internal helpers
// ============================================================================

// noise_composition_density: mean of 1/n_genes over components with non-trivial
// composition mass, matching Julia's noise_composition_density.
template<int N>
static double noise_composition_density(const BmmData<N>& data) {
    if (data.components.empty()) return 0.0;
    double acc = 0.0;
    double n_comps = 0.0;
    for (const auto& comp : data.components) {
        if (comp.composition_params.sum_counts <= 1e-3) continue;
        acc += 1.0 / std::max(comp.composition_params.n_genes, 1);
        n_comps += 1.0;
    }
    return acc / std::max(n_comps, 1.0);
}

// ============================================================================
// maximize (M-step)
// ============================================================================

template<int N>
void maximize(BmmData<N>& data, bool freeze_composition, bool freeze_position) {
    int nc = data.n_components();
    int n  = data.n_molecules();

    // Group molecule indices by component (1-based → 0-based component index)
    auto ids_by_comp = split_ids(data.assignment, nc, /*drop_zero=*/true);

    // Resize cluster_per_cell if needed
    if (!data.cluster_per_molecule.empty()) {
        data.cluster_per_cell.assign(nc, 0);
    }

    // Small chunks of a few cells: per-cell work varies widely with the cell
    // size, and preemption on shared hosts is the dominant source of imbalance.
    parallel_for(0, nc, 2, [&](int ci) {
        const auto& mol_ids = ids_by_comp[ci];
        int np = static_cast<int>(mol_ids.size());

        const std::vector<double>* nuc_probs =
            data.nuclei_prob_per_molecule.empty() ? nullptr : &data.nuclei_prob_per_molecule;

        data.components[ci].maximize_indexed(
            data.position_data,
            data.composition_data,
            mol_ids.data(), np,
            nuc_probs,
            data.min_nuclei_frac,
            freeze_position,
            freeze_composition
        );

        // cluster_per_cell: mode of cluster_per_molecule for this component
        if (!data.cluster_per_molecule.empty() && np > 0) {
            std::unordered_map<int, int> cnt;
            for (int k = 0; k < np; ++k) {
                cnt[data.cluster_per_molecule[mol_ids[k]]]++;
            }
            int best = 0, best_cnt = 0;
            for (auto& [cl, c] : cnt) {
                if (c > best_cnt) { best_cnt = c; best = cl; }
            }
            data.cluster_per_cell[ci] = best;
        }
    });

    data.noise_density = data.noise_position_density * noise_composition_density(data);
    if (std::isinf(data.noise_density) || std::isnan(data.noise_density)) {
        spdlog::warn("Unexpected noise density: {}", data.noise_density);
        data.noise_density = 0.0;
    }
}

// ============================================================================
// adjust_densities_by_prior_segmentation (E-step helper)
// Port of Julia lines 81-102 in bmm_algorithm.jl
// ============================================================================

template<int N>
static void adjust_densities_by_prior_segmentation(
    std::vector<double>& denses,
    const std::vector<int>& adj_classes,
    int segment_id,
    int largest_cell_id,      // 1-based component ID
    const BmmData<N>& data
) {
    double psc = data.prior_seg_confidence;
    double seg_prior_pow = psc * std::exp(3.0 * psc);
    int seg_size = data.n_molecules_per_segment[segment_id - 1];  // segment_id is 1-based

    // largest_cell_size: min(mols of largest_cell in this segment + 1, seg_size)
    int lc_idx = largest_cell_id - 1;
    int lc_in_seg = 0;
    {
        auto it = data.components[lc_idx].n_molecules_per_segment.find(segment_id);
        if (it != data.components[lc_idx].n_molecules_per_segment.end()) {
            lc_in_seg = it->second;
        }
    }
    int largest_cell_size = std::min(lc_in_seg + 1, seg_size);

    for (int j = 0; j < static_cast<int>(adj_classes.size()); ++j) {
        int c_adj = adj_classes[j];  // 1-based
        if (c_adj == largest_cell_id) continue;

        int c_idx = c_adj - 1;
        int main_seg = (c_idx < static_cast<int>(data.main_segment_per_cell.size()))
                       ? data.main_segment_per_cell[c_idx] : 0;

        int n_cell_mols = 0;
        {
            auto it = data.components[c_idx].n_molecules_per_segment.find(segment_id);
            if (it != data.components[c_idx].n_molecules_per_segment.end()) {
                n_cell_mols = it->second;
            }
        }
        int n_cell_mols_per_seg = std::min(n_cell_mols + 1, seg_size);

        if (main_seg == segment_id || main_seg == 0) {
            // Part against over-segmentation
            denses[j] *= std::pow(1.0 - psc, 0.5)
                         * std::pow(static_cast<double>(n_cell_mols_per_seg)
                                    / largest_cell_size, psc);
        } else {
            // Part against under-segmentation and overlap
            denses[j] *= std::pow(1.0 - psc, 0.5)
                         * std::pow(1.0 - static_cast<double>(n_cell_mols_per_seg)
                                    / seg_size, seg_prior_pow);
        }
    }
}

// ============================================================================
// expect_dirichlet_spatial (E-step)
// ============================================================================

// Seed for the per-chunk RNG stream of the multi-threaded E-step. Streams are
// keyed by (iteration, fixed-size chunk index), so the outcome is independent
// of scheduling and of the thread count, and differs between iterations.
static std::uint64_t estep_stream_seed(std::uint64_t rng_salt, std::int64_t chunk_idx) {
    auto mix = [](std::uint64_t x) {
        x += 0x9E3779B97F4A7C15ULL;
        x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
        x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
        return x ^ (x >> 31);
    };
    return mix(mix(rng_salt + 1) ^ mix(static_cast<std::uint64_t>(chunk_idx) + 0x100000001B3ULL));
}

// Fixed chunk size of the E-step loop. Chunk boundaries must not depend on
// the thread count: the multi-threaded RNG stream is keyed by chunk index.
// 256 balances per-chunk RNG setup against load balance on shared hosts.
static constexpr std::int64_t kEstepChunkSize = 256;

template<int N>
EstepStats expect_dirichlet_spatial(BmmData<N>& data, bool stochastic, std::uint64_t rng_salt) {
    int n = data.n_molecules();
    bool has_segments = !data.segment_per_molecule.empty();
    bool has_clusters  = !data.cluster_per_molecule.empty();

    // Jacobi snapshot: read from old_assignment during parallel phase
    std::vector<int> old_assignment = data.assignment;
    std::vector<int> new_assignment(n, 0);

    int n_threads = thread_pool_size();

    // Per-thread buffers
    std::vector<JuliaIntDoubleDict>  component_weights_buf(n_threads);
    std::vector<std::vector<int>>    adj_classes_buf(n_threads);
    std::vector<std::vector<double>> adj_weights_buf(n_threads);
    std::vector<std::vector<double>> denses_buf(n_threads);
    constexpr size_t reserve_hint = 64;
    for (int t = 0; t < n_threads; ++t) {
        adj_classes_buf[t].reserve(reserve_hint);
        adj_weights_buf[t].reserve(reserve_hint);
        denses_buf[t].reserve(reserve_hint);
    }

    // For single-thread parity, continue the same RNG stream used by earlier
    // preprocessing steps such as duplicate-point jitter in normalize_points.
    // Multi-threaded runs draw from per-chunk streams keyed by
    // (rng_salt, chunk index), so results do not depend on scheduling or on
    // the number of threads.
    const bool single_threaded = (n_threads <= 1);

    run_parallel_chunks(0, n, kEstepChunkSize, Scheduling::Dynamic,
        [&](std::int64_t chunk_begin, std::int64_t chunk_end, int ti) {
        std::optional<Xoshiro256pp> chunk_rng;
        if (!single_threaded) {
            chunk_rng.emplace(estep_stream_seed(rng_salt, chunk_begin / kEstepChunkSize));
        }
        auto& comp_weights = component_weights_buf[ti];
        auto& adj_classes  = adj_classes_buf[ti];
        auto& adj_weights  = adj_weights_buf[ti];
        auto& denses       = denses_buf[ti];

        for (int mol_id = static_cast<int>(chunk_begin); mol_id < static_cast<int>(chunk_end); ++mol_id) {

        // ---- aggregate_adjacent_component_weights ----
        comp_weights.clear();
        adj_classes.clear();
        adj_weights.clear();

        int  nc_adj = data.adj_list.neighbor_count(mol_id);
        const int32_t* nb_ids  = data.adj_list.neighbor_ids(mol_id);
        const double*  nb_wts  = data.adj_list.neighbor_weights(mol_id);

        double bg_comp_weight = 0.0;
        for (int ai = 0; ai < nc_adj; ++ai) {
            int nb   = nb_ids[ai];
            int c_id = old_assignment[nb];
            double cw = nb_wts[ai];
            if (c_id == 0) {
                bg_comp_weight += cw;
            } else {
                comp_weights.add(c_id, cw);
            }
        }
        comp_weights.for_each([&](int c_id, double cw) {
            adj_classes.push_back(c_id);
            adj_weights.push_back(cw);
        });

        int n_adj = static_cast<int>(adj_classes.size());
        if (n_adj == 0 && data.confidence[mol_id] >= 1.0) {
            // No adjacent cells and high confidence → stays noise
            new_assignment[mol_id] = 0;
            continue;
        }

        // ---- expect_density_for_molecule ----
        const double* x = data.position_data.col(mol_id).data();
        int gene         = data.composition_data[mol_id];
        double conf      = data.confidence[mol_id];
        int mol_cluster  = has_clusters ? data.cluster_per_molecule[mol_id] : -1;
        int segment_id   = has_segments ? data.segment_per_molecule[mol_id] : 0;

        denses.resize(n_adj);
        int largest_cell_id   = 0;
        int largest_cell_size = 0;

        for (int j = 0; j < n_adj; ++j) {
            int c_adj = adj_classes[j];    // 1-based
            int c_idx = c_adj - 1;         // 0-based
            const auto& comp = data.components[c_idx];

            double c_dens = conf
                * std::exp(data.mrf_strength * adj_weights[j])
                * comp.pdf(x, gene, data.use_gene_smoothing);

            // TODO(parity): Julia currently uses `< length(cluster_per_cell)`
            // here, which skips the last component from this penalty path.
            // That looks like an indexing quirk rather than intended model
            // behavior, but we keep it for parity for now.
            // Cluster penalty
            if (has_clusters
                && c_adj > 0
                && c_adj < static_cast<int>(data.cluster_per_cell.size())
                && data.cluster_per_cell[c_idx] != mol_cluster) {
                c_dens *= data.cluster_penalty_mult;
            }

            // Track largest cell for prior segmentation
            if (has_segments && segment_id > 0) {
                int main_seg = (c_idx < static_cast<int>(data.main_segment_per_cell.size()))
                               ? data.main_segment_per_cell[c_idx] : 0;
                if (main_seg == segment_id || main_seg == 0) {
                    int cur_mols = 0;
                    auto it = comp.n_molecules_per_segment.find(segment_id);
                    if (it != comp.n_molecules_per_segment.end()) cur_mols = it->second;
                    int seg_size = (segment_id <= static_cast<int>(data.n_molecules_per_segment.size()))
                                   ? data.n_molecules_per_segment[segment_id - 1] : 1;
                    int cur_mols_per_seg = std::min(cur_mols + 1, seg_size);

                    if (cur_mols_per_seg > largest_cell_size
                        || (cur_mols_per_seg == largest_cell_size
                            && comp.n_samples > (largest_cell_id > 0
                                ? data.components[largest_cell_id-1].n_samples : 0))) {
                        largest_cell_size = cur_mols_per_seg;
                        largest_cell_id   = c_adj;
                    }
                }
            }

            denses[j] = c_dens;
        }

        // Prior segmentation adjustment
        if (has_segments && segment_id > 0 && largest_cell_id > 0) {
            adjust_densities_by_prior_segmentation<N>(
                denses, adj_classes, segment_id, largest_cell_id, data);
        }

        // Noise term: only added when confidence < 1.0
        if (conf < 1.0) {
            denses.push_back(
                (1.0 - conf)
                * std::exp(data.mrf_strength * bg_comp_weight)
                * data.noise_density);
            adj_classes.push_back(0);
        }

        // ---- estimate_molecule_cell_assignment ----
        int n_total = static_cast<int>(adj_classes.size());
        double sum_d = 0.0;
        for (double d : denses) sum_d += d;

        if (sum_d < 1e-100) {
            new_assignment[mol_id] = 0;
        } else if (!stochastic) {
            int best = 0;
            double best_d = -1.0;
            for (int j = 0; j < n_total; ++j) {
                if (denses[j] > best_d) { best_d = denses[j]; best = j; }
            }
            new_assignment[mol_id] = adj_classes[best];
        } else {
            if (single_threaded) {
                new_assignment[mol_id] =
                    fsample(adj_classes.data(), denses.data(), n_total, global_xoshiro_rng());
            } else {
                new_assignment[mol_id] = fsample(adj_classes.data(), denses.data(), n_total, *chunk_rng);
            }
        }
        }
    });

    // Apply new assignments sequentially (preserves segment bookkeeping correctness)
    std::int64_t n_changed = 0;
    for (int mol_id = 0; mol_id < n; ++mol_id) {
        if (data.assignment[mol_id] != new_assignment[mol_id]) ++n_changed;
        data.assign(mol_id, new_assignment[mol_id]);
    }
    return {n_changed};
}

// ============================================================================
// drop_unused_components
// ============================================================================

template<int N>
void drop_unused_components(BmmData<N>& data, int min_n_samples) {
    int nc = data.n_components();
    auto n_mols = data.num_molecules_per_cell();  // 0-indexed: n_mols[i] for component i

    // Build compact id_map: old 1-based → new 1-based (0 = dropped)
    std::vector<int> id_map(nc, 0);
    int new_idx = 0;
    for (int i = 0; i < nc; ++i) {
        if (n_mols[i] >= min_n_samples) {
            id_map[i] = ++new_idx;
        }
    }

    if (new_idx == nc) return;  // nothing to drop

    // Remap assignment
    for (int& a : data.assignment) {
        if (a > 0) a = id_map[a - 1];
    }

    // Compact components vector
    int wi = 0;
    for (int i = 0; i < nc; ++i) {
        if (id_map[i] > 0) {
            if (wi != i) data.components[wi] = std::move(data.components[i]);
            ++wi;
        }
    }
    data.components.erase(data.components.begin() + wi, data.components.end());

    // Resize bookkeeping arrays
    if (!data.main_segment_per_cell.empty()) {
        data.main_segment_per_cell.resize(wi);
    }
    if (!data.cluster_per_cell.empty()) {
        data.cluster_per_cell.resize(wi);
    }
}

// ============================================================================
// split_cells_by_connected_components
// ============================================================================

template<int N>
void split_cells_by_connected_components(BmmData<N>& data) {
    int nc = data.n_components();
    if (nc == 0) return;

    const std::vector<int> assignment_snapshot = data.assignment;
    auto ids_per_cell = split_ids(data.assignment, nc, /*drop_zero=*/true);

    struct SplitScratch {
        std::unordered_map<int,int> mol_pos;
        std::vector<int> label;
        std::vector<int> queue;
        std::vector<int> cc_size;
    };
    std::vector<SplitScratch> scratch_buf(thread_pool_size());

    parallel_for(0, nc, 1, [&](int cell_id_0, int ti) {
        auto& scratch = scratch_buf[ti];
        auto& mol_pos = scratch.mol_pos;
        auto& label = scratch.label;
        auto& queue = scratch.queue;
        auto& cc_size = scratch.cc_size;

        const auto& mol_ids = ids_per_cell[cell_id_0];
            if (mol_ids.size() <= 1) return;

            int cell_id_1 = cell_id_0 + 1;  // 1-based

            // Build fast lookup: mol_id → position in mol_ids
            mol_pos.clear();
            mol_pos.reserve(mol_ids.size());
            for (int k = 0; k < static_cast<int>(mol_ids.size()); ++k) {
                mol_pos[mol_ids[k]] = k;
            }

            // BFS to find connected components among mol_ids
            int nm = static_cast<int>(mol_ids.size());
            label.assign(nm, -1);
            int n_cc = 0;

            for (int start = 0; start < nm; ++start) {
                if (label[start] >= 0) continue;
                label[start] = n_cc;
                queue.clear();
                queue.push_back(start);

                for (size_t head = 0; head < queue.size(); ++head) {
                    int k = queue[head];
                    int mol = mol_ids[k];
                    int nc_adj = data.adj_list.neighbor_count(mol);
                    const int32_t* nb_ids = data.adj_list.neighbor_ids(mol);
                    for (int ai = 0; ai < nc_adj; ++ai) {
                        int nb = nb_ids[ai];
                        if (assignment_snapshot[nb] != cell_id_1) continue;
                        auto it = mol_pos.find(nb);
                        if (it == mol_pos.end()) continue;
                        int nk = it->second;
                        if (label[nk] >= 0) continue;
                        label[nk] = n_cc;
                        queue.push_back(nk);
                    }
                }
                ++n_cc;
            }

            if (n_cc <= 1) return;

            // Find largest connected component
            cc_size.assign(n_cc, 0);
            for (int l : label) cc_size[l]++;
            int largest_cc = static_cast<int>(
                std::max_element(cc_size.begin(), cc_size.end()) - cc_size.begin());

            // Reassign non-largest CC molecules to noise. Each molecule appears
            // in only one cell group, so these writes are disjoint across threads.
            for (int k = 0; k < nm; ++k) {
                if (label[k] != largest_cc) {
                    data.assignment[mol_ids[k]] = 0;
                }
            }
    });
}

// ============================================================================
// estimate_assignment_by_history
// ============================================================================

template<int N>
std::pair<std::vector<int>, std::vector<double>>
estimate_assignment_by_history(const BmmData<N>& data) {
    int n = data.n_molecules();

    if (data.assignment_history.empty()) {
        // Fallback: return current assignment with 0.5 confidence
        std::vector<double> conf(n, 0.5);
        return {data.assignment, conf};
    }

    // Build guid → local 1-based ID map
    std::unordered_map<int,int> guid_map;
    for (int i = 0; i < data.n_components(); ++i) {
        guid_map[data.components[i].guid] = i + 1;
    }
    // current_guids includes 0 (noise)
    std::unordered_set<int> current_guids;
    for (auto& [g, _] : guid_map) current_guids.insert(g);
    current_guids.insert(0);

    int n_hist = static_cast<int>(data.assignment_history.size());

    std::vector<int> reassignment(n, 0);
    std::vector<double> match_frac(n, 0.0);

    for (int i = 0; i < n; ++i) {
        // Count frequency of each GUID across history (restricted to current_guids)
        std::unordered_map<int,int> freq;
        int valid = 0;
        for (int t = 0; t < n_hist; ++t) {
            int g = data.assignment_history[t][i];
            if (current_guids.count(g)) {
                freq[g]++;
                valid++;
            }
        }

        int best_guid = 0, best_cnt = 0;
        for (auto& [g, c] : freq) {
            if (c > best_cnt) { best_cnt = c; best_guid = g; }
        }

        reassignment[i] = (guid_map.count(best_guid)) ? guid_map[best_guid] : 0;
        match_frac[i] = (valid > 0) ? static_cast<double>(best_cnt) / valid : 0.0;
    }

    return {reassignment, match_frac};
}

// ============================================================================
// bmm — main EM loop
// ============================================================================

template<int N>
void bmm(BmmData<N>& data,
         int min_molecules_drop, int n_iters,
         int assignment_history_depth, bool verbose,
         int component_split_step, bool refine,
         bool freeze_composition, bool freeze_position, bool freeze_components,
         double tol, int min_molecules_display)
{
    // Display threshold: matches Julia's min_molecules_per_cell for the progress bar.
    // Drop threshold: matches Julia's hardcoded min_n_samples=2 in drop_unused_components!
    const int disp_thresh = (min_molecules_display > 0) ? min_molecules_display : min_molecules_drop;

    // Multi-level diagnostic thresholds (mirrors Julia's trace_nums):
    // [max(round(t * min_molecules_per_cell), 1) for t in [0.5, 1.0, 2.0, 5.0]]
    // but we use the display threshold as the "1.0" level.
    // Simplified: report at >=1, >=drop_thresh, >=disp_thresh (omit duplicates).
    // This gives comparable output to Julia's tracer thresholds.
    auto build_diag_str = [&]() -> std::string {
        auto n_mols_vec = data.num_molecules_per_cell();
        int n1 = 0, nd = 0, ndisp = 0;
        for (int m : n_mols_vec) {
            if (m >= 1)           ++n1;
            if (m >= min_molecules_drop)  ++nd;
            if (m >= disp_thresh) ++ndisp;
        }
        int n_noise = 0;
        for (int a : data.assignment) if (a == 0) ++n_noise;
        double noise_pct = 100.0 * n_noise / std::max(data.n_molecules(), 1);

        // Format: "noise=X%, total=N, >=drop=N, >=disp=N"
        // Omit redundant columns when thresholds coincide.
        std::string s = fmt::format("noise={:.1f}%", noise_pct);
        s += fmt::format(", total={}", n1);
        if (min_molecules_drop > 1)
            s += fmt::format(", >={}_segs={}", min_molecules_drop, nd);
        if (disp_thresh != min_molecules_drop)
            s += fmt::format(", >={}_cells={}", disp_thresh, ndisp);
        return s;
    };

    constexpr int n_iters_without_update = 20;
    std::vector<double> change_fracs;
    change_fracs.reserve(n_iters);

    // Initial trace + maximize to warm-start parameters
    trace_n_components(data, disp_thresh);
    maximize(data, freeze_composition, freeze_position);

    for (int iter = 1; iter <= n_iters; ++iter) {
        // Update prior probabilities: each component's prior = n_samples
        for (auto& comp : data.components) {
            comp.prior_probability = static_cast<double>(comp.n_samples);
        }

        data.update_n_mols_per_segment();

        // E-step — track assignment changes for convergence.
        // rng_salt = iteration index: multi-threaded draws come from per-chunk
        // streams keyed by (iteration, chunk), so different every iteration.
        EstepStats estep_stats = expect_dirichlet_spatial(data, /*stochastic=*/true,
                                                          /*rng_salt=*/static_cast<std::uint64_t>(iter));

        // Compute fraction of changed assignments
        if (tol > 0.0) {
            int n_total   = data.n_molecules();
            change_fracs.push_back(
                static_cast<double>(estep_stats.n_changed) / std::max(n_total, 1));
        }

        // Periodic connected component splitting
        if ((iter % component_split_step == 0) || (iter == n_iters)) {
            split_cells_by_connected_components(data);
        }

        // Drop components with fewer than min_molecules_drop molecules.
        // Matches Julia's drop_unused_components!(data) which has hardcoded default min_n_samples=2.
        if (!freeze_components) {
            drop_unused_components(data, min_molecules_drop);
        }

        // M-step
        maximize(data, freeze_composition, freeze_position);

        // Tracing
        trace_n_components(data, disp_thresh);
        trace_assignment_history(data, assignment_history_depth);

        if (verbose) {
            spdlog::info("Iter {:4d}/{}: {}", iter, n_iters, build_diag_str());
        }

        // Convergence check
        if (tol > 0.0 && iter >= n_iters_without_update) {
            int look_back = std::min(n_iters_without_update,
                                     static_cast<int>(change_fracs.size()));
            double worst = 0.0;
            for (int t = static_cast<int>(change_fracs.size()) - look_back;
                 t < static_cast<int>(change_fracs.size()); ++t) {
                if (change_fracs[t] > worst) worst = change_fracs[t];
            }
            if (worst < tol) {
                spdlog::info("Algorithm converged after {} iterations. Max change frac: {:.5f}",
                             iter, worst);
                break;
            }
        }
    }

    // Refinement phase
    if (refine) {
        if (!data.assignment_history.empty()) {
            auto [new_assign, conf] = estimate_assignment_by_history(data);
            data.assignment = std::move(new_assign);
            data.assignment_confidence = std::move(conf);
            maximize(data);
        }

        if (!freeze_components) {
            drop_unused_components(data, 1);  // keep any cell with >= 1 molecule (matches Julia)
        }
        maximize(data);

        if (verbose) {
            spdlog::info("Post-refine: {}", build_diag_str());
        }
    }
}

// Explicit instantiations
template void bmm<2>(BmmData<2>&, int, int, int, bool, int, bool, bool, bool, bool, double, int);
template void bmm<3>(BmmData<3>&, int, int, int, bool, int, bool, bool, bool, bool, double, int);
template EstepStats expect_dirichlet_spatial<2>(BmmData<2>&, bool, std::uint64_t);
template EstepStats expect_dirichlet_spatial<3>(BmmData<3>&, bool, std::uint64_t);
template void maximize<2>(BmmData<2>&, bool, bool);
template void maximize<3>(BmmData<3>&, bool, bool);
template void drop_unused_components<2>(BmmData<2>&, int);
template void drop_unused_components<3>(BmmData<3>&, int);
template void split_cells_by_connected_components<2>(BmmData<2>&);
template void split_cells_by_connected_components<3>(BmmData<3>&);
template std::pair<std::vector<int>, std::vector<double>> estimate_assignment_by_history<2>(const BmmData<2>&);
template std::pair<std::vector<int>, std::vector<double>> estimate_assignment_by_history<3>(const BmmData<3>&);

} // namespace baysor
