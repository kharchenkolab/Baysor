#include "baysor/processing/bmm_algorithm/bmm_algorithm.h"
#include "baysor/processing/bmm_algorithm/tracing.h"
#include "baysor/utils/general.h"
#include "baysor/utils/julia_int_dict.h"
#include "baysor/utils/thread_pool.h"

#include <spdlog/spdlog.h>
#include <spdlog/fmt/fmt.h>

#include <algorithm>
#include <cstddef>
#include <memory_resource>
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

// Group molecule ids by component into a reusable CSR (ids ascending inside
// each group, i.e. the same groups as split_ids(assignment, nc,
// drop_zero=true)) instead of nc freshly allocated vectors per call.
// Parallel counting sort: per-worker histograms over static blocks, prefix
// sums in (component, block) order, then each block scatters its molecules
// in ascending order. The result does not depend on the number of workers.
// Work-shared.
static void group_phase(ParallelRegion& region, const std::vector<int>& assignment, int nc,
                        IdsByComponent& out, std::vector<int>& hist) {
    const int n = static_cast<int>(assignment.size());
    const int n_workers = region.n_workers();
    region.single([&]() {
        hist.assign(static_cast<size_t>(n_workers) * nc, 0);
        out.offsets.resize(nc + 1);
    });
    region.for_chunks(0, n, 0, Scheduling::Static,
        [&](std::int64_t b, std::int64_t e, int w) {
        int* h = hist.data() + static_cast<size_t>(w) * nc;
        for (std::int64_t i = b; i < e; ++i) {
            const int f = assignment[i];
            if (f > 0 && f <= nc) h[f - 1]++;
        }
    });
    region.single([&]() {
        // hist[w][c] becomes the first output position of block w in group c.
        int pos = 0;
        for (int c = 0; c < nc; ++c) {
            out.offsets[c] = pos;
            for (int w = 0; w < n_workers; ++w) {
                int& h = hist[static_cast<size_t>(w) * nc + c];
                const int cnt = h;
                h = pos;
                pos += cnt;
            }
        }
        out.offsets[nc] = pos;
        out.ids.resize(pos);
    });
    region.for_chunks(0, n, 0, Scheduling::Static,
        [&](std::int64_t b, std::int64_t e, int w) {
        int* h = hist.data() + static_cast<size_t>(w) * nc;
        for (std::int64_t i = b; i < e; ++i) {
            const int f = assignment[i];
            if (f > 0 && f <= nc) out.ids[h[f - 1]++] = static_cast<int>(i);
        }
    });
}

// ============================================================================
// maximize (M-step)
// ============================================================================

// M-step from `ids_by_comp`, the grouping of the current assignment.
// Work-shared.
template<int N>
static void maximize_phase(ParallelRegion& region, BmmData<N>& data,
                           const IdsByComponent& ids_by_comp,
                           bool freeze_composition, bool freeze_position) {
    int nc = data.n_components();
    BmmWorkspace& ws = data.workspace;

    region.single([&]() {
        // Resize cluster_per_cell if needed
        if (!data.cluster_per_molecule.empty()) {
            data.cluster_per_cell.assign(nc, 0);
        }
        // Per-worker arena for the cluster-mode maps below
        if (static_cast<int>(ws.cluster_mode_arena.size()) != region.n_workers()) {
            ws.cluster_mode_arena.assign(region.n_workers(), std::vector<std::byte>(16384));
        }
    });

    // Small chunks of a few cells: per-cell work varies widely with the cell
    // size, and preemption on shared hosts is the dominant source of imbalance.
    region.for_each(0, nc, 2, [&](int ci, int w) {
        const int* mol_ids = ids_by_comp.begin(ci);
        int np = ids_by_comp.size(ci);

        const std::vector<double>* nuc_probs =
            data.nuclei_prob_per_molecule.empty() ? nullptr : &data.nuclei_prob_per_molecule;

        data.components[ci].maximize_indexed(
            data.position_data,
            data.composition_data,
            mol_ids, np,
            nuc_probs,
            data.min_nuclei_frac,
            freeze_position,
            freeze_composition
        );

        // cluster_per_cell: mode of cluster_per_molecule for this component
        if (!data.cluster_per_molecule.empty() && np > 0) {
            // A std::unordered_map on a per-worker arena instead of the heap:
            // libstdc++'s hashing, bucket policy and iteration order do not
            // depend on the allocator, so ties resolve exactly as before.
            auto& arena = ws.cluster_mode_arena[w];
            std::pmr::monotonic_buffer_resource res(arena.data(), arena.size());
            std::pmr::unordered_map<int, int> cnt(&res);
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

    region.single([&]() {
        data.noise_density = data.noise_position_density * noise_composition_density(data);
        if (std::isinf(data.noise_density) || std::isnan(data.noise_density)) {
            spdlog::warn("Unexpected noise density: {}", data.noise_density);
            data.noise_density = 0.0;
        }
    });
}

template<int N>
void maximize(BmmData<N>& data, bool freeze_composition, bool freeze_position) {
    BmmWorkspace& ws = data.workspace;
    parallel_region([&](ParallelRegion& region) {
        // Group molecule indices by component (1-based → 0-based component index)
        group_phase(region, data.assignment, data.n_components(), ws.ids_by_comp, ws.group_hist);
        maximize_phase(region, data, ws.ids_by_comp, freeze_composition, freeze_position);
    });
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
    const BmmData<N>& data,
    double seg_prior_pow,     // psc * exp(3 psc), loop-invariant
    double sqrt_one_m_psc     // pow(1 - psc, 0.5), loop-invariant
) {
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
            denses[j] *= sqrt_one_m_psc
                         * std::pow(static_cast<double>(n_cell_mols_per_seg)
                                    / largest_cell_size, data.prior_seg_confidence);
        } else {
            // Part against under-segmentation and overlap
            denses[j] *= sqrt_one_m_psc
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

// E-step over all molecules into workspace.new_assignment. Work-shared:
// every participant of the region must call it.
template<int N>
static void estep_phase(ParallelRegion& region, BmmData<N>& data, bool stochastic,
                        std::uint64_t rng_salt) {
    int n = data.n_molecules();
    bool has_segments = !data.segment_per_molecule.empty();
    bool has_clusters  = !data.cluster_per_molecule.empty();

    // Jacobi update: data.assignment is not written during the parallel
    // phase (results go to new_assignment), so it is read in place.
    const std::vector<int>& old_assignment = data.assignment;
    BmmWorkspace& ws = data.workspace;
    std::vector<int>& new_assignment = ws.new_assignment;

    // Per-worker buffers, kept across iterations
    auto& component_weights_buf = ws.component_weights;
    auto& adj_classes_buf = ws.adj_classes;
    auto& adj_weights_buf = ws.adj_weights;
    auto& denses_buf = ws.denses;
    region.single([&]() {
        // Every molecule is written below, so no initialization is needed.
        new_assignment.resize(n);
        const int n_workers = region.n_workers();
        if (static_cast<int>(adj_classes_buf.size()) != n_workers) {
            constexpr size_t reserve_hint = 64;
            adj_classes_buf.assign(n_workers, {});
            adj_weights_buf.assign(n_workers, {});
            denses_buf.assign(n_workers, {});
            for (int t = 0; t < n_workers; ++t) {
                adj_classes_buf[t].reserve(reserve_hint);
                adj_weights_buf[t].reserve(reserve_hint);
                denses_buf[t].reserve(reserve_hint);
            }
        }
        // Fresh 16-slot dictionaries per call: the table size persists after a
        // grow() and determines the slot order, so it must not leak from one
        // E-step call into the next.
        component_weights_buf.assign(n_workers, JuliaIntDoubleDict());
    });

    // Loop-invariant factors of adjust_densities_by_prior_segmentation
    const double psc = data.prior_seg_confidence;
    const double seg_prior_pow = psc * std::exp(3.0 * psc);
    const double sqrt_one_m_psc = std::pow(1.0 - psc, 0.5);

    // For single-thread parity, continue the same RNG stream used by earlier
    // preprocessing steps such as duplicate-point jitter in normalize_points.
    // Multi-threaded runs draw from per-chunk streams keyed by
    // (rng_salt, chunk index), so results do not depend on scheduling or on
    // the number of threads.
    const bool single_threaded = (thread_pool_size() <= 1);

    region.for_chunks(0, n, kEstepChunkSize, Scheduling::Dynamic,
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
        // Accumulate the neighbour weights per distinct class in first-seen
        // order (typically 1-3 classes, linear search). The Julia Dict slot
        // order depends only on the insertion order of distinct keys, and each
        // key's sum is formed in the same neighbour order, so feeding the
        // per-class sums into the dict in first-seen order reproduces the
        // dict's slot order and values bitwise. Single-class molecules skip
        // the dict entirely.
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
                const int m = static_cast<int>(adj_classes.size());
                int j = 0;
                while (j < m && adj_classes[j] != c_id) ++j;
                if (j == m) {
                    adj_classes.push_back(c_id);
                    adj_weights.push_back(cw);
                } else {
                    adj_weights[j] += cw;
                }
            }
        }
        if (adj_classes.size() > 1) {
            comp_weights.clear();
            for (size_t j = 0; j < adj_classes.size(); ++j) {
                comp_weights.add(adj_classes[j], adj_weights[j]);
            }
            adj_classes.clear();
            adj_weights.clear();
            comp_weights.for_each([&](int c_id, double cw) {
                adj_classes.push_back(c_id);
                adj_weights.push_back(cw);
            });
        }

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
                denses, adj_classes, segment_id, largest_cell_id, data,
                seg_prior_pow, sqrt_one_m_psc);
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
}

// Apply workspace.new_assignment; returns the number of changed molecules
// (the same value on every participant). Work-shared.
template<int N>
static std::int64_t apply_phase(ParallelRegion& region, BmmData<N>& data) {
    BmmWorkspace& ws = data.workspace;
    const int n = data.n_molecules();
    std::int64_t n_changed = 0;

    if (!data.segment_per_molecule.empty()) {
        // assign() updates the per-segment maps of the old and new component;
        // call it only for changed molecules (unchanged ones are no-ops in
        // assign), in ascending order, so every map receives the same
        // operation sequence as before.
        region.single([&]() {
            std::int64_t cnt = 0;
            for (int mol_id = 0; mol_id < n; ++mol_id) {
                const int a = ws.new_assignment[mol_id];
                if (data.assignment[mol_id] != a) {
                    ++cnt;
                    data.assign(mol_id, a);
                }
            }
            ws.n_changed = cnt;
        });
        return ws.n_changed;
    }

    // Without prior segments assign() only stores the value: count the
    // changes in parallel and swap the vectors.
    region.single([&]() { ws.worker_count.assign(region.n_workers(), 0); });
    region.for_chunks(0, n, 0, Scheduling::Static,
        [&](std::int64_t b, std::int64_t e, int w) {
        std::int64_t cnt = 0;
        for (std::int64_t i = b; i < e; ++i) {
            cnt += (data.assignment[i] != ws.new_assignment[i]);
        }
        ws.worker_count[w] = cnt;
    });
    region.single([&]() {
        std::int64_t cnt = 0;
        for (std::int64_t c : ws.worker_count) cnt += c;
        ws.n_changed = cnt;
        data.assignment.swap(ws.new_assignment);
    });
    n_changed = ws.n_changed;
    return n_changed;
}

template<int N>
EstepStats expect_dirichlet_spatial(BmmData<N>& data, bool stochastic, std::uint64_t rng_salt) {
    std::int64_t n_changed = 0;
    parallel_region([&](ParallelRegion& region) {
        estep_phase(region, data, stochastic, rng_salt);
        const std::int64_t c = apply_phase(region, data);
        if (region.is_master()) n_changed = c;
    });
    return {n_changed};
}

// ============================================================================
// drop_unused_components
// ============================================================================

// Drop components with fewer than min_n_samples molecules. `ids_by_comp` must
// be the grouping of the current assignment; it is rebuilt when components
// are dropped, so on return it is again the grouping of the (remapped)
// assignment. Work-shared.
template<int N>
static void drop_phase(ParallelRegion& region, BmmData<N>& data, int min_n_samples,
                       IdsByComponent& ids_by_comp) {
    BmmWorkspace& ws = data.workspace;
    const int nc = data.n_components();
    region.single([&]() {
        // Build compact id_map: old 1-based → new 1-based (0 = dropped)
        ws.id_map.assign(nc, 0);
        int new_idx = 0;
        for (int i = 0; i < nc; ++i) {
            if (ids_by_comp.size(i) >= min_n_samples) {
                ws.id_map[i] = ++new_idx;
            }
        }
        ws.n_kept = new_idx;
    });
    if (ws.n_kept == nc) return;  // nothing to drop (read after the barrier)

    // Remap assignment
    region.for_chunks(0, data.n_molecules(), 0, Scheduling::Static,
        [&](std::int64_t b, std::int64_t e, int) {
        for (std::int64_t i = b; i < e; ++i) {
            int& a = data.assignment[i];
            if (a > 0) a = ws.id_map[a - 1];
        }
    });

    region.single([&]() {
        // Compact components vector
        int wi = 0;
        for (int i = 0; i < nc; ++i) {
            if (ws.id_map[i] > 0) {
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
    });

    group_phase(region, data.assignment, data.n_components(), ids_by_comp, ws.group_hist);
}

template<int N>
void drop_unused_components(BmmData<N>& data, int min_n_samples) {
    BmmWorkspace& ws = data.workspace;
    parallel_region([&](ParallelRegion& region) {
        group_phase(region, data.assignment, data.n_components(), ws.ids_by_comp, ws.group_hist);
        drop_phase(region, data, min_n_samples, ws.ids_by_comp);
    });
}

// ============================================================================
// split_cells_by_connected_components
// ============================================================================

// Reassign all but the largest connected component of every cell to noise.
// `ids_per_cell` must be the grouping of the current assignment; on return
// it is stale. Work-shared.
template<int N>
static void split_phase(ParallelRegion& region, BmmData<N>& data,
                        const IdsByComponent& ids_per_cell) {
    const int nc = data.n_components();
    BmmWorkspace& ws = data.workspace;

    // Position of each molecule inside its cell's id list. Every molecule
    // belongs to one cell, so the writes of different cells are disjoint and
    // one shared array replaces a per-cell hash map. A neighbour passes the
    // assignment check below only if it is in the same cell, i.e. in
    // mol_ids, so the lookup always hits.
    std::vector<int>& mol_pos = ws.mol_pos;
    region.single([&]() {
        mol_pos.resize(data.n_molecules());
        if (static_cast<int>(ws.split_scratch.size()) != region.n_workers()) {
            ws.split_scratch.assign(region.n_workers(), {});
        }
        for (auto& sc : ws.split_scratch) sc.dropped.clear();
    });

    // data.assignment is only read in this loop: the molecules to drop are
    // collected per worker and reset after it.
    region.for_each(0, nc, 1, [&](int cell_id_0, int ti) {
        auto& scratch = ws.split_scratch[ti];
        auto& label = scratch.label;
        auto& queue = scratch.queue;
        auto& cc_size = scratch.cc_size;

        const int* mol_ids = ids_per_cell.begin(cell_id_0);
        const int nm = ids_per_cell.size(cell_id_0);
        if (nm <= 1) return;

        int cell_id_1 = cell_id_0 + 1;  // 1-based

        // Fast lookup: mol_id → position in mol_ids
        for (int k = 0; k < nm; ++k) {
            mol_pos[mol_ids[k]] = k;
        }

        // BFS to find connected components among mol_ids
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
                    if (data.assignment[nb] != cell_id_1) continue;
                    int nk = mol_pos[nb];
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

        // Molecules of the non-largest components go to noise
        for (int k = 0; k < nm; ++k) {
            if (label[k] != largest_cc) {
                scratch.dropped.push_back(mol_ids[k]);
            }
        }
    });

    region.for_each(0, region.n_workers(), 1, [&](int w) {
        for (int mol : ws.split_scratch[w].dropped) data.assignment[mol] = 0;
    });
}

template<int N>
void split_cells_by_connected_components(BmmData<N>& data) {
    if (data.n_components() == 0) return;
    BmmWorkspace& ws = data.workspace;
    parallel_region([&](ParallelRegion& region) {
        group_phase(region, data.assignment, data.n_components(), ws.ids_by_comp, ws.group_hist);
        split_phase(region, data, ws.ids_by_comp);
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

    const AssignmentHistory& history = data.assignment_history;
    const int n_hist = static_cast<int>(history.size());

    std::vector<int> reassignment(n, 0);
    std::vector<double> match_frac(n, 0.0);

    // Parallel over fixed blocks of molecules. Each block reconstructs its
    // molecules' rows from the newest entry and the backward deltas (sorted
    // by molecule, so one cursor per delta advances monotonically). The vote
    // per molecule is unchanged: a fresh hash map fed in history order (here
    // on a per-block arena; libstdc++'s iteration order does not depend on
    // the allocator), so ties resolve as before.
    constexpr std::int64_t kBlock = 4096;
    parallel_for(0, (static_cast<std::int64_t>(n) + kBlock - 1) / kBlock, 1, [&](std::int64_t blk) {
        const int b = static_cast<int>(blk * kBlock);
        const int e = static_cast<int>(std::min<std::int64_t>(n, b + kBlock));

        std::vector<size_t> cursor(std::max(n_hist - 1, 0));
        for (int t = 0; t + 1 < n_hist; ++t) {
            const auto& d = history.changes_after(t);
            cursor[t] = static_cast<size_t>(std::lower_bound(d.begin(), d.end(), b,
                [](const AssignmentHistory::Change& c, int mol) { return c.mol < mol; }) - d.begin());
        }
        std::vector<int> row(n_hist);
        std::vector<std::byte> arena(4096);

        for (int i = b; i < e; ++i) {
            row[n_hist - 1] = history.newest()[i];
            for (int t = n_hist - 2; t >= 0; --t) {
                const auto& d = history.changes_after(t);
                size_t& c = cursor[t];
                if (c < d.size() && d[c].mol == i) {
                    row[t] = d[c].value;
                    ++c;
                } else {
                    row[t] = row[t + 1];
                }
            }

            // Count frequency of each GUID across history (restricted to current_guids)
            std::pmr::monotonic_buffer_resource res(arena.data(), arena.size());
            std::pmr::unordered_map<int, int> freq(&res);
            int valid = 0;
            for (int t = 0; t < n_hist; ++t) {
                int g = row[t];
                if (current_guids.count(g)) {
                    freq[g]++;
                    valid++;
                }
            }

            int best_guid = 0, best_cnt = 0;
            for (auto& [g, c] : freq) {
                if (c > best_cnt) { best_cnt = c; best_guid = g; }
            }

            auto it = guid_map.find(best_guid);
            reassignment[i] = (it != guid_map.end()) ? it->second : 0;
            match_frac[i] = (valid > 0) ? static_cast<double>(best_cnt) / valid : 0.0;
        }
    });

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
    // Called right after an M-step: the per-cell counts come from its grouping.
    const IdsByComponent& groups = data.workspace.ids_by_comp;
    auto build_diag_str = [&]() -> std::string {
        int n1 = 0, nd = 0, ndisp = 0;
        for (int c = 0; c < data.n_components(); ++c) {
            const int m = groups.size(c);
            if (m >= 1)           ++n1;
            if (m >= min_molecules_drop)  ++nd;
            if (m >= disp_thresh) ++ndisp;
        }
        const int n_noise = data.n_molecules() - static_cast<int>(groups.ids.size());
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

        // The grouping of the last M-step is still the grouping of the
        // current assignment.
        data.update_n_mols_per_segment(groups);

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
        trace_n_components(data, disp_thresh, groups);
        // With tol == 0 all n_iters iterations run, so only the last
        // assignment_history_depth entries can survive the trimming.
        if (tol > 0.0 || iter > n_iters - assignment_history_depth) {
            trace_assignment_history(data, assignment_history_depth);
        }

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
