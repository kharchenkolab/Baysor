// The parallel E-step apply with prior segments must leave every
// component's per-segment map exactly as the serial BmmData::assign() loop
// does (same contents and same iteration order, which drives the
// main-segment tie-breaks).
#include <gtest/gtest.h>

#include "baysor/processing/bmm_algorithm/bmm_algorithm.h"
#include "baysor/utils/general.h"
#include "baysor/utils/thread_pool.h"

#include <Eigen/Dense>
#include <random>
#include <utility>
#include <vector>

namespace {

using baysor::BmmData;

class PoolSizeGuard {
public:
    explicit PoolSizeGuard(int n) : old_(baysor::thread_pool_size()) {
        baysor::set_thread_pool_size(n);
    }
    ~PoolSizeGuard() { baysor::set_thread_pool_size(old_); }

private:
    int old_;
};

// Molecules on a jittered grid, a ring graph plus grid-neighbour edges,
// cells on a coarse grid, prior segments that disagree with the cells, and
// low confidence: many molecules change cell (or go to noise) per E-step.
BmmData<2> noisy_data_with_segments() {
    const int side = 60, n = side * side, n_genes = 5;
    std::mt19937 rng(17);
    std::uniform_real_distribution<double> jitter(-0.3, 0.3);
    std::uniform_int_distribution<int> gene(0, n_genes - 1);

    BmmData<2> data;
    data.position_data.resize(2, n);
    data.composition_data.resize(n);
    data.confidence.assign(n, 0.7);
    data.segment_per_molecule.resize(n);
    for (int i = 0; i < n; ++i) {
        const int x = i % side, y = i / side;
        data.position_data(0, i) = x + jitter(rng);
        data.position_data(1, i) = y + jitter(rng);
        data.composition_data[i] = gene(rng);
        // Segments: diagonal stripes, some molecules without a segment.
        data.segment_per_molecule[i] = ((x + y) % 7 == 0) ? 0 : 1 + (x + 2 * y) / 9 % 13;
    }

    std::vector<int> src, dst;
    std::vector<double> wt;
    for (int i = 0; i < n; ++i) {
        const int x = i % side;
        if (x + 1 < side) { src.push_back(i); dst.push_back(i + 1); wt.push_back(1.0); }
        if (i + side < n) { src.push_back(i); dst.push_back(i + side); wt.push_back(0.8); }
    }
    data.adj_list = baysor::AdjList::from_edge_list(src.data(), dst.data(), wt.data(),
                                                     static_cast<int>(src.size()), n);

    const int cells_side = 6;
    const double step = static_cast<double>(side) / cells_side;
    for (int cy = 0; cy < cells_side; ++cy) {
        for (int cx = 0; cx < cells_side; ++cx) {
            Eigen::Vector2d center((cx + 0.5) * step, (cy + 0.5) * step);
            baysor::MvNormal<2> pos(center, Eigen::Matrix2d::Identity() * step * step / 4.0);
            baysor::CategoricalSmoothed comp(n_genes, 1.0);
            comp.set_uniform_counts(1.0f);
            data.components.emplace_back(pos, comp, std::nullopt,
                                         static_cast<int>(data.components.size()) + 1);
        }
    }
    data.assignment.resize(n);
    for (int i = 0; i < n; ++i) {
        const int cx = std::min(cells_side - 1, static_cast<int>((i % side) / step));
        const int cy = std::min(cells_side - 1, static_cast<int>((i / side) / step));
        data.assignment[i] = 1 + cy * cells_side + cx;
    }
    for (auto& c : data.components) c.n_samples = 100;
    data.max_component_guid = static_cast<int>(data.components.size());

    int max_seg = 0;
    for (int s : data.segment_per_molecule) max_seg = std::max(max_seg, s);
    data.n_molecules_per_segment = baysor::count_array(data.segment_per_molecule, max_seg, true);
    data.noise_position_density = 1e-3;
    data.noise_density = 1e-3;
    data.prior_seg_confidence = 0.5;
    data.update_n_mols_per_segment();
    return data;
}

std::vector<std::vector<std::pair<int, int>>> maps_in_iteration_order(const BmmData<2>& data) {
    std::vector<std::vector<std::pair<int, int>>> out;
    for (const auto& c : data.components) {
        out.emplace_back(c.n_molecules_per_segment.begin(), c.n_molecules_per_segment.end());
    }
    return out;
}

} // namespace

TEST(BmmParallelApply, SegmentMapsMatchSerialAssignLoop) {
    for (int n_threads : {1, 3, 8}) {
        PoolSizeGuard guard(n_threads);
        for (std::uint64_t salt : {1u, 2u, 3u}) {
            auto data = noisy_data_with_segments();
            auto reference = data;  // deep copy, maps included

            auto stats = baysor::expect_dirichlet_spatial(data, /*stochastic=*/true, salt);

            // Replay the result through the serial assign() loop.
            std::int64_t n_changed = 0;
            for (int i = 0; i < data.n_molecules(); ++i) {
                if (reference.assignment[i] != data.assignment[i]) ++n_changed;
                reference.assign(i, data.assignment[i]);
            }
            ASSERT_GT(n_changed, data.n_molecules() / 50) << "fixture must produce many changes";
            EXPECT_EQ(stats.n_changed, n_changed);
            EXPECT_EQ(maps_in_iteration_order(data), maps_in_iteration_order(reference))
                << "threads " << n_threads << " salt " << salt;
        }
    }
}
