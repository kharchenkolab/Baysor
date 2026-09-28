#pragma once

// Test-only seam for boundary-estimation helpers that normally live in an
// anonymous namespace inside boundary_estimation.cpp. NOT public API: it
// exists solely so unit tests can reach the parameterized convergence guards
// (max_iters / max_border_len) whose warning branches production callers
// never trigger; it may change or disappear with the internals it wraps.
//
// Production code calls the original implementations directly with the
// constants below, which are also the wrappers' defaults. Forwarding preserves
// identical behavior; no production behavior changes.

#include <array>
#include <utility>
#include <vector>
#include <Eigen/Dense>

namespace baysor {
namespace internal {

/// Iteration cap of the admixture filter in find_border_without_admixture.
inline constexpr int kBorderFilterMaxIters = 100;
/// Maximal number of edges walked when tracing one polygon border.
inline constexpr int kMaxBorderLength = 10000;

/// Forwards to the anonymous-namespace implementation used by
/// build_polygons_for_cells.
std::vector<std::pair<int, int>> find_border_without_admixture(
    const std::vector<std::array<int, 3>>& triangles,
    const Eigen::MatrixXd& pos_data,
    const Eigen::MatrixXd& non_cell_pos,
    int max_iters = kBorderFilterMaxIters);

/// Forwards to the anonymous-namespace implementation used by
/// build_polygons_for_cells and boundary_polygons_from_grid.
std::vector<int> border_edges_to_poly(
    const std::vector<std::pair<int, int>>& border_edges,
    int max_border_len = kMaxBorderLength);

} // namespace internal
} // namespace baysor
