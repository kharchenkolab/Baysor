#pragma once

// Internals of boundary_estimation.cpp, exposed for unit tests only.

#include <Eigen/Dense>
#include <array>
#include <utility>
#include <vector>

namespace baysor {
namespace internal {

/// Border edges of a triangulated cell after excluding, for up to max_iters
/// passes, the border triangles that contain admixture points.
std::vector<std::pair<int, int>> find_border_without_admixture(
    const std::vector<std::array<int, 3>>& triangles,
    const Eigen::MatrixXd& pos_data,
    const Eigen::MatrixXd& non_cell_pos,
    int max_iters = 100);

/// Vertex cycle walked along the border edges; {} if they do not close into
/// one within max_border_len steps.
std::vector<int> border_edges_to_poly(
    const std::vector<std::pair<int, int>>& border_edges,
    int max_border_len = 10000);

} // namespace internal
} // namespace baysor
