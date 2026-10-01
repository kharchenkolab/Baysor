#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <Eigen/Dense>

namespace baysor {

/// Polygon collection: cell name/ID -> Nx2 matrix of boundary vertices
using PolygonCollection = std::unordered_map<std::string, Eigen::MatrixXd>;
using PolygonStack = std::vector<std::pair<std::string, PolygonCollection>>;

/// Estimate boundary polygons for all cells from molecule positions and assignments
PolygonCollection boundary_polygons(
    const Eigen::MatrixXd& pos_data,       // dims x n_molecules
    const std::vector<int>& cell_labels,
    const std::vector<std::string>* cell_names = nullptr,
    double offset_rel = 0.01
);

/// Estimate polygons from a grid-based segmentation mask
std::vector<Eigen::MatrixXd> boundary_polygons_from_grid(
    const Eigen::Matrix<uint32_t, Eigen::Dynamic, Eigen::Dynamic>& grid_labels,
    int grid_step = 1
);

/// Auto boundary estimation with optional per-z-slice polygons for 3D data.
/// When the data has more distinct z values than `max_z_slices`, the z-stack
/// is binned into `max_z_slices` layers (`max_z_slices` must be >= 1).
std::pair<PolygonCollection, PolygonStack> boundary_polygons_auto(
    const Eigen::MatrixXd& pos_data,
    const std::vector<int>& assignment,
    bool estimate_per_z = false,
    const std::vector<std::string>* cell_names = nullptr,
    bool verbose = true,
    int max_z_slices = 10
);

} // namespace baysor
