#pragma once

#include <Eigen/Dense>

namespace baysor {

/// Compute 2D convex hull of a set of points.
/// Returns a matrix where each column is a hull vertex in clockwise order
/// (signed shoelace area is negative), matching Julia's `convex_hull.jl`.
Eigen::MatrixXd convex_hull(const Eigen::MatrixXd& points);

/// Compute area of a 2D polygon (columns = vertices)
double polygon_area(const Eigen::MatrixXd& polygon);

} // namespace baysor
