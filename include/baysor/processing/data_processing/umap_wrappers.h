#pragma once

#include <Eigen/Dense>

namespace baysor {

/// Run UMAP on a column-major data matrix.
/// @param data       ndim x nobs matrix (column per observation)
/// @param ndim_out   Output dimensionality (2 or 3)
/// @param n_neighbors  Number of nearest neighbors for the UMAP graph
/// @param n_epochs   Number of optimization epochs
/// @param seed       Random seed
/// @returns ndim_out x nobs embedding matrix
Eigen::MatrixXd umap_embed(
    const Eigen::MatrixXd& data,
    int ndim_out,
    int n_neighbors = 15,
    int n_epochs    = 200,
    int seed        = 42,
    double spread   = 1.0,
    double min_dist = 0.1
);

/// Run UMAP on a precomputed symmetric distance matrix.
/// Extracts k-nearest neighbors per row and runs UMAP.
/// @param dist_mat   n x n symmetric distance matrix (diagonal = 0)
/// @param ndim_out   Output dimensionality
/// @param n_neighbors  Neighbors to extract from dist_mat
/// @param n_epochs   Optimization epochs
/// @param seed       Random seed
/// @returns ndim_out x n embedding matrix
Eigen::MatrixXd umap_embed_precomputed(
    const Eigen::MatrixXd& dist_mat,
    int ndim_out,
    int n_neighbors = 15,
    int n_epochs    = 500,
    int seed        = 42
);

} // namespace baysor
