using UMAP
using NearestNeighbors
using Statistics

import MultivariateStats
import MultivariateStats.transform
import MultivariateStats.fit

struct UmapFit
    embedding::Array;
    nn_tree;
    pca_transform::MultivariateStats.PCA;
    nn_interpolate::Int;
end

function fit(::Type{UmapFit}, x::Array{Float64, 2}; n_pcs::Int=100, n_components::Int=2, nn_interpolate::Int=5, kwargs...)::UmapFit
    pca = fit(MultivariateStats.PCA, x, maxoutdim=n_pcs)
    nn_tree = KDTree(transform(pca, x))
    embedding = umap(x, n_components, kwargs...)

    return UmapFit(embedding, nn_tree, pca, nn_interpolate)
end

function transform(transformation::UmapFit, x::Array{Float64, 2}; dist_offset::Float64=1e-10)::Array{Float64, 2}
    x_pc = transform(transformation.pca_transform, x)
    indices, distances = knn(transformation.nn_tree, x_pc, transformation.nn_interpolate)

    res = [sum(transformation.embedding[:,ids]' ./ (dists .+ dist_offset), dims=1) ./ sum(1 ./ (dists .+ dist_offset))
        for (ids, dists) in zip(indices, distances)]

    return copy(vcat(res...)')
end