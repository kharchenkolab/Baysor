using UMAP
using NearestNeighbors
using Statistics

# import Distances: SemiMetric, Euclidean
import MultivariateStats
import MultivariateStats.transform
import MultivariateStats.fit
# import NearestNeighborDescent: DescentGraph, search

struct UmapFit
    embedding::Array;
    nn_tree;
    pca_transform::MultivariateStats.PCA;
    nn_interpolate::Int;
end

function fit(::Type{UmapFit}, x::Array{Float64, 2}; n_pcs::Int=15, n_components::Int=2, nn_interpolate::Int=5, kwargs...)::UmapFit # , metric::SemiMetric=Euclidean()
    pca = fit(MultivariateStats.PCA, x, maxoutdim=n_pcs)
    x_pc = transform(pca, x)
    # nn_tree = DescentGraph(x_pc, nn_interpolate, metric)
    # indices, distances = search(transformation.nn_tree, x_pc, transformation.nn_interpolate)
    # embedding = umap(x_pc, n_components; metric=metric, kwargs...)
    nn_tree = KDTree(x_pc)
    embedding = umap(x_pc, n_components; kwargs...)


    return UmapFit(embedding, nn_tree, pca, nn_interpolate)
end

function transform(transformation::UmapFit, x::Array{Float64, 2}; dist_offset::Float64=1e-10)::Array{Float64, 2}
    x_pc = transform(transformation.pca_transform, x)
    # indices, distances = search(transformation.nn_tree, x_pc, transformation.nn_interpolate)
    indices, distances = knn(transformation.nn_tree, x_pc, transformation.nn_interpolate)

    res = [sum(transformation.embedding[:,ids]' ./ (dists .+ dist_offset), dims=1) ./ sum(1 ./ (dists .+ dist_offset))
        for (ids, dists) in zip(indices, distances)]

    return copy(vcat(res...)')
end