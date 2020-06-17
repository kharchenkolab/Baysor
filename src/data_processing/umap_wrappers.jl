using Distances
using UMAP
using NearestNeighbors
using Statistics

import MultivariateStats
import MultivariateStats.transform
import MultivariateStats.fit

struct UmapFit{TT <: NNTree}
    embedding::Matrix{<:Real};
    nn_tree::TT;
    pca_transform::MultivariateStats.PCA;
    nn_interpolate::Int;
end

function fit(::Type{UmapFit}, x::Array{Float64, 2}, pca::MultivariateStats.PCA; n_components::Int=2, nn_interpolate::Int=5,
        spread=2.0, min_dist=0.1, metric::MT where MT <: NearestNeighbors.MinkowskiMetric = Euclidean(), kwargs...)::UmapFit
    x_pc = transform(pca, x)
    nn_tree = KDTree(x_pc, metric)
    embedding = umap(x_pc, n_components; spread=spread, min_dist=min_dist, metric=metric, kwargs...)

    return UmapFit(embedding, nn_tree, pca, nn_interpolate)
end

fit(::Type{UmapFit}, x::Array{Float64, 2}; n_pcs::Int=15, kwargs...)::UmapFit =
    fit(UmapFit, x, fit(MultivariateStats.PCA, x, maxoutdim=n_pcs); kwargs...)

function transform(transformation::UmapFit, x::Array{Float64, 2}; dist_offset::Float64=1e-10)::Array{Float64, 2}
    x_pc = transform(transformation.pca_transform, x)

    # This transformation is much faster than the implementation from UMAP.jl, even for input dimensionality = 50.
    # Probably, because of the slow NearestNeighborDescent package.
    indices, distances = knn(transformation.nn_tree, x_pc, transformation.nn_interpolate)
    res = [mapslices(x -> wmean(x, 1 ./ (dists .+ dist_offset)), transformation.embedding[:, ids], dims=2) for (ids, dists) in zip(indices, distances)]
    return hcat(res...)
end