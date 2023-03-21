using Distances
using UMAP
using NearestNeighbors
using Statistics
using Base.Threads

import MultivariateStats
import MultivariateStats.transform
import MultivariateStats.fit

struct UmapFit{TT <: NNTree}
    embedding::Matrix{<:Real};
    nn_tree::TT;
    nn_interpolate::Int;
end

function fit(
        ::Type{UmapFit}, x::AbstractMatrix{<:Real}; n_components::Int=2, nn_interpolate::Int=5,
        spread=2.0, min_dist=0.1, metric::MT where MT <: NearestNeighbors.MinkowskiMetric = Euclidean(), kwargs...
    )
    nn_tree = KDTree(x, metric)
    embedding = umap(x, n_components; spread=spread, min_dist=min_dist, metric=metric, kwargs...)

    return UmapFit(embedding, nn_tree, nn_interpolate)
end

function transform(transformation::UmapFit, x::AbstractMatrix{<:Real}; dist_offset::Float64=1e-10)
    # This transformation is much faster than the implementation from UMAP.jl, even for input dimensionality = 50.
    # Probably, because of the slow NearestNeighborDescent package.
    indices, distances = knn_parallel(transformation.nn_tree, x, transformation.nn_interpolate)
    res = zeros(eltype(transformation.embedding), size(transformation.embedding, 1), size(indices, 1))
    @threads for oi in axes(res, 2)
        for di in axes(res, 1)
            res[di,oi] = wmean(transformation.embedding[di, indices[oi]], 1 ./ (distances[oi] .+ dist_offset))
        end
    end

    return res
end