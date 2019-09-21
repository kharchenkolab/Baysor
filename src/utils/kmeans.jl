import Clustering
import Random

using Distances
using Statistics

function calinski_harabasz_score(points::Matrix, labels::Vector{Int})
    n_labels = length(unique(labels))

    extra_disp, intra_disp = 0., 0.
    c_mean = mean(points, dims=2)
    for ids in split_ids(labels)
        cluster_k = points[:, ids]
        mean_k = mean(cluster_k, dims=2)
        extra_disp += length(ids) * sum((mean_k .- c_mean) .^ 2)
        intra_disp += sum((cluster_k .- mean_k) .^ 2)
    end

    return intra_disp == 0. ? 1. : extra_disp * (size(points, 2) - n_labels) / (intra_disp * (n_labels - 1.))
end

assign_to_centers(points::Matrix{T} where T <: Real, means::Matrix{Float64}, distance::D where D <: SemiMetric) =
    vec(mapslices(x -> findmin(x)[2], pairwise(distance, points, means; dims=2), dims=2))

function kmeans!(points::Matrix{T} where T <: Real, means::Matrix{Float64}; min_cluster_size::Int=0,
        max_iters::Int=10000, distance::D where D <: SemiMetric = Euclidean(), robust::Bool=true)
    min_cluster_size = min(size(points, 2), min_cluster_size)
    avg_func = robust ? median : mean

    assignment = zeros(Int, size(points, 2))
    prev_assignment = assignment

    for i in 1:max_iters
        prev_assignment = assignment
        assignment = assign_to_centers(points, means, distance)
        if all(prev_assignment .== assignment)
            break
        end

        means = hcat([length(ids) > 0 ? avg_func(points[:, ids], dims=2) : points[:, rand(1:size(points, 2))] for ids in split_ids(assignment)]...)
    end

    if any(prev_assignment .!= assignment)
        @warn "K-means didn't converge for $max_iters iterations"
    end

    n_points_per_clust = count_array(assignment)
    if any(n_points_per_clust .< min_cluster_size)
        means = means[:, n_points_per_clust .>= min_cluster_size]
        kmeans!(points, means; min_cluster_size=min_cluster_size, max_iters=max_iters, distance=distance, robust=robust)
    end

    return means, assignment
end

kmeans(points::Matrix{T} where T <: Real, k::Int; init_alg::Symbol=:kmpp, kwargs...) =
    kmeans!(points, points[:, Clustering.initseeds(init_alg, points, k)]; kwargs...)

function kmeans_stable(points::Matrix{T} where T <: Real, args...; n_inits::Int=10, kwargs...)
    res = [kmeans(points, args...; kwargs...) for i in 1:n_inits]
    scores = [calinski_harabasz_score(points, r[2]) for r in res]

    centers, labels = res[findmax(scores)[2]]
    centers = centers[:, sortperm(1:size(centers, 2), lt= (i, j) -> centers[1, i] < centers[1, j])]
    return centers, labels
end