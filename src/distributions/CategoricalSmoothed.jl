using Distributions
using LinearAlgebra
using StatsBase

mutable struct CategoricalSmoothed{CT <: Real} <: Distributions.Distribution{Distributions.Multivariate,Distributions.Discrete}
    counts::Vector{CT};
    smooth::Float64;
    sum_counts::CT;

    CategoricalSmoothed(counts::Vector{IT}; smooth::Real=1.0, sum_counts::IT=sum(counts)) where IT <: Real = new{IT}(counts, Float64(smooth), sum_counts)
end

function pdf(dist::CategoricalSmoothed, x::Int; use_smoothing::Bool=true)::Float64
    cnt = dist.counts[x]
    if !use_smoothing
        return cnt / dist.sum_counts
    end

    # Can't simply add smoothing to each vector component because of sparsity
    return fmax(Float64(cnt), smooth) / (dist.sum_counts + smooth)
end

function maximize!(dist::CategoricalSmoothed, x::T where T <: AbstractArray{Int, 1}, confidences::T2 where T2 <: AbstractVector{Float64})
    count_array!(dist.counts, x, confidences)
    dist.sum_counts = sum(confidences)
    return dist
end
