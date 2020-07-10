using Distributions
using LinearAlgebra
using StatsBase

mutable struct CategoricalSmoothed{CT <: Real} <: Distributions.Distribution{Distributions.Multivariate,Distributions.Discrete}
    counts::Vector{CT};
    smooth::Float64;
    n_samples::CT;

    CategoricalSmoothed(counts::Vector{IT}; smooth::Real=1.0, n_samples::IT=sum(counts)) where IT <: Real = new{IT}(counts, Float64(smooth), n_samples)
end

pdf(dist::CategoricalSmoothed, x::Int; use_smoothing::Bool=true)::Float64 =
    pdf(counts(dist), dist.n_samples, x, dist.smooth; use_smoothing=use_smoothing)

function pdf(cnt::Vector{CT}, cnt_sum::CT, x::Int, smooth::Float64; use_smoothing::Bool=true)::Float64 where CT <: Real
    cnt = cnt[x]
    if !use_smoothing
        return cnt / cnt_sum
    end

    # Can't simply add smoothing to each vector component because of sparsity
    return fmax(Float64(cnt), smooth) / (cnt_sum + smooth)
end

function pdf(dist::CategoricalSmoothed, x::Int, prior::Array{Int, 1}, prior_sum::Int; use_smoothing::Bool=true)::Float64
    if prior_sum == 0 # For speed
        return pdf(dist, x; use_smoothing=use_smoothing)
    end

    cnt, cnt_sum = counts(dist)[x], dist.n_samples
    cnt_adj = cnt + prior[x]
    cnt_adj_sum = cnt_sum + prior_sum

    if !use_smoothing
        return cnt_adj / cnt_adj_sum
    end

    if cnt_adj < dist.smooth
        cnt_adj_sum += (dist.smooth - cnt_adj)
        cnt_adj = dist.smooth
    end

    # In case we have very few molecules per cell, prior is incorrect, but it suppresses smoothing
    if cnt < dist.smooth
        cnt_sum += (dist.smooth - cnt)
        return max(dist.smooth / cnt_sum, cnt_adj / cnt_adj_sum) # If cnt_sum is large enough, dist.smooth / cnt_sum is negligible
    end

    return cnt_adj / cnt_adj_sum
end

function maximize!(dist::CategoricalSmoothed, x::T where T <: AbstractArray{Int, 1}, confidences::T2 where T2 <: AbstractVector{Float64})
    count_array!(dist.counts, x, confidences)
    dist.n_samples = sum(confidences)
    return dist
end

counts(d::CategoricalSmoothed) = d.counts