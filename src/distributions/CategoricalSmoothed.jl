using Distributions
using LinearAlgebra
using StatsBase

mutable struct CategoricalSmoothed <: Distributions.Distribution{Distributions.Multivariate,Distributions.Discrete}
    counts::Array{Int, 1};
    smooth::Float64;
    n_samples::Int;

    CategoricalSmoothed(counts::Array{Int, 1}; smooth::Number=1.0, n_samples::Int=sum(counts)) = new(counts, smooth, n_samples)
end

# function pdf(dist::CategoricalSmoothed, x::Int; use_smoothing::Bool=true)::Float64
#     cnt, cnt_sum = counts(dist)[x], dist.n_samples

#     if !use_smoothing || (cnt >= dist.smooth)
#         return cnt / cnt_sum
#     end

#     # Can't simply add smoothing to each vector component because of sparsity
#     cnt_sum += (dist.smooth - cnt)
#     return dist.smooth / cnt_sum
# end

pdf(dist::CategoricalSmoothed, x::Int; use_smoothing::Bool=true)::Float64 =
    pdf(counts(dist), dist.n_samples, x, dist.smooth; use_smoothing=use_smoothing)

function pdf(cnt::Vector{Int}, cnt_sum::Int, x::Int, smooth::Float64; use_smoothing::Bool=true)::Float64
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

# maximize(dist::CategoricalSmoothed, x::Array{Int, 1}) =
#     CategoricalSmoothed(count_array(x, max_value=length(dist.counts)), smooth=dist.smooth, n_samples=length(x))

function maximize!(dist::CategoricalSmoothed, x::T where T <: AbstractArray{Int, 1})
    count_array!(dist.counts, x)
    dist.n_samples = length(x)
    return dist
end

counts(d::CategoricalSmoothed) = d.counts