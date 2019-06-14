using Distributions
using LinearAlgebra
using StatsBase

mutable struct SingleTrialMultinomial <: Distributions.Distribution{Distributions.Multivariate,Distributions.Discrete}
    counts::Array{Int, 1};
    smooth::Float64;
    n_samples::Int;

    SingleTrialMultinomial(counts::Array{Int, 1}; smooth::Number=1.0, n_samples::Int=0) = new(counts, smooth, n_samples)
end

# returns first element if it's presented in the data, but also correctly process case with huge smooth
function pdf(dist::SingleTrialMultinomial, x::Union{Int, Array{Int, 1}}; use_smoothing::Bool=true, prior::Union{Nothing, Array{Int, 1}}=nothing,
             prior_sum::Int=0)::Float64
    cnt, cnt_sum = counts(dist)[x], dist.n_samples
    if prior !== nothing
        cnt += prior[x]
        cnt_sum += prior_sum
    end

    if !use_smoothing
        return cnt / cnt_sum
    end

    return max(cnt, dist.smooth) / (cnt_sum + dist.smooth)
end

maximize(dist::SingleTrialMultinomial, x::Array{Int, 1}) =
    SingleTrialMultinomial(count_array(x, max_value=length(dist.counts)), smooth=dist.smooth, n_samples=length(x))

function maximize!(dist::SingleTrialMultinomial, x::Array{Int, 1})
    count_array!(dist.counts, x)
    dist.n_samples = length(x)
    return dist
end

counts(d::SingleTrialMultinomial) = d.counts