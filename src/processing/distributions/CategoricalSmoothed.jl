using Distributions
using LinearAlgebra

mutable struct CategoricalSmoothed{CT <: Real} <: Distributions.Distribution{Distributions.Multivariate,Distributions.Discrete}
    counts::Vector{CT};
    smooth::Float64;
    sum_counts::CT;
    n_genes::Int;

    CategoricalSmoothed(counts::Vector{IT}; smooth::Real=1.0, sum_counts::IT=sum(counts)) where IT <: Real =
        new{IT}(counts, Float64(smooth), sum_counts, sum(counts .> 0))
end

function pdf(dist::CategoricalSmoothed, x::Int; use_smoothing::Bool=true)::Float64
    cnt = dist.counts[x]
    if !use_smoothing
        return cnt / dist.sum_counts
    end

    # Can't simply add smoothing to each vector component because of sparsity
    return fmax(Float64(cnt), dist.smooth) / (dist.sum_counts + dist.smooth)
end

function maximize!(dist::CategoricalSmoothed, x::T where T <: Union{AbstractArray{Int, 1}, AbstractArray{Union{Missing, Int}, 1}}, confidences::T2 where T2 <: AbstractVector{Float64})
    dist.counts .= 0
    dist.sum_counts = 0.0
    dist.n_genes = 0
    for i in 1:length(x)
        c = confidences[i]
        v = x[i]
        ismissing(v) && continue

        if dist.counts[v] ≈ 0
            dist.n_genes += 1
        end
        dist.counts[v] += c
        dist.sum_counts += c
    end

    return dist
end

function maximize!(dist::CategoricalSmoothed, x::T where T <: Union{AbstractArray{Int, 1}, AbstractArray{Union{Missing, Int}, 1}})
    dist.counts .= 0
    dist.sum_counts = 0.0
    dist.n_genes = 0
    for v in x
        ismissing(v) && continue
        if dist.counts[v] ≈ 0
            dist.n_genes += 1
        end

        dist.counts[v] += 1
        dist.sum_counts += 1
    end

    return dist
end
