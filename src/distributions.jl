using Distributions
using LinearAlgebra
using PDMats
using StatsBase

# MvNormal

function log_pdf(μx::Float64, μy::Float64, Σ::PDMats.PDiagMat{Float64,Array{Float64,1}}, x::Float64, y::Float64)::Float64
    std_x = sqrt(Σ.diag[1])
    std_y = sqrt(Σ.diag[2])

    ndx = (x - μx) / std_x
    ndy = (y - μy) / std_y

    log_div2 = log(2 * π * std_x * std_y)

    return -(ndx * ndx + ndy * ndy) / 2 - log_div2
end

function log_pdf(μx::Float64, μy::Float64, Σ::Array{Float64,2}, x::Float64, y::Float64)::Float64
    std_x = sqrt(Σ[1, 1])
    std_y = sqrt(Σ[2, 2])
    ρ = max(min(1 - 1e-2, Σ[1, 2] / (std_x * std_y)), -1 + 1e-2)

    ndx = (x - μx) / std_x
    ndy = (y - μy) / std_y;

    div1 = 2 * (1 - ρ^2)
    log_div2 = log(2 * π * std_x * std_y * sqrt(1 - ρ^2))

    return -(ndx * ndx + ndy * ndy - 2 * ρ * ndx * ndy) / div1 - log_div2
end

log_pdf(μx::Float64, μy::Float64, Σ::PDMats.PDMat{Float64,Array{Float64,2}}, x::Float64, y::Float64) = log_pdf(μx, μy, Σ.mat, x, y)
log_pdf(d::MvNormal, x::Float64, y::Float64)::Float64 = log_pdf(d.μ[1], d.μ[2], d.Σ, x, y)

pdf(d::MvNormal, x::Float64, y::Float64) = exp(log_pdf(d, x, y))
shape(d::MvNormal) = Matrix(d.Σ)

function robust_cov(x::Array{Float64, 2}; prop::Float64=0.1)
    v1 = mad(x[:,1], normalize=true)^2
    v2 = mad(x[:,2], normalize=true)^2

    dev1 = x[:,1] .- trim_mean(x[:,1]; prop=prop)
    dev2 = x[:,2] .- trim_mean(x[:,2]; prop=prop)
    covar = mean(winsor(dev1 .* dev2; prop=prop))

    return [v1 covar; covar v2]
end

function maximize(dist::Distributions.MvNormal, x::Array{Float64, 2})::Distributions.MvNormal
    if size(x, 2) == 0
        return dist
    end

    x = copy(x')
    μ = [trim_mean(x[:,1], prop=0.2), trim_mean(x[:,2], prop=0.2)]
    if size(x, 1) <= 2
        return dist
    end

    Σ = adjust_cov_matrix(robust_cov(x; prop=0.1))

    if det(Σ) ≈ 0
        @show Σ
        throw(ErrorException("Singular covariance matrix"))
    end

    return MvNormal(μ, Σ)
end

function adjust_cov_matrix(Σ::Array{Float64, 2}; max_iter::Int=100,
                           cov_modifier::Array{Float64, 2}=Matrix(Diagonal(1e-4 .* ones(size(Σ, 2)))), tol=1e-10)::Array{Float64, 2}
    if isposdef(Σ) && det(Σ) > tol
        return Σ
    end

    res = copy(Σ)
    for it in 1:max_iter
        if isposdef(res) && det(res) > tol
            break
        end

        res += cov_modifier
        cov_modifier *= 1.5
    end

    return res
end

# SingleTrialMultinomial

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

# Scaled Inverse Chi Square

struct ScaledInverseChisq <: Distributions.Distribution{Distributions.Univariate,Distributions.Continuous}
    ν::Int
    τ::Float64
end

rand(d::ScaledInverseChisq, args::Vararg{Int}) = d.τ^2 .* d.ν ./ rand(Chisq(d.ν), args...)
pdf(d::ScaledInverseChisq, x) = pdf(Chisq(d.ν), d.τ^2 * d.ν / x)

# Scaled T

struct TScaled <: Distributions.Distribution{Distributions.Univariate,Distributions.Continuous}
    μ::Float64;
    σ::Float64;
    ν::Float64;

    t_dist::TDist;

    TScaled(μ::Real, σ::Real, ν::Real) = new(μ, σ, ν, TDist(ν))
end

pdf(d::TScaled, x::Float64) = pdf(d.t_dist, (x - d.μ) / σ)
rand(d::TScaled) = rand(d.t_dist) * d.σ + d.μ

# Normal-Gamma

"""
    NormalGamma(μ, κ, α, β)

Mean was estimated from κ observations with sample mean μ , and precision was estimated from 2α observations with sample mean μ and sum of squared deviations 2β.

# Arguments
- `μ::Float64`: mean of the distribution
- `κ::Float64`: number of observations, used to estimate mean
- `α::Float64`: 0.5 * number of observations used to estimate precision
- `β::Float64`: 0.5 * sum of squared deviations
"""
struct NormalGamma <: Distributions.Distribution{Distributions.Multivariate,Distributions.Continuous}
    μ::Float64;
    κ::Float64;
    α::Float64;
    β::Float64;

    t_scaled::TScaled;
    Γ_dist::Gamma;

    NormalGamma(μ::Real, κ::Real, α::Real, β::Real) =
        new(μ, κ, α, β, TScaled(μ, sqrt(β / (κ * α)), 2 * α), Gamma(α, 1 / β)) # Gamma here uses different parametrization
end

pdf(d::NormalGamma, μ::Float64, λ::Float64) = pdf(Normal(d.μ, 1 / sqrt(d.κ * λ)), μ) * pdf(d.Γ_dist, λ)
rand(d::NormalGamma) = (rand(d.t_scaled), rand(d.Γ_dist))
# var(d::NormalGamma) = d.β / d.α
