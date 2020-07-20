using Distributions
using LinearAlgebra
using StatsBase
using StaticArrays

import LinearAlgebra.isposdef

CovMat = MMatrix{2, 2, Float64, 4}
MeanVec = MVector{2, Float64}

struct MvNormalF
    μ::MeanVec;
    Σ::CovMat;
    MvNormalF(μ::MeanVec, Σ::CovMat=CovMat(diagm(0 => ones(2)))) = new(μ, Σ)
    MvNormalF(μ::T where T <: AbstractArray{Float64, 1}) = MvNormalF(MeanVec(μ))
    MvNormalF(μ::T1 where T1 <: AbstractArray{Float64, 1}, Σ::T2 where T2 <: AbstractArray{Float64, 2}) =
        MvNormalF(MeanVec(μ), CovMat(Σ))
end

function wmean!(μ::MeanVec, values::T where T <: AbstractMatrix{<:Real}, weights::T1 where T1 <: AbstractVector{<:Real})
    @assert length(μ) == size(values, 1)
    μ .= 0.0
    sum_w = 0.0
    for ci in 1:size(values, 2)
        w = fmax(weights[ci], 0.01)
        for ri in 1:length(μ)
            μ[ri] += values[ri, ci] * w
        end
        sum_w += w
    end

    μ ./= sum_w
    return μ
end

@inline @inbounds function log_pdf(d::MvNormalF, x::Float64, y::Float64)::Float64
    std_x = sqrt(d.Σ[1])
    std_y = sqrt(d.Σ[4])
    ρ = fmax(fmin(1 - 1e-2, d.Σ[3] / (std_x * std_y)), -1 + 1e-2)

    ndx = (x - d.μ[1]) / std_x
    ndy = (y - d.μ[2]) / std_y;

    div1 = 2 * (1 - ρ^2)
    log_div2 = log(2 * π * std_x * std_y * sqrt(1 - ρ^2))

    return -(ndx * ndx + ndy * ndy - 2 * ρ * ndx * ndy) / div1 - log_div2
end

pdf(d::MvNormalF, x::Float64, y::Float64) = exp(log_pdf(d, x, y))
shape(d::MvNormalF) = Matrix(d.Σ)

function robust_cov(x::Array{Float64, 2}; prop::Float64=0.1)
    v1 = mad(x[:,1], normalize=true)^2
    v2 = mad(x[:,2], normalize=true)^2

    dev1 = x[:,1] .- trim_mean(x[:,1]; prop=prop)
    dev2 = x[:,2] .- trim_mean(x[:,2]; prop=prop)
    covar = mean(winsor(dev1 .* dev2; prop=prop))

    return [v1 covar; covar v2]
end

estimate_sample_cov(args...; kwargs...) =
    estimate_sample_cov!(CovMat(zeros(2, 2)), args...; kwargs...)

function estimate_sample_cov!(Σ::CovMat, x::T where T <: AbstractMatrix{Float64}, weights::T2 where T2 <: AbstractVector{<:Real}; μ::MeanVec)
    # https://en.wikipedia.org/wiki/Estimation_of_covariance_matrices#Intrinsic_expectation
    Σ .= 0.0
    sum_w = 0.0
    for i in 1:size(x, 2)
        w = fmax(weights[i], 0.01)
        for c in 1:size(Σ, 1)
            for r in 1:size(Σ, 1)
                Σ[r, c] += w * (x[c, i] - μ[c]) * (x[r, i] - μ[r])
            end
        end
        sum_w += w
    end

    Σ ./= sum_w
    return Σ
end

function maximize!(dist::MvNormalF, x::T, confidences::T2 where T2 <: AbstractVector{<:Real})::MvNormalF where T <: AbstractMatrix{Float64}
    if size(x, 2) == 0
        return dist
    end

    wmean!(dist.μ, x, confidences)
    if size(x, 2) <= 2
        return dist
    end

    estimate_sample_cov!(dist.Σ, x, confidences; μ=dist.μ)
    return dist
end

isposdef(A::StaticMatrix) = LinearAlgebra.isposdef(cholesky(A))

function adjust_cov_matrix!(Σ::T; max_iter::Int=100, cov_modifier::Float64=1e-4, tol=1e-10)::T where T <: Union{CovMat, Matrix{Float64}}
    if !issymmetric(Σ)
        error("Non-symmetric matrix: $Σ")
    end

    if isposdef(Σ) && det(Σ) > tol
        return Σ
    end

    for it in 1:max_iter
        if isposdef(Σ) && det(Σ) > tol
            break
        end

        @inbounds Σ[diagind(Σ)] .+= cov_modifier
        cov_modifier *= 1.5
    end

    return Σ
end