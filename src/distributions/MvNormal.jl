using Distributions
using LinearAlgebra
using StatsBase
using StaticArrays

import LinearAlgebra.isposdef

CovMat = (MMatrix{T, T, Float64} where T)
MeanVec = (MVector{T, Float64} where T)

struct MvNormalF{L}
    μ::MeanVec{L};
    Σ::CovMat{L};
    MvNormalF(μ::MeanVec{N}, Σ::CovMat{N}=CovMat{N}(diagm(0 => ones(N)))) where N = new{N}(μ, Σ)
    function MvNormalF(μ::T where T <: AbstractVector{Float64})
        if length(μ) == 2
            return MvNormalF(MeanVec{2}(μ))
        end

        return MvNormalF(MeanVec{3}(μ))
    end
        
    function MvNormalF(μ::T1 where T1 <: AbstractVector{Float64}, Σ::T2 where T2 <: AbstractMatrix{Float64})
        if length(μ) == 2
            return MvNormalF(MeanVec{2}(μ), CovMat{2}(Σ))
        end

        return MvNormalF(MeanVec{3}(μ), CovMat{3}(Σ))
    end
end

function wmean!(μ::MV where MV <: MeanVec, values::T where T <: AbstractMatrix{<:Real}, weights::T1 where T1 <: AbstractVector{<:Real})
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

@inline @inbounds function log_pdf(d::MvNormalF{2}, x::SVector{2, Float64})::Float64
    std_x = sqrt(d.Σ[1])
    std_y = sqrt(d.Σ[4])
    ρ = fmax(fmin(1 - 1e-2, d.Σ[3] / (std_x * std_y)), -1 + 1e-2)

    ndx = (x[1] - d.μ[1]) / std_x
    ndy = (x[2] - d.μ[2]) / std_y

    div1 = 2 * (1 - ρ^2)
    log_div2 = log(2 * π * std_x * std_y * sqrt(1 - ρ^2))

    return -(ndx * ndx + ndy * ndy - 2 * ρ * ndx * ndy) / div1 - log_div2
end

@inline function log_pdf(d::MvNormalF{3}, x::SVector{3, Float64})::Float64 # TODO: pre-estimate divider and  inv(d.Σ)
    divider = (2π)^3 * det(d.Σ)
    dx = x .- d.μ
    return -0.5 * (dx') * inv(d.Σ) * dx - 0.5 * log(divider)
end

pdf(d::MvNormalF, x::SVector{N, Float64} where N) = exp(log_pdf(d, x))
shape(d::MvNormalF) = Matrix(d.Σ)

function estimate_sample_cov!(Σ::CovMat{N}, x::T where T <: AbstractMatrix{Float64}, weights::Nothing=nothing; μ::MeanVec{N}) where N
    # https://en.wikipedia.org/wiki/Estimation_of_covariance_matrices#Intrinsic_expectation
    Σ .= 0.0
    for i in 1:size(x, 2)
        for c in 1:size(Σ, 1)
            for r in 1:size(Σ, 1)
                Σ[r, c] += (x[c, i] - μ[c]) * (x[r, i] - μ[r])
            end
        end
    end

    Σ ./= size(x, 2)
    return Σ
end

function estimate_sample_cov!(Σ::CovMat{N}, x::T where T <: AbstractMatrix{Float64}, weights::T2 where T2 <: AbstractVector{<:Real}; μ::MeanVec{N}) where N
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

function maximize!(dist::MvNormalF, x::T, confidences::T2 where T2 <: Union{AbstractVector{<:Real}, Nothing}=nothing)::MvNormalF where T <: AbstractMatrix{Float64}
    if size(x, 2) == 0
        return dist
    end

    if confidences === nothing
        mean!(dist.μ, x)
    else
        wmean!(dist.μ, x, confidences)
    end
    if size(x, 2) <= length(dist.μ)
        return dist
    end

    estimate_sample_cov!(dist.Σ, x, confidences; μ=dist.μ)
    return dist
end

isposdef(A::StaticMatrix) = LinearAlgebra.isposdef(cholesky(A))

function adjust_cov_matrix!(Σ::T; max_iter::Int=100, cov_modifier::Float64=1e-4, tol=1e-10)::T where T <: Union{CovMat{2}, CovMat{3}, Matrix{Float64}}
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