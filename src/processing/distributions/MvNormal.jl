using Distributions
using LinearAlgebra
using StaticArrays

import LinearAlgebra.isposdef

CovMat = (MMatrix{T, T, Float64, T2} where T2 where T)
MeanVec = (MVector{T, Float64} where T)

mutable struct MvNormalF{L, L2}
    μ::MeanVec{L};
    Σ::CovMat{L, L2};
    pdf_divider::Float64;
    Σ_inv::CovMat{L, L2};
end

MvNormalF(μ::MeanVec{N}, Σ::CovMat{N}=CovMat{N}(diagm(0 => ones(N)))) where N =
    MvNormalF(μ, Σ, norm_pdf_divider(Σ), inv(Σ))

function MvNormalF(μ::T where T <: AbstractVector{Float64})
    if length(μ) == 2
        return MvNormalF(MeanVec{2}(μ))
    end

    return MvNormalF(MeanVec{3}(μ))
end

function MvNormalF(μ::T1 where T1 <: AbstractVector{Float64}, Σ::T2 where T2 <: AbstractMatrix{Float64})
    if length(μ) == 2
        return MvNormalF(MeanVec{2}(μ), CovMat{2, 4}(Σ))
    end

    return MvNormalF(MeanVec{3}(μ), CovMat{3, 9}(Σ))
end

norm_pdf_divider(Σ::CovMat{T}) where T = 0.5 * log((2π)^3 * det(Σ))

function wmean!(μ::MV where MV <: MeanVec, values::T where T <: AbstractMatrix{<:Real}, weights::T1 where T1 <: AbstractVector{<:Real})
    @assert length(μ) == size(values, 1)
    μ .= 0.0
    sum_w = 0.0
    for ci in 1:size(values, 2)
        w = fmax(weights[ci], 0.01)
        for ri in eachindex(μ)
            μ[ri] += values[ri, ci] * w
        end
        sum_w += w
    end

    μ ./= sum_w
    return μ
end

@inline @inbounds vXvt(a::Float64, b::Float64, m::CovMat{2, 4}) =
    a*a*m[1, 1] + b*b*m[2, 2] + 2*a*b*m[1, 2]

@inline @inbounds vXvt(a::Float64, b::Float64, c::Float64, m::CovMat{3, 9}) =
    a*a*m[1, 1] + b*b*m[2, 2] + c*c*m[3, 3] + 2*a*b*m[1, 2] + 2*a*c*m[1, 3] + 2*b*c*m[2, 3]

log_pdf(d::MvNormalF{2, 4}, x::AbstractVector{Float64})::Float64 =
    -0.5 * vXvt(x[1] - d.μ[1], x[2] - d.μ[2], d.Σ_inv) - d.pdf_divider

log_pdf(d::MvNormalF{3, 9}, x::AbstractVector{Float64})::Float64 =
    -0.5 * vXvt(x[1] - d.μ[1], x[2] - d.μ[2], x[3] - d.μ[3], d.Σ_inv) - d.pdf_divider

pdf(d::MvNormalF, x::AbstractVector{Float64}) = exp(log_pdf(d, x))

shape(d::MvNormalF) = Matrix(d.Σ)

function estimate_sample_cov!(Σ::CovMat{N, N2} where N2, x::T where T <: AbstractMatrix{Float64}; μ::MeanVec{N}) where N
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

function maximize!(dist::MvNormalF, x::T; center_probs::T2 where T2 <: Union{AbstractVector{<:Real}, Nothing}=nothing)::MvNormalF where T <: AbstractMatrix{Float64}
    if size(x, 2) == 0
        return dist
    end

    if center_probs === nothing
        mean!(dist.μ, x)
    else
        wmean!(dist.μ, x, center_probs)
    end
    if size(x, 2) <= length(dist.μ)
        return dist
    end

    estimate_sample_cov!(dist.Σ, x; μ=dist.μ)

    dist.Σ_inv .= inv(dist.Σ) # TODO: make it consistent with mutability of MvNormal
    dist.pdf_divider = norm_pdf_divider(dist.Σ)

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