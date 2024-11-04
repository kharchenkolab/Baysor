using Distributions
using LinearAlgebra
using StaticArrays

import LinearAlgebra.isposdef

CovMat = (MMatrix{T, T, Float64, T2} where T2 where T)
MeanVec = (MVector{T, Float64} where T)

struct ShapePrior{N}
    std_values::MeanVec{N};
    std_value_stds::MeanVec{N};
    n_samples::Int;
end

distributions(prior::ShapePrior) = Normal.(prior.std_values, prior.std_value_stds)
rand(prior::ShapePrior) = rand.(distributions(prior))

var_posterior(prior::ShapePrior{N}, eigen_values::T where T <: Union{MeanVec{N}, StaticArray{Tuple{N},Float64,1}}; n_samples::TR where TR <: Real) where N =
    var_posterior.(prior.std_values, prior.std_value_stds, eigen_values; n_samples=n_samples, prior_n_samples=prior.n_samples)

# f(dx) = sign(dx) * sqrt(|dx/std|) * std
function var_posterior(prior_std::Float64, prior_std_std::Float64, eigen_value::Float64; n_samples::TR1 where TR1 <: Real, prior_n_samples::TR2 where TR2 <: Real)
    d_std = (sqrt(eigen_value) - prior_std)
    std_adj = (prior_std + sign(d_std) * (sqrt(abs(d_std) / prior_std_std + 1) - 1) * prior_std_std)

    # The formula below only corrects variance estimates for cells with low number of molecules
    # Cells with high number of molecules still get their correction above
    return ((prior_n_samples * prior_std + n_samples * std_adj) / (prior_n_samples + n_samples)) ^ 2
end

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

function estimate_sample_cov(x::AbstractMatrix{<:Real})
    CM, MV = (size(x, 1) == 2) ? (CovMat{2}, MeanVec{2}) : (CovMat{3}, MeanVec{3})
    return estimate_sample_cov!(zeros(CM), x; μ=MV(mean(x, dims=2)[:]))
end

function estimate_sample_cov!(Σ::CovMat{N, N2} where N2, x::T where T <: AbstractMatrix{<:Real}; μ::MeanVec{N}) where N
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

function update_cache!(dist::MvNormalF)
    dist.Σ_inv .= inv(dist.Σ)
    dist.pdf_divider = norm_pdf_divider(dist.Σ)
end

function maximize!(
        dist::MvNormalF{N}, x::T;
        center_probs::T2 where T2 <: Union{AbstractVector{<:Real}, Nothing}=nothing,
        shape_prior::Union{Nothing, ShapePrior{N}}=nothing,
        n_samples::Int=-1
    ) where T <: AbstractMatrix{<:Real} where N
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
    if shape_prior !== nothing
        @assert n_samples >= 0
        adjust_cov_by_prior!(dist.Σ, shape_prior; n_samples=n_samples)
    end

    update_cache!(dist)
    return dist
end

# Copied from LinearAlgebra implementation
isposdef(A::StaticMatrix) =
    LinearAlgebra.ishermitian(A) && LinearAlgebra.isposdef(LinearAlgebra.cholesky(LinearAlgebra.Hermitian(A); check = false))

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

function adjust_cov_by_prior!(Σ::CovMat{N}, prior::ShapePrior{N}; n_samples::TR where TR <: Real) where N
    if (Σ[2, 1] / max(Σ[1, 1], Σ[2, 2])) < 1e-5 # temporary fix untill https://github.com/JuliaArrays/StaticArrays.jl/pull/694 is merged
        Σ[1, 2] = Σ[2, 1] = 0.0
    end

    fact = eigen(Σ)
    eigen_values_posterior = var_posterior(prior, fmax.(fact.values, 0.0); n_samples=n_samples)
    Σ .= fact.vectors * SMatrix{N, N, Float64}(diagm(0 => eigen_values_posterior)) * inv(fact.vectors)
    Σ .= fmax.(Σ, Σ')

    return adjust_cov_matrix!(Σ)
end