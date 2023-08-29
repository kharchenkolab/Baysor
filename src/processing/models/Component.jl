using Distributions
using StaticArrays

mutable struct Component{N, CT}
    position_params::MvNormalF{N};
    composition_params::CT;
    n_samples::Int;
    prior_probability::Float64;
    confidence::Float64;

    shape_prior::Union{Nothing, ShapePrior{N}};
    n_molecules_per_segment::Dict{Int, Int};

    guid::Int;

    Component(
        position_params::MvNormalF{L}, composition_params::CT;
        n_samples::Int=0, confidence::Float64=1.0, shape_prior::Union{Nothing, ShapePrior{L}}=nothing, guid::Int=-1
    ) where {L, CT} =
        new{L, CT}(position_params, composition_params, n_samples, 1.0, confidence, shape_prior, Dict{Int, Int}(), guid)
end

function maximize!(
        c::Component{N, CT} where {N, CT}, pos_data::T1 where T1 <: AbstractMatrix{Float64},
        comp_data::T2 where T2 <: Union{AbstractVector{Int}, AbstractVector{Union{Int, Missing}}, AbstractMatrix{Float64}};
        nuclei_probs::Union{<:AbstractVector{Float64}, Nothing}=nothing, min_nuclei_frac::Float64=0.1
    )
    c.n_samples = size(pos_data, 2)
    maximize!(c.composition_params, comp_data)
    maximize!(c.position_params, pos_data; center_probs=nuclei_probs, c.shape_prior, c.n_samples)
    if (nuclei_probs !== nothing) && (length(nuclei_probs) > 1)
        c.confidence = quantile(nuclei_probs, 1 - min_nuclei_frac)
    end

    return c
end

pdf(comp::Component{N, CT} where {N, CT}, x::AbstractVector{Float64}, gene::Missing; use_smoothing::Bool=true) =
    comp.prior_probability * comp.confidence * pdf(comp.position_params, x)

pdf(comp::Component{N, CategoricalSmoothed{FT}} where {N, FT}, x::AbstractVector{Float64}, gene::Int64; use_smoothing::Bool=true) =
    comp.prior_probability * comp.confidence * pdf(comp.position_params, x) * pdf(comp.composition_params, gene; use_smoothing=use_smoothing)

function pdf(comp::Component{T, MvNormalF{M, N}} where {T, M, N}, x::AbstractVector{Float64}, gene::AbstractVector{<:Real}; use_smoothing::Bool=true)
    # comp.prior_probability * comp.confidence * pdf(comp.position_params, x) * pdf(comp.composition_params, gene)
    p_pos = pdf(comp.position_params, x)
    p_comp = pdf(comp.composition_params, gene)

    # weighted geometric mean
    wp = 1.0
    # wc = 0.25
    wc = 0.5
    return comp.prior_probability * comp.confidence * ((p_pos^wp) * (p_comp^wc)) ^ (1 / (wp + wc))
    # return comp.prior_probability, comp.confidence, p_pos, p_comp, (wp * p_pos + wc * p_comp) ^ (1 / (wp + wc))
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
