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
        nuclei_probs::Union{<:AbstractVector{Float64}, Nothing}=nothing, min_nuclei_frac::Float64=0.1,
        freeze_position::Bool=false, freeze_composition::Bool=false
    )
    c.n_samples = size(pos_data, 2)

    if !freeze_composition
        maximize!(c.composition_params, comp_data)
    end

    if !freeze_position
        maximize!(c.position_params, pos_data; center_probs=nuclei_probs, c.shape_prior, c.n_samples)
    end

    if (nuclei_probs !== nothing) && (length(nuclei_probs) > 1)
        c.confidence = quantile(nuclei_probs, 1 - min_nuclei_frac)
    end

    return c
end

pdf(comp::Component{N, CT} where {N, CT}, x::AbstractVector{Float64}, gene::Missing; use_smoothing::Bool=true) =
    comp.prior_probability * comp.confidence * pdf(comp.position_params, x)

pdf(comp::Component{N, CategoricalSmoothed{FT}} where {N, FT}, x::AbstractVector{Float64}, gene::Int64; use_smoothing::Bool=true) =
    comp.prior_probability * comp.confidence * pdf(comp.position_params, x) * pdf(comp.composition_params, gene; use_smoothing=use_smoothing)

pdf(comp::Component{T, CT} where {T, CT}, x::AbstractVector{Float64}, gene::AbstractVector{<:Real}; use_smoothing::Bool=true) =
    comp.prior_probability * comp.confidence * pdf(comp.position_params, x) * pdf(comp.composition_params, gene)
