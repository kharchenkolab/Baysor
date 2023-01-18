using Distributions
using StaticArrays

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

    return ((prior_n_samples * prior_std + n_samples * std_adj) / (prior_n_samples + n_samples)) ^ 2
end

function sample_var(d::Normal)
    @assert d.μ > 0
    while true
        v = rand(d)
        if v > 0
            return v.^2
        end
    end
end
sample_var(prior::ShapePrior) = sample_var.(distributions(prior))

mutable struct Component{N}
    position_params::MvNormalF{N};
    composition_params::CategoricalSmoothed;
    n_samples::Int;
    prior_probability::Float64;
    confidence::Float64;

    shape_prior::Union{Nothing, ShapePrior{N}};
    n_molecules_per_segment::Dict{Int, Int};

    guid::Int;

    Component(position_params::MvNormalF{L}, composition_params::CategoricalSmoothed;
              n_samples::Int=0, confidence::Float64=1.0, shape_prior::Union{Nothing, ShapePrior{L}}=nothing, guid::Int=-1) where L =
        new{L}(position_params, composition_params, n_samples, 1.0, confidence, shape_prior, Dict{Int, Int}(), guid)
end

function maximize!(c::Component{N} where N, pos_data::T1 where T1 <: AbstractMatrix{Float64}, comp_data::T2 where T2 <: Union{AbstractVector{Int}, AbstractVector{Union{Int, Missing}}};
        nuclei_probs::Union{<:AbstractVector{Float64}, Nothing}=nothing, min_nuclei_frac::Float64=0.1)
    c.n_samples = size(pos_data, 2)
    maximize!(c.composition_params, comp_data)
    if nuclei_probs === nothing
        maximize!(c.position_params, pos_data)
    else
        maximize!(c.position_params, pos_data; center_probs=nuclei_probs)
        if length(nuclei_probs) > 1
            c.confidence = quantile(nuclei_probs, 1 - min_nuclei_frac)
        end
    end

    if c.shape_prior !== nothing
        adjust_cov_by_prior!(c.position_params.Σ, c.shape_prior; n_samples=c.n_samples)
    end

    c.position_params = MvNormalF(c.position_params.μ, c.position_params.Σ) # TODO: make it consistent with immutability

    return c
end

pdf(comp::Component{N}, x::AbstractVector{Float64}, gene::Missing; use_smoothing::Bool=true) where N =
    comp.prior_probability * comp.confidence * pdf(comp.position_params, x)

pdf(comp::Component{N}, x::AbstractVector{Float64}, gene::Int64; use_smoothing::Bool=true) where N =
    comp.prior_probability * comp.confidence * pdf(comp.position_params, x) * pdf(comp.composition_params, gene; use_smoothing=use_smoothing)

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
