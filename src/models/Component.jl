using DataFrames
using Distributions
using NearestNeighbors
using StatsBase

struct CellCenter
    μ::Array{Float64, 1};
    Σ::Array{Float64, 2};
    n_degrees_of_freedom::Int;
    # confidence::Float64; # TODO: uncomment and account for these parameters
    # distance_to_others::Array{Float64, 1}; # Probability of the center to coexist with the other center

    # CellCenter(μ::Array{Float64, 1}, Σ::Array{Float64, 2}, n_degrees_of_freedom::Int) = new(μ, Σ, n_degrees_of_freedom, 1.0, Float64[])
end

struct ShapePrior
    std_values::Array{Float64, 1};
    std_value_stds::Array{Float64, 1};
    n_samples::Int;
end

distributions(prior::ShapePrior) = Normal.(prior.std_values, prior.std_value_stds)
rand(prior::ShapePrior) = rand.(distributions(prior))

var_posterior(prior::ShapePrior, eigen_values::Array{Float64, 1}; n_samples::Int) =
    var_posterior.(prior.std_values, prior.std_value_stds, eigen_values; n_samples=n_samples, prior_n_samples=prior.n_samples)

# f(dx) = sign(dx) * sqrt(|dx/std|) * std
function var_posterior(prior_std::Float64, prior_std_std::Float64, eigen_value::Float64; n_samples::Int, prior_n_samples::Int)
    d_std = (sqrt(eigen_value) - prior_std)
    std_adj = (prior_std + sign(d_std) * (sqrt(abs(d_std) / prior_std_std + 1) - 1) * prior_std_std)

    return ((prior_n_samples * prior_std + n_samples * std_adj) / (prior_n_samples + n_samples)) ^ 2
end

function sample_var(d::Normal)
    @assert d.μ > 0
    while true
        v = rand(d)
        if v > 0
            return v
        end
    end
end
sample_var(prior::ShapePrior) = sample_var.(distributions(prior))

# pdf(prior::ShapePrior, Σ::Array{Float64, 2}) = pdf.(distributions(prior), eigen(Σ).values)
# logpdf(prior::ShapePrior, Σ::Array{Float64, 2}) = logpdf.(distributions(prior), eigen(Σ).values)

# struct GeneCountPrior
#     vectors_per_gene_type::Array{Array{Int, 2}, 1};
#     expected_prob_per_gene_type::Array{Float64, 1};

#     function GeneCountPrior(vectors_per_gene_type::Array{Array{Int, 2}, 1})
#         gp = new()
#     end
# end

mutable struct Component
    position_params::MvNormal;
    composition_params::SingleTrialMultinomial;
    n_samples::Int;
    prior_weight::Float64;
    prior_probability::Float64;
    can_be_dropped::Bool;

    center_prior::Union{Nothing, CellCenter};
    shape_prior::Union{Nothing, ShapePrior};

    gene_count_prior::Array{Int, 1};
    gene_count_prior_sum::Int;

    guid::Int;

    Component(position_params::MvNormal, composition_params::SingleTrialMultinomial; prior_weight::Float64, can_be_dropped::Bool,
              n_samples::Int=0, center_prior::Union{Nothing, CellCenter}=nothing, shape_prior::Union{Nothing, ShapePrior}=nothing,
              gene_count_prior::Array{Int, 1}=zeros(Int, length(counts(composition_params))), guid::Int=-1) =
        new(position_params, composition_params, n_samples, prior_weight, 1.0, can_be_dropped, center_prior, shape_prior, gene_count_prior, sum(gene_count_prior), guid)
end

function maximize!(c::Component, pos_data::Array{Float64, 2}, comp_data::Array{Int, 1})
    maximize!(c.composition_params, comp_data);

    pos_new = maximize(c.position_params, pos_data);
    μ, Σ = pos_new.μ, Matrix(pos_new.Σ)

    if c.center_prior !== nothing
        μ, Σ = normal_posterior(μ, c.center_prior.μ, Σ, c.center_prior.Σ, n=size(pos_data, 2), n_prior=c.center_prior.n_degrees_of_freedom)
        # μ = normal_posterior(μ, c.center_prior.μ, Σ, c.center_prior.Σ, size(pos_data, 2))[1]
    end

    if c.shape_prior !== nothing
        adjust_cov_by_prior!(Σ, c.shape_prior; n_samples=size(pos_data, 2))
    end

    try
        c.position_params = MvNormal(μ, Σ)
    catch
        error("Can't maximize position params. μ: '$μ', Σ: '$Σ', n_samples: $(size(pos_data ,1)).")
    end

    return c
end

# pdf(params::Component, x::Float64, y::Float64, gene::Int64; use_smoothing::Bool=true)::Float64 =
#     pdf(params.position_params, x, y) *
#     pdf(params.composition_params, gene, params.gene_count_prior, params.gene_count_prior_sum; use_smoothing=use_smoothing)

function pdf(params::Component, x::Float64, y::Float64, gene::Int64; use_smoothing::Bool=true)::Float64
    pos_pdf = pdf(params.position_params, x, y)

    if params.gene_count_prior_sum == 0
        return pos_pdf * pdf(params.composition_params, gene; use_smoothing=use_smoothing)
    else
        return pos_pdf * pdf(params.gene_count_prior, params.gene_count_prior_sum, gene, params.composition_params.smooth; use_smoothing=use_smoothing)
    end
end


function adjust_cov_by_prior!(Σ::Array{Float64, 2}, prior::ShapePrior; n_samples::Int)
    fact = eigen(Σ)
    eigen_values_posterior = var_posterior(prior, fact.values; n_samples=n_samples)
    Σ .= fact.vectors * diagm(0 => eigen_values_posterior) * inv(fact.vectors)
    Σ[1, 2] = Σ[2, 1]

    return adjust_cov_matrix!(Σ)
end

function normal_posterior(μ::AbstractArray{Float64, 1}, μ_prior::AbstractArray{Float64, 1}, Σ::AbstractArray{Float64, 2}, Σ_prior::AbstractArray{Float64, 2};
                          n::Int, n_prior::Int)
    μ_post = (n .* μ .+ n_prior .* μ_prior) ./ (n + n_prior)
    Σ_post = (n .* Σ .+ n_prior .* Σ_prior .+ (n * n_prior / (n + n_prior)) .* ((μ .- μ_prior) * (μ .- μ_prior)')) ./ (n + n_prior)

    return μ_post, adjust_cov_matrix!(Σ_post)
end

position(c::Component) = c.position_params.μ

# set_shape_prior!(c::Component, var_arr::Array{Float64, 1}) = begin c.shape_prior = ShapePrior(c.shape_prior.n_samples_var, var_arr) end
set_shape_prior!(c::Component, var_arr::Array{Float64, 1}) = error("Not implemented")

eigen_values(c::Component) = eigen(shape(c.position_params)).values