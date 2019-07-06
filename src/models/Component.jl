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
    n_samples_var::Int;
    eigen_values::Array{Float64, 1};
end

distributions(prior::ShapePrior) = Gamma.(Ref(prior.n_samples_var), prior.eigen_values ./ prior.n_samples_var)

# argmax of posterior of Inverse Chi-squared
var_posterior(prior::ShapePrior, eigen_values::Array{Float64, 1}, n::Int) =
    [(prior.n_samples_var * pv + n * ev) / (prior.n_samples_var + n) for (pv, ev) in zip(prior.eigen_values, eigen_values)]

sample_var(prior::ShapePrior, args...) = [rand(d, args...) for d in distributions(prior)]
# pdf(prior::ShapePrior, Σ::Array{Float64, 2}) = pdf.(distributions(prior), eigen(Σ).values)
logpdf(prior::ShapePrior, Σ::Array{Float64, 2}) = logpdf.(distributions(prior), eigen(Σ).values)

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

    Component(position_params::MvNormal, composition_params::SingleTrialMultinomial; prior_weight::Float64, can_be_dropped::Bool,
              n_samples::Int=0, center_prior::Union{Nothing, CellCenter}=nothing, shape_prior::Union{Nothing, ShapePrior}=nothing,
              gene_count_prior::Array{Int, 1}=zeros(Int, length(counts(composition_params)))) =
        new(position_params, composition_params, n_samples, prior_weight, 1.0, can_be_dropped, center_prior, shape_prior, gene_count_prior, sum(gene_count_prior))
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
        Σ = adjust_cov_by_prior(Σ, size(pos_data, 2), c.shape_prior)
    end

    try
        c.position_params = MvNormal(μ, Σ)
    catch
        @show Σ
    end

    return c
end

pdf(params::Component, x::Float64, y::Float64, gene::Int64; use_smoothing::Bool=true)::Float64 =
    pdf(params.position_params, x, y) *
    pdf(params.composition_params, gene, use_smoothing=use_smoothing, prior=params.gene_count_prior, prior_sum=params.gene_count_prior_sum)

function adjust_cov_by_prior(Σ::Array{Float64, 2}, n_samples::Int, prior::ShapePrior)
    fact = eigen(Σ)
    eigen_values_posterior = var_posterior(prior, fact.values, n_samples)
    d = fact.vectors * diagm(0 => eigen_values_posterior) * inv(fact.vectors)
    d[1, 2] = d[2, 1]

    return adjust_cov_matrix(d)
end

function normal_posterior(μ::AbstractArray{Float64, 1}, μ_prior::AbstractArray{Float64, 1}, Σ::AbstractArray{Float64, 2}, Σ_prior::AbstractArray{Float64, 2};
                          n::Int, n_prior::Int)
    μ_post = (n .* μ .+ n_prior .* μ_prior) ./ (n + n_prior)
    Σ_post = (n .* Σ .+ n_prior .* Σ_prior .+ (n * n_prior / (n + n_prior)) .* ((μ .- μ_prior) * (μ .- μ_prior)')) ./ (n + n_prior)

    return μ_post, adjust_cov_matrix(Σ_post)
end

position(c::Component) = c.position_params.μ

set_shape_prior!(c::Component, var_arr::Array{Float64, 1}) = begin c.shape_prior = ShapePrior(c.shape_prior.n_samples_var, var_arr) end

eigen_values(c::Component) = eigen(shape(c.position_params)).values