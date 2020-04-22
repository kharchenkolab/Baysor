using Random

function sample_distribution!(data::BmmData, shape_prior::ShapePrior; center_prior::Union{CellCenter, Nothing}=nothing)::MvNormalF
    μ = (center_prior === nothing) ? sample_center!(data) : rand(MvNormal(Vector(center_prior.μ), Matrix(center_prior.Σ))) # TODO: should we stor it as non-static arrays from the beginning?
    Σ = Array(Diagonal(shuffle(sample_var(shape_prior))))

    return MvNormalF(μ, Σ)
end

function sample_center!(data::BmmData; cache_size::Int=10000)
    if length(data.center_sample_cache) == 0
        data.center_sample_cache = sample(1:size(data.x, 1), Weights(confidence(data)), cache_size)
    end
    μ = position_data(data)[:, pop!(data.center_sample_cache)]
    return μ
end

# DEPRECATED?
function sample_n_samples(data::BmmData; n_trials_max::Int=5)
    for i in 1:n_trials_max # Try sample first from all components, as filtration takes too long
        n_samples = sample(data.components).n_samples;
        if n_samples != 0
            return n_samples
        end
    end

    non_zero_comps = filter(c -> c.n_samples > 0, data.components)
    if length(non_zero_comps) == 0
        error("All components have n_samples == 0")
    end

    return sample(non_zero_comps).n_samples;
end

function sample_composition_params(data::BmmData)
    # n_samples = sample_n_samples(data)
    # gene_probs = data.gene_probs_given_single_transcript[:,sample(1:size(data.gene_probs_given_single_transcript, 2))];
    # gene_counts = rand(Multinomial(n_samples, gene_probs));

    # TODO: fix this
    samp_comp = sample(data.components);
    gene_counts = samp_comp.composition_params.counts .+ samp_comp.gene_count_prior
    return CategoricalSmoothed(gene_counts, smooth=data.distribution_sampler.composition_params.smooth, n_samples=sum(gene_counts));
end

function maximize_from_prior!(comp::Component, data::BmmData)
    sampler = data.distribution_sampler
    shape_prior = comp.shape_prior === nothing ? sampler.shape_prior : comp.shape_prior
    center_prior = comp.center_prior === nothing ? sampler.center_prior : comp.center_prior
    comp.position_params = sample_distribution!(data, shape_prior, center_prior=center_prior);
    comp.composition_params = sample_composition_params(data);

    return comp;
end

function sample_distribution(data::BmmData; guid::Int)
    sampler = data.distribution_sampler
    position_params = sample_distribution!(data, sampler.shape_prior, center_prior=sampler.center_prior);

    if :cluster_centers in keys(data.misc)
        gene_count_prior = data.misc[:cluster_centers][rand(1:size(data.misc[:cluster_centers], 1)),:]
        return Component(position_params, deepcopy(sample(data.components).composition_params);
                     prior_weight=sampler.prior_weight, can_be_dropped=true, gene_count_prior=round.(Int, gene_count_prior .* 100000),
                     center_prior=deepcopy(sampler.center_prior), shape_prior=deepcopy(sampler.shape_prior), guid=guid);
    end

    composition_params = sample_composition_params(data);

    return Component(position_params, composition_params; prior_weight=sampler.prior_weight, can_be_dropped=true,
                     center_prior=deepcopy(sampler.center_prior), shape_prior=deepcopy(sampler.shape_prior), guid=guid);
end