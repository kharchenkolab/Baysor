using Random

function sample_distribution!(data::BmmData, shape_prior::ShapePrior)::MvNormalF
    μ = sample_center!(data) # TODO: should we stor it as non-static arrays from the beginning?
    Σ = Array(Diagonal(shuffle(sample_var(shape_prior))))

    return MvNormalF(μ, Σ)
end

function sample_center!(data::BmmData; cache_size::Int=10000)
    if length(data.center_sample_cache) == 0
        data.center_sample_cache = sample(1:size(data.x, 1), Weights(confidence(data)), cache_size)
    end

    return position_data(data)[:, pop!(data.center_sample_cache)]
end

function sample_composition_params(data::BmmData)
    # TODO: fix this
    gene_counts = sample(data.components).composition_params.counts
    return CategoricalSmoothed(gene_counts, smooth=data.distribution_sampler.composition_params.smooth, sum_counts=sum(gene_counts));
end

function maximize_from_prior!(comp::Component, data::BmmData)
    shape_prior = something(comp.shape_prior, data.distribution_sampler.shape_prior)
    comp.position_params = sample_distribution!(data, shape_prior);
    comp.composition_params = sample_composition_params(data);

    return comp;
end

function sample_distribution!(data::BmmData; guid::Int)
    sampler = data.distribution_sampler
    position_params = sample_distribution!(data, sampler.shape_prior);
    composition_params = sample_composition_params(data);

    return Component(position_params, composition_params; can_be_dropped=true, shape_prior=deepcopy(sampler.shape_prior), guid=guid);
end