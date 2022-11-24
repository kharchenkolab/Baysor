using Random

function sample_position_params!(data::BmmData, shape_prior::ShapePrior{N})::MvNormalF{N} where N
    μ = sample_center!(data) # TODO: should we store it as non-static arrays from the beginning?
    Σ = Array(Diagonal(shuffle(sample_var(shape_prior))))

    return MvNormalF(μ, Σ)
end

function sample_center!(data::BmmData; cache_size::Int=10000)
    if length(data.center_sample_cache) == 0 # otherwise, sampling takes too long
        data.center_sample_cache = sample(1:size(data.x, 1), Weights(confidence(data)), cache_size)
    end

    return position_data(data)[:, pop!(data.center_sample_cache)]
end

function sample_composition_params(data::BmmData)
    gene_counts = shuffle!(deepcopy(sample(data.components).composition_params.counts))
    return CategoricalSmoothed(gene_counts, smooth=data.distribution_sampler.composition_params.smooth, sum_counts=sum(gene_counts));
end

function sample_distribution!(data::BmmData)
    data.max_component_guid += 1
    shape_prior = data.distribution_sampler.shape_prior
    position_params = sample_position_params!(data, shape_prior);
    composition_params = sample_composition_params(data);

    return Component(position_params, composition_params; shape_prior=deepcopy(shape_prior), guid=data.max_component_guid);
end