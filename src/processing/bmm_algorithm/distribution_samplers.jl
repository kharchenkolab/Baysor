using Random

function sample_position_params(μ::AbstractVector{<:Real}, shape_prior::ShapePrior{N})::MvNormalF{N} where N
    return MvNormalF(μ, Array(Diagonal(shuffle(sample_var(shape_prior)))))
end

function sample_centers!(data::BmmData, n::Int=1; cache_size::Int=20000)
    if length(data.center_sample_cache) < n # otherwise, sampling takes too long
        cache_size = max(cache_size, n)
        data.center_sample_cache = sample(1:size(data.x, 1), Weights(confidence(data)), cache_size)
    end

    ids = data.center_sample_cache[(end-n+1):end]
    data.center_sample_cache = data.center_sample_cache[1:(end-n)]
    return position_data(data)[:, ids]
end

function sample_composition_params(data::BmmData)
    return CategoricalSmoothed(
        shuffle!(deepcopy(sample(data.components).composition_params.counts)),
        smooth=data.distribution_sampler.composition_params.smooth
    )
end

function sample_distribution(data::BmmData, μ::AbstractVector{<:Real}, new_guid::Int)
    shape_prior = data.distribution_sampler.shape_prior
    position_params = sample_position_params(μ, shape_prior);
    composition_params = sample_composition_params(data);

    return Component(position_params, composition_params; shape_prior=deepcopy(shape_prior), guid=new_guid);
end

function sample_distributions!(data::BmmData, n::Int)
    μ = sample_centers!(data, n)
    dists = map(1:n) do i
        @spawn sample_distribution(data, μ[:,i], data.max_component_guid + i)
    end

    dists = fetch.(dists)
    data.max_component_guid += n
    return dists
end

function sample_distribution!(data::BmmData)
    data.max_component_guid += 1
    μ = sample_centers!(data, 1)[:,1]
    return sample_distribution(data, μ, data.max_component_guid)
end
