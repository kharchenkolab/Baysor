using Random

function sample_distribution(data::BmmData, shape_prior::ShapePrior; center_prior::Union{CellCenter, Nothing}=nothing)::MvNormal
    μ = (center_prior === nothing) ? position_data(data)[:, sample(1:size(data.x, 1))] : rand(MvNormal(center_prior.μ, center_prior.Σ))
    Σ = Array(Diagonal(shuffle(sample_var(shape_prior))))

    return MvNormal(μ, Σ)
end

function sample_n_samples(data::BmmData; n_trials_max::Int=500)
    for i in 1:n_trials_max
        n_samples = sample(data.components).n_samples;
        if n_samples != 0
            return n_samples
        end
    end

    error("Too many components have n_samples == 0")
end

function sample_composition_params(data::BmmData)
    n_samples = sample_n_samples(data)
    # gene_probs = data.gene_probs_given_single_transcript[:,sample(1:size(data.gene_probs_given_single_transcript, 2))];
    # gene_counts = rand(Multinomial(n_samples, gene_probs));
    samp_comp = sample(data.components);
    gene_counts = samp_comp.composition_params.counts .+ samp_comp.gene_count_prior
    return SingleTrialMultinomial(gene_counts, smooth=data.distribution_sampler.composition_params.smooth, n_samples=n_samples);
end

function maximize_from_prior!(comp::Component, data::BmmData)
    sampler = data.distribution_sampler
    comp.position_params = sample_distribution(data, sampler.shape_prior, center_prior=sampler.center_prior);
    comp.composition_params = sample_composition_params(data);

    return comp;
end

function sample_distribution(data::BmmData)
    sampler = data.distribution_sampler
    position_params = sample_distribution(data, sampler.shape_prior, center_prior=sampler.center_prior);
    composition_params = sample_composition_params(data);

    return Component(position_params, composition_params; prior_weight=sampler.prior_weight, can_be_dropped=true,
                     center_prior=deepcopy(sampler.center_prior), shape_prior=deepcopy(sampler.shape_prior));
end