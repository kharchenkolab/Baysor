import Random
import MultivariateStats

### Main functions

function score_doublets(data::BmmData; kwargs...)
    cm = extract_gene_matrix_from_distributions2(data.components)
    if :cluster_centers in keys(data.misc)
        return score_doublets(cm, data.misc[:cluster_centers]; kwargs...)
    end

    return score_doublets(cm; kwargs...)
end

function score_doublets(cm::Matrix{TR} where TR <: Real, n_expression_clusters::Int; distance::T where T <: Distances.SemiMetric, min_molecules_per_cell::Int, kwargs...)
    real_cell_inds = findall(vec(sum(cm, dims=1)) .>= min_molecules_per_cell);

    cm_norm = cm ./ max.(sum(cm, dims=1), 1e-50);
    feature_mtx, feature_space_centers = estimate_expression_clusters(cm_norm[:, real_cell_inds], n_expression_clusters; distance=distance, kwargs...)
    clust_per_cell = assign_to_centers(feature_mtx, feature_space_centers, distance);
    cluster_centers = hcat([mean(cm_norm[:, ids], dims=2) for ids in split(real_cell_inds, clust_per_cell)]...);

    return score_doublets(cm, cluster_centers; min_molecules_per_cell=min_molecules_per_cell, kwargs...)
end

"""
kwargs are passed to estimate_expression_clusters. Important parameters: min_cluster_size, n_pcs
"""
function score_doublets(cm::Matrix{TR} where TR <: Real, cluster_centers::Matrix{Float64}; min_molecules_per_cell::Int,
        na_value::Union{Float64, Missing, Nothing}=missing, min_impact_level::Float64=0.1, kwargs...)
    real_cell_inds = findall(vec(sum(cm, dims=1)) .>= min_molecules_per_cell);

    comb_improvements, mixing_fractions = factorize_doublet_expression_ll(Float64.(cm[:, real_cell_inds]), cluster_centers)[3:4]

    comb_improvements[comb_improvements .< min_impact_level] .= 0.;

    all_scores = repeat(Union{Float64, typeof(na_value)}[na_value], size(cm, 2))
    all_fractions = repeat(Union{Float64, typeof(na_value)}[na_value], size(cm, 2))
    all_pvals = repeat(Union{Float64, typeof(na_value)}[na_value], size(cm, 2))

    all_scores[real_cell_inds] .= comb_improvements
    all_fractions[real_cell_inds] .= 1 .- mixing_fractions

    return all_scores, all_fractions
end

function score_doublets_by_local_clusters(cell_assignment::Vector{Int}, cluster_assignment::Vector; na_value::Union{Float64, Missing, Nothing}=missing)
    mol_clusts_per_cell = split(denserank(cluster_assignment), cell_assignment .+ 1)[2:end];
    non_empty_cells = findall(length.(mol_clusts_per_cell) .> 0);
    all_scores = repeat(Union{Float64, typeof(na_value)}[na_value], length(mol_clusts_per_cell));
    all_scores[non_empty_cells] .= 1 .- maximum.(prob_array.(mol_clusts_per_cell[non_empty_cells]));
    return all_scores
end

"""
kwargs are passed to estimate_expression_clusters. Important parameters: min_cluster_size, n_pcs
"""
function split_components_by_expression!(data::BmmData, n_splitted_clusters::Int; n_expression_clusters::Int, min_molecules_per_cell::Int,
    distance::T where T <: Distances.SemiMetric = Distances.CosineDist(), pvalue_threshold::Float64=0.01, improvement_threshold::Float64=0.1, kwargs...)
    cm = extract_gene_matrix_from_distributions2(data.components);
    real_cell_inds = findall(vec(sum(cm, dims=1)) .>= min_molecules_per_cell);
    cluster_centers = nothing

    if :cluster_centers in keys(data.misc)
        cluster_centers = data.misc[:cluster_centers]
    else
        cm_norm = cm ./ max.(sum(cm, dims=1), 1e-50);
        feature_mtx, feature_space_centers = estimate_expression_clusters(cm_norm[:, real_cell_inds], n_expression_clusters; distance=distance, kwargs...)
        clust_per_cell = assign_to_centers(feature_mtx, feature_space_centers, distance);
        cluster_centers = hcat([mean(cm_norm[:, ids], dims=2) for ids in split(real_cell_inds, clust_per_cell)]...);
    end

    base1_ids, base2_ids, comb_improvements = factorize_doublet_expression_ll(cm[:, real_cell_inds], cluster_centers)[1:3]
    imp_cells = findall(comb_improvements .>= improvement_threshold)
    perm_pvals = estimate_spatial_separation_pvals(data, cluster_centers, base1_ids[imp_cells], base2_ids[imp_cells], real_cell_inds[imp_cells])

    for id in real_cell_inds[imp_cells[perm_pvals .< pvalue_threshold]]
        split_component!(data, id, n_splitted_clusters);
    end

    maximize!(data, min_molecules_per_cell)
end

### Algorithm

function split_component!(data::BmmData, component_id::Int, n_clusters::Int)
    mol_ids = findall(data.assignment .== component_id)
    if length(mol_ids) == 0
        @warn "Component $component_id is already empty"
        return data
    end

    init_params = cell_centers_with_clustering(position_data(data)[:, mol_ids], n_clusters; scale=data.distribution_sampler.shape_prior.std_values[1]);

    parent_comp = data.components[component_id];
    comps, prior, assignment = initial_distributions(data.x[mol_ids,:], init_params; size_prior=parent_comp.shape_prior,
        new_component_weight=data.distribution_sampler.prior_weight, gene_smooth=parent_comp.composition_params.smooth,
        gene_num=maximum(data.composition_data));

    parent_comp.n_samples = 0;
    data.assignment[mol_ids] .= length(data.components) .+ assignment;
    append!(data.components, comps);

    return data
end

function estimate_expression_clusters(feature_mtx::Matrix{Float64}, n_clusters::Int, real_cell_inds::T where T<: AbstractArray{Int, 1} = 1:size(feature_mtx, 2); n_pcs::Int=0, kwargs...)
    if n_pcs > 0
        pca = MultivariateStats.fit(MultivariateStats.PCA, feature_mtx[:, real_cell_inds]; maxoutdim=n_pcs);
        feature_mtx = MultivariateStats.transform(pca, feature_mtx);
    end

    k_centers = kmeans_stable(feature_mtx[:, real_cell_inds], n_clusters; kwargs...)[1];

    return feature_mtx, k_centers
end

function estimate_spatial_separation_pvals(data::BmmData, cm_cluster_centers::Matrix{Float64}, base1_id_per_cell::Vector{Int}, base2_id_per_cell::Vector{Int}, comp_id_per_cell::Vector{Int})
    perm_pvals = Float64[]
    probs_per_imp_cell = []

    for (b1i, b2i, ai) in zip(base1_id_per_cell, base2_id_per_cell, comp_id_per_cell)
        mol_ids = findall(data.assignment .== ai);
        probs = cm_cluster_centers[b2i, data.composition_data[mol_ids]] ./ vec(sum(cm_cluster_centers[[b2i, b1i], data.composition_data[mol_ids]], dims=1));
        push!(perm_pvals, mass_center_permutation_test(copy(data.position_data[:, mol_ids]'), probs, n_iters=2000))
        push!(probs_per_imp_cell, probs)
    end

    return perm_pvals
end

### Likelihood-based factorizaton

@inline log_ll(expression::Vector{T} where T<: Real, gene_probs::Vector{Float64}; prob_pseudocount::Float64=1e-10)::Float64 =
    -mean(expression .* log.(gene_probs .+ prob_pseudocount)) # With mean instead sum, likelihood scale is more meaningful, which is required for improvement_threshold

likelihood_dist(counts::Matrix{T} where T <: Real, gene_probs::Matrix{Float64})::Matrix{Float64} =
    hcat([[log_ll(counts[:, j], gene_probs[i,:]) for j in 1:size(counts, 2)] for i in 1:size(gene_probs, 1)]...)

function factorize_doublet_expression_ll(counts::Matrix{Float64}, cluster_centers::Matrix{Float64})
    nn_dists = likelihood_dist(counts, cluster_centers);
    base1_ids = vec(mapslices(x -> findmin(x)[2], nn_dists, dims=2));

    opt_comps = [[linsearch_gs(w -> log_ll(counts[:, ci], cluster_centers[base1_ids[ci], :] .* w + cluster_centers[i, :] .* (1 - w)), 0.0, 1.0)
        for i in 1:size(cluster_centers, 1)] for ci in 1:size(counts, 2)];

    opt_dists = hcat([[v[2] for v in x] for x in opt_comps]...);

    base2_ids = vec(mapslices(x -> findmin(x)[2], opt_dists, dims=1))
    # comb_improvements = vec(mapslices(x -> (maximum(x) - minimum(x)) / maximum(x), opt_dists, dims=1));
    comb_improvements = vec(mapslices(x -> 1 - exp(minimum(x) - maximum(x)), opt_dists, dims=1)); # (max(likelihood) - min(likelihood)) / max(likelihood)
    mixing_fractions = hcat([[v[1] for v in x] for x in opt_comps]...)[CartesianIndex.(base2_ids, 1:length(base2_ids))];

    return base1_ids, base2_ids, comb_improvements, mixing_fractions
end

### Utils

# function center_dist(pos_data::Matrix{Float64}, probs::Vector{Float64})
#     cl1_cent = [median(pos_data[:, i], ProbabilityWeights(probs)) for i in 1:2]
#     cl2_cent = [median(pos_data[:, i], ProbabilityWeights(1 .- probs)) for i in 1:2]
#     return sum((cl1_cent .- cl2_cent) .^ 2)
# end

function distance_between_mass_centers(pos_data::Matrix{Float64}, probs::Vector{Float64})
    cl1_cent = sum(pos_data .* probs, dims=1) ./ sum(probs);
    cl2_cent = sum(pos_data .* (1 .- probs), dims=1) ./ sum(1 .- probs);
    return sum((cl1_cent .- cl2_cent) .^ 2)
end

function mass_center_permutation_test(pos_data::Matrix{Float64}, probs::Vector{Float64}; n_iters::Int=1000)
    obs_dist = distance_between_mass_centers(pos_data, probs)
    perm_dists = [distance_between_mass_centers(pos_data, Random.shuffle(probs)) for i in 1:n_iters]
    return mean(obs_dist .< perm_dists)
end