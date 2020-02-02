import Random
import MultivariateStats

### Main functions

"""
kwargs are passed to estimate_expression_clusters. Important parameters: min_cluster_size, n_pcs
"""
function score_doublets(data::BmmData; n_expression_clusters::Int, distance::T where T <: Distances.SemiMetric,
        min_molecules_per_cell::Int, kwargs...)
    cm = extract_gene_matrix_from_distributions2(data.components);
    real_cell_inds = findall(vec(sum(cm, dims=1)) .>= min_molecules_per_cell);

    cm_norm = cm ./ max.(sum(cm, dims=1), 1e-50);
    feature_mtx, cluster_centers = estimate_expression_clusters(cm_norm[:, real_cell_inds], n_expression_clusters; distance=distance, kwargs...)
    comb_improvements, mixing_fractions = factorize_doublet_expression(feature_mtx, cluster_centers; distance=distance)[3:4]

    all_scores = repeat(Union{Float64, Missing}[missing], length(data.assignment))
    all_fractions = repeat(Union{Float64, Missing}[missing], length(data.assignment))
    all_scores[real_cell_inds] .= comb_improvements
    all_fractions[real_cell_inds] .= 1 .- mixing_fractions

    return all_scores, all_fractions
end

"""
kwargs are passed to estimate_expression_clusters. Important parameters: min_cluster_size, n_pcs
"""
function split_components_by_expression!(data::BmmData, n_splitted_clusters::Int; n_expression_clusters::Int, distance::T where T <: Distances.SemiMetric,
        min_molecules_per_cell::Int, pvalue_threshold::Float64=0.01, improvement_threshold::Float64=0.1, kwargs...)
    cm = extract_gene_matrix_from_distributions2(data.components);
    real_cell_inds = findall(vec(sum(cm, dims=1)) .>= min_molecules_per_cell);

    cm_norm = cm ./ max.(sum(cm, dims=1), 1e-50);
    feature_mtx, cluster_centers = estimate_expression_clusters(cm_norm[:, real_cell_inds], n_expression_clusters; distance=distance, kwargs...)
    clust_per_cell = assign_to_centers(feature_mtx, cluster_centers, distance);
    clust_cm_centers = hcat([mean(cm_norm[:, real_cell_inds[ids]], dims=2) for ids in split_ids(clust_per_cell)]...);

    # Here can be update of gene priors

    base1_ids, base2_ids, comb_improvements = factorize_doublet_expression(feature_mtx, cluster_centers; distance=distance)[1:3]
    imp_cells = findall(comb_improvements .>= improvement_threshold)
    perm_pvals = estimate_spatial_separation_pvals(data, clust_cm_centers, base1_ids[imp_cells], base2_ids[imp_cells], real_cell_inds[imp_cells])

    for id in real_cell_inds[imp_cells[perm_pvals .< pvalue_threshold]]
        split_component!(data, id, n_splitted_clusters);
    end
end

### Algorithm

function split_component!(data::BmmData, component_id::Int, n_clusters::Int)
    mol_ids = findall(data.assignment .== component_id)
    if length(mol_ids) == 0
        @warn "Component $component_id is already empty"
        return data
    end

    init_params = cell_centers_with_clustering(data.x[mol_ids,:], n_clusters; scale=nothing); # TODO: real scale
    parent_comp = data.components[component_id];
    comps, prior, assignment = initial_distributions(data.x[mol_ids,:], init_params; size_prior=parent_comp.shape_prior,
        new_component_weight=data.distribution_sampler.prior_weight, gene_smooth=parent_comp.composition_params.smooth,
        gene_num=maximum(data.composition_data));

    parent_comp.n_samples = 0;
    data.assignment[mol_ids] .= length(data.components) .+ assignment;
    append!(data.components, comps);

    return data
end

function factorize_doublet_expression(feature_mtx::Matrix{Float64}, cluster_centers::Matrix{Float64}; distance::T where T <: Distances.SemiMetric)
    nn_dists = Distances.pairwise(Distances.CosineDist(), feature_mtx, cluster_centers; dims=2);
    base1_ids = vec(mapslices(x -> findmin(x)[2], nn_dists, dims=2));

    opt_comps = [[linsearch_gs(w -> distance(cluster_centers[:, base1_ids[ci]] .* w + cluster_centers[:, i] .* (1 - w), feature_mtx[:, ci]), 0.0, 1.0)
        for i in 1:size(cluster_centers, 2)] for ci in 1:size(feature_mtx, 2)];

    opt_dists = hcat([[v[2] for v in x] for x in opt_comps]...);

    base2_ids = vec(mapslices(x -> findmin(x)[2], opt_dists, dims=1))
    comb_improvements = vec(mapslices(x -> (maximum(x) - minimum(x)) / maximum(x), opt_dists, dims=1));
    mixing_fractions = hcat([[v[1] for v in x] for x in opt_comps]...)[CartesianIndex.(base2_ids, 1:length(base2_ids))];

    return base1_ids, base2_ids, comb_improvements, mixing_fractions
end

function estimate_expression_clusters(cm_norm::Matrix{Float64}, n_clusters::Int; n_pcs::Int=0, kwargs...)
    if n_pcs > 0
        pca = MultivariateStats.fit(MultivariateStats.PCA, cm_norm; maxoutdim=n_pcs);
        cm_norm = MultivariateStats.transform(pca, cm_norm);
    end

    k_centers = kmeans_stable(cm_norm, n_clusters; kwargs...)[1];

    return cm_norm, k_centers
end

function estimate_spatial_separation_pvals(data::BmmData, cm_cluster_centers::Matrix{Float64}, base1_id_per_cell::Vector{Int}, base2_id_per_cell::Vector{Int}, comp_id_per_cell::Vector{Int})
    perm_pvals = Float64[]
    probs_per_imp_cell = []

    for (b1i, b2i, ai) in zip(base1_id_per_cell, base2_id_per_cell, comp_id_per_cell)
        mol_ids = findall(data.assignment .== ai);
        probs = cm_cluster_centers[data.composition_data[mol_ids], b2i] ./ vec(sum(cm_cluster_centers[data.composition_data[mol_ids], [b2i, b1i]], dims=2));
        push!(perm_pvals, mass_center_permutation_test(copy(data.position_data[:, mol_ids]'), probs, n_iters=2000))
        push!(probs_per_imp_cell, probs)
    end

    return perm_pvals
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