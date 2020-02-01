using NearestNeighbors
using DataFrames
using ProgressMeter

import Distances
import Distributions
import StatsBase
import MultivariateStats
import NearestNeighborDescent

"""
    estimate_gene_probs_given_single_transcript(cm, n_molecules_per_cell)

Compute probabilities of gene expression given expression of a transcript: p(g_i | t_k) = P[i, k]
"""
function estimate_gene_probs_given_single_transcript(cm::Array{Float64, 2}, n_molecules_per_cell::Array{Int, 1})::Array{Float64, 2}
    prior_cell_probs = n_molecules_per_cell ./ sum(n_molecules_per_cell);
    probs = Array{Float64, 1}[]

    for t in 1:size(cm, 1)
        prior_transcript_probs = cm[t,:] .* prior_cell_probs;
        push!(probs, [sum(cm[g,:] .* prior_transcript_probs) for g in 1:size(cm, 1)] ./ sum(prior_transcript_probs))
    end

    return hcat(probs...)
end

function extract_gene_matrix_from_distributions(components::Array{Component, 1}, n_genes::Int)::Array{Float64, 2}
    counts_per_cell = [counts(c.composition_params) for c in components];
    for i in 1:size(counts_per_cell, 1)
        n_genes_cur = size(counts_per_cell[i], 1)
        if n_genes_cur < n_genes
            append!(counts_per_cell[i], zeros(n_genes - n_genes_cur))
        end
    end

    count_matrix = hcat(counts_per_cell...);
    return count_matrix ./ sum(count_matrix, dims=1);
end

function knn_by_expression(count_matrix::Array{T, 2} where T <: Real, n_molecules_per_cell::Array{Int, 1};
                           k::Int=15, min_molecules_per_cell::Int=10, n_prin_comps::Int=0)::Array{Array{Int, 1}, 1}
    count_matrix_norm = count_matrix ./ sum(count_matrix, dims=1)
    neighborhood_matrix = count_matrix_norm
    if n_prin_comps > 0
        pca = MultivariateStats.fit(MultivariateStats.PCA, count_matrix_norm; maxoutdim=n_prin_comps);
        neighborhood_matrix = MultivariateStats.transform(pca, count_matrix_norm)
    end

    if maximum(n_molecules_per_cell) < min_molecules_per_cell
        @warn "No cells pass min_molecules threshold ($min_molecules_per_cell). Resetting it to 1"
        min_molecules_per_cell = 1
    end

    real_cell_inds = findall(n_molecules_per_cell .>= min_molecules_per_cell);

    if length(real_cell_inds) < k
        @warn "Number of large cells ($(length(real_cell_inds))) is lower than the requested number of nearest neighbors ($k)"
        k = length(real_cell_inds)
    end

    # neighb_inds = knn(KDTree(neighborhood_matrix[:, real_cell_inds]), neighborhood_matrix, k, true)[1]
    kd_tree = KDTree(neighborhood_matrix[:, real_cell_inds]);
    neighb_inds = vcat(pmap(inds -> knn(kd_tree, neighborhood_matrix[:,inds], k, true)[1], split(1:length(real_cell_inds), n_parts=nprocs()))...)
    # neighb_inds = Array{Array{Array{Int, 1}, 1}, 1}(undef, Threads.nthreads());
    # splitted_inds = split(1:size(neighborhood_matrix, 2), n_parts=length(neighb_inds));
    # @threads for i in 1:length(splitted_inds)
    #     neighb_inds[i] = knn(kd_tree, neighborhood_matrix[:,splitted_inds[i]], k, true)[1]
    # end
    # neighb_inds = vcat(neighb_inds...);

    return [real_cell_inds[inds] for inds in neighb_inds]
end

function update_gene_count_priors!(components::Array{Component, 1}, neighb_inds::Array{Array{Int, 1}, 1})
    for (comp_dst, inds) in zip(components, neighb_inds)
        comp_dst.gene_count_prior .= 0
        for comp_src in components[inds]
            comp_dst.gene_count_prior .+= counts(comp_src.composition_params)
        end

        comp_dst.gene_count_prior_sum = sum(comp_dst.gene_count_prior)
    end
end

function smooth_size_prior_knn!(components::Array{Component, 1}, neighb_inds::Array{Array{Int, 1}, 1})
    sizes_per_cell = [hcat(eigen_values.(components[inds])...) for inds in neighb_inds];
    mean_sizes_per_cell = [vec(mapslices(trim_mean, sizes, dims=2)) for sizes in sizes_per_cell];

    for (prior_means, c) in zip(mean_sizes_per_cell, components)
        set_shape_prior!(c, prior_means)
    end
end

function smooth_size_prior_global!(bm_data_arr::Array{BmmData, 1}; set_individual_priors::Bool=false)
    sizes_per_cell = hcat(vcat([eigen_values(bm_data.components) for bm_data in bm_data_arr]...)...)
    n_mols_per_cell = vcat(num_of_molecules_per_cell.(bm_data_arr)...)

    n_min, n_max = minimum(n_mols_per_cell), maximum(n_mols_per_cell)
    threshold = 0.01 * (n_max - n_min)
    mean_prior = vec(median(sizes_per_cell[:, (n_mols_per_cell .>= (threshold + n_min)) .& (n_mols_per_cell .<= (threshold + n_max))], dims=2))

    for bm_data in bm_data_arr
        set_shape_prior!(bm_data.distribution_sampler, mean_prior)

        for c in vcat([bmd.components for bmd in bm_data_arr]...)
            if set_individual_priors || c.n_samples == 0
                set_shape_prior!(c, mean_prior)
            end
        end
    end

    return bm_data_arr
end

function update_gene_prior!(bmm_data_arr::Array{BmmData,1}, count_matrix::Array{Float64, 2}, n_molecules_per_cell::Array{Int, 1})
    probs = estimate_gene_probs_given_single_transcript(count_matrix, n_molecules_per_cell)
    for bm_data in bmm_data_arr
        bm_data.gene_probs_given_single_transcript = probs
    end

    return bmm_data_arr
end

function update_priors!(bmm_data_arr::Array{BmmData,1}; use_cell_type_size_prior::Bool, use_global_size_prior::Bool, smooth_expression::Bool,
                        min_molecules_per_cell::Int, n_prin_comps::Int)
    n_molecules_per_cell = vcat(num_of_molecules_per_cell.(bmm_data_arr)...);
    components = vcat([bm.components for bm in bmm_data_arr]...);

    components = components[n_molecules_per_cell .> 0]
    n_molecules_per_cell = n_molecules_per_cell[n_molecules_per_cell .> 0]

    count_matrix = extract_gene_matrix_from_distributions(components, maximum([maximum(ed.x[:gene]) for ed in bmm_data_arr]));

    if use_cell_type_size_prior || smooth_expression
        neighb_inds = knn_by_expression(count_matrix, n_molecules_per_cell, min_molecules_per_cell=min_molecules_per_cell, n_prin_comps=n_prin_comps)

        if use_cell_type_size_prior
            smooth_size_prior_knn!(components, neighb_inds);
        end

        if smooth_expression
            update_gene_count_priors!(components, neighb_inds);
        end
    end

    smooth_size_prior_global!(bmm_data_arr, set_individual_priors=(!use_cell_type_size_prior && use_global_size_prior))
    update_gene_prior!(bmm_data_arr, count_matrix, n_molecules_per_cell)

    return bmm_data_arr
end

## NEW

function extract_gene_matrix_from_distributions2(components::Array{Component, 1}, n_genes::Int=length(components[1].composition_params.counts))::Array{Float64, 2}
    counts_per_cell = [counts(c.composition_params) for c in components];
    for i in 1:size(counts_per_cell, 1)
        n_genes_cur = size(counts_per_cell[i], 1)
        if n_genes_cur < n_genes
            append!(counts_per_cell[i], zeros(n_genes - n_genes_cur))
        end
    end

    res_mat = hcat(counts_per_cell...);
    res_mat[:, [c.n_samples == 0 for c in components]] .= 0

    return res_mat;
end

function knn_by_expression2(count_matrix::Array{T, 2} where T <: Real,
                            n_molecules_per_cell::Array{Int, 1}=round.(Int, vec(sum(count_matrix, dims=1)));
                            k::Int=15, min_molecules_per_cell::Int=10, n_prin_comps::Int=0, refine::Bool=true,
                            refine_mult::Real = 2, use_correlations::Bool=false)::Array{Array{Int, 1}, 1}
    count_matrix_norm = count_matrix ./ max.(sum(count_matrix, dims=1), 1e-50)
    # count_matrix_norm = log2.(1e-6 .+ count_matrix ./ sum(count_matrix, dims=1))
    neighborhood_matrix = count_matrix_norm

    real_cell_inds = findall(n_molecules_per_cell .>= min_molecules_per_cell);

    if n_prin_comps > 0
        pca = MultivariateStats.fit(MultivariateStats.PCA, count_matrix_norm[:, real_cell_inds]; maxoutdim=n_prin_comps);
        neighborhood_matrix = MultivariateStats.transform(pca, count_matrix_norm)
    end

    if maximum(n_molecules_per_cell) < min_molecules_per_cell
        @warn "No cells pass min_molecules threshold ($min_molecules_per_cell). Resetting it to 1"
        min_molecules_per_cell = 1
    end

    if length(real_cell_inds) < k
        @warn "Number of large cells ($(length(real_cell_inds))) is lower than the requested number of nearest neighbors ($k)"
        k = length(real_cell_inds)
    end

    neighb_inds = Int[]
    k0 = (refine ? refine_mult * k : k)
    if use_correlations
        p_dists = Distances.pairwise(Distances.CosineDist(), neighborhood_matrix, neighborhood_matrix[:, real_cell_inds], dims=2);
        res = [real_cell_inds[inds] for inds in eachrow(mapslices(sortperm, p_dists, dims=2)[:,1:(k0 + 1)])]
        neighb_inds = [v[v .!= i][1:k] for (i, v) in enumerate(res)]
    else
        kd_tree = KDTree(neighborhood_matrix[:, real_cell_inds]);
        neighb_inds = [real_cell_inds[ids[2:end]] for ids in knn(kd_tree, neighborhood_matrix, k0 + 1, true)[1]]
    end

    if !refine
        return neighb_inds
    end

    return refine_knn_by_expression.(Ref(neighborhood_matrix), neighb_inds, k=k)
end

"""
    This function is very sensitive to proper selection of metric and reduction for NN search
"""
function knn_by_expression_robust(count_matrix::Array{T, 2} where T <: Real,
                                  n_molecules_per_cell::Array{Int, 1}=round.(Int, vec(sum(count_matrix, dims=1)));
                                  k::Int=10, min_molecules_per_cell::Int=10, distance::Union{T2, Nothing} where T2 <: Distances.SemiMetric = nothing,
                                  n_prin_comps::Union{Int, Nothing}=nothing, approximate::Bool=false, kwargs...) #::Array{Array{Int, 1}, 1}
    neighborhood_matrix = count_matrix ./ max.(sum(count_matrix, dims=1), 1e-50)
    real_cell_inds = findall(n_molecules_per_cell .>= min_molecules_per_cell);

    if maximum(n_molecules_per_cell) < min_molecules_per_cell
        @warn "No cells pass min_molecules threshold ($min_molecules_per_cell). Resetting it to 1"
        min_molecules_per_cell = 1
    end

    if length(real_cell_inds) < k
        @warn "Number of large cells ($(length(real_cell_inds))) is lower than the requested number of nearest neighbors ($k)"
        k = length(real_cell_inds)
    end

    if n_prin_comps === nothing
        n_prin_comps = max(min(ceil(Int, 0.75 * median(sum(count_matrix .> 0, dims=1))), 100), 3)
    end

    if n_prin_comps > 0
        pca = MultivariateStats.fit(MultivariateStats.PCA, neighborhood_matrix[:, real_cell_inds]; maxoutdim=n_prin_comps);
        neighborhood_matrix = MultivariateStats.transform(pca, neighborhood_matrix)
    end

    if distance === nothing
        distance = (size(neighborhood_matrix, 1) > 20) ? Distances.CosineDist() : Distances.Euclidean()
    end

    if approximate
        graph = NearestNeighborDescent.DescentGraph(neighborhood_matrix[:, real_cell_inds], k, distance; kwargs...);
        g_ids, g_dists = NearestNeighborDescent.search(graph, neighborhood_matrix, k);
        nn_per_cell = [real_cell_inds[ids[sortperm(dists)]] for (ids, dists) in zip(eachcol(g_ids), eachcol(g_dists))];
    else
        if typeof(distance) <: Distances.Metric
            kd_tree = KDTree(neighborhood_matrix[:, real_cell_inds]);
            nn_per_cell = [real_cell_inds[ids] for ids in knn(kd_tree, neighborhood_matrix, k + 1, true)[1]]
        else
            @warn "Using bruteforce NN computation. It can take long..."
            p_dists = Distances.pairwise(distance, neighborhood_matrix, neighborhood_matrix[:, real_cell_inds], dims=2);
            nn_per_cell = [real_cell_inds[inds] for inds in eachrow(mapslices(sortperm, p_dists, dims=2)[:,1:(k + 1)])]
        end

        nn_per_cell = [v[v .!= i][1:k] for (i, v) in enumerate(nn_per_cell)]
    end

    # return nn_per_cell

    nn_variance_per_cell = [quantile(vec(Distances.pairwise(distance, neighborhood_matrix[:, ids], dims=2)), 0.95) for ids in nn_per_cell];
    non_doublet_cells = unique([ids[findmin(nn_variance_per_cell[ids])[2]] for ids in nn_per_cell]);
    is_cell_non_doublet = falses(length(nn_variance_per_cell));
    is_cell_non_doublet[non_doublet_cells] .= true;

    non_doublet_nn_per_cell = [ids[findfirst(is_cell_non_doublet[ids])] for ids in nn_per_cell];
    # return nn_per_cell #, is_cell_non_doublet
    # return [ids[is_cell_non_doublet[ids]] for ids in nn_per_cell];
    # return [[ids[findfirst(is_cell_non_doublet[ids])]] for ids in nn_per_cell]
    return nn_per_cell[non_doublet_nn_per_cell] #, nn_per_cell, [ids[is_cell_non_doublet[ids]] for ids in nn_per_cell];
end

function refine_knn_by_expression(neighborhood_matrix::Matrix{Float64}, neighb_inds::Vector{Int}; k::Int, max_iter::Int = 100)
    kd_tree = KDTree(neighborhood_matrix[:, neighb_inds]);
    adj_cells = neighb_inds[1:k];

    for i in 1:max_iter
        adj_cells_prev = adj_cells;
        adj_cells = neighb_inds[knn(kd_tree, vec(median(neighborhood_matrix[:, adj_cells], dims=2)), k, true)[1]];

        if all(sort(adj_cells) .== sort(adj_cells_prev))
            break
        end
    end

    return adj_cells
end

"""
    k-means based
"""
function update_gene_count_priors!(components::Vector{Component}; min_molecules_per_cell::Int, n_pcs::Int, n_clusters::Int, min_cluster_size::Int,
                                   distance=Distances.Euclidean())
    cm = extract_gene_matrix_from_distributions2(components);
    neighborhood_matrix = cm ./ max.(sum(cm, dims=1), 1e-50);

    real_cell_inds = findall(vec(sum(cm, dims=1)) .>= min_molecules_per_cell);

    pca = MultivariateStats.fit(MultivariateStats.PCA, neighborhood_matrix[:, real_cell_inds]; maxoutdim=n_pcs);
    neighborhood_matrix = MultivariateStats.transform(pca, neighborhood_matrix);

    k_centers = kmeans_stable(neighborhood_matrix[:, real_cell_inds], n_clusters, n_inits=10, min_cluster_size=min_cluster_size)[1];
    clust_per_cell = assign_to_centers(neighborhood_matrix, k_centers, distance);

    var_to_mean_ratios = hcat([vec(var(cm[:, ids], dims=2) ./ max.(mean(cm[:, ids], dims=2), 1e-20)) for ids in split_ids(clust_per_cell)]...)
    prior_rate_per_clust = 1 ./ max.(var_to_mean_ratios, 1.0);
    clust_centers = hcat([median(cm[:, ids], dims=2) for ids in Baysor.split_ids(clust_per_cell)]...);

    for (i, c) in enumerate(components)
        c_clust = clust_per_cell[i]
        cur_prior_rate = prior_rate_per_clust[:, c_clust]
        c.gene_count_prior = round.(Int, clust_centers[:, c_clust] .* cur_prior_rate .+ c.composition_params.counts .* (1 .- cur_prior_rate))
        c.gene_count_prior_sum = sum(c.gene_count_prior)
    end
end