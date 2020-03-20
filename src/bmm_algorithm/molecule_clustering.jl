import Clustering
using ProgressMeter
using DataFrames
using Statistics
using StatsBase
using Random

## Full EM. Probably will be needed when scRNA-seq is used as prior

function maximize_molecule_clusters!(cell_type_exprs::Matrix{Float64}, cell_type_exprs_norm::Matrix{Float64}, genes::Vector{Int},
        assignment_probs::Matrix{Float64})
    cell_type_exprs .= 0.0;
    for i in 1:length(genes)
        t_gene = genes[i];

        for j in 1:size(cell_type_exprs, 1)
            cell_type_exprs[j, t_gene] += assignment_probs[j, i];
        end
    end

    # TODO: here should be adjustment based on prior from scRNA-seq using adj_value_norm

    cell_type_exprs .= (cell_type_exprs .+ 1) ./ (sum(cell_type_exprs, dims=2) .+ 1);
    cell_type_exprs_norm .= cell_type_exprs ./ sum(cell_type_exprs, dims=2);
end

function expect_molecule_clusters!(assignment_probs::Matrix{Float64}, cell_type_exprs::Matrix{Float64}, cell_type_exprs_norm::Matrix{Float64}, genes::Vector{Int},
        confidence::Vector{Float64}, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}; new_prob::Float64=0.05)
    total_ll = 0.0
    for i in 1:length(genes)
        gene = genes[i]
        cur_weights = adjacent_weights[i]
        cur_points = adjacent_points[i]

        assignment_probs[:, i] .= 0.0
        dense_sum = 0.0
        for ri in 1:size(assignment_probs, 1)
            c_d = 0.0
            adj_prob = 1.0 # probability that there are no neighbors from this component
            for j in 1:length(cur_points)
                a_p = assignment_probs[ri, cur_points[j]]
                c_d += confidence[cur_points[j]] * cur_weights[j] * a_p
                adj_prob *= 1 - a_p
            end

            assignment_probs[ri, i] = cell_type_exprs[ri, gene] * exp(c_d) * (1 - adj_prob)
            dense_sum += assignment_probs[ri, i]
        end

        total_ll += log10(dense_sum)
        for ri in 1:size(assignment_probs, 1)
            assignment_probs[ri, i] = (1 - new_prob) * assignment_probs[ri, i] / dense_sum + new_prob * cell_type_exprs_norm[ri, gene]
        end
    end

    return total_ll
end

function cluster_molecules_on_mrf(df_spatial::DataFrame, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}};
        n_clusters::Int, nn_num::Int=30, n_mols_init::Int=min(max(5000, round(Int, sqrt(size(df_spatial, 1)) * 30)), 50000), confidence_threshold::Float64=0.95, kwargs...)
    c_mol_ids = select_ids_uniformly(df_spatial.x::Vector{Float64}, df_spatial.y::Vector{Float64}, df_spatial.confidence::Vector{Float64},
        n_mols_init, confidence_threshold=confidence_threshold);
    # gene_ids_per_adj = getindex.(Ref(df_spatial.gene), adjacent_points)[c_mol_ids];

    # TODO: better initialization algorithm faster than O(n^2)
    t_pd = Baysor.position_data(df_spatial)
    gene_ids_per_adj = getindex.(Ref(df_spatial.gene), knn(KDTree(t_pd[:, df_spatial.confidence .> confidence_threshold]), t_pd[:, c_mol_ids], nn_num)[1])

    dist_mat = pairwise_jaccard(gene_ids_per_adj);
    t_clusts = Clustering.cutree(Clustering.hclust(dist_mat, linkage=:ward), k=n_clusters);
    ct_exprs_init = copy(hcat([prob_array(vcat(x...), max_value=maximum(df_spatial.gene)) for x in split(gene_ids_per_adj, t_clusts)]...)');

    return cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights; cell_type_exprs=ct_exprs_init, kwargs...)
end

function cluster_molecules_on_mrf(genes::Vector{Int}, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}},
        confidence::Vector{Float64}=ones(length(genes)); k::Int=1, new_prob::Float64=0.05, tol::Float64=0.01, do_maximize::Bool=true,
        min_iters::Int=div(length(genes), 300), max_iters::Int=3*min_iters,
        cell_type_exprs::Union{Matrix{Float64}, Nothing}=nothing, verbose::Bool=true, progress::Union{Progress, Nothing}=nothing)
    if cell_type_exprs === nothing
        if k <= 1
            error("Either k or cell_type_exprs must be specified")
        end

        cell_type_exprs = copy(hcat(prob_array.(split(genes, rand(1:k, length(genes))), max_value=maximum(genes))...)')
    end

    # TODO: filter molecules with low confidence to prevent wasting clusters on them

    cell_type_exprs = (cell_type_exprs .+ 1) ./ (sum(cell_type_exprs, dims=2) .+ 1)
    cell_type_exprs_norm = cell_type_exprs ./ sum(cell_type_exprs, dims=2)

    assignment_probs = cell_type_exprs_norm[:, genes];
    assignment_probs_prev = deepcopy(assignment_probs)
    max_diffs = Float64[]

    if verbose && progress === nothing
        progress = Progress(max_iters, 0.3)
    end

    n_iters = 0
    for i in 1:max_iters
        n_iters = i
        assignment_probs_prev .= assignment_probs
        expect_molecule_clusters!(assignment_probs, cell_type_exprs, cell_type_exprs_norm, genes, confidence, adjacent_points, adjacent_weights, new_prob=new_prob)
        if do_maximize
            maximize_molecule_clusters!(cell_type_exprs, cell_type_exprs_norm, genes, assignment_probs)
        end

        push!(max_diffs, estimate_difference_l0(assignment_probs, assignment_probs_prev))

        if verbose
            next!(progress)
        end

        if (max_diffs[end] < 0.05) && (new_prob > 1e-5)
            new_prob = 0.0
            assignment = vec(mapslices(x -> findmax(x)[2], assignment_probs, dims=1));
            cell_type_exprs = remove_unused_molecule_clusters!(assignment, cell_type_exprs, genes, confidence);
            continue
        end

        if (i >= min_iters) && (max_diffs[end] < tol) # TODO: need a better criterion. Differences are not monotone, so this threshold doesn't meen much
            if verbose
                finish!(progress)
            end
            break
        end
    end

    if verbose
        println("Algorithm stopped after $n_iters iterations. Error: $(round(max_diffs[end], sigdigits=3)). Converged: $(max_diffs[end] <= tol).")
    end

    if do_maximize
        maximize_molecule_clusters!(cell_type_exprs, cell_type_exprs_norm, genes, assignment_probs)
    end

    assignment = vec(mapslices(x -> findmax(x)[2], assignment_probs, dims=1));

    return cell_type_exprs_norm, assignment, max_diffs, assignment_probs
end

## Sampling of assignments. Converges faster for clustering.

function filter_correlated_molecule_clusters!(cell_type_exprs::Matrix{Float64}, assignment::Vector{Int}; correlation_threshold::Float64=0.95)
    cors = cor(cell_type_exprs');
    cors[diagind(cors)] .= 0
    triu!(cors);
    max_cor_ids = vec(mapslices(findmax, cors, dims=1));

    was_filtering = false
    for i1 in 1:length(max_cor_ids)
        c, i2 = max_cor_ids[i1]
        if c < correlation_threshold
            continue
        end

        was_filtering = true
        cell_type_exprs[i1, :] .= 0
        assignment[assignment .== i1] .= i2
    end

    return was_filtering
end

function remove_unused_molecule_clusters!(assignment::Vector{Int}, cell_type_exprs::Matrix{Float64}, genes::Vector{Int}, confidence::Union{Vector{Float64}, Nothing}=nothing;
        frac_of_expected::Float64=0.05)
    mol_fracs_per_type = (confidence === nothing) ? prob_array(assignment) : (sum.(split(confidence, assignment)) ./ length(assignment))
    frac_threshold = frac_of_expected / maximum(assignment)
    real_type_ids = findall(mol_fracs_per_type .>= frac_threshold)
    if length(real_type_ids) == length(mol_fracs_per_type)
        return cell_type_exprs
    end

    id_map = zeros(Int, size(cell_type_exprs, 1))
    id_map[real_type_ids] .= 1:length(real_type_ids)

    cell_type_exprs = cell_type_exprs[real_type_ids,:]

    cell_type_exprs_norm = cell_type_exprs ./ sum(cell_type_exprs, dims=2)
    for i in 1:length(assignment)
        if mol_fracs_per_type[assignment[i]] < frac_threshold
            assignment[i] = findmax(cell_type_exprs_norm[:, genes[i]])[2]
        else
            assignment[i] = id_map[assignment[i]]
        end
    end

    return cell_type_exprs
end


function maximize_molecule_clusters!(cell_type_exprs::Matrix{Float64}, genes::Vector{Int}, assignment::Vector{Int})
    cell_type_exprs[unique(assignment), :] .= 0.0

    for i in 1:length(assignment)
        cell_type_exprs[assignment[i], genes[i]] += 1.0
    end
end

function expect_molecule_clusters!(assignment::Vector{Int}, cell_type_exprs::Matrix{Float64}, genes::Vector{Int}, confidence::Vector{Float64},
        adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}; new_prob::Float64=0.001)
    denses = zeros(size(cell_type_exprs, 1))
    cell_type_exprs = (cell_type_exprs .+ 1) ./ (sum(cell_type_exprs, dims=2) .+ 1) # It doesn't sum to 1, but it simulates sparse pseudocounts
    for i in 1:length(genes)
        gene = genes[i]
        cur_weights = adjacent_weights[i]
        cur_points = adjacent_points[i]

        denses .= 0.0

        if rand() < new_prob # allowing to link new component improves global convergence
            assignment[i] = fsample(cell_type_exprs[:, gene])
        else
            for j in 1:length(cur_points)
                ca = assignment[cur_points[j]]
                denses[ca] = cell_type_exprs[ca, gene] * confidence[cur_points[j]]
            end

            for j in 1:length(cur_weights)
                denses[assignment[cur_points[j]]] *= cur_weights[j]
            end

            assignment[i] = fsample(denses)
        end
    end
end

function cluster_molecules_on_mrf_sampling(genes::Vector{Int}, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}},
        confidence::Vector{Float64}=ones(length(genes));
        k::Int=1, n_iters::Int=1000, history_depth::Int=200, new_prob::Float64=0.05, min_mols_per_type::Int=-1, correlation_threshold::Float64=0.95,
        cell_type_exprs::Union{Matrix{Float64}, Nothing}=nothing, assignment::Union{Vector{Int}, Nothing}=nothing,
        verbose::Bool=true, progress::Union{Progress, Nothing}=nothing)
    assignment_history = Vector{Int}[]
    adjacent_weights = [exp.(aw) for aw in adjacent_weights]

    if cell_type_exprs === nothing
        if k <= 1
            error("Either k or cell_type_exprs must be specified")
        end

        cell_type_exprs = zeros(k, maximum(genes));

        assignment = rand(1:k, length(genes));
        maximize_molecule_clusters!(cell_type_exprs, genes, assignment)
    else
        cell_type_exprs = deepcopy(cell_type_exprs)
        if assignment === nothing
            assignment = [fsample(cell_type_exprs[:, g]) for g in genes]
        else
            assignment = deepcopy(assignment)
        end
    end

    if min_mols_per_type < 0
        min_mols_per_type = round(Int, 0.05 * length(genes) / size(cell_type_exprs, 1))
    end

    if verbose && progress === nothing
        progress = Progress(n_iters, 0.3)
    end

    for i in 1:n_iters
        expect_molecule_clusters!(assignment, cell_type_exprs, genes, confidence, adjacent_points, adjacent_weights, new_prob=new_prob)
        maximize_molecule_clusters!(cell_type_exprs, genes, assignment)
        if (n_iters - i) == (history_depth - 1)
            cell_type_exprs = remove_unused_molecule_clusters!(assignment, cell_type_exprs, genes; min_mols_per_type=min_mols_per_type)
            new_prob = 0.0
        end

        if (n_iters - i) < history_depth
            push!(assignment_history, deepcopy(assignment))
        else
            if filter_correlated_molecule_clusters!(cell_type_exprs, assignment, correlation_threshold=correlation_threshold)
                maximize_molecule_clusters!(cell_type_exprs, genes, assignment)
            end
        end

        if verbose
            next!(progress)
        end
    end

    assignment_cons = vec(mapslices(mode, hcat(assignment_history...), dims=2));

    return cell_type_exprs ./ sum(cell_type_exprs, dims=2), assignment_cons, assignment_history
end

## Utils

function build_molecule_graph(df_spatial::DataFrame; kwargs...)
    edge_list, adjacent_dists = adjacency_list(df_spatial; kwargs...);

    real_edge_length = quantile(adjacent_dists, 0.3);
    adjacent_weights = real_edge_length ./ max.(adjacent_dists, real_edge_length);

    return convert_edge_list_to_adj_list(edge_list, adjacent_weights; n_verts=size(df_spatial, 1));
end

"""
    Adjust value based on prior. Doesn't penalize values < σ, penalize sub-linearly values in [σ; 3σ], and penalize super-linarly all >3σ
"""
@inline function adj_value_norm(x::Float64, μ::Float64, σ::Float64)::Float64
    dx = x - μ
    z = abs(dx) / σ
    adj_mult = 1 / (1 + sqrt(2/3) - 2/3)
    if z < adj_mult
        return x
    end
    return μ + sign(dx) * sqrt(z) * σ * adj_mult
end

function select_ids_uniformly(xs::Vector{Float64}, ys::Vector{Float64}, confidence::Vector{Float64}, n::Int; confidence_threshold::Float64=0.95)::Vector{Int}
    dense_ids = findall(confidence .>= confidence_threshold)
    if length(dense_ids) < n
        return dense_ids
    end

    xs = xs[dense_ids] .- minimum(xs[dense_ids])
    ys = ys[dense_ids] .- minimum(ys[dense_ids])

    n_y_bins = round(Int, sqrt(n) * maximum(ys) / maximum(xs))
    n_x_bins = div(n, n_y_bins)

    index_vals = round.(Int, xs .* (n_x_bins / maximum(xs))) .* n_y_bins .+ round.(Int, ys .* (n_y_bins / maximum(ys)));
    return dense_ids[sortperm(index_vals)[unique(round.(Int, range(1, length(dense_ids), length=n)))]];
end