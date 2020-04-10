import Clustering
using ProgressMeter
using DataFrames
using Statistics
using StatsBase
using Random

## Full EM. Probably will be needed when scRNA-seq is used as prior

function maximize_molecule_clusters!(cell_type_exprs::Matrix{Float64}, cell_type_exprs_norm::Matrix{Float64}, genes::Vector{Int},
        confidence::Vector{Float64}, assignment_probs::Matrix{Float64})
    cell_type_exprs .= 0.0;
    for i in 1:length(genes)
        t_gene = genes[i];
        t_conf = confidence[i];

        for j in 1:size(cell_type_exprs, 1)
            cell_type_exprs[j, t_gene] += t_conf * assignment_probs[j, i];
        end
    end

    # TODO: here should be adjustment based on prior from scRNA-seq using adj_value_norm

    cell_type_exprs .= (cell_type_exprs .+ 1) ./ (sum(cell_type_exprs, dims=2) .+ 1);
    cell_type_exprs_norm .= cell_type_exprs ./ sum(cell_type_exprs, dims=2);
end

"""
Params:
- adjacent_weights: must be multiplied by confidence of the corresponding adjacent_point
"""
function expect_molecule_clusters!(assignment_probs::Matrix{Float64}, cell_type_exprs::Matrix{Float64}, cell_type_exprs_norm::Matrix{Float64}, genes::Vector{Int},
        adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}; new_prob::Float64=0.05, comp_weights::Vector{Float64}=ones(size(assignment_probs, 1)))
    total_ll = 0.0
    for i in 1:length(genes)
        gene = genes[i]
        cur_weights = adjacent_weights[i]
        cur_points = adjacent_points[i]

        dense_sum = 0.0
        for ri in 1:size(assignment_probs, 1)
            c_d = 0.0
            adj_prob = 1.0 # probability that there are no neighbors from this component
            for j in 1:length(cur_points) # TODO: can try to use sparsity to optimize it. Can store BitMatrix with info about a_p > 1e-10
                a_p = assignment_probs[ri, cur_points[j]]
                c_d += cur_weights[j] * a_p
                adj_prob *= 1 - a_p
            end

            assignment_probs[ri, i] = cell_type_exprs[ri, gene] * exp(c_d) * (1 - adj_prob) * comp_weights[ri]
            dense_sum += assignment_probs[ri, i]
        end

        total_ll += log10(dense_sum)
        # TODO: for initial test, add lock array here and check whether taking it takes time
        for ri in 1:size(assignment_probs, 1)
            assignment_probs[ri, i] = (1 - new_prob) * assignment_probs[ri, i] / dense_sum + new_prob * cell_type_exprs_norm[ri, gene]
        end
    end

    return total_ll
end

function cluster_molecules_on_mrf(df_spatial::DataFrame, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}};
        n_clusters::Int, confidence_threshold::Float64=0.95, kwargs...)

    cor_mat = pairwise_gene_spatial_cor(df_spatial.gene, df_spatial.confidence, adjacent_points, adjacent_weights; confidence_threshold=confidence_threshold);
    ica_fit = fit(MultivariateStats.ICA, cor_mat, n_clusters, maxiter=100000);
    ct_exprs_init = copy((abs.(ica_fit.W) ./ sum(abs.(ica_fit.W), dims=1))')

    return cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights; cell_type_exprs=ct_exprs_init, kwargs...)
end

function cluster_molecules_on_mrf(genes::Vector{Int}, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}, confidence::Vector{Float64};
        n_clusters::Int=1, new_prob::Float64=0.05, tol::Float64=0.01, do_maximize::Bool=true, max_iters::Int=div(length(genes), 200), n_iters_without_update::Int=20,
        min_mols_per_cell::Int=0,
        cell_type_exprs::Union{Matrix{Float64}, Nothing}=nothing, assignment::Union{Vector{Int}, Nothing}=nothing,
        verbose::Bool=true, progress::Union{Progress, Nothing}=nothing, weights_pre_adjusted::Bool=false)
    if cell_type_exprs === nothing
        if n_clusters <= 1
            error("Either n_clusters or cell_type_exprs must be specified")
        end

        cell_type_exprs = copy(hcat(prob_array.(split(genes, rand(1:n_clusters, length(genes))), max_value=maximum(genes))...)')
    end

    if !weights_pre_adjusted
        adjacent_weights = [adjacent_weights[i] .* confidence[adjacent_points[i]] for i in 1:length(adjacent_weights)] # instead of multiplying each time in expect
    end

    cell_type_exprs = (cell_type_exprs .+ 1) ./ (sum(cell_type_exprs, dims=2) .+ 1)
    cell_type_exprs_norm = cell_type_exprs ./ sum(cell_type_exprs, dims=2)

    assignment_probs = zeros(size(cell_type_exprs, 1), length(genes));

    if assignment === nothing
        assignment_probs .= cell_type_exprs_norm[:, genes];
    else
        assignment_probs[CartesianIndex.(assignment, 1:length(assignment))] .= 1.0;
    end

    assignment_probs_prev = deepcopy(assignment_probs)
    max_diffs = Float64[]
    change_fracs = Float64[]

    if verbose && progress === nothing
        progress = Progress(max_iters, 0.3)
    end

    comp_weights = ones(size(assignment_probs, 1))
    n_iters = 0
    for i in 1:max_iters
        n_iters = i
        assignment_probs_prev .= assignment_probs

        expect_molecule_clusters!(assignment_probs, cell_type_exprs, cell_type_exprs_norm, genes, adjacent_points, adjacent_weights,
            new_prob=new_prob, comp_weights=comp_weights)

        if (min_mols_per_cell > 1) && (i > n_iters_without_update) && (i % 10 == 0)
            assignment_probs, n_mols_per_comp_per_clust, real_clust_ids = filter_small_molecule_clusters( # filter empty clusters
                genes, confidence, adjacent_points, assignment_probs, cell_type_exprs; min_mols_per_cell=1)

            if length(real_clust_ids) != length(comp_weights)
                assignment_probs_prev = assignment_probs_prev[real_clust_ids,:]
                cell_type_exprs = cell_type_exprs[real_clust_ids, :]
                cell_type_exprs_norm = cell_type_exprs_norm[real_clust_ids, :]
            end
            comp_weights = [max(sum(x[x .>= min_mols_per_cell]) / sum(x), 0.001) for x in n_mols_per_comp_per_clust]
        end

        if do_maximize
            maximize_molecule_clusters!(cell_type_exprs, cell_type_exprs_norm, genes, confidence, assignment_probs)
        end

        md, cf = estimate_difference_l0(assignment_probs, assignment_probs_prev, col_weights=confidence)
        push!(max_diffs, md)
        push!(change_fracs, cf)

        if verbose
            next!(progress)
        end

        if (max_diffs[end] < 2 * tol) && (new_prob > 1e-5)
            new_prob = 0.0
            assignment = vec(mapslices(x -> findmax(x)[2], assignment_probs, dims=1));
            continue
        end

        if (i > n_iters_without_update) && (maximum(max_diffs[(end - n_iters_without_update):end]) < tol)
            if verbose
                finish!(progress)
            end
            break
        end
    end

    if verbose
        @info "Algorithm stopped after $n_iters iterations. Error: $(round(max_diffs[end], sigdigits=3)). Converged: $(max_diffs[end] <= tol)."
    end

    if do_maximize
        maximize_molecule_clusters!(cell_type_exprs, cell_type_exprs_norm, genes, confidence, assignment_probs)
    end

    assignment = vec(mapslices(x -> findmax(x)[2], assignment_probs, dims=1));

    return (exprs=cell_type_exprs_norm, assignment=assignment, diffs=max_diffs, assignment_probs=assignment_probs, change_fracs=change_fracs)
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

function filter_small_molecule_clusters(genes::Vector{Int}, confidence::Vector{Float64}, adjacent_points::Vector{Vector{Int}},
        assignment_probs::Matrix{Float64}, cell_type_exprs::Matrix{Float64}; min_mols_per_cell::Int, confidence_threshold::Float64=0.95)

    assignment = vec(mapslices(x -> findmax(x)[2], assignment_probs, dims=1));
    conn_comps_per_clust = get_connected_component_per_label(assignment, adjacent_points, 1;
        confidence=confidence, confidence_threshold=confidence_threshold)[1];
    n_mols_per_comp_per_clust = [length.(c) for c in conn_comps_per_clust];

    real_clust_ids = findall(maximum.(n_mols_per_comp_per_clust) .>= min_mols_per_cell)
    # real_clust_ids = Int[]
    # if method == :size
    #     real_clust_ids = findall(maximum.(n_mols_per_comp_per_clust) .>= min_mols_per_cell)
    # elseif method == :frac
    #     real_clust_ids = findall([sum(x[x .>= 30]) / sum(x) for x in n_mols_per_comp_per_clust] .> 0.5) # TODO: remove this method
    # else
    #     error("Unknown method: $method")
    # end

    if length(real_clust_ids) == size(assignment_probs, 1)
        return assignment_probs, n_mols_per_comp_per_clust, real_clust_ids
    end

    assignment_probs = assignment_probs[real_clust_ids,:];
    cell_type_exprs = cell_type_exprs[real_clust_ids,:];

    for i in findall(vec(sum(assignment_probs, dims=1) .< 1e-10))
        assignment_probs[:, i] .= cell_type_exprs[:, genes[i]]
    end

    assignment_probs ./= sum(assignment_probs, dims=1);

    return assignment_probs, n_mols_per_comp_per_clust[real_clust_ids], real_clust_ids
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

function cluster_molecules_on_mrf_sampling(genes::Vector{Int}, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}, confidence::Vector{Float64};
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

function pairwise_gene_spatial_cor(genes::Vector{Int}, confidence::Vector{Float64}, adjacent_points::Array{Vector{Int}, 1}, adjacent_weights::Array{Vector{Float64}, 1};
        confidence_threshold::Float64=0.95)::Matrix{Float64}
    gene_cors = zeros(maximum(genes), maximum(genes))
    sum_weight_per_gene = zeros(maximum(genes))
    for gi in 1:length(genes)
        cur_adj_points = adjacent_points[gi]
        cur_adj_weights = adjacent_weights[gi]
        g2 = genes[gi]
        if confidence[gi] < confidence_threshold
            continue
        end

        for ai in 1:length(cur_adj_points)
            if confidence[cur_adj_points[ai]] < confidence_threshold
                continue
            end

            g1 = genes[cur_adj_points[ai]]
            cw = cur_adj_weights[ai]
            gene_cors[g2, g1] += cw
            sum_weight_per_gene[g1] += cw
            sum_weight_per_gene[g2] += cw
        end
    end

    for ci in 1:length(sum_weight_per_gene)
        for ri in 1:length(sum_weight_per_gene)
            gene_cors[ri, ci] /= fmax(sqrt(sum_weight_per_gene[ri] * sum_weight_per_gene[ci]), 0.1)
        end
    end

    return gene_cors
end