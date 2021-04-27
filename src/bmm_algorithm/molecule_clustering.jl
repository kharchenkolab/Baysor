using ProgressMeter
using DataFrames
using Statistics
using StatsBase
using Random

## Full EM. Probably will be needed when scRNA-seq is used as prior

function maximize_molecule_clusters!(cell_type_exprs::Matrix{Float64}, genes::Vector{Int}, confidence::Vector{Float64}, assignment_probs::Matrix{Float64};
        prior_exprs::Union{Matrix{Float64}, Nothing}=nothing, prior_stds::Union{Matrix{Float64}, Nothing}=nothing, add_pseudocount::Bool=false)
    cell_type_exprs .= 0.0;
    for i in 1:length(genes)
        t_gene = genes[i];
        t_conf = confidence[i];

        for j in 1:size(cell_type_exprs, 1)
            cell_type_exprs[j, t_gene] += t_conf * assignment_probs[j, i];
        end
    end

    if prior_exprs !== nothing
        mult = sum(cell_type_exprs, dims=2)
        cell_type_exprs .= adj_value_norm.(cell_type_exprs, prior_exprs .* mult, prior_stds .* mult)
    end

    if add_pseudocount
        cell_type_exprs .= (cell_type_exprs .+ 1) ./ (sum(cell_type_exprs, dims=2) .+ 1);
    else
        cell_type_exprs ./= sum(cell_type_exprs, dims=2);
    end
end

"""
Params:
- adjacent_weights: must be multiplied by confidence of the corresponding adjacent_point
"""
function expect_molecule_clusters!(assignment_probs::Matrix{Float64}, cell_type_exprs::Matrix{Float64}, genes::Vector{Int},
        adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}; new_prob::Float64=0.05)
    total_ll = 0.0
    for i in 1:length(genes)
        gene = genes[i]
        cur_weights = adjacent_weights[i]
        cur_points = adjacent_points[i]

        dense_sum = 0.0
        ct_dense_sum = 0.0 # Because of pseudocounts, cell_type_exprs aren't normalized
        for ri in 1:size(assignment_probs, 1)
            c_d = 0.0
            for j in 1:length(cur_points) # TODO: can try to use sparsity to optimize it. Can store BitMatrix with info about a_p > 1e-10
                a_p = assignment_probs[ri, cur_points[j]]
                if a_p < 1e-5
                    continue
                end

                c_d += cur_weights[j] * a_p
            end

            ctp = cell_type_exprs[ri, gene]
            assignment_probs[ri, i] = ctp * exp(c_d)
            dense_sum += assignment_probs[ri, i]
            ct_dense_sum += ctp
        end

        total_ll += log10(dense_sum)
        # TODO: for initial test, add lock array here and check whether taking it takes time
        for ri in 1:size(assignment_probs, 1)
            assignment_probs[ri, i] = (1 - new_prob) * assignment_probs[ri, i] / dense_sum + new_prob * cell_type_exprs[ri, gene] / ct_dense_sum
        end
    end

    return total_ll
end


function cluster_molecules_on_mrf(df_spatial::DataFrame, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}};
        n_clusters::Int, confidence_threshold::Float64=0.95, kwargs...)

    cor_mat = pairwise_gene_spatial_cor(df_spatial.gene, df_spatial.confidence, adjacent_points, adjacent_weights; confidence_threshold=confidence_threshold);
    ica_fit = fit(MultivariateStats.ICA, cor_mat, n_clusters, maxiter=10000);
    ct_exprs_init = copy((abs.(ica_fit.W) ./ sum(abs.(ica_fit.W), dims=1))')

    return cluster_molecules_on_mrf(df_spatial.gene, adjacent_points, adjacent_weights, df_spatial.confidence; cell_type_exprs=ct_exprs_init, kwargs...)
end


function cluster_molecules_on_mrf(genes::Vector{Int}, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}, confidence::Vector{Float64};
        n_clusters::Int=1, new_prob::Float64=0.05, tol::Float64=0.01, do_maximize::Bool=true, max_iters::Int=max(10000, div(length(genes), 200)), n_iters_without_update::Int=20,
        cell_type_exprs::Union{<:AbstractMatrix{Float64}, Nothing}=nothing, assignment::Union{Vector{Int}, Nothing}=nothing, assignment_probs::Union{Matrix{Float64}, Nothing}=nothing,
        verbose::Bool=true, progress::Union{Progress, Nothing}=nothing, weights_pre_adjusted::Bool=false, weight_mult::Float64=1.0, init_mod::Int=10000, kwargs...)
    # Initialization
    if cell_type_exprs === nothing
        if n_clusters <= 1
            if assignment !== nothing
                n_clusters = maximum(assignment)
            else
                error("Either n_clusters, assignment or cell_type_exprs must be specified")
            end
        end

        if init_mod < 0
            cell_type_exprs = copy(hcat(prob_array.(split(genes, rand(1:n_clusters, length(genes))), max_value=maximum(genes))...)')
        else
            gene_probs = prob_array(genes)
            cell_type_exprs = gene_probs' .* (0.95 .+ (hcat([[hash(x1 * x2^2) for x2 in 1:n_clusters] for x1 in 1:length(gene_probs)]...) .% init_mod) ./ 100000) # determenistic way of adding pseudo-random noise
        end
        # cell_type_exprs ./= sum(cell_type_exprs, dims=2)
        cell_type_exprs = (cell_type_exprs .+ 1) ./ (sum(cell_type_exprs, dims=2) .+ 1)
    else
        cell_type_exprs = deepcopy(Matrix(cell_type_exprs))
    end
    ## On the first iteration we don't have a pseudocount, but it shouldn't be a problem

    if !weights_pre_adjusted
        adjacent_weights = [weight_mult .* adjacent_weights[i] .* confidence[adjacent_points[i]] for i in 1:length(adjacent_weights)] # instead of multiplying each time in expect
    end

    if assignment_probs === nothing
        assignment_probs = zeros(size(cell_type_exprs, 1), length(genes));

        if assignment === nothing
            assignment_probs .= cell_type_exprs[:, genes];
            assignment_probs ./= sum(assignment_probs, dims=1)
        else
            assignment_probs[CartesianIndex.(assignment, 1:length(assignment))] .= 1.0;
        end
    else
        assignment_probs = deepcopy(assignment_probs)
    end

    assignment_probs_prev = deepcopy(assignment_probs)
    max_diffs = Float64[]
    change_fracs = Float64[]

    if verbose && progress === nothing
        progress = Progress(max_iters, 0.3)
    end

    # EM iterations
    penalize_small_clusters = false
    n_iters = 0
    for i in 1:max_iters
        n_iters = i
        assignment_probs_prev .= assignment_probs

        expect_molecule_clusters!(assignment_probs, cell_type_exprs, genes, adjacent_points, adjacent_weights, new_prob=new_prob)

        if do_maximize
            maximize_molecule_clusters!(cell_type_exprs, genes, confidence, assignment_probs; add_pseudocount=true, kwargs...)
        end

        md, cf = estimate_difference_l0(assignment_probs, assignment_probs_prev, col_weights=confidence)
        push!(max_diffs, md)
        push!(change_fracs, cf)

        prog_vals = [("Iteration", i), ("Max. difference", md), ("Fraction of assignment changed", cf)]
        if verbose
            next!(progress, showvalues=prog_vals)
        end

        if (max_diffs[end] < tol) & ((new_prob > 1e-10) | !penalize_small_clusters)
            new_prob = 0.0
            penalize_small_clusters = true
            continue
        end

        if (i > n_iters_without_update) && (maximum(max_diffs[(end - n_iters_without_update):end]) < tol)
            if verbose
                finish!(progress, showvalues=prog_vals)
            end
            break
        end
    end

    if verbose
        @info "Algorithm stopped after $n_iters iterations. Error: $(round(max_diffs[end], sigdigits=3)). Converged: $(max_diffs[end] <= tol)."
    end

    if do_maximize
        maximize_molecule_clusters!(cell_type_exprs, genes, confidence, assignment_probs, add_pseudocount=false)
    end

    assignment = vec(mapslices(x -> findmax(x)[2], assignment_probs, dims=1));

    return (exprs=cell_type_exprs, assignment=assignment, diffs=max_diffs, assignment_probs=assignment_probs, change_fracs=change_fracs)
end

function filter_small_molecule_clusters(genes::Vector{Int}, confidence::Vector{Float64}, adjacent_points::Vector{Vector{Int}},
        assignment_probs::Matrix{Float64}, cell_type_exprs::Matrix{Float64}; min_mols_per_cell::Int, confidence_threshold::Float64=0.95)

    assignment = vec(mapslices(x -> findmax(x)[2], assignment_probs, dims=1));
    conn_comps_per_clust = get_connected_components_per_label(assignment, adjacent_points, 1;
        confidence=confidence, confidence_threshold=confidence_threshold)[1];
    n_mols_per_comp_per_clust = [length.(c) for c in conn_comps_per_clust];

    real_clust_ids = findall(maximum.(n_mols_per_comp_per_clust) .>= min_mols_per_cell)

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

## Utils

"""
    Adjust value based on prior. Doesn't penalize values < σ, penalize linearly values in [σ; 3σ], and super-linarly all >= 3σ
"""
@inline function adj_value_norm(x::Float64, μ::Float64, σ::Float64)::Float64
    dx = x - μ
    z = abs(dx) / σ
    if z < 1
        return x
    end

    if z < 3
        return μ + sign(dx) * (1 + (z - 1) / 4) * σ
    end

    return μ + sign(dx) * (sqrt(z) + 1.5 - sqrt(3)) * σ
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