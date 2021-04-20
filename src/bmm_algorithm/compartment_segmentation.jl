function expect_molecule_compartments!(assignment_probs::Matrix{Float64}, adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}, is_locked::BitVector)
    for i in 1:length(adjacent_points)
        if is_locked[i]
            continue
        end
        cur_weights = adjacent_weights[i]
        cur_points = adjacent_points[i]

        dense_sum = 0.0
        for ri in 1:size(assignment_probs, 1)
            c_d = 0.0
            for j in 1:length(cur_points)
                c_d += cur_weights[j] * assignment_probs[ri, cur_points[j]]
            end

            dense_sum += assignment_probs[ri, i] = exp(c_d)
        end

        assignment_probs[:, i] ./= dense_sum
    end
end

init_molecule_compartments(pos_data::Matrix{Float64}, genes::Vector{Int}, comp_genes::Vector{Vector{String}}; kwargs...) =
    init_molecule_compartments(pos_data, genes, [Dict(g => 1.0 for g in gs) for gs in comp_genes]; kwargs...)

function init_molecule_compartments(pos_data::Matrix{Float64}, genes::Vector{Int}, comp_genes::Vector{Dict{String, Float64}}; 
        gene_names::Vector{String}, nn_id::Int)
    id_per_gene = Dict(g => i for (i,g) in enumerate(gene_names))
    comp_genes = [Dict(id_per_gene[g] => v for (g,v) in gsd) for gsd in comp_genes]

    assignment_probs = zeros(length(comp_genes) + 1, length(genes))
    for (i,gs) in enumerate(comp_genes)
        for (g,p) in gs
            assignment_probs[i, genes .== g] .= p
        end
    end
    
    nn_ids = knn(KDTree(pos_data), pos_data, nn_id + 1)[1];
    is_assigned = sum(assignment_probs, dims=1)[:] .> 1e-5;
    assignment_probs[end, [!any(is_assigned[ids]) for ids in nn_ids]] .= 1.0;
    
    is_locked = sum(assignment_probs, dims=1)[:] .> 1e-5;
    assignment_probs[:, .!is_locked] .= 1. / size(assignment_probs, 1);
    return assignment_probs, is_locked
end

function segment_molecule_compartments(pos_data::Matrix{Float64}, genes::Vector{Int}, 
        adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}, confidence::Vector{Float64}; 
        comp_genes::Vector{<:Union{Vector{String}, Dict{String, Float64}}}, gene_names::Vector{String}, nn_id::Int, n_iters_without_update=20,
        weights_pre_adjusted::Bool=false, weight_mult::Float64=1.0, tol::Float64=0.01, max_iter::Int=500, verbose::Bool=true)

    if !weights_pre_adjusted
        adjacent_weights = [weight_mult .* adjacent_weights[i] .* confidence[adjacent_points[i]] for i in 1:length(adjacent_weights)] # instead of multiplying each time in expect
    end

    # Init
    assignment_probs, is_locked = init_molecule_compartments(pos_data, genes, comp_genes; gene_names=gene_names, nn_id=nn_id);

    max_diffs, change_fracs = Float64[], Float64[]
    assignment_probs_prev = deepcopy(assignment_probs)
    progress = verbose ? Progress(max_iter) : nothing
    for i in 1:max_iter
        assignment_probs_prev .= assignment_probs
        # Iteration
        expect_molecule_compartments!(assignment_probs, adjacent_points, adjacent_weights, is_locked)
        
        # Stop criterion
        md, cf = estimate_difference_l0(assignment_probs, assignment_probs_prev, col_weights=confidence)
        push!(max_diffs, md)
        push!(change_fracs, cf)

        if progress !== nothing
            next!(progress, showvalues = [(:iteration, i), (:max_diff, md), (:change_frac, cf)])
        end
        if (i > n_iters_without_update) && (maximum(max_diffs[(end - n_iters_without_update):end]) < tol)
            if verbose
                finish!(progress)
            end
            break
        end
    end

    # Wrap results
    assignment = vec(mapslices(x -> findmax(x)[2], assignment_probs, dims=1));
    return (assignment=assignment, diffs=max_diffs, assignment_probs=assignment_probs, change_fracs=change_fracs)
end

adjust_mrf_with_compartments(adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}, args...; kwargs...) =
    adjust_mrf_with_compartments!(deepcopy(adjacent_points), deepcopy(adjacent_weights), args...; kwargs...)

function adjust_mrf_with_compartments!(adjacent_points::Vector{Vector{Int}}, adjacent_weights::Vector{Vector{Float64}}, nuclei_probs::Vector{Float64}, cyto_probs::Vector{Float64}; min_nuc_prob::Float64=0.25, min_weight::Float64=0.01)
    for (m1,ids) in enumerate(adjacent_points)
        nuc_prob = nuclei_probs[m1]
        (min_nuc_prob > 0.25) || continue
        weights = adjacent_weights[m1]
        for (i2,m2) in enumerate(ids)
            weights[i2] *= 1 - cyto_probs[m2] * nuc_prob
        end

        mask = (weights .> min_weight)
        any(mask) || continue

        adjacent_points[m1] = adjacent_points[m1][mask]
        adjacent_weights[m1] = adjacent_weights[m1][mask]
    end
    return adjacent_points, adjacent_weights
end