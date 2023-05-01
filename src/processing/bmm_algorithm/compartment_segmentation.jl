function expect_molecule_compartments!(assignment_probs::Matrix{Float64}, adj_list::AdjList, is_locked::BitVector)
    for i in eachindex(adj_list.ids)
        if is_locked[i]
            continue
        end
        cur_points = adj_list.ids[i]
        cur_weights = adj_list.weights[i]

        if length(cur_points) == 0
            continue
        end

        dense_sum = 0.0
        for ri in 1:size(assignment_probs, 1)
            c_d = 0.0
            for j in eachindex(cur_points)
                c_d += cur_weights[j] * assignment_probs[ri, cur_points[j]]
            end

            dense_sum += assignment_probs[ri, i] = exp(c_d)
        end

        assignment_probs[:, i] ./= dense_sum
    end
end

init_nuclei_cyto_compartments(args...; nuclei_genes::T, cyto_genes::T, kwargs...) where T<:Union{Vector{String}, Dict{String, Float64}} =
    init_molecule_compartments(args..., [nuclei_genes, cyto_genes]; kwargs...)

init_molecule_compartments(pos_data::Matrix{Float64}, genes::Vector{Int}, comp_genes::Vector{Vector{String}}; kwargs...) =
    init_molecule_compartments(pos_data, genes, [Dict(g => 1.0 for g in gs) for gs in comp_genes]; kwargs...)

function init_molecule_compartments(pos_data::Matrix{Float64}, genes::Vector{Int}, comp_genes::Vector{Dict{String, Float64}};
        gene_names::Vector{String}, scale::Float64, required_comps::Union{Vector{Int}, Nothing}=nothing)

    id_per_gene = Dict(g => i for (i,g) in enumerate(gene_names))
    comp_genes = [Dict(id_per_gene[g] => v for (g,v) in gsd) for gsd in comp_genes]

    default_probs = vcat(ones(length(comp_genes)) ./ length(comp_genes), 0.0)
    assignment_probs = zeros(length(comp_genes) + 1, length(genes))
    assignment_probs .= default_probs

    is_assigned = [falses(size(pos_data, 2)) for c in comp_genes]
    for (i,gs) in enumerate(comp_genes)
        for (g,p) in gs
            mask = (genes .== g)
            gene_vec = setindex!(zeros(length(comp_genes) + 1), p, i)
            assignment_probs[:, mask] .= gene_vec .+ (1 - p) .* default_probs
            if (required_comps === nothing) || (i in required_comps)
                is_assigned[i] .|= mask
            end
        end
    end

    in_compartment_probs = Float64[]
    is_locked = any(hcat(is_assigned...), dims=2)[:]
    dists = zeros(size(pos_data, 2))

    for ia in is_assigned
        dists = max.(dists, maximum.(knn(KDTree(pos_data[:,ia]), pos_data, 2)[2]))
    end
    in_compartment_probs = [1.0 - max(min((d / scale - 1.) / 2., 1.), 0.) for d in dists]

    is_locked .&= (in_compartment_probs .> 0.99)
    is_locked .|= (in_compartment_probs .< 0.01)

    assignment_probs[end, :] .= 1 .- in_compartment_probs;
    assignment_probs[1:(end-1), :] .*= in_compartment_probs';

    return assignment_probs, is_locked
end

function segment_molecule_compartments(
        assignment_probs::Matrix{Float64}, is_locked::BitVector, adj_list::AdjList, confidence::Vector{Float64};
        n_iters_without_update=20, weight_mult::Float64=1.0, tol::Float64=0.01, max_iter::Int=500, verbose::Bool=true
    )

    adj_weights = [weight_mult .* adj_list.weights[i] .* confidence[adj_list.ids[i]] for i in eachindex(adj_list.ids)] # instead of multiplying each time in expect
    adj_list = AdjList(adj_list.ids, adj_weights)

    max_diffs, change_fracs = Float64[], Float64[]
    assignment_probs_prev = deepcopy(assignment_probs)
    progress = verbose ? Progress(max_iter) : nothing

    # Run
    for i in 1:max_iter
        assignment_probs_prev .= assignment_probs
        # Iteration
        expect_molecule_compartments!(assignment_probs, adj_list, is_locked)

        # Stop criterion
        md, cf = estimate_difference_l0(assignment_probs, assignment_probs_prev, col_weights=confidence)
        push!(max_diffs, md)
        push!(change_fracs, cf)

        prog_vals = [("Iteration", i), ("Max. difference", md), ("Fraction of probs changed", cf)]
        if progress !== nothing
            next!(progress, showvalues=prog_vals)
        end
        if (i > n_iters_without_update) && (maximum(max_diffs[(end - n_iters_without_update):end]) < tol)
            if verbose
                finish!(progress, showvalues=prog_vals)
            end
            break
        end
    end

    # Wrap results
    assignment = vec(mapslices(x -> findmax(x)[2], assignment_probs, dims=1));
    return (assignment=assignment, diffs=max_diffs, assignment_probs=assignment_probs, change_fracs=change_fracs)
end

adjust_mrf_with_compartments(adj_list::AdjList, args...; kwargs...) =
    adjust_mrf_with_compartments!(deepcopy(adj_list), args...; kwargs...)

function adjust_mrf_with_compartments!(
        adj_list::AdjList, nuclei_probs::Vector{Float64}, cyto_probs::Vector{Float64};
        min_nuc_prob::Float64=0.25, min_weight::Float64=0.01
    )
    for (m1,ids) in enumerate(adj_list.ids)
        nuc_prob = nuclei_probs[m1]
        (min_nuc_prob > 0.25) || continue
        weights = adj_list.weights[m1]
        for (i2,m2) in enumerate(ids)
            weights[i2] *= 1 - cyto_probs[m2] * nuc_prob
        end

        mask = (weights .> min_weight)
        any(mask) || continue

        adj_list.ids[m1] = adj_list.ids[m1][mask]
        adj_list.weights[m1] = adj_list.weights[m1][mask]
    end
    return adj_list
end

## Wrappers

function estimate_molecule_compartments(df_spatial::DataFrame, gene_names::Vector{String}; nuclei_genes::String, cyto_genes::String, scale::Float64)
    @info "Estimating compartment regions..."
    comp_genes = Dict{String, Vector{String}}()

    # Parse genes from CLI
    nuclei_genes, cyto_genes = map(("nuclei" => nuclei_genes, "cyto" => cyto_genes)) do (k,g)
        c_genes = split_string_list(g)

        missing_genes = c_genes[.!in.(c_genes, Ref(Set(gene_names)))]
        (length(missing_genes) == 0) || @warn "Genes $(join(missing_genes, ',')) are missing from the data"
        c_genes = intersect(c_genes, gene_names)
        length(c_genes) > 0 || error("No genes left in $k after filtration")
        c_genes
    end

    # Run segmentation
    adj_list = build_molecule_graph(df_spatial, filter=false);

    init_probs, is_locked = init_nuclei_cyto_compartments(
        position_data(df_spatial), df_spatial.gene; gene_names=gene_names, scale=scale,
        nuclei_genes=nuclei_genes, cyto_genes=cyto_genes
    );

    comp_segs = segment_molecule_compartments(init_probs, is_locked, adj_list, df_spatial.confidence);
    # TODO: comp_genes is always empty now

    @info "Done"

    id_per_gene = Dict(g => i for (i,g) in enumerate(gene_names))
    comp_genes = [[id_per_gene[g] for g in gs] for gs in values(comp_genes)]

    compartment = ["Nuclei", "Cyto", "Unknown"][comp_segs.assignment]

    return comp_segs, comp_genes, compartment
end
