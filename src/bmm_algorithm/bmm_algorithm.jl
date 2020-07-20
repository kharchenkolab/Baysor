using DataFrames
using Distributions
using NearestNeighbors
using StatsBase

using Base.Threads
using Dates: now, DateTime

import DistributedArrays
import LightGraphs

function drop_unused_components!(data::BmmData; min_n_samples::Int=2)
    non_noise_ids = findall(data.assignment .> 0)
    n_mols_per_cell = num_of_molecules_per_cell(data) # Can't use comp.n_samples because it's updated only on the maximization step
    existed_ids = findall(n_mols_per_cell .>= min_n_samples)
    id_map = zeros(Int, length(data.components))
    id_map[existed_ids] .= 1:length(existed_ids)

    data.assignment[non_noise_ids] .= id_map[data.assignment[non_noise_ids]]
    data.components = data.components[existed_ids]
end

adjacent_component_ids(assignment::Array{Int, 1}, adjacent_points::Array{Int, 1})::Array{Int, 1} =
    [x for x in Set(assignment[adjacent_points]) if x > 0]

"""
...
# Arguments
- `assignment::Array{Int, 1}`:
- `adjacent_points::Array{Int, 1}`:
- `adjacent_weights::Array{Float64, 1}`: weights here mean `1 / distance` between two adjacent points
"""
@inline function adjacent_component_weights!(comp_weights::Vector{Float64}, comp_ids::Vector{Int}, component_weights::Dict{Int, Float64},
        assignment::Vector{Int}, adjacent_points::Vector{Int}, adjacent_weights::Vector{Float64}, confidences::Vector{Float64})
    empty!(component_weights)
    empty!(comp_weights)
    empty!(comp_ids)
    zero_comp_weight = 0.0

    @inbounds @simd for i in 1:length(adjacent_weights)
        c_point = adjacent_points[i]
        c_id = assignment[c_point]
        c_conf = confidences[c_point]
        cw = adjacent_weights[i]
        if c_id == 0
            zero_comp_weight += cw * (1 - c_conf)
        else
            component_weights[c_id] = get(component_weights, c_id, 0.0) + cw * c_conf
        end
    end

    for (k, v) in component_weights
        push!(comp_ids, k)
        push!(comp_weights, v)
    end

    return zero_comp_weight
end

function adjust_densities_by_prior_segmentation!(denses::Vector{Float64}, segment_id::Int, largest_cell_id::Int,
        data::BmmData, adj_classes::Vector{Int}, adj_weights::Vector{Float64})
    seg_prior_pow = data.prior_seg_confidence * exp(3*data.prior_seg_confidence)
    seg_size = data.n_molecules_per_segment[segment_id]
    largest_cell_size = min(get(data.components[largest_cell_id].n_molecules_per_segment, segment_id, 0) + 1, seg_size)

    for j in eachindex(adj_weights)
        c_adj = adj_classes[j]
        if c_adj == largest_cell_id
            continue
        end

        main_seg = data.main_segment_per_cell[c_adj]
        n_cell_mols_per_seg = min(get(data.components[c_adj].n_molecules_per_segment, segment_id, 0) + 1, seg_size)

        if (segment_id == main_seg) || (main_seg == 0) # Part against over-segmentation
            denses[j] *= (1. - data.prior_seg_confidence)^0.5 * (n_cell_mols_per_seg / largest_cell_size)^data.prior_seg_confidence
        else # Part against under-segmentation and overlap
            denses[j] *= (1. - data.prior_seg_confidence)^0.5 * (1. - n_cell_mols_per_seg / seg_size)^seg_prior_pow
        end
    end
end

function expect_dirichlet_spatial!(data::BmmData; stochastic::Bool=true, noise_density_threshold::Float64=1e-100)
    component_weights = Dict{Int, Float64}();
    adj_classes = Int[]
    adj_weights = Float64[]
    denses = Float64[]

    has_seg_prior = !isempty(data.segment_per_molecule)
    adj_classes_global = get_global_adjacent_classes(data)

    for i in 1:size(data.x, 1)
        x::Float64 = position_data(data)[1,i]
        y::Float64 = position_data(data)[2,i]
        gene::Int = composition_data(data)[i]
        confidence::Float64 = data.confidence[i]
        mol_cluster::Int = isempty(data.cluster_per_molecule) ? 0 : data.cluster_per_molecule[i]
        segment_id = has_seg_prior ? data.segment_per_molecule[i] : 0

        # Looks like it's impossible to optimize further, even with vectorization. It means that creating vectorized version of expect_dirichlet_spatial makes few sense
        zero_comp_weight = adjacent_component_weights!(adj_weights, adj_classes, component_weights, data.assignment,
            data.adjacent_points[i], data.adjacent_weights[i], data.confidence)

        if i in keys(adj_classes_global)
            n1 = length(adj_classes)
            append!(adj_classes, adj_classes_global[i])
            append!(adj_weights, ones(length(adj_classes) - n1) .* data.real_edge_weight)
        end

        empty!(denses)
        largest_cell_id = 0
        largest_cell_size = 0
        for j in eachindex(adj_weights)
            c_adj = adj_classes[j]
            cc = data.components[c_adj]
            c_dens = confidence * adj_weights[j] * cc.prior_probability * pdf(cc, x, y, gene, use_smoothing=data.use_gene_smoothing)
            if (c_adj > 0) && (c_adj < length(data.cluster_per_cell)) && (data.cluster_per_cell[c_adj] != mol_cluster)
                c_dens *= data.cluster_penalty_mult
            end

            if (segment_id > 0) && (data.main_segment_per_cell[c_adj] in (segment_id, 0))
                # Find an adjacent cell with the largest number of molecules per cell after assigning the current point to it
                cur_mols_per_seg = min(get(cc.n_molecules_per_segment, segment_id, 0) + 1, data.n_molecules_per_segment[segment_id])
                if (cur_mols_per_seg > largest_cell_size) || ((cur_mols_per_seg == largest_cell_size) && (cc.n_samples > data.components[largest_cell_id].n_samples))
                    largest_cell_size = cur_mols_per_seg
                    largest_cell_id = c_adj
                end
            end

            push!(denses, c_dens)
        end

        if (segment_id > 0) & (largest_cell_id > 0)
            adjust_densities_by_prior_segmentation!(denses, segment_id, largest_cell_id, data, adj_classes, adj_weights)
        end

        if sum(denses) < noise_density_threshold
            assign!(data, i, 0) # Noise class
            continue
        end

        if confidence < 1.0
            push!(denses, (1 - confidence) * fmax(data.real_edge_weight, zero_comp_weight) * data.noise_density)
            push!(adj_classes, 0)
        end

        if !stochastic
            assign!(data, i, adj_classes[findmax(denses)[2]])
            continue
        end

        assign!(data, i, fsample(adj_classes, denses))
    end
end

function update_prior_probabilities!(components::Array{Component, 1}, new_component_weight::Float64)
    c_weights = [max(c.n_samples, new_component_weight) for c in components]
    prior_probs = rand(Distributions.Dirichlet(c_weights))
    # t_mmc = 150.0;
    # t_dist = Normal(t_mmc, t_mmc * 3)
    # prior_probs = pdf.(t_dist, c_weights)
    # prior_probs[c_weights .< t_mmc] .= c_weights[c_weights .< t_mmc] ./ t_mmc .* pdf(t_dist, t_mmc)

    for (c, p) in zip(components, prior_probs)
        c.prior_probability = p
    end
end

function maximize!(data::BmmData)
    ids_by_assignment = split_ids(data.assignment .+ 1, max_factor=length(data.components)+1)[2:end]

    @inbounds @views for i in 1:length(data.components)
        p_ids = ids_by_assignment[i]
        maximize!(data.components[i], position_data(data)[:, p_ids], composition_data(data)[p_ids], confidence(data)[p_ids])
    end

    data.noise_density = estimate_noise_density_level(data)
    if isinf(data.noise_density)
        error("Infinity noise density")
    end

    if length(data.cluster_per_molecule) > 0
        data.cluster_per_cell = [isempty(ids) ? 0 : mode(data.cluster_per_molecule[ids]) for ids in ids_by_assignment]
    end
end

function noise_composition_density(data::BmmData)::Float64
#     mean([mean(c.composition_params.counts[c.composition_params.counts .> 0] ./ c.composition_params.sum_counts) for c in data.components]);
    acc = 0.0 # Equivalent to the above commented expression if sum_counts == sum(counts)
    n_comps = 0.0
    for (ci, c) in enumerate(data.components)
        if c.composition_params.sum_counts < 1e-3
            continue
        end

        inn_acc = 0
        for v in c.composition_params.counts
            inn_acc += (v > 0)
        end

        acc += 1.0 / inn_acc
        n_comps += 1.0
    end

    return acc / fmax(n_comps, 1.0)
end

function estimate_noise_density_level(data::BmmData)::Float64
    composition_density = noise_composition_density(data)

    std_vals = data.distribution_sampler.shape_prior.std_values;
    position_density = pdf(MultivariateNormal([0.0, 0.0], diagm(0 => std_vals.^2)), 3 .* std_vals)

    return position_density * composition_density
end

append_empty_component!(data::BmmData) =
    push!(data.components, sample_distribution!(data))[end]

function append_empty_components!(data::BmmData, new_component_frac::Float64)
    for i in 1:round(Int, new_component_frac * length(data.components))
        append_empty_component!(data)
    end
end

function get_global_adjacent_classes(data::BmmData)::Dict{Int, Array{Int, 1}}
    adj_classes_global = Dict{Int, Array{Int, 1}}()
    for (cur_id, comp) in enumerate(data.components)
        if comp.n_samples > 0
            continue
        end

        nearest_id = knn(data.position_knn_tree, comp.position_params.μ, 1)[1][1]
        for t_id in data.adjacent_points[nearest_id]
            if t_id in keys(adj_classes_global)
                push!(adj_classes_global[t_id], cur_id)
            else
                adj_classes_global[t_id] = [cur_id]
            end
        end

        if nearest_id in keys(adj_classes_global)
            push!(adj_classes_global[nearest_id], cur_id)
        else
            adj_classes_global[nearest_id] = [cur_id]
        end
    end

    return adj_classes_global
end

function build_cell_graph(assignment::Vector{Int}, adjacent_points::Array{Vector{Int}, 1}, mol_ids::Vector{Int}, cell_id::Int; confidence::Union{Vector{Float64}, Nothing}=nothing, confidence_threshold::Float64=0.5)
    cur_adj_point = Vector{Int}[]
    if confidence === nothing
        cur_adj_point = filter.(p -> assignment[p] == cell_id, adjacent_points[mol_ids])
    else
        mol_ids = mol_ids[confidence[mol_ids] .> confidence_threshold]
        cur_adj_point = filter.(p -> (assignment[p] == cell_id) & (confidence[p] .> confidence_threshold), view(adjacent_points, mol_ids))
    end

    vert_id_per_mol_id = Dict(mi => vi for (vi, mi) in enumerate(mol_ids))

    for points in cur_adj_point
        for j in 1:length(points)
            points[j] = vert_id_per_mol_id[points[j]]
        end
    end

    neg = LightGraphs.SimpleGraphs.cleanupedges!(cur_adj_point)
    return LightGraphs.SimpleGraph(neg, cur_adj_point), mol_ids
end

function get_connected_components_per_label(assignment::Vector{Int}, adjacent_points::Array{Vector{Int}, 1}, min_molecules_per_cell::Int; kwargs...)
    mol_ids_per_cell = split(1:length(assignment), assignment .+ 1)[2:end]
    real_cell_ids = findall(length.(mol_ids_per_cell) .>= min_molecules_per_cell)
    graph_per_cell = [build_cell_graph(assignment, adjacent_points, mol_ids_per_cell[ci], ci; kwargs...)[1] for ci in real_cell_ids];
    return LightGraphs.connected_components.(graph_per_cell), real_cell_ids, mol_ids_per_cell
end

function split_cells_by_connected_components!(data::BmmData; add_new_components::Bool, min_molecules_per_cell::Int)
    conn_comps_per_cell, real_cell_ids, mol_ids_per_cell = get_connected_components_per_label(data.assignment, data.adjacent_points,
        min_molecules_per_cell; confidence=data.confidence)

    for (cell_id, conn_comps) in zip(real_cell_ids, conn_comps_per_cell)
        if length(conn_comps) < 2
            continue
        end

        largest_cc_id = findmax(length.(conn_comps))[2]
        for (ci, c_ids) in enumerate(conn_comps)
            if ci == largest_cc_id
                continue
            end

            mol_ids = mol_ids_per_cell[cell_id][c_ids]

            if add_new_components
                append_empty_component!(data)
                data.assignment[mol_ids] .= length(data.components)
            else
                data.assignment[mol_ids] .= 0
            end
        end
    end
end

log_em_state(data::BmmData, iter_num::Int, time_start::DateTime) =
    @info "BMM part done for $(now() - time_start) in $iter_num iterations. #Components: $(sum(num_of_molecules_per_cell(data) .> 0)). " *
        "Noise level: $(round(mean(data.assignment .== 0) * 100, digits=3))%"

track_progress!(progress::Nothing) = nothing
track_progress!(progress::RemoteChannel) = put!(progress, true)
track_progress!(progress::Progress) = next!(progress)

function bmm!(data::BmmData; min_molecules_per_cell::Int, n_iters::Int=1000, log_step::Int=4, verbose=true,
              new_component_frac::Float64=0.05, new_component_weight::Float64=0.3,
              assignment_history_depth::Int=0, trace_components::Bool=false, progress::Union{Progress, RemoteChannel, Nothing}=nothing,
              component_split_step::Int=3)
    time_start = now()

    if verbose
        println("BMM fit started")
    end

    if (assignment_history_depth > 0) && !(:assignment_history in keys(data.tracer))
        data.tracer[:assignment_history] = Vector{Int}[]
    end

    if trace_components && !(:component_history in keys(data.tracer))
        data.tracer[:component_history] = Vector{Component}[]
    end

    it_num = 0
    trace_n_components!(data, min_molecules_per_cell);

    maximize!(data)

    for i in 1:n_iters
        append_empty_components!(data, new_component_frac)
        update_prior_probabilities!(data.components, new_component_weight)
        update_n_mols_per_segment!(data)

        expect_dirichlet_spatial!(data)

        if (i % component_split_step == 0) || (i == n_iters)
            split_cells_by_connected_components!(data; add_new_components=(new_component_frac > 1e-10), min_molecules_per_cell=(i == n_iters ? 0 : min_molecules_per_cell))
        end

        drop_unused_components!(data)
        maximize!(data)

        it_num = i
        if verbose && i % log_step == 0
            log_em_state(data, it_num, time_start)
        end

        trace_n_components!(data, min_molecules_per_cell);
        trace_assignment_history!(data, assignment_history_depth; use_guids=!trace_components)
        if trace_components
            trace_component_history!(data)
        end

        track_progress!(progress)
    end

    if verbose
        log_em_state(data, it_num, time_start)
    end

    return data
end

run_bmm_parallel(bm_data_arr, args...; kwargs...) =
    run_bmm_parallel!(deepcopy(bm_data_arr), args...; kwargs...)

function run_bmm_parallel!(bm_data_arr::Array{BmmData, 1}, n_iters::Int; min_molecules_per_cell::Int, kwargs...)::BmmData
    @info "Algorithm start"
    p = Progress(n_iters * length(bm_data_arr))
    @threads for i in 1:length(bm_data_arr)
        bmm!(bm_data_arr[i]; n_iters=n_iters, min_molecules_per_cell=min_molecules_per_cell, verbose=false, progress=p, kwargs...)
    end

    bm_data_merged = merge_bm_data(bm_data_arr)

    if :assignment_history in keys(bm_data_merged.tracer)
        bm_data_merged.assignment = estimate_assignment_by_history(bm_data_merged)[1];
        maximize!(bm_data_merged)
    end

    drop_unused_components!(bm_data_merged; min_n_samples=1)

    @info "Done!"

    return bm_data_merged
end