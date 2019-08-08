using DataFrames
using Distributions
using NearestNeighbors
using StatsBase

using Base.Threads
using Dates: now, DateTime

import DistributedArrays

function drop_unused_components!(data::BmmData; min_n_samples::Int=2, force::Bool=false)
    non_noise_ids = findall(data.assignment .> 0)
    existed_ids = findall([(c.n_samples >= min_n_samples) || (!c.can_be_dropped && !force) for c in data.components])
    id_map = zeros(Int, length(data.components))
    id_map[existed_ids] = 1:length(existed_ids)

    data.assignment[non_noise_ids] = id_map[data.assignment[non_noise_ids]]
    data.components = data.components[existed_ids]

    @assert all(num_of_molecules_per_cell(data) .== [c.n_samples for c in data.components])
end

# point is adjacent to itself
adjacent_component_ids(assignment::Array{Int, 1}, adjacent_points::Array{Int, 1})::Array{Int, 1} =
    [x for x in Set(assignment[adjacent_points]) if x > 0]

"""
...
# Arguments
- `assignment::Array{Int, 1}`:
- `adjacent_points::Array{Int, 1}`:
- `adjacent_weights::Array{Float64, 1}`: weights here mean `1 / distance` between two adjacent points
"""
@inline function adjacent_component_weights(assignment::Array{Int, 1}, adjacent_points::Array{Int, 1}, adjacent_weights::Array{Float64, 1})
    component_weights = Dict{Int, Float64}()
    zero_comp_weight = 0.0
    for (c_id, cw) in zip(assignment[adjacent_points], adjacent_weights)
        if c_id == 0
            zero_comp_weight += cw
            continue
        end

        component_weights[c_id] = get(component_weights, c_id, 0.0) + cw
    end

    return collect(keys(component_weights)), collect(values(component_weights)), zero_comp_weight
end

function expect_dirichlet_spatial!(data::BmmData, adj_classes_global::Dict{Int, Array{Int, 1}}=Dict{Int, Array{Int, 1}}(); stochastic::Bool=true, noise_density_threshold::Float64=1e-100)
    for i in 1:size(data.x, 1)
        x, y = position_data(data)[:,i]
        gene = composition_data(data)[i]
        confidence = data.confidence[i]

        adj_classes, adj_weights, zero_comp_weight = adjacent_component_weights(data.assignment, data.adjacent_points[i], data.adjacent_weights[i])

        adj_global = get(adj_classes_global, i, Int[]);
        if length(adj_global) > 0
            append!(adj_classes, adj_global)
            append!(adj_weights, ones(length(adj_global)) .* data.real_edge_weight)
        end

        denses = confidence .* adj_weights .* [c.prior_probability * pdf(c, x, y, gene) for c in data.components[adj_classes]]

        if sum(denses) < noise_density_threshold
            assign!(data, i, 0) # Noise class
            continue
        end

        if confidence < 1.0
            append!(denses, (1 - confidence) * max(data.real_edge_weight, zero_comp_weight) * data.noise_density)
            append!(adj_classes, 0)
        end

        if !stochastic
            assign!(data, i, adj_classes[findmax(denses)[2]])
            continue
        end

        assign!(data, i, sample(adj_classes, Weights(denses)))
    end
end

function update_prior_probabilities!(components)
    c_weights = [max(c.n_samples, c.prior_weight) for c in components]
    v = [rand(Beta(1 + n, n_sum)) for (n, n_sum) in zip(c_weights[1:end-1], cumsum(c_weights[end:-1:2])[end:-1:1])];

    prior_probs = vcat(1, cumprod(1 .- v[1:end-1])) .* v;
    push!(prior_probs, 1 - sum(prior_probs));

    for (c, p) in zip(components, prior_probs)
        c.prior_probability = p
    end
end

function maximize!(c::Component, pos_data::Array{Float64, 2}, comp_data::Array{Int, 1}, data::BmmData)
    c.n_samples = size(pos_data, 2)

    if size(pos_data, 2) == 0
        return maximize_from_prior!(c, data)
    end

    return maximize!(c, pos_data, comp_data)
end

function maximize_prior!(data::BmmData, min_molecules_per_cell::Int)
    if data.update_priors == :no
        return
    end

    components = (data.update_priors == :all) ? data.components : filter(c -> !c.can_be_dropped, data.components)
    components = filter(c -> c.n_samples >= min_molecules_per_cell, components)

    if length(components) < 2
        return
    end

    mean_shape = vec(median(hcat(eigen_values.(components)...), dims=2))
    set_shape_prior!(data.distribution_sampler, mean_shape)

    for c in data.components
        set_shape_prior!(c, mean_shape)
    end
end

function maximize!(data::BmmData, min_molecules_per_cell::Int; do_maximize_prior::Bool=true)
    ids_by_assignment = split(1:length(data.assignment), data.assignment .+ 1)[2:end]
    append!(ids_by_assignment, [Int[] for i in length(ids_by_assignment):length(data.components)])

    pos_data_by_assignment = [position_data(data)[:, p_ids] for p_ids in ids_by_assignment]
    comp_data_by_assignment = [composition_data(data)[p_ids] for p_ids in ids_by_assignment]

    # @threads for i in 1:length(data.components)
    for i in 1:length(data.components)
        # data.components[i] = maximize(data.components[i], pos_data_by_assignment[i], comp_data_by_assignment[i])
        maximize!(data.components[i], pos_data_by_assignment[i], comp_data_by_assignment[i], data)
    end
    # data.components = pmap(v -> maximize(v...), zip(data.components, pos_data_by_assignment, comp_data_by_assignment))

    if do_maximize_prior
        maximize_prior!(data, min_molecules_per_cell)
    end

    data.noise_density = estimate_noise_density_level(data)
end

function estimate_noise_density_level(data::BmmData)
    composition_density = mean([mean(c.composition_params.counts[c.composition_params.counts .> 0] ./ c.composition_params.n_samples) for c in data.components])

    std_vals = data.distribution_sampler.shape_prior.std_values;
    position_density = pdf(MultivariateNormal([0.0, 0.0], diagm(0 => std_vals.^2)), 3 .* std_vals)

    return position_density * composition_density
end

function append_empty_components!(data::BmmData, new_component_frac::Float64)
    for i in 1:round(Int, new_component_frac * length(data.components))
        data.max_component_guid += 1
        push!(data.components, sample_distribution(data; guid=data.max_component_guid))
    end
end

function get_global_adjacent_classes(data::BmmData)::Dict{Int, Array{Int, 1}}
    adj_classes_global = Dict{Int, Array{Int, 1}}()
    for (cur_id, comp) in enumerate(data.components)
        if comp.n_samples > 0
            continue
        end

        for t_id in data.adjacent_points[knn(data.position_knn_tree, comp.position_params.μ, 1)[1][1]]
            if t_id in keys(adj_classes_global)
                push!(adj_classes_global[t_id], cur_id)
            else
                adj_classes_global[t_id] = [cur_id]
            end
        end
    end

    return adj_classes_global
end

log_em_state(data::BmmData, iter_num::Int, time_start::DateTime) =
    @info "EM part done for $(now() - time_start) in $iter_num iterations. #Components: $(sum(num_of_molecules_per_cell(data) .> 0)). " *
        "Noise level: $(round(mean(data.assignment .== 0) * 100, digits=3))%"

function bmm!(data::BmmData; min_molecules_per_cell::Int, n_iters::Int=1000, log_step::Int=4, verbose=true, new_component_frac::Float64=0.05,
              assignment_history_depth::Int=0, channel::Union{RemoteChannel, Nothing}=nothing)
    time_start = now()

    if verbose
        println("EM fit started")
    end

    it_num = 0

    if !("assignment_history" in keys(data.tracer))
        data.tracer["assignment_history"] = Vector{Int}[]
    end
    adj_classes_global = Dict{Int, Array{Int, 1}}()

    trace_prior_shape!(data);
    trace_n_components!(data, min_molecules_per_cell);

    maximize!(data, min_molecules_per_cell; do_maximize_prior=false)

    for i in 1:n_iters
        expect_dirichlet_spatial!(data, adj_classes_global)
        maximize!(data, min_molecules_per_cell)

        trace_prior_shape!(data);
        trace_n_components!(data, min_molecules_per_cell);

        drop_unused_components!(data)

        append_empty_components!(data, new_component_frac)
        update_prior_probabilities!(data.components)
        adj_classes_global = get_global_adjacent_classes(data)

        it_num = i
        if verbose && i % log_step == 0
            log_em_state(data, it_num, time_start)
        end

        trace_assignment_history!(data, assignment_history_depth)

        if channel !== nothing
            put!(channel, true)
        end
    end

    if verbose
        log_em_state(data, it_num, time_start)
    end

    return data
end

fetch_bm_func(bd::BmmData) = [bd]

function pmap_progress(func, arr::DistributedArrays.DArray, n_steps::Int, args...; kwargs...)
    p = Progress(n_steps)
    channel = RemoteChannel(()->Channel{Bool}(n_steps), 1)

    exception = nothing
    @sync begin # Parallel progress tracing
        @async while take!(channel)
            next!(p)
        end

        @async begin
            try
                futures = [remotecall(func, procs()[i], arr[i], args...; channel=channel, kwargs...) for i in 1:length(arr)];
                arr = fetch.(wait.(futures))
            catch ex
                exception = ex
            finally
                put!(channel, false)
            end
        end
    end

    if exception !== nothing
        throw(exception)
    end

    return arr
end

function push_data_to_workers(bm_data_arr::Array{BmmData, 1})::DistributedArrays.DArray
    if length(bm_data_arr) > nprocs()
        n_new_workers = length(bm_data_arr) - nprocs()
        @info "Creating $n_new_workers new workers. Total $(nprocs() + n_new_workers) workers."
        addprocs(n_new_workers)
        eval(:(@everywhere using Baysor))
    end

    futures = [remotecall(fetch_bm_func, procs()[i], bm_data_arr[i]) for i in 1:length(bm_data_arr)]

    try
        return DistributedArrays.DArray(futures);
    catch ex
        bds = fetch.(wait.(futures))
        err_idx = findfirst(typeof.(bm_data_arr) .!== Baysor.BmmData)
        if err_idx === nothing
            throw(ex)
        end

        error(bds[err_idx])
    end
end

run_bmm_parallel(bm_data_arr, args...; kwargs...) = 
    run_bmm_parallel!(deepcopy(bm_data_arr), args...; kwargs...)

function run_bmm_parallel!(bm_data_arr::Array{BmmData, 1}, n_iters::Int; min_molecules_per_cell::Int, kwargs...)::BmmData
    bm_data_arr = deepcopy(bm_data_arr)
    @info "Pushing data to workers"
    da = push_data_to_workers(bm_data_arr)

    @info "Algorithm start"
    bm_data_arr = pmap_progress(bmm!, da, n_iters * length(da); n_iters=n_iters, min_molecules_per_cell=min_molecules_per_cell, verbose=false, kwargs...)
    bm_data_merged = merge_bm_data(bm_data_arr)

    @info "Done!"

    return bm_data_merged
end