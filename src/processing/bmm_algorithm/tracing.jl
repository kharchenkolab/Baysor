function trace_n_components!(data::BmmData, min_molecules_per_cell::Int)
    if !(:n_components in keys(data.tracer))
        data.tracer[:n_components] = Dict{Int, Int}[]
    end

    trace_nums = unique(max.(round.(Int, [0.0, 0.5, 1.0, 2.0, 5.0] .* min_molecules_per_cell), 1))
    n_mols_per_cell = num_of_molecules_per_cell(data)
    push!(data.tracer[:n_components], Dict(tn => sum(n_mols_per_cell .>= tn) for tn in trace_nums));
end

function trace_assignment_history!(data::BmmData, assignment_history_depth::Int)
    (assignment_history_depth > 0) || return

    push!(data.tracer[:assignment_history], global_assignment_ids(data))
    if length(data.tracer[:assignment_history]) > assignment_history_depth
        data.tracer[:assignment_history] = data.tracer[:assignment_history][(end - assignment_history_depth + 1):end]
    end
end

function merge_tracers(tracers::Array{Dict{Symbol, Any}, 1})::Dict{Symbol, Any}
    res = Dict{Symbol, Any}();

    tracers = filter(tr -> :n_components in keys(tr), tracers)
    if length(tracers) == 0
        return res
    end

    n_components_per_tracer = [tr[:n_components] for tr in tracers]
    n_components_res = deepcopy(n_components_per_tracer[1])

    for it in 1:length(n_components_res)
        for k in keys(n_components_res[it])
            for nc in n_components_per_tracer[2:end]
                n_components_res[it][k] += nc[it][k]
            end
        end
    end
    res[:n_components] = n_components_res

    assignments = get.(tracers, Ref(:assignment_history), Ref([]))
    if any(length.(assignments) .> 0)
        if !all(length.(assignments) .== length(assignments[1]))
            @warn "Tracers have different length of assignment_history. Can't merge it."
        else
            res[:assignment_history] = vcat.(assignments...)
        end
    end

    return res
end