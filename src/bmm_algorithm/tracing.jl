# function trace_n_components!(data::BmmData, min_molecules_per_cell::Int, iter_num::Int)
function trace_n_components!(data::BmmData, min_molecules_per_cell::Int)
    if !(:n_components in keys(data.tracer))
        data.tracer[:n_components] = Dict{Int, Int}[]
    end

    trace_nums = unique(max.(round.(Int, [0.0, 0.5, 1.0, 2.0, 5.0] .* min_molecules_per_cell), 1))
    n_mols_per_cell = num_of_molecules_per_cell(data)
    push!(data.tracer[:n_components], Dict(tn => sum(n_mols_per_cell .>= tn) for tn in trace_nums));
end

function trace_prior_shape!(data::BmmData)
    if !(:prior_shape in keys(data.tracer))
        data.tracer[:prior_shape] = Array{Float64, 1}[]
    end

    push!(data.tracer[:prior_shape], data.distribution_sampler.shape_prior.std_values);
end

function trace_assignment_history!(data::BmmData, assignment_history_depth::Int; use_guids::Bool=true)
    if assignment_history_depth == 0
        return
    end

    if use_guids
        push!(data.tracer[:assignment_history], global_assignment_ids(data))
    else
        push!(data.tracer[:assignment_history], deepcopy(data.assignment))
    end
    if length(data.tracer[:assignment_history]) > assignment_history_depth
        data.tracer[:assignment_history] = data.tracer[:assignment_history][(end - assignment_history_depth + 1):end]
    end
end

function trace_component_history!(data::BmmData) # TODO: remove it
    push!(data.tracer[:component_history], deepcopy(data.components))
end

function merge_tracers(tracers::Array{Dict{Symbol, Any}, 1})::Dict{Symbol, Any}
    # @warn "`merge_tracers` doesn't merge 'prior_shape', as it can't be aggregated over frames"

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