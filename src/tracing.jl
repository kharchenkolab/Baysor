# function trace_n_components!(data::BmmData, min_molecules_per_cell::Int, iter_num::Int)
function trace_n_components!(data, min_molecules_per_cell::Int)
    if !("n_components" in keys(data.tracer))
        data.tracer["n_components"] = Dict{Int, Int}[]
    end

    trace_nums = unique(round.(Int, [0.0, 0.5, 1.0, 2.0, 5.0] .* min_molecules_per_cell))
    push!(data.tracer["n_components"], Dict(tn => sum(num_of_molecules_per_cell(data) .>= tn) for tn in trace_nums));
end