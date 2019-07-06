# function trace_n_components!(data::BmmData, min_molecules_per_cell::Int, iter_num::Int)
function trace_n_components!(data::BmmData, min_molecules_per_cell::Int)
    if !("n_components" in keys(data.tracer))
        data.tracer["n_components"] = Dict{Int, Int}[]
    end

    trace_nums = unique(max.(round.(Int, [0.0, 0.5, 1.0, 2.0, 5.0] .* min_molecules_per_cell), 1))
    push!(data.tracer["n_components"], Dict(tn => sum(num_of_molecules_per_cell(data) .>= tn) for tn in trace_nums));
end

function trace_prior_shape!(data::BmmData)
    if !("prior_shape" in keys(data.tracer))
        data.tracer["prior_shape"] = Array{Float64, 1}[]
    end

    push!(data.tracer["prior_shape"], data.distribution_sampler.shape_prior.eigen_values);
end

function merge_tracers(tracers::Array{Dict{String, Any}, 1})::Dict{String, Any}
    @warn "`merge_tracers` doesn't merge 'prior_shape', as it can't be aggregated over frames" # TODO: fix it

    tracers = filter(tr -> "n_components" in keys(tr), tracers)
    if length(tracers) == 0
        return Dict{String, Any}()
    end

    n_components_per_tracer = [tr["n_components"] for tr in tracers]
    n_components_res = deepcopy(n_components_per_tracer[1])

    for it in 1:length(n_components_res)
        for k in keys(n_components_res[it])
            for nc in n_components_per_tracer[2:end]
                n_components_res[it][k] += nc[it][k]
            end
        end
    end

    return Dict{String, Any}("n_components" => n_components_res)
end