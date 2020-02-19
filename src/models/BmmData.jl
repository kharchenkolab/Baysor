using DataFrames
using DataFramesMeta
using LinearAlgebra
using NearestNeighbors
using Statistics
using StatsBase: countmap

mutable struct BmmData

    # Static data
    x::DataFrame;
    position_data::Array{Float64, 2};
    composition_data::Array{Int, 1};
    confidence::Array{Float64, 1}

    adjacent_points::Array{Array{Int, 1}, 1};
    adjacent_weights::Array{Array{Float64, 1}, 1};
    real_edge_weight::Float64;
    position_knn_tree::KDTree;
    knn_neighbors::Array{Array{Int, 1}, 1};

    # Distribution-related
    components::Array{Component, 1};
    distribution_sampler::Component;
    assignment::Array{Int, 1};
    max_component_guid::Int;

    noise_density::Float64;

    gene_probs_given_single_transcript::Array{Float64, 2};

    # Utils
    tracer::Dict{String, Any};

    # Parameters
    update_priors::Symbol;

    """
    ...
    # Arguments
    - `components::Array{Component, 1}`:
    - `x::DataFrame`:
    - `adjacent_points::Array{Array{Int, 1}, 1}`:
    - `adjacent_weights::Array{Array{Float64, 1}, 1}`: edge weights, used for smoothness penalty
    - `real_edge_weight::Float64`: weight of an edge for "average" real point
    - `distribution_sampler::Component`:
    - `assignment::Array{Int, 1}`:
    - `k_neighbors::Int=20`:
    - `update_priors::Symbol=:no`: method of prior updates. Possible values: `:no` (no update), `:all` (use all distribitions) and `:centers` (only use distributions, based on prior centers)
    """
    function BmmData(components::Array{Component, 1}, x::DataFrame, adjacent_points::Array{Array{Int, 1}, 1},
                     adjacent_weights::Array{Array{Float64, 1}, 1}, real_edge_weight::Float64, distribution_sampler::Component,
                     assignment::Array{Int, 1}; k_neighbors::Int=20, update_priors::Symbol=:no)
        @assert maximum(assignment) <= length(components)
        @assert minimum(assignment) >= 0
        @assert length(assignment) == size(x, 1)

        if !all(s in names(x) for s in [:x, :y, :gene])
            error("`x` data frame must have columns 'x', 'y' and 'gene'")
        end

        p_data = position_data(x)
        position_knn_tree = KDTree(p_data)
        knn_neighbors = knn(position_knn_tree, p_data, k_neighbors)[1]

        n_genes = maximum(composition_data(x))

        if update_priors == :centers && all(c.can_be_dropped for c in components)
            update_priors = :no
        end

        x = deepcopy(x)
        if !(:confidence in names(x))
            x[:confidence] = 0.95
        end

        self = new(x, p_data, composition_data(x), confidence(x), adjacent_points, adjacent_weights, real_edge_weight,
                   position_knn_tree, knn_neighbors, components, deepcopy(distribution_sampler), assignment, length(components),
                   0.0, ones(n_genes, n_genes) ./ n_genes, Dict{String, Any}(), update_priors)

        for c in self.components
            c.n_samples = 0
        end

        for c_id in assignment[assignment .> 0]
            self.components[c_id].n_samples += 1
        end

        # Resolve component guids
        guids = [c.guid for c in self.components]

        if maximum(guids) <= 0
            self.max_component_guid = length(self.components)
            for (i,c) in enumerate(self.components)
                c.guid = i
            end
        else
            if minimum(guids) <= 0
                error("Either all or no guids can be <= 0")
            end

            self.max_component_guid = maximum(guids)
        end

        return self
    end
end

position_data(df::AbstractDataFrame)::Array{Float64, 2} = Matrix{Float64}(df[:, [:x, :y]])'
position_data(data::BmmData)::Array{Float64, 2} = data.position_data
composition_data(df::AbstractDataFrame)::Array{Int, 1} = df[!, :gene]
composition_data(data::BmmData)::Array{Int, 1} = data.composition_data
confidence(df::AbstractDataFrame)::Array{Float64, 1} = df[!, :confidence]
confidence(data::BmmData)::Array{Float64, 1} = data.confidence

num_of_molecules_per_cell(data::BmmData) = count_array(data.assignment .+ 1, max_value=length(data.components) + 1)[2:end]

function assign!(data::BmmData, point_ind::Int, component_id::Int)
    @assert component_id <= length(data.components) "Too large component id: $component_id, maximum available: $(length(data.components))"
    data.assignment[point_ind] = component_id
end

function merge_bm_data(bmm_data_arr::Array{BmmData, 1}; reestimate_triangulation::Bool=false)
    @assert length(bmm_data_arr) > 0

    # Spatail DataFrame
    x = vcat([deepcopy(bd.x) for bd in bmm_data_arr]...)

    # Components
    components = [deepcopy(bd.components) for bd in bmm_data_arr]
    tracers = [deepcopy(bd.tracer) for bd in bmm_data_arr]

    ## Update GUIDs
    max_guid = 0
    for (bd, comps, tr) in zip(bmm_data_arr, components, tracers)
        for c in comps
            c.guid += max_guid
        end

        if "assignment_history" in keys(tr)
            for ah in tr["assignment_history"]
                ah[ah .> 0] .+= max_guid
            end
        end

        max_guid += bd.max_component_guid
    end

    components = vcat(components...)

    # Adjacency lists
    adjacent_points = Array{Int64,1}[]
    adjacent_weights = Array{Float64,1}[]

    if reestimate_triangulation
        adjacent_points, adjacent_weights = adjacency_list(x)
    else
        ap_offset = 0;
        for bd in bmm_data_arr
            append!(adjacent_points, [deepcopy(ap) .+ ap_offset for ap in bd.adjacent_points])
            append!(adjacent_weights, deepcopy(bd.adjacent_weights))
            ap_offset += size(bd.x, 1)
        end
    end

    # Assignments
    assignments = Array{Int64,1}[]
    assignment_offset = 0;
    for bd in bmm_data_arr
        cur_assignment = deepcopy(bd.assignment) .+ assignment_offset
        cur_assignment[bd.assignment .== 0] .= 0
        push!(assignments, cur_assignment)
        assignment_offset += length(bd.components)
    end

    k_neighbors=length(bmm_data_arr[1].knn_neighbors[1])

    res = BmmData(components, x, adjacent_points, adjacent_weights, bmm_data_arr[1].real_edge_weight,
        deepcopy(bmm_data_arr[1].distribution_sampler), vcat(assignments...); k_neighbors=k_neighbors,
        update_priors=bmm_data_arr[1].update_priors)

    res.tracer = merge_tracers(tracers)

    return res
end

function estimate_assignment_by_history(data::BmmData)
    if !("assignment_history" in keys(data.tracer)) || (length(data.tracer["assignment_history"]) == 0)
        @warn "Data has no saved history of assignments. Fall back to the basic assignment"
        return data.assignment, ones(length(data.assignment)) / 2
    end

    guid_map = Dict(c.guid => i for (i,c) in enumerate(data.components))
    current_guids = Set(vcat(collect(keys(guid_map)), [0]));
    assignment_mat = hcat(data.tracer["assignment_history"]...);

    reassignment = mapslices(assignment_mat, dims=2) do row
        c_row = row[in.(row, Ref(current_guids))]
        if length(c_row) == 0
            return 0
        end

        c_counts = countmap(c_row);
        count_vals = collect(values(c_counts))
        return maximum(collect(keys(c_counts))[count_vals .== maximum(count_vals)])
    end

    return get.(Ref(guid_map), vec(reassignment), 0), vec(mean(assignment_mat .== reassignment, dims=2))
end

function get_cell_stat_df(data::BmmData, segmented_df::Union{DataFrame, Nothing}=nothing; add_qc::Bool=true)
    df = DataFrame(:cell => 1:length(data.components))

    centers = hcat([c.position_params.μ for c in data.components]...)

    for s in [:x, :y]
        df[!,s] = mean.(split(data.x[!,s], data.assignment .+ 1, max_factor=length(data.components) + 1)[2:end])
    end

    df[!,:has_center] = [c.center_prior !== nothing for c in data.components]
    df[!,:x_prior], df[!,:y_prior] = [[(c.center_prior === nothing) ? NaN : c.center_prior.μ[i] for c in data.components] for i in 1:2]

    if add_qc
        if segmented_df === nothing
            segmented_df = get_segmentation_df(data);
        end
        seg_df_per_cell = split(segmented_df, segmented_df.cell .+ 1; max_factor=length(data.components)+1)[2:end];
        pos_data_per_cell = position_data.(seg_df_per_cell);

        df[!,:n_transcripts] = size.(pos_data_per_cell, 2);
        large_cell_mask = (df.n_transcripts .> 2)

        df[!,:density] = fill(NaN, size(df, 1))
        df[!,:elongation] = fill(NaN, size(df, 1))
        df[!,:area] = fill(NaN, size(df, 1))

        df.area[large_cell_mask] = area.(convex_hull.(pos_data_per_cell[large_cell_mask]));
        df.density[large_cell_mask] = df.n_transcripts[large_cell_mask] ./ df.area[large_cell_mask];
        df.elongation[large_cell_mask] = [x[2] / x[1] for x in eigvals.(cov.(transpose.(pos_data_per_cell[large_cell_mask])))];
        if :confidence in names(segmented_df)
            df[!,:avg_confidence] = [mean(df.confidence) for df in seg_df_per_cell]
        end
    end

    return df[num_of_molecules_per_cell(data) .> 0,:]
end

function get_segmentation_df(data::BmmData, gene_names::Union{Nothing, Array{String, 1}}=nothing; use_assignment_history::Bool=true)
    df = deepcopy(data.x)
    df[!,:cell] = deepcopy(data.assignment);

    if use_assignment_history && ("assignment_history" in keys(data.tracer)) && (length(data.tracer["assignment_history"]) > 1)
        df[!,:cell], df[!,:assignment_confidence] = estimate_assignment_by_history(data)
    end

    df[!,:is_noise] = (df.cell .== 0);

    if gene_names !== nothing
        df[!,:gene] = gene_names[df[!,:gene]]
    end

    return df
end

function global_assignment_ids(data::BmmData)::Vector{Int}
    cur_guids = [c.guid for c in data.components]
    res = deepcopy(data.assignment)
    non_noise_mask = (res .> 0)
    res[non_noise_mask] .= cur_guids[res[non_noise_mask]]

    return res
end

function convert_segmentation_df_to_cm(segmentation_df::DataFrame; noise_id::Int=0) # TODO: move to other file
    coll_df = deepcopy(segmentation_df[:, [:cell, :gene]])
    coll_df[!, :val] .= 1
    cm = unstack(@by(coll_df, [:cell, :gene], val=sum(:val)), :cell, :gene, :val);
    return coalesce.(@where(cm, :cell .!= noise_id), 0);
end