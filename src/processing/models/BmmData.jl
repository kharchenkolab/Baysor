using LinearAlgebra
using NearestNeighbors
using Statistics
using StatsBase: countmap

mutable struct BmmData{L, CT}
    # Static data
    ## Segmentation data
    x::DataFrame;
    position_data::Matrix{Float64};
    composition_data::Vector{<:Union{Int, Missing}};
    confidence::Vector{Float64}

    cluster_per_molecule::Vector{Int}
    segment_per_molecule::Vector{Int}
    nuclei_prob_per_molecule::Vector{Float64}

    ## MRF
    adj_list::AdjList;
    real_edge_weight::Float64;
    position_knn_tree::KDTree;
    knn_neighbors::Array{Vector{Int}, 1};

    # Distribution-related
    components::Array{Component{L, CT}, 1};
    assignment::Vector{Int};
    max_component_guid::Int;

    noise_position_density::Float64;
    noise_density::Float64;

    cluster_per_cell::Vector{Int};

    # Prior segmentation

    n_molecules_per_segment::Vector{Int}
    main_segment_per_cell::Vector{Int}

    # Utils
    tracer::Dict{Symbol, Any}
    misc::Dict{Symbol, Any}

    # Parameters
    prior_seg_confidence::Float64
    cluster_penalty_mult::Float64
    use_gene_smoothing::Bool
    min_nuclei_frac::Float64
    mrf_strength::Float64
end

"""
...
# Arguments
- `components::Array{Component, 1}`:
- `x::DataFrame`:
- `adj_list::AdjList`:
- `adjacent_weights::Array{Array{Float64, 1}, 1}`: edge weights, used for smoothness penalty
- `real_edge_weight::Float64`: weight of an edge for "average" real point
- `assignment::Array{Int, 1}`:
"""
function BmmData(
        components::Array{Component{N, CT}, 1}, x::DataFrame, adj_list::AdjList, assignment::Vector{Int};
        real_edge_weight::Float64=1.0, k_neighbors::Int=20, noise_position_density::Float64=0.0,
        cluster_penalty_mult::Float64=0.25, use_gene_smoothing::Bool=true, prior_seg_confidence::Float64=0.5,
        min_nuclei_frac::Float64=0.1, mrf_strength::Float64=0.1, na_genes::Vector{Int}=Int[]
    ) where {N, CT}
    @assert maximum(assignment) <= length(components)
    @assert minimum(assignment) >= 0
    @assert length(assignment) == size(x, 1)

    metadata_used = Symbol[];

    if !all(s in propertynames(x) for s in [:x, :y, :gene])
        error("`x` data frame must have columns 'x', 'y' and 'gene'")
    end

    p_data = position_data(x)
    @assert size(p_data, 1) == N

    position_knn_tree = KDTree(p_data)
    knn_neighbors = knn(position_knn_tree, p_data, k_neighbors)[1]

    comp_data = [ifelse(g in na_genes, missing, g) for g in composition_data(x)]

    x = deepcopy(x)
    if !(:confidence in propertynames(x))
        x[!, :confidence] .= 0.95
    end

    nuclei_probs = :nuclei_probs in propertynames(x) ? x[!, :nuclei_probs] : Float64[]
    cluster_per_molecule = :cluster in propertynames(x) ? x.cluster : Int[]
    self = BmmData{N, CT}(
        x, p_data, comp_data, confidence(x), cluster_per_molecule, Int[], nuclei_probs,
        adj_list, real_edge_weight, position_knn_tree, knn_neighbors,
        components, assignment, length(components), noise_position_density, 0.0,
        Int[],
        Int[], Int[], # prior segmentation info
        Dict{Symbol, Any}(), Dict{Symbol, Any}(),
        prior_seg_confidence, cluster_penalty_mult, use_gene_smoothing, min_nuclei_frac, mrf_strength
    )

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

    if :prior_segmentation in propertynames(x)
        self.segment_per_molecule = deepcopy(x.prior_segmentation);
        self.n_molecules_per_segment = count_array(self.segment_per_molecule, drop_zero=true);
        update_n_mols_per_segment!(self);
    end

    for c in [:confidence, :nuclei_probs, :cluster, :prior_segmentation]
        if c in propertynames(x)
            push!(metadata_used, c)
        end
    end

    @info "Using the following additional information about molecules: $(metadata_used)"

    return self
end

function position_data(df::AbstractDataFrame)::Matrix{Float64}
    if (:z in propertynames(df))
        return copy(Matrix{Float64}(df[:, [:x, :y, :z]])')
    end

    return copy(Matrix{Float64}(df[:, [:x, :y]])')
end
@inline position_data(data::BmmData)::Matrix{Float64} = data.position_data

@inline composition_data(df::AbstractDataFrame)::Union{Vector{Int}, Vector{Union{Missing, Int}}} = df.gene
@inline composition_data(data::BmmData{T, CategoricalSmoothed{FT}} where {T, FT}) = data.composition_data

@inline composition_data(data::T where T <: Union{AbstractDataFrame, BmmData{BT, CategoricalSmoothed{FT}}} where {BT, FT}, ids::AbstractVector{Int}) =
    view(composition_data(data), ids)
@inline composition_data(data::BmmData{T, CategoricalSmoothed{FN}} where {T, FN}, id::Int) =
    composition_data(data)[id]

@inline confidence(df::AbstractDataFrame)::Vector{Float64} = df.confidence
@inline confidence(data::BmmData)::Vector{Float64} = data.confidence

num_of_molecules_per_cell(data::BmmData) = count_array(data.assignment, max_value=length(data.components), drop_zero=true)

function assign!(data::BmmData, point_ind::Int, component_id::Int)
    old_id = data.assignment[point_ind]
    if old_id == component_id
        return
    end

    @assert component_id <= length(data.components) "Too large component id: $component_id, maximum available: $(length(data.components))"

    segment_id = isempty(data.segment_per_molecule) ? 0 : data.segment_per_molecule[point_ind]
    if segment_id > 0
        if component_id > 0
            data.components[component_id].n_molecules_per_segment[segment_id] =
                get(data.components[component_id].n_molecules_per_segment, segment_id, 0) + 1
        end

        if old_id > 0
            data.components[old_id].n_molecules_per_segment[segment_id] -= 1
        end
    end

    data.assignment[point_ind] = component_id
end

function estimate_assignment_by_history(data::BmmData)
    """
    Estimate the assignment by the history of assignments. Returns the new assignment (local ids) and its confidence.
    Only the components present in the final assignment are considered.
    """
    # TODO: it doesn't guarantee connectivity. Can try to run deterministic EM, or use some better estimate here
    if !(:assignment_history in keys(data.tracer)) || (length(data.tracer[:assignment_history]) == 0)
        @warn "Data has no saved history of assignments. Fall back to the basic assignment"
        return data.assignment, ones(length(data.assignment)) ./ 2
    end

    assignment_history::Vector{Vector{eltype(data.assignment)}} = data.tracer[:assignment_history]
    guid_map = get_component_guid_map(data)
    current_guids = Set(vcat(collect(keys(guid_map)), [0]));
    assignment_mat::Matrix{eltype(data.assignment)} = hcat(assignment_history...);

    reassignment = Vector{eltype(data.assignment)}(undef, size(assignment_mat, 1))
    for i in 1:size(assignment_mat, 1)
        # Find the most frequent value in the row, but only among the current components
        c_row = filter!(x -> (x in current_guids), assignment_mat[i, :])
        reassignment[i] = (length(c_row) == 0) ? 0 : mode(c_row)
    end

    return get.(Ref(guid_map), reassignment, 0), vec(mean(assignment_mat .== reassignment, dims=2))
end

function get_cell_qc_df(
        segmented_df::DataFrame, cell_assignment::Vector{Int};
        component_lifespan::Union{Nothing, Dict{Int, Int}}=nothing,
        sigdigits::Int=4, max_cell::Int=maximum(cell_assignment),
        dapi_arr::Union{Matrix{<:Real}, Nothing}=nothing # TODO: use it
    )
    seg_df_per_cell = split(segmented_df, cell_assignment; max_factor=max_cell, drop_zero=true);
    pos_data_per_cell = [position_data(df)[1:2,:] for df in seg_df_per_cell];

    df = DataFrame(:n_transcripts => size.(pos_data_per_cell, 2));
    large_cell_mask = (df.n_transcripts .> 2)

    df[!,:density] = fill(NaN, size(df, 1))
    df[!,:elongation] = fill(NaN, size(df, 1))
    df[!,:area] = fill(NaN, size(df, 1))

    df.area[large_cell_mask] = round.(area.(convex_hull.(pos_data_per_cell[large_cell_mask])), sigdigits=sigdigits);
    df.density[large_cell_mask] = round.(df.n_transcripts[large_cell_mask] ./ df.area[large_cell_mask], sigdigits=sigdigits);
    df.elongation[large_cell_mask] = [round(x[2] / x[1], sigdigits=sigdigits) for x in eigvals.(cov.(transpose.(pos_data_per_cell[large_cell_mask])))];
    if :confidence in propertynames(segmented_df)
        df[!,:avg_confidence] = round.([mean(df.confidence) for df in seg_df_per_cell], sigdigits=sigdigits)
    end

    if :assignment_confidence in propertynames(segmented_df)
        df[!,:avg_assignment_confidence] = round.([mean(df.assignment_confidence) for df in seg_df_per_cell], sigdigits=sigdigits)
    end

    if :cluster in propertynames(segmented_df)
        df[!,:max_cluster_frac] = [isempty(df) ? NaN : (maximum(values(countmap(df.cluster))) / length(df.cluster)) for df in seg_df_per_cell]
    end

    if component_lifespan !== nothing
        df[!, :lifespan] = get.(Ref(component_lifespan), 1:max_cell, -1)
    end

    return df
end

function get_cell_stat_df(
        data::BmmData, segmented_df::Union{DataFrame, Nothing}=nothing, assignment::Vector{<:Integer}=data.assignment;
        add_qc::Bool=true, sigdigits::Int=4, run_id::String=""
    )
    df = DataFrame(:cell => get_cell_name.(1:length(data.components); type=:cell, run_id))

    for s in intersect([:x, :y, :z], propertynames(data.x))
        df[!,s] = mean.(split(data.x[!,s], assignment, max_factor=length(data.components), drop_zero=true))
    end

    if !isempty(data.cluster_per_cell)
        df[!, :cluster] = data.cluster_per_cell
    end

    if add_qc
        if segmented_df === nothing
            segmented_df = get_segmentation_df(data; run_id);
        end

        component_lifespan = :assignment_history in keys(data.tracer) ? estimate_component_lifespan(data) : nothing
        df = hcat(df, get_cell_qc_df(
            segmented_df, assignment;
            component_lifespan=component_lifespan,
            sigdigits=sigdigits, max_cell=length(data.components)
        ))
    end

    return df[num_of_molecules_per_cell(data) .> 0,:]
end

function get_segmentation_df(data::BmmData, gene_names::Union{Nothing, Array{String, 1}}=nothing; run_id::String="")
    df = deepcopy(data.x)
    df[!,:cell] = get_cell_name.(data.assignment; type=:cell, run_id);

    if (:assignment_history in keys(data.tracer)) && (length(data.tracer[:assignment_history]) > 1)
        # data.assignment is already adjusted based on the history if `refine=true` in `bmm` (default)
        assignment_mat::Matrix{eltype(data.assignment)} = hcat(data.tracer[:assignment_history]...)
        df[!,:assignment_confidence] = round.(mean(assignment_mat .== global_assignment_ids(data), dims=2)[:], digits=5)
    end

    if :confidence in propertynames(df)
        df.confidence = round.(df.confidence, digits=5)
    end

    df[!,:is_noise] = (data.assignment .== 0);

    if gene_names !== nothing
        df[!,:gene] = gene_names[df[!,:gene]]
    end

    if !isempty(data.cluster_per_molecule)
        df[!, :cluster] = data.cluster_per_molecule
    end

    return df
end

function convert_segmentation_to_counts(
        genes::Vector{Int}, cell_assignment::Vector{Int}, gene_names::Nothing=nothing;
        drop_empty_labels::Bool=false, n_genes::Int=maximum(genes)
    )
    if drop_empty_labels
        if minimum(cell_assignment) == 0
            cell_assignment = denserank(cell_assignment) .- 1
        else
            cell_assignment = denserank(cell_assignment)
        end
    end

    cm = sphstack(
        count_array_sparse.(split(genes, cell_assignment, drop_zero=true), total=n_genes)
    )

    return cm
end

function convert_segmentation_to_counts(
        genes::Vector{Int}, cell_assignment::Vector{Int}, gene_names::Vector{String};
        run_id::String="", kwargs...
    )
    # DEPRECATED?
    cm = convert_segmentation_to_counts(genes, cell_assignment; n_genes=length(gene_names), kwargs...) |> Matrix

    cm = DataFrame(cm, get_cell_name.(1:size(cm, 2); type=:cell, run_id))
    cm[!, :gene] = gene_names
    cm = cm[:, vcat(end, 1:end-1)]

    return cm
end

function get_segmentation_results(data::BmmData, gene_names::Union{Vector{String}, Nothing}=nothing; run_id::String="", add_qc::Bool=true)
    segmented_df = get_segmentation_df(data, gene_names; run_id)
    cell_stat_df = get_cell_stat_df(data, segmented_df; add_qc=add_qc, run_id)
    cm = convert_segmentation_to_counts(data.x.gene, data.assignment) # We don't pass gene_names here to keep the matrix sparse

    return segmented_df, cell_stat_df, cm
end

function global_assignment_ids(data::BmmData)::Vector{Int}
    cur_guids = [c.guid for c in data.components]
    res = deepcopy(data.assignment)
    non_noise_mask = (res .> 0)
    res[non_noise_mask] .= cur_guids[res[non_noise_mask]]

    return res
end

function update_n_mols_per_segment!(bm_data::BmmData)
    if isempty(bm_data.segment_per_molecule)
        return
    end

    # Estimate number of molecules per segment per component
    for comp in bm_data.components
        empty!(comp.n_molecules_per_segment)
    end

    for i in 1:length(bm_data.assignment)
        c_cell = bm_data.assignment[i]
        c_seg = bm_data.segment_per_molecule[i]
        if (c_cell == 0) || (c_seg == 0)
            continue
        end

        bm_data.components[c_cell].n_molecules_per_segment[c_seg] = get(bm_data.components[c_cell].n_molecules_per_segment, c_seg, 0) + 1
    end

    # Estimate the main segment per cell
    resize!(bm_data.main_segment_per_cell, length(bm_data.components))
    bm_data.main_segment_per_cell .= 0

    for ci in 1:length(bm_data.components)
        if isempty(bm_data.components[ci].n_molecules_per_segment)
            continue
        end

        i_max = 0
        f_max = 0.0
        for (si,nms) in bm_data.components[ci].n_molecules_per_segment
            seg_size = bm_data.n_molecules_per_segment[si]
            if ((nms / seg_size) > (f_max + 1e-10)) || ((nms == seg_size) && (seg_size > bm_data.n_molecules_per_segment[i_max]))
                f_max = nms / seg_size
                i_max = si
            end
        end
        bm_data.main_segment_per_cell[ci] = i_max
    end
end

function estimate_component_lifespan(data::BmmData; use_guids::Bool=false)
    comp_lifespans = estimate_component_lifespan(data.tracer[:assignment_history])
    use_guids && return comp_lifespans

    guid_map = get_component_guid_map(data)
    comp_lifespans = Dict(guid_map[k] => v for (k,v) in comp_lifespans if k in keys(guid_map));
    return comp_lifespans
end

function estimate_component_lifespan(assignment_history::Vector{Vector{Int}})
    comp_lifespans = Dict{eltype(assignment_history[1]), Int}()
    for it in reverse(eachindex(assignment_history))
        for cmp in assignment_history[it]
            (cmp > 0) || continue
            comp_lifespans[cmp] = length(assignment_history) - it + 1
        end
    end

    return comp_lifespans
end

get_component_guid_map(data::BmmData) =
    Dict(c.guid => i for (i,c) in enumerate(data.components))