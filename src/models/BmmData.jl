using DataFrames
using DataFramesMeta
using LinearAlgebra
using NearestNeighbors
using Statistics
using StatsBase: countmap

mutable struct BmmData

    # Static data
    x::DataFrame;
    position_data::Matrix{Float64};
    composition_data::Vector{Int};
    confidence::Vector{Float64}

    adjacent_points::Array{Vector{Int}, 1};
    adjacent_weights::Array{Vector{Float64}, 1};
    real_edge_weight::Float64;
    position_knn_tree::KDTree;
    knn_neighbors::Array{Vector{Int}, 1};

    # Distribution-related
    components::Array{Component, 1};
    distribution_sampler::Component;
    assignment::Vector{Int};
    max_component_guid::Int;

    noise_density::Float64;

    cluster_per_molecule::Vector{Int};
    cluster_per_cell::Vector{Int};
    gene_probs_given_single_transcript::Matrix{Float64}; # DEPRECATED?

    center_sample_cache::Vector{Int}

    # Utils
    tracer::Dict{Symbol, Any};
    misc::Dict{Symbol, Any};

    # Parameters
    update_priors::Symbol; # DEPRECATED?
    cluster_penalty_mult::Float64;
    use_gene_smoothing::Bool;

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
    function BmmData(components::Array{Component, 1}, x::DataFrame, adjacent_points::Array{Array{Int, 1}, 1}, adjacent_weights::Array{Array{Float64, 1}, 1},
                     real_edge_weight::Float64, distribution_sampler::Component, assignment::Array{Int, 1};
                     k_neighbors::Int=20, update_priors::Symbol=:no, cluster_per_molecule::Union{Symbol, Vector{Int}}=:cluster, cluster_penalty_mult::Float64=0.25,
                     cluster_per_cell::Vector{Int}=Vector{Int}(), use_gene_smoothing::Bool=true)
        @assert maximum(assignment) <= length(components)
        @assert minimum(assignment) >= 0
        @assert length(assignment) == size(x, 1)

        if !all(s in names(x) for s in [:x, :y, :gene])
            error("`x` data frame must have columns 'x', 'y' and 'gene'")
        end

        if isa(cluster_per_molecule, Symbol)
            if cluster_per_molecule in names(x)
                cluster_per_molecule = x[:, cluster_per_molecule]
            else
                cluster_per_molecule = Vector{Int}()
            end
        elseif (length(cluster_per_molecule) > 0) && (length(cluster_per_molecule) != size(x, 1))
            error("cluster_per_molecule has length $(length(cluster_per_molecule)), but $(size(x, 1)) is expected")
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
            x[!, :confidence] .= 0.95
        end

        self = new(x, p_data, composition_data(x), confidence(x), adjacent_points, adjacent_weights, real_edge_weight,
                   position_knn_tree, knn_neighbors, components, deepcopy(distribution_sampler), assignment, length(components),
                   0.0, cluster_per_molecule, deepcopy(cluster_per_cell), ones(n_genes, n_genes) ./ n_genes, Int[],
                   Dict{Symbol, Any}(), Dict{Symbol, Any}(), update_priors, cluster_penalty_mult, use_gene_smoothing)

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

position_data(df::AbstractDataFrame)::Matrix{Float64} = Matrix{Float64}(df[:, [:x, :y]])'
position_data(data::BmmData)::Matrix{Float64} = data.position_data
composition_data(df::AbstractDataFrame)::Vector{Int} = df.gene
composition_data(data::BmmData)::Vector{Int} = data.composition_data
confidence(df::AbstractDataFrame)::Vector{Float64} = df.confidence
confidence(data::BmmData)::Vector{Float64} = data.confidence

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

        if :assignment_history in keys(tr)
            for ah in tr[:assignment_history]
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

    cluster_per_molecule = vcat([bmd.cluster_per_molecule for bmd in bmm_data_arr]...)
    cluster_per_cell = vcat([bmd.cluster_per_cell for bmd in bmm_data_arr]...)

    res = BmmData(components, x, adjacent_points, adjacent_weights, bmm_data_arr[1].real_edge_weight,
        deepcopy(bmm_data_arr[1].distribution_sampler), vcat(assignments...); k_neighbors=k_neighbors,
        update_priors=bmm_data_arr[1].update_priors, cluster_per_molecule=cluster_per_molecule,
        cluster_penalty_mult=bmm_data_arr[1].cluster_penalty_mult, cluster_per_cell=cluster_per_cell,
        use_gene_smoothing=bmm_data_arr[1].use_gene_smoothing)

    res.tracer = merge_tracers(tracers)

    return res
end

function estimate_assignment_by_history(data::BmmData)
    if !(:assignment_history in keys(data.tracer)) || (length(data.tracer[:assignment_history]) == 0)
        @warn "Data has no saved history of assignments. Fall back to the basic assignment"
        return data.assignment, ones(length(data.assignment)) / 2
    end

    guid_map = Dict(c.guid => i for (i,c) in enumerate(data.components))
    current_guids = Set(vcat(collect(keys(guid_map)), [0]));
    assignment_mat = hcat(data.tracer[:assignment_history]...);

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

function get_cell_stat_df(data::BmmData, segmented_df::Union{DataFrame, Nothing}=nothing; add_qc::Bool=true, do_score_doublets::Bool=true, min_molecules_per_cell::Int=3, sigdigits::Int=4)
    df = DataFrame(:cell => 1:length(data.components))

    centers = hcat([Vector(c.position_params.μ) for c in data.components]...)

    for s in [:x, :y]
        df[!,s] = mean.(split(data.x[!,s], data.assignment .+ 1, max_factor=length(data.components) + 1)[2:end])
    end

    df[!,:has_center] = [c.center_prior !== nothing for c in data.components]

    if any([c.center_prior !== nothing for c in data.components])
        df[!,:x_prior], df[!,:y_prior] = [[(c.center_prior === nothing) ? NaN : c.center_prior.μ[i] for c in data.components] for i in 1:2]
    end

    if !isempty(data.cluster_per_cell)
        df[!, :cluster] = data.cluster_per_cell
    end

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

        df.area[large_cell_mask] = round.(area.(convex_hull.(pos_data_per_cell[large_cell_mask])), sigdigits=sigdigits);
        df.density[large_cell_mask] = round.(df.n_transcripts[large_cell_mask] ./ df.area[large_cell_mask], sigdigits=sigdigits);
        df.elongation[large_cell_mask] = [round(x[2] / x[1], sigdigits=sigdigits) for x in eigvals.(cov.(transpose.(pos_data_per_cell[large_cell_mask])))];
        if :confidence in names(segmented_df)
            df[!,:avg_confidence] = round.([mean(df.confidence) for df in seg_df_per_cell], sigdigits=sigdigits)
        end

        if do_score_doublets & !isempty(data.cluster_per_molecule)
            cluster_centers = convert_segmentation_to_counts(composition_data(data), data.cluster_per_molecule);
            cluster_centers = cluster_centers' ./ sum(cluster_centers', dims=2);

            cm = convert_segmentation_to_counts(composition_data(data), data.assignment);
            doublet_fractions = score_doublets(cm, cluster_centers; min_molecules_per_cell=min_molecules_per_cell)[2];
            df[!, :doublet_score] = round.(min.(doublet_fractions, 0.5) .* 2, sigdigits=sigdigits)
        end
    end

    return df[num_of_molecules_per_cell(data) .> 0,:]
end

function get_segmentation_df(data::BmmData, gene_names::Union{Nothing, Array{String, 1}}=nothing; use_assignment_history::Bool=true)
    df = deepcopy(data.x)
    df[!,:cell] = deepcopy(data.assignment);

    if use_assignment_history && (:assignment_history in keys(data.tracer)) && (length(data.tracer[:assignment_history]) > 1)
        df[!,:cell], df[!,:assignment_confidence] = estimate_assignment_by_history(data)
        df.assignment_confidence .= round.(df.assignment_confidence, digits=5)
    end

    if :confidence in names(df)
        df.confidence = round.(df.confidence, digits=5)
    end

    df[!,:is_noise] = (df.cell .== 0);

    if gene_names !== nothing
        df[!,:gene] = gene_names[df[!,:gene]]
    end

    if !isempty(data.cluster_per_molecule)
        df[!, :cluster] = data.cluster_per_molecule
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
