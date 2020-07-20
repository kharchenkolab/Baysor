using DataFrames
import MultivariateStats
using NearestNeighbors
using StatsBase: countmap, denserank
using Statistics

using Base.Threads

import CSV
import Distances

# Parse parameters

function default_param_value(param::Symbol, min_molecules_per_cell::Union{Int, Nothing};
                             n_molecules::Union{Int, Nothing}=nothing, n_genes::Union{Int, Nothing}=nothing)
    if min_molecules_per_cell === nothing
        error("Either `$param` or `min_molecules_per_cell` must be provided")
    end

    min_molecules_per_cell = max(min_molecules_per_cell, 3)

    if param == :min_transcripts_per_segment
        return max(round(Int, min_molecules_per_cell / 4), 2)
    end

    if param == :confidence_nn_id
        return max(div(min_molecules_per_cell, 2) + 1, 5)
    end

    if param == :composition_neighborhood
        if n_genes === nothing
            return max(min_molecules_per_cell, 3)
        end

        return max(div(n_genes, 10), min_molecules_per_cell, 3)
    end

    if param == :n_gene_pcs
        if n_genes === nothing
            error("Either `$param` or `n_genes` must be provided")
        end

        return min(max(div(n_genes, 3), 30), 100, n_genes)
    end

    if param == :n_cells_init
        if n_molecules === nothing
            error("Either `$param` or `n_molecules` must be provided")
        end

        return div(n_molecules, min_molecules_per_cell) * 2
    end
end

default_if_not_provided(value::Union{<:Real, Nothing}, param_name::Symbol, args...; kwargs...) =
    (value === nothing) ? default_param_value(param_name, args...; kwargs...) : value

parse_scale_std(scale_std::Float64, ::Real) = scale_std
parse_scale_std(scale_std::Nothing, scale::Real) = 0.25 * scale
function parse_scale_std(scale_std::String, scale::Real)
    scale_std = strip(scale_std)
    if scale_std[end] == '%'
        return parse(Float64, scale_std[1:end-1]) * scale / 100.0
    end

    return parse(Float64, scale_std)
end

# Initialize cell positions

function position_data_by_assignment(pos_data::T where T<: AbstractMatrix{<:Real}, assignment::Vector{<:Integer})
    filt_ids = findall(assignment .> 0)
    return [pos_data[:, ids] for ids in split(filt_ids, assignment[filt_ids])]
end

function covs_from_assignment(pos_data::T where T<: AbstractMatrix{<:Real}, assignment::Vector{<:Integer}; min_size::Int=1)::Array{CovMat, 1}
    pos_data_by_assignment = position_data_by_assignment(pos_data, assignment)
    covs = estimate_sample_cov.(pos_data_by_assignment);
    mean_covs = MeanVec(vec(median(hcat(Vector.(eigvals.(covs[size.(pos_data_by_assignment, 2) .>= min_size]))...), dims=2)))

    for i in findall(size.(pos_data_by_assignment, 2) .<= min_size)
        covs[i] = CovMat(diagm(0 => deepcopy(mean_covs)))
    end

    return covs
end

cell_centers_uniformly(spatial_df::DataFrame, args...; kwargs...) =
    cell_centers_uniformly(position_data(spatial_df), args..., (:confidence in names(spatial_df)) ? spatial_df.confidence : nothing; kwargs...)

function cell_centers_uniformly(pos_data::T where T<: AbstractMatrix{<:Real}, n_clusters::Int,
        confidences::Union{Vector{Float64}, Nothing}=nothing; scale::Union{<:Real, Nothing})
    n_samples = (confidences === nothing) ? size(pos_data, 2) : ceil(Int, sum(confidences))
    n_clusters = min(n_clusters, n_samples)

    cluster_centers = pos_data[:, select_ids_uniformly(pos_data', confidences; n=n_clusters, confidence_threshold=0.25)]
    cluster_labels = vcat(knn(KDTree(cluster_centers), pos_data, 1)[1]...)

    covs = (scale === nothing) ? covs_from_assignment(pos_data, cluster_labels) : Float64(scale) ^ 2
    return InitialParams(copy(cluster_centers'), covs, cluster_labels)
end

# Build MRF

function convert_edge_list_to_adj_list(edge_list::Matrix{Int}, edge_weights::Union{Vector{Float64}, Nothing}=nothing; n_verts::Int=maximum(edge_list))
    res_ids = [vcat(v...) for v in zip(split(edge_list[2,:], edge_list[1,:], max_factor=n_verts),
                                       split(edge_list[1,:], edge_list[2,:], max_factor=n_verts))];

    if edge_weights === nothing
        return res_ids
    end

    res_weights = [vcat(v...) for v in zip(split(edge_weights, edge_list[1,:], max_factor=n_verts),
                                           split(edge_weights, edge_list[2,:], max_factor=n_verts))];

    return res_ids, res_weights
end

function estimate_local_composition_similarities(df_spatial::DataFrame, edge_list::Matrix{Int}; composition_neighborhood::Int, n_gene_pcs::Int, min_similarity::Float64=0.01)
    neighb_cm = neighborhood_count_matrix(df_spatial, composition_neighborhood);
    neighb_cm = MultivariateStats.transform(MultivariateStats.fit(MultivariateStats.PCA, neighb_cm, maxoutdim=n_gene_pcs), neighb_cm);

    return max.(map(i -> cor(neighb_cm[:, edge_list[1,i]], neighb_cm[:, edge_list[2,i]]), 1:size(edge_list, 2)), min_similarity);
end

function build_molecule_graph(df_spatial::DataFrame; min_edge_quant::Float64=0.3, use_local_gene_similarities::Bool=false, n_gene_pcs::Int=0, composition_neighborhood::Int=0, kwargs...)
    if use_local_gene_similarities && ((composition_neighborhood == 0) || (n_gene_pcs == 0))
        @warn "composition_neighborhood=$composition_neighborhood and n_gene_pcs=$n_gene_pcs, while use_local_gene_similarities=true. Force use_local_gene_similarities=false."
        use_local_gene_similarities = false
    end

    edge_list, adjacent_dists = adjacency_list(df_spatial; kwargs...);

    min_edge_length = quantile(adjacent_dists, min_edge_quant);

    adjacent_weights = min_edge_length ./ max.(adjacent_dists, min_edge_length)
    if use_local_gene_similarities
        adjacent_weights .*= estimate_local_composition_similarities(df_spatial, edge_list; composition_neighborhood=composition_neighborhood, n_gene_pcs=n_gene_pcs)
    end

    adjacent_points, adjacent_weights = convert_edge_list_to_adj_list(edge_list, adjacent_weights; n_verts=size(df_spatial, 1));
    return adjacent_points, adjacent_weights, adjacent_dists
end

function load_df(data_path; x_col::Symbol=:x, y_col::Symbol=:y, gene_col::Symbol=:gene, min_molecules_per_gene::Int=0, kwargs...)
    df_spatial = read_spatial_df(data_path; x_col=x_col, y_col=y_col, gene_col=gene_col, kwargs...)

    gene_counts = StatsBase.countmap(df_spatial[!, :gene]);
    large_genes = Set{String}(collect(keys(gene_counts))[collect(values(gene_counts)) .> min_molecules_per_gene]);
    df_spatial = df_spatial[in.(df_spatial[!, :gene], Ref(large_genes)),:];

    df_spatial[!, :x] = Array{Float64, 1}(df_spatial[!, :x])
    df_spatial[!, :y] = Array{Float64, 1}(df_spatial[!, :y])
    df_spatial[!, :gene], gene_names = encode_genes(df_spatial[!, :gene]);
    return df_spatial, gene_names;
end

# Initialize BmmData

"""
    main function for initialization of bm_data
"""
function initial_distribution_arr(df_spatial::DataFrame; n_frames::Int, n_frames_return::Int=0, n_cells_init::Union{Int, Nothing}=nothing,
                                  scale::T where T<: Real, scale_std::Union{<:Real, String, Nothing}=nothing,
                                  confidence_nn_id::Union{Int, Nothing}=nothing, min_molecules_per_cell::Union{<:Integer, Nothing}=nothing, composition_neighborhood::Union{Int, Nothing}=nothing,
                                  n_gene_pcs::Union{Int, Nothing}=nothing, prior_seg_confidence::Float64=0.5, kwargs...)::Array{BmmData, 1}
    df_spatial = deepcopy(df_spatial)

    ## Parse parameters
    confidence_nn_id = default_if_not_provided(confidence_nn_id, :confidence_nn_id, min_molecules_per_cell)
    n_cells_init = default_if_not_provided(n_cells_init, :n_cells_init, min_molecules_per_cell, n_molecules=size(df_spatial, 1))
    n_cells_init = max(div(n_cells_init, n_frames), 1)

    composition_neighborhood = default_if_not_provided(composition_neighborhood, :composition_neighborhood, min_molecules_per_cell, n_genes=maximum(df_spatial.gene))
    n_gene_pcs = default_if_not_provided(n_gene_pcs, :n_gene_pcs, min_molecules_per_cell, n_genes=maximum(df_spatial.gene))

    scale_std = parse_scale_std(scale_std, scale)

    ## Estimate confidence
    if confidence_nn_id > 0
        @info "Estimate confidence per molecule"
        prior_segmentation = (:prior_segmentation in names(df_spatial)) ? df_spatial.prior_segmentation : nothing
        append_confidence!(df_spatial, prior_segmentation; nn_id=confidence_nn_id, prior_confidence=prior_seg_confidence)
        @info "Done"
    end

    ## Split data
    dfs_spatial = n_frames > 1 ? split_spatial_data(df_spatial, n_frames) : [df_spatial]
    @info "#frames: $(length(dfs_spatial)); mean number of molecules per frame: $(median(size.(dfs_spatial, 1)))."

    if (n_frames_return > 0) && (n_frames_return < n_frames)
        dfs_spatial = dfs_spatial[1:n_frames_return]
    end

    ## Initialize BmmData array
    @info "Initializing algorithm. Scale: $scale, scale std: $scale_std, initial #clusters: $n_cells_init."
    size_prior = ShapePrior(Float64[scale, scale], Float64[scale_std, scale_std], min_molecules_per_cell);

    bm_datas_res = Array{BmmData, 1}(undef, length(dfs_spatial))
    for i in 1:length(dfs_spatial)
        init_params = cell_centers_uniformly(dfs_spatial[i], n_cells_init; scale=scale)
        bm_datas_res[i] = initialize_bmm_data(dfs_spatial[i], init_params; size_prior=size_prior,
            composition_neighborhood=composition_neighborhood, n_gene_pcs=n_gene_pcs, prior_seg_confidence=prior_seg_confidence, kwargs...)
    end

    return bm_datas_res;
end

function initialize_bmm_data(df_spatial::DataFrame, args...; composition_neighborhood::Int=0, n_gene_pcs::Int=0,
        use_local_gene_similarities::Bool=true, adjacency_type::Symbol=:triangulation, prior_seg_confidence::Float64=0.5, kwargs...)::BmmData
    adjacent_points, adjacent_weights, adjacent_dists = build_molecule_graph(df_spatial; use_local_gene_similarities=use_local_gene_similarities,
        n_gene_pcs=n_gene_pcs, composition_neighborhood=composition_neighborhood, adjacency_type=adjacency_type)

    components, sampler, assignment = initial_distributions(df_spatial, args...; kwargs...)

    return BmmData(components, df_spatial, adjacent_points, adjacent_weights, 1.0, sampler, assignment, prior_seg_confidence=prior_seg_confidence)
end

function initial_distributions(df_spatial::DataFrame, initial_params::InitialParams; size_prior::ShapePrior,
                               gene_smooth::Real=1.0, gene_num::Int=maximum(df_spatial[!,:gene]))
    position_distrubutions = [MvNormalF(initial_params.centers[i,:], initial_params.covs[i]) for i in 1:size(initial_params.centers, 1)]
    gene_distributions = [CategoricalSmoothed(ones(Float64, gene_num), smooth=Float64(gene_smooth)) for i in 1:length(position_distrubutions)]

    components = [Component(pd, gd, shape_prior=deepcopy(size_prior)) for (pd, gd) in zip(position_distrubutions, gene_distributions)]

    gene_sampler = CategoricalSmoothed(ones(Float64, gene_num))
    sampler = Component(MvNormalF(zeros(2)), gene_sampler, shape_prior=deepcopy(size_prior)) # position_params are never used

    return components, sampler, initial_params.assignment
end

# Split data by frames

function split_spatial_data(df::DataFrame, n::Int, key::Symbol)::Array{DataFrame, 1}
    factor = vec(sum(hcat([df[!,key] .<= quantile(df[!,key], q) for q in range(1 / n, stop=1.0, length=n)]...), dims=2))
    return split(df, factor)
end

split_spatial_data(df::DataFrame, n_hor::Int, n_ver::Int) = vcat(split_spatial_data.(split_spatial_data(df, n_ver, :y), n_hor, :x)...)

function split_spatial_data(df::DataFrame, n::Int) # TODO: very approximate separation. Example: n=3.
    df_sizes = Dict(s => maximum(df[!,s]) - minimum(df[!,s]) for s in [:x, :y]);
    x_elongation = df_sizes[:x] / sum(values(df_sizes))
    a = round(Int, sqrt(x_elongation * n / (1 - x_elongation))) # solution of "a * b = n; a / (a + b) = x_elongation"

    return split_spatial_data(df, a, max(floor(Int, n / a), 1))
end

split_spatial_data(df::DataFrame; mean_mols_per_frame::Int) = split_spatial_data(df, round(Int, size(df, 1) / mean_mols_per_frame))

# Utils

function staining_value_per_transcript(df_spatial::DataFrame, staining::MT where MT <: AbstractMatrix{T}; quiet::Bool=false)::Vector{T} where T <: Real
    x_vals = round.(Int, df_spatial.x)
    y_vals = round.(Int, df_spatial.y)
    if !quiet && ((maximum(x_vals) > size(staining, 2)) || (maximum(y_vals) > size(staining, 1)))
        @warn "Maximum transcript coordinates are $((maximum(y_vals), maximum(x_vals))), which is larger than the DAPI size: $(size(staining)). Filling it with 0."
    end

    if !quiet && ((minimum(x_vals) < 1) || (minimum(y_vals) < 1))
        @warn "Minimum transcript coordinates are < 1: $((minimum(y_vals), minimum(x_vals))). Filling it with 0."
    end

    if !quiet && ((maximum(x_vals) < 0.5 * size(staining, 2)) || (maximum(y_vals) < 0.5 * size(staining, 1)))
        @warn "Maximum transcript coordinates are $((maximum(y_vals), maximum(x_vals))), which is much smaller than the DAPI size: $(size(staining)). May be result of an error."
    end

    inds = findall((x_vals .> 0) .& (x_vals .< size(staining, 2)) .& (y_vals .> 0) .& (y_vals .< size(staining, 1)))
    staining_vals = zeros(T, length(x_vals))
    for i in inds
        staining_vals[i] = staining[y_vals[i], x_vals[i]];
    end

    return staining_vals
end

function encode_genes(gene_list)
    gene_names = unique(gene_list);
    gene_ids = Dict(zip(gene_names, 1:length(gene_names)))
    return [gene_ids[g] for g in gene_list], gene_names
end

# This function is a determenistic analogue of sampling. It picks points in a manner that preserves the distributions across x and y.
function select_ids_uniformly(vals::Union{Vector{<:Real}, <:AbstractMatrix{<:Real}}, confidence::Union{Vector{Float64}, Nothing}=nothing; n::Int, confidence_threshold::Float64=0.95)::Vector{Int}
    if n <= 1
        error("n must be > 1")
    end

    high_conf_ids = (confidence===nothing) ? (1:size(vals, 1)) : findall(confidence .>= confidence_threshold)
    if length(high_conf_ids) < n
        @warn "n=$n, which is > length(high_conf_ids) ($(length(high_conf_ids)))"
        return high_conf_ids
    end

    vals = sum(vals, dims=2)[high_conf_ids]
    return high_conf_ids[sortperm(vals)[unique(round.(Int, range(1, length(high_conf_ids), length=n)))]];
end