import MultivariateStats
using NearestNeighbors
using StatsBase: countmap, denserank
using Statistics

using Base.Threads

import Distances

# Parse parameters

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

position_data_by_assignment(pos_data::Matrix{Float64}, assignment::Vector{Int}) =
    [pos_data[:, ids] for ids in split_ids(assignment, drop_zero=true)]

function covs_from_assignment(pos_data::Matrix{Float64}, assignment::Vector{Int}; min_size::Int=1)::Array{<:CovMat, 1}
    pos_data_by_assignment = position_data_by_assignment(pos_data, assignment)
    CM, MV = (size(pos_data, 1) == 2) ? (CovMat{2}, MeanVec{2}) : (CovMat{3}, MeanVec{3})

    covs = [estimate_sample_cov!(zeros(CM), pd; μ=MV(mean(pd, dims=2)[:])) for pd in pos_data_by_assignment];
    mean_covs = MV(vec(median(hcat(Vector.(eigvals.(covs[size.(pos_data_by_assignment, 2) .>= min_size]))...), dims=2)))

    for i in findall(size.(pos_data_by_assignment, 2) .<= min_size)
        covs[i] = CM(diagm(0 => deepcopy(mean_covs)))
    end

    return covs
end

cell_centers_from_assignment(pos_data::Matrix{Float64}, assignment::Vector{Int}; min_size::Int=1)::Array{<:MeanVec, 1} =
    [MeanVec{length(x)}(vec(x)) for x in mean.(position_data_by_assignment(pos_data, assignment), dims=2)]

cell_centers_uniformly(spatial_df::DataFrame, args...; kwargs...) =
    cell_centers_uniformly(position_data(spatial_df), args..., (:confidence in propertynames(spatial_df)) ? spatial_df.confidence : nothing; kwargs...)

function cell_centers_uniformly(
        pos_data::Matrix{Float64}, n_clusters::Int,
        confidences::Union{Vector{Float64}, Nothing}=nothing; scale::Union{Float64, Nothing}
    )
    n_clusters = min(n_clusters, size(pos_data, 2))

    cluster_centers = pos_data[:, select_ids_uniformly(pos_data', confidences; n=n_clusters, confidence_threshold=0.25)]
    cluster_labels = vcat(knn(KDTree(cluster_centers), pos_data, 1)[1]...)

    covs = nothing
    if scale === nothing
        covs = covs_from_assignment(pos_data, cluster_labels)
    elseif size(pos_data, 1) == 2
        scale = scale ^ 2
        covs = [CovMat{2, 4}(diagm(0 => (ones(2) .* scale))) for _ in 1:size(cluster_centers, 2)]
    elseif size(pos_data, 1) == 3
        scale = scale ^ 2
        covs = [CovMat{3, 9}(diagm(0 => (ones(3) .* scale))) for _ in 1:size(cluster_centers, 2)]
    else
        error("Unexpected number of dimensions: $(size(pos_data, 1))")
    end

    return InitialParams(copy(cluster_centers'), covs, cluster_labels)
end

# Build MRF

function convert_edge_list_to_adj_list(edge_list::Matrix{Int}, edge_weights::Union{Vector{Float64}, Nothing}=nothing; n_verts::Int=maximum(edge_list))
    res_ids = map(vcat,
        split(edge_list[2,:], edge_list[1,:], max_factor=n_verts),
        split(edge_list[1,:], edge_list[2,:], max_factor=n_verts)
    )

    if edge_weights === nothing
        return res_ids
    end

    res_weights = map(vcat,
        split(edge_weights, edge_list[1,:], max_factor=n_verts),
        split(edge_weights, edge_list[2,:], max_factor=n_verts)
    );

    return res_ids, res_weights
end

function estimate_local_composition_similarities(df_spatial::DataFrame, edge_list::Matrix{Int}; composition_neighborhood::Int, n_gene_pcs::Int, min_similarity::Float64=0.01)
    neighb_cm = neighborhood_count_matrix(df_spatial, composition_neighborhood);
    neighb_cm = MultivariateStats.transform(MultivariateStats.fit(MultivariateStats.PCA, neighb_cm, maxoutdim=n_gene_pcs), neighb_cm);

    return max.(map(i -> cor(neighb_cm[:, edge_list[1,i]], neighb_cm[:, edge_list[2,i]]), 1:size(edge_list, 2)), min_similarity);
end

function build_molecule_graph(
        df_spatial::DataFrame; min_edge_quant::Float64=0.3, use_local_gene_similarities::Bool=false,
        n_gene_pcs::Int=0, composition_neighborhood::Int=0, kwargs...
    )
    if use_local_gene_similarities && ((composition_neighborhood == 0) || (n_gene_pcs == 0))
        @warn "composition_neighborhood=$composition_neighborhood and n_gene_pcs=$n_gene_pcs, while use_local_gene_similarities=true. Force use_local_gene_similarities=false."
        use_local_gene_similarities = false
    end

    if use_local_gene_similarities && any(ismissing.(composition_data(df_spatial)))
        @warn "use_local_gene_similarities is not supported with missing gene values. Setting to false."
        use_local_gene_similarities = false
    end

    edge_list, adjacent_dists = adjacency_list(df_spatial; kwargs...);

    min_edge_length = quantile(adjacent_dists, min_edge_quant);

    adjacent_weights = min_edge_length ./ max.(adjacent_dists, min_edge_length)
    if use_local_gene_similarities
        adjacent_weights .*= estimate_local_composition_similarities(df_spatial, edge_list; composition_neighborhood=composition_neighborhood, n_gene_pcs=n_gene_pcs)
    end

    adjacent_points, adjacent_weights = convert_edge_list_to_adj_list(edge_list, adjacent_weights; n_verts=size(df_spatial, 1));
    return AdjList(adjacent_points, adjacent_weights)
end

# Initialize BmmData

function initialize_shape_prior(scale::Float64, scale_std::Float64; min_molecules_per_cell::Int, is_3d::Bool)
    is_3d || return ShapePrior{2}(Float64[scale, scale], Float64[scale_std, scale_std], min_molecules_per_cell)
    return ShapePrior{3}(Float64[scale, scale, scale], Float64[scale_std, scale_std, scale_std], min_molecules_per_cell)
end

"""
    main function for initialization of bm_data
"""
function initialize_bmm_data(
        df_spatial::DataFrame; min_molecules_per_cell::Int, adj_list::AdjList,
        n_cells_init::Int, scale::T where T<: Real, scale_std::Union{<:Real, String, Nothing}=nothing,
        prior_seg_confidence::Float64=0.5, verbose::Bool=true, kwargs...
    )::BmmData
    df_spatial = deepcopy(df_spatial)
    scale_std = parse_scale_std(scale_std, scale)

    ## Initialize BmmData array
    verbose && @info "Initializing algorithm. Scale: $scale, scale std: $scale_std, initial #components: $n_cells_init, #molecules: $(size(df_spatial, 1))."
    size_prior = initialize_shape_prior(Float64(scale), scale_std; min_molecules_per_cell, is_3d=(:z in propertynames(df_spatial)))
    init_params = cell_centers_uniformly(df_spatial, n_cells_init; scale=scale)

    components, assignment = initial_distributions(df_spatial, init_params; size_prior=size_prior)

    std_vals = size_prior.std_values
    noise_position_density = pdf(MultivariateNormal(zeros(length(std_vals)), diagm(0 => std_vals.^2)), 3 .* std_vals)

    return BmmData(
        components, df_spatial, adj_list, assignment;
        prior_seg_confidence, noise_position_density, kwargs...
    )
end

function initial_distributions(df_spatial::DataFrame, initial_params::InitialParams; size_prior::ShapePrior{N},
                               gene_smooth::Real=1.0, gene_num::Int=maximum(skipmissing(composition_data(df_spatial)))) where N
    position_distrubutions = [MvNormalF(initial_params.centers[i,:], initial_params.covs[i]) for i in 1:size(initial_params.centers, 1)]
    gene_distributions = [CategoricalSmoothed(ones(Float64, gene_num), smooth=Float64(gene_smooth)) for i in 1:length(position_distrubutions)]

    components = [Component(pd, gd, shape_prior=deepcopy(size_prior)) for (pd, gd) in zip(position_distrubutions, gene_distributions)]

    gene_sampler = CategoricalSmoothed(ones(Float64, gene_num))

    return components, initial_params.assignment
end

function initialize_position_params_from_assignment(pos_data::Matrix{Float64}, cell_assignment::Vector{Int})
    centers = cell_centers_from_assignment(pos_data, cell_assignment);
    covs = covs_from_assignment(pos_data, cell_assignment);
    adjust_cov_matrix!.(covs);

    return MvNormalF.(centers, covs);
end

# Utils

# This function is a determenistic analogue of sampling. It picks points in a manner that preserves the distributions across x and y.
function select_ids_uniformly(
        vals::Union{Vector{<:Real}, <:AbstractMatrix{<:Real}},
        confidence::Union{Vector{Float64}, Nothing}=nothing;
        n::Int, confidence_threshold::Float64=0.95
    )::Vector{Int}
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