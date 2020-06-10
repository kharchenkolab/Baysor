using DataFrames
import MultivariateStats
using NearestNeighbors
using StatsBase: countmap, denserank
using Statistics

using Base.Threads

import CSV
import Distances

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

        return (min_molecules_per_cell > 10) ? max(div(n_genes, 10), min_molecules_per_cell) : max(min_molecules_per_cell, 3)
    end

    if param == :n_gene_pcs
        if n_genes === nothing
            error("Either `$param` or `n_genes` must be provided")
        end

        return min(max(div(n_genes, 3), 10), 100)
    end

    if param == :n_cells_init
        if n_molecules === nothing
            error("Either `$param` or `n_molecules` must be provided")
        end

        return div(n_molecules, min_molecules_per_cell) * 2
    end
end

parse_scale_std(scale_std::Float64, ::Real) = scale_std
parse_scale_std(scale_std::Nothing, scale::Real) = 0.25 * scale
function parse_scale_std(scale_std::String, scale::Real)
    scale_std = strip(scale_std)
    if scale_std[end] == '%'
        return parse(Float64, scale_std[1:end-1]) * scale / 100.0
    end

    return parse(Float64, scale_std)
end

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

function adjust_field_weights_by_dapi!(bm_data::BmmData, dapi::Matrix{Float64}; min_weight::Float64=0.01)
    if !(:dapi_brightness in names(bm_data.x))
        staining_vals = staining_value_per_transcript(bm_data.x, dapi)
        if maximum(staining_vals) < 1e-20
            @warn "Maximum value of dapi brightness per molecule is < 1e-20. Random field can't be adjusted with it."
            return
        end
        bm_data.x[!, :dapi_brightness] = max.(staining_vals, 0.01 * maximum(staining_vals));
    end

    for p1 in 1:length(bm_data.adjacent_points)
        x1, y1 = round(Int, bm_data.x.x[p1]), round(Int, bm_data.x.y[p1])
        for (i, p2) in enumerate(bm_data.adjacent_points[p1])
            bm_data.adjacent_weights[p1][i] *= bm_data.x.dapi_brightness[p2] / bm_data.x.dapi_brightness[p1];

            x2, y2 = round(Int, bm_data.x.x[p2]), round(Int, bm_data.x.y[p2])
            if (x1 < 1) || (x2 < 1) || (y1 < 1) || (y2 < 1) || (x1 > size(dapi, 2)) || (x2 > size(dapi, 2)) || (y1 > size(dapi, 1)) || (y2 > size(dapi, 1))
                continue
            end

            min_br_on_line = trace_values_along_line(dapi, x1, y1, x2, y2) |> minimum
            min_br_on_ends = min(bm_data.x.dapi_brightness[p2], bm_data.x.dapi_brightness[p1])
            bm_data.adjacent_weights[p1][i] *= max(min_br_on_line / min_br_on_ends, min_weight)
        end
    end
end

function encode_genes(gene_list)
    gene_names = unique(gene_list);
    gene_ids = Dict(zip(gene_names, 1:length(gene_names)))
    return [gene_ids[g] for g in gene_list], gene_names
end

function position_data_by_assignment(pos_data::T where T<: AbstractArray{Float64, 2}, assignment::Array{Int, 1})
    filt_ids = findall(assignment .> 0)
    return [pos_data[:, ids] for ids in split(filt_ids, assignment[filt_ids])]
end

function initial_params_from_assignment(spatial_df::DataFrame, assignment::Array{Int, 1})
    assignment = deepcopy(assignment)
    assignment[assignment .> 0] .= denserank(assignment[assignment .> 0])
    center_data = center_data_from_assignment(spatial_df, assignment, cov_mult=1.0)
    return InitialParams(copy(position_data(center_data.centers)'), CovMat.(center_data.center_covs), assignment)
end

# DEPRECATED?
center_data_from_assignment(spatial_df::DataFrame, assignment_col::Symbol; kwargs...)  =
    center_data_from_assignment(spatial_df, Array(spatial_df[!, assignment_col]); kwargs...)

# DEPRECATED?
function center_data_from_assignment(spatial_df::DataFrame, assignment::Array{Int, 1}; cov_mult::Float64=0.5)
    pos_data = position_data(spatial_df)
    cluster_centers = hcat([vec(mean(pos_data, dims=2)) for pos_data in position_data_by_assignment(pos_data, assignment)]...)
    covs = covs_from_assignment(pos_data, assignment)
    scale = estimate_scale_from_centers(mean.(eigvals.(covs)) .^ 0.5)
    return CenterData(DataFrame(cluster_centers', [:x, :y]), [cov_mult .* m for m in covs], scale...)
end

function covs_from_assignment(pos_data::T where T<: AbstractArray{Float64, 2}, assignment::Array{Int, 1}; min_size::Int=1)::Array{CovMat, 1}
    pos_data_by_assignment = position_data_by_assignment(pos_data, assignment)
    covs = sample_cov.(pos_data_by_assignment);
    mean_stds = MeanVec(vec(median(hcat(Vector.(eigvals.(covs[size.(pos_data_by_assignment, 2) .> min_size]))...), dims=2)))

    for i in findall(size.(pos_data_by_assignment, 2) .<= min_size)
        covs[i] = CovMat(diagm(0 => deepcopy(mean_stds)))
    end

    # return adjust_cov_matrix!.(covs)
    return covs
end

cell_centers_with_clustering(spatial_df::DataFrame, args...; kwargs...) =
    cell_centers_with_clustering(position_data(spatial_df), args...; kwargs...)

function cell_centers_with_clustering(pos_data::T where T<: AbstractArray{Float64, 2}, n_clusters::Int; scale::Union{Real, Nothing})
    n_clusters = min(n_clusters, size(pos_data, 2))

    cluster_centers = pos_data[:, select_ids_uniformly(pos_data[1,:], pos_data[2,:]; n=n_clusters)]
    cluster_labels = kshiftlabels(pos_data, cluster_centers);

    covs = (scale === nothing) ? covs_from_assignment(pos_data, cluster_labels) : Float64(scale) ^ 2
    return InitialParams(copy(cluster_centers'), covs, cluster_labels)
end

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
    neighb_cm = MultivariateStats.projection(MultivariateStats.fit(MultivariateStats.PCA, neighb_cm', maxoutdim=n_gene_pcs));

    return max.(map(x -> cor(x...), zip(eachrow.([neighb_cm[edge_list[1,:], :], neighb_cm[edge_list[2,:], :]])...)), min_similarity);
end

function initialize_bmm_data(df_spatial::DataFrame, args...; composition_neighborhood::Int, n_gene_pcs::Int, update_priors::Symbol=:no, real_edge_quant::Float64=0.3,
                             use_local_gene_similarities::Bool=true, adjacency_type::Symbol=:triangulation, prior_seg_confidence::Float64=0.5, kwargs...)::BmmData
    edge_list, adjacent_dists = adjacency_list(df_spatial, adjacency_type=adjacency_type)
    real_edge_length = quantile(adjacent_dists, real_edge_quant)
    min_edge_length = quantile(adjacent_dists, 0.1)

    adjacent_weights = nothing
    if use_local_gene_similarities
        gene_comp_sims = estimate_local_composition_similarities(df_spatial, edge_list; composition_neighborhood=composition_neighborhood, n_gene_pcs=n_gene_pcs)
        adjacent_weights = gene_comp_sims ./ max.(adjacent_dists, min_edge_length)
    else
        adjacent_weights = 1 ./ max.(adjacent_dists, min_edge_length)
    end

    adjacent_points, adjacent_weights = convert_edge_list_to_adj_list(edge_list, adjacent_weights; n_verts=size(df_spatial, 1))

    max_weight = quantile(vcat(adjacent_weights...), 0.975) # increases robustnes
    adjacent_weights = [min.(w, max_weight) for w in adjacent_weights]

    components, sampler, assignment = initial_distributions(df_spatial, args...; kwargs...)
    return BmmData(components, df_spatial, adjacent_points, adjacent_weights, 1 ./ real_edge_length, sampler, assignment, update_priors=update_priors, prior_seg_confidence=prior_seg_confidence)
end

"""
    main function for initialization of bm_data
"""
function initial_distribution_arr(df_spatial::DataFrame, args...; n_frames::Int, n_frames_return::Int=0, n_cells_init::Union{Int, Nothing}=nothing,
                                  confidence_nn_id::Union{Int, Nothing}=nothing, min_molecules_per_cell::Union{Int, Nothing}=nothing, composition_neighborhood::Union{Int, Nothing}=nothing,
                                  n_gene_pcs::Union{Int, Nothing}=nothing, prior_seg_confidence::Float64=0.5, kwargs...)::Array{BmmData, 1}
    df_spatial = deepcopy(df_spatial)

    confidence_nn_id = something(confidence_nn_id, default_param_value(:confidence_nn_id, min_molecules_per_cell))

    if confidence_nn_id > 0
        @info "Estimate confidence per molecule"
        prior_segmentation = (:prior_segmentation in names(df_spatial)) ? df_spatial.prior_segmentation : nothing
        append_confidence!(df_spatial, prior_segmentation; nn_id=confidence_nn_id, prior_confidence=prior_seg_confidence)
        @info "Done"
    end

    dfs_spatial = n_frames > 1 ? split_spatial_data(df_spatial, n_frames) : [df_spatial]
    @info "#frames: $(length(dfs_spatial)); mean number of molecules per frame: $(median(size.(dfs_spatial, 1)))."

    if (n_frames_return > 0) && (n_frames_return < n_frames)
        dfs_spatial = dfs_spatial[1:n_frames_return]
    end

    n_cells_init = something(n_cells_init, default_param_value(:n_cells_init, min_molecules_per_cell, n_molecules=size(df_spatial, 1)))
    composition_neighborhood = something(composition_neighborhood, default_param_value(:composition_neighborhood, min_molecules_per_cell, n_genes=maximum(df_spatial.gene)))
    n_gene_pcs = something(n_gene_pcs, default_param_value(:n_gene_pcs, min_molecules_per_cell, n_genes=maximum(df_spatial.gene)))

    return initial_distribution_arr(dfs_spatial, args...; n_cells_init=n_cells_init, min_molecules_per_cell=min_molecules_per_cell,
        composition_neighborhood=composition_neighborhood, n_gene_pcs=n_gene_pcs, prior_seg_confidence=prior_seg_confidence, kwargs...)
end

# DEPRECATED?
function initial_distribution_arr(dfs_spatial::Array{DataFrame, 1}, centers::CenterData; n_cells_init::Int,
                                  scale::Union{Number, Nothing}=nothing, scale_std::Union{Float64, String, Nothing}=nothing,
                                  new_component_weight::Number=0.2, center_component_weight::Number=1.0, n_degrees_of_freedom_center::Union{Int, Nothing}=nothing,
                                  update_priors::Symbol=:no, min_molecules_per_cell::Union{Int, Nothing}=nothing, kwargs...)::Array{BmmData, 1}
    scale = something(scale, centers.scale_estimate)
    scale_std = parse_scale_std(something(scale_std, centers.scale_std_estimate), scale)

    @info "Initializing algorithm. Scale: $scale, scale std: $scale_std, initial #clusters: $n_cells_init."

    n_cells_init = max(div(n_cells_init, length(dfs_spatial)), 1)

    n_degrees_of_freedom_center = something(n_degrees_of_freedom_center, default_param_value(:n_degrees_of_freedom_center, min_molecules_per_cell))
    centers.center_covs = something(centers.center_covs, [diagm(0 => [scale / 2, scale / 2] .^ 2) for i in 1:size(centers.centers, 1)])

    centers_per_frame = subset_by_coords.(Ref(centers), dfs_spatial);

    size_prior = ShapePrior(Float64[scale, scale], Float64[scale_std, scale_std], min_molecules_per_cell);

    if any([size(c.centers, 1) == 0 for c in centers_per_frame])
        if (n_cells_init == 0)
            error("Some frames don't contain cell centers. Try to reduce number of frames or provide better segmentation or increase n-cells-init.")
        else
            @warn "Some frames don't contain cell centers. Possibly, coordinates of DAPI and transcripts are not aligned or DAPI should be transposed."
        end
    end

    bm_datas_res = Array{BmmData, 1}(undef, length(dfs_spatial))
    @threads for i in 1:length(dfs_spatial)
        bm_datas_res[i] = initialize_bmm_data(dfs_spatial[i], centers_per_frame[i]; size_prior=size_prior, new_component_weight=new_component_weight,
            prior_component_weight=center_component_weight, default_std=scale, n_degrees_of_freedom_center=n_degrees_of_freedom_center,
            update_priors=update_priors, n_cells_init=n_cells_init, kwargs...)
    end
    return bm_datas_res;
end

function initial_distribution_arr(dfs_spatial::Array{DataFrame, 1}; scale::Number, scale_std::Union{Float64, String, Nothing}=nothing, n_cells_init::Int,
                                  new_component_weight::Number=0.2, min_molecules_per_cell::Union{Int, Nothing}=nothing, kwargs...)::Array{BmmData, 1}
    scale_std = parse_scale_std(scale_std, scale)

    @info "Initializing algorithm. Scale: $scale, scale std: $scale_std, initial #clusters: $n_cells_init."

    size_prior = ShapePrior(Float64[scale, scale], Float64[scale_std, scale_std], min_molecules_per_cell);
    n_cells_init = max(div(n_cells_init, length(dfs_spatial)), 1)

    bm_datas_res = Array{BmmData, 1}(undef, length(dfs_spatial))
    for i in 1:length(dfs_spatial)
        init_params = cell_centers_with_clustering(dfs_spatial[i], n_cells_init; scale=scale)
        bm_datas_res[i] = initialize_bmm_data(dfs_spatial[i], init_params; size_prior=size_prior, new_component_weight=new_component_weight, kwargs...)
    end
    return bm_datas_res;
end

# DEPRECATED?
function initial_distribution_data(df_spatial::DataFrame, prior_centers::CenterData; n_degrees_of_freedom_center::Int, default_std::Union{Real, Nothing}=nothing,
           gene_num::Int=maximum(df_spatial[!,:gene]), noise_confidence_threshold::Float64=0.1, n_cells_init::Int=0)
    mtx_centers = position_data(prior_centers.centers)
    pos_data = position_data(df_spatial)
    can_be_dropped = falses(size(mtx_centers, 2))
    n_added_clusters = 0

    if n_cells_init > size(mtx_centers, 2)
        # TODO: use select_ids_uniformly instead if decide to keep this function
        cluster_centers = kshiftmedoids(pos_data, n_cells_init)[1];

        if size(mtx_centers, 2) > 0
            dists_to_prior = getindex.(knn(KDTree(mtx_centers), cluster_centers, 1)[2], 1);
            cluster_centers = cluster_centers[:, dists_to_prior .> sort(dists_to_prior)[size(mtx_centers, 2)]];
        end

        n_added_clusters = size(cluster_centers, 2)
        mtx_centers = hcat(mtx_centers, cluster_centers);
        can_be_dropped = vcat(can_be_dropped, trues(n_added_clusters));
    end

    assignment = kshiftlabels(pos_data, mtx_centers);
    if :confidence in names(df_spatial)
        assignment[df_spatial[!,:confidence] .< noise_confidence_threshold] .= 0
    end

    covs = (default_std === nothing) ?
        covs_from_assignment(pos_data, assignment) :
        [Float64[default_std 0; 0 default_std] .^ 2 for i in 1:size(mtx_centers, 2)]

    pos_distributions = [MvNormalF(vec(mtx_centers[:, i]), cov) for (i, cov) in enumerate(covs)]

    center_priors = [CellCenter(pd.μ, deepcopy(cov), n_degrees_of_freedom_center)
        for (pd,cov) in zip(pos_distributions[1:length(prior_centers.center_covs)], prior_centers.center_covs)]
    center_priors = vcat(center_priors, [nothing for x in 1:n_added_clusters])

    return assignment, pos_distributions, center_priors, can_be_dropped
end

# DEPRECATED?
"""
    Creates `distributions with `prior_centers` centers and `default_std` standard deviation

    # Arguments
    - `df_spatial::DataFrame`: DataFrame with columns `:x`, `:y` and `:gene`
    - `prior_centers::CenterData`: info on centers, extracted from DAPIs. Must have `center_covs` filled.
    - `size_prior::Union{ShapePrior, Nothing}`: shape prior for position_params
    - `new_component_weight::Float64`:
    - `prior_component_weight::Float64`:
    - `n_degrees_of_freedom_center::Int`:
    - `default_std::Union{Real, Nothing}=nothing`: initial std for position_params in components. Estimated from assignment if `nothing` is passed.
    - `gene_num::Int=maximum(df_spatial[:gene])`: total number of genes in the dataset
    - `noise_confidence_threshold::Float64=0.1`: all molecules with confidence < `noise_confidence_threshold` are assigned to noise class
    - `n_cells_init::Int=0`: number of clusters for initialization. Ignored if less than number of centers. Otherwise corresponding clusters are added with KShift clustering
"""
function initial_distributions(df_spatial::DataFrame, prior_centers::CenterData; size_prior::Union{ShapePrior, Nothing}, new_component_weight::Float64,
                               prior_component_weight::Float64, n_degrees_of_freedom_center::Int, default_std::Union{Real, Nothing}=nothing,
                               gene_num::Int=maximum(df_spatial[!,:gene]), noise_confidence_threshold::Float64=0.1,
                               n_cells_init::Int=0)
    if new_component_weight < 1e-10
        n_cells_init = 0
    end

    assignment, pos_distributions, center_priors, can_be_dropped = initial_distribution_data(df_spatial, prior_centers;
        n_degrees_of_freedom_center=n_degrees_of_freedom_center, default_std=default_std, gene_num=gene_num,
        noise_confidence_threshold=noise_confidence_threshold, n_cells_init=n_cells_init)

    gene_prior = CategoricalSmoothed(ones(Int, gene_num));
    n_mols_per_center = count_array(assignment .+ 1, max_value=length(pos_distributions) + 1)[2:end];
    components = [Component(pd, deepcopy(gene_prior), shape_prior=deepcopy(size_prior), center_prior=cp,
                            n_samples=n, prior_weight=prior_component_weight, can_be_dropped=dr)
                    for (pd, cp, n, dr) in zip(pos_distributions, center_priors, n_mols_per_center, can_be_dropped)];

    ids_per_comp = split(collect(1:length(assignment)), assignment .+ 1)[2:end]
    for (ids, comp) in zip(ids_per_comp, components)
        if comp.n_samples > 0
            maximize!(comp.position_params, position_data(df_spatial[ids,:]))
            maximize!(comp.composition_params, composition_data(df_spatial[ids,:]))
        end
    end

    gene_sampler = CategoricalSmoothed(ones(Int, gene_num))
    sampler = Component(MvNormalF(zeros(2)), gene_sampler, shape_prior=deepcopy(size_prior), prior_weight=new_component_weight, can_be_dropped=false) # position_params are never used

    return components, sampler, assignment
end

function initial_distributions(df_spatial::DataFrame, initial_params::InitialParams; size_prior::ShapePrior, new_component_weight::Float64,
                               gene_smooth::Real=1.0, gene_num::Int=maximum(df_spatial[!,:gene]))
    position_distrubutions = [MvNormalF(initial_params.centers[i,:], initial_params.covs[i]) for i in 1:size(initial_params.centers, 1)]
    gene_distributions = [CategoricalSmoothed(ones(Int, gene_num), smooth=Float64(gene_smooth)) for i in 1:length(position_distrubutions)]

    components = [Component(pd, gd, shape_prior=deepcopy(size_prior), prior_weight=new_component_weight, can_be_dropped=true)
                    for (pd, gd) in zip(position_distrubutions, gene_distributions)]

    gene_sampler = CategoricalSmoothed(ones(Int, gene_num))
    sampler = Component(MvNormalF(zeros(2)), gene_sampler, shape_prior=deepcopy(size_prior), prior_weight=new_component_weight, can_be_dropped=false) # position_params are never used

    return components, sampler, initial_params.assignment
end

# DEPRECATED?
function filter_small_components(c_components::Array{Array{Int, 1}, 1}, adjacent_points::Array{Array{Int, 1}, 1}, df_spatial::DataFrame;
                                 min_molecules_per_cell::Int=10)
    c_components = c_components[length.(c_components) .> min_molecules_per_cell];

    presented_ids = sort(vcat(c_components...))
    map_dict = Dict(presented_ids .=> 1:length(presented_ids));

    c_components = [get.(Ref(map_dict), comp, 0) for comp in c_components];
    adjacent_points = [get.(Ref(map_dict), ids, 0) for ids in adjacent_points[presented_ids]];
    adjacent_points = [v[v .!= 0] for v in adjacent_points];

    return c_components, adjacent_points, df_spatial[presented_ids, :]
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
    # return filter_background(df_spatial), gene_names;
end

function split_spatial_data(df::DataFrame, n::Int, key::Symbol)::Array{DataFrame, 1}
    factor = vec(sum(hcat([df[!,key] .<= quantile(df[!,key], q) for q in range(1 / n, stop=1.0, length=n)]...), dims=2))
    return split(df, factor)
end

split_spatial_data(df::DataFrame, n_hor::Int, n_ver::Int) = vcat(split_spatial_data.(split_spatial_data(df, n_ver, :y), n_hor, :x)...)

function split_spatial_data(df::DataFrame, n::Int) # TODO: very approximate separation. Example: n=3.
    df_sizes = Dict(s => maximum(df[!,s]) - minimum(df[!,s]) for s in [:x, :y]);
    x_elongation = df_sizes[:x] / sum(values(df_sizes))
    a = round(Int, sqrt(x_elongation * n / (1 - x_elongation))) # solution of "a * b = n; a / (a + b) = x_elongation"

    return split_spatial_data(df, a, round(Int, n / a))
end

split_spatial_data(df::DataFrame; mean_mols_per_frame::Int) = split_spatial_data(df, round(Int, size(df, 1) / mean_mols_per_frame))
