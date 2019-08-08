using DataFrames
using NearestNeighbors
using StatsBase: countmap, denserank
using Statistics

import CSV
import Distances

function default_param_value(param::Symbol, min_molecules_per_cell::Union{Int, Nothing}; n_molecules::Union{Int, Nothing}=nothing)
    if min_molecules_per_cell === nothing
        error("Either `$param` or `min_molecules_per_cell` must be provided")
    end

    min_molecules_per_cell = max(min_molecules_per_cell, 3)

    if param == :n_degrees_of_freedom_center
        return min_molecules_per_cell
    end

    if param == :confidence_nn_id
        return max(div(min_molecules_per_cell, 2) + 1, 3)
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

function append_confidence!(df_spatial::DataFrame; nn_id::Int, border_quantiles::Tuple{Float64, Float64}=(0.3, 0.975))
    pos_data = position_data(df_spatial);
    mean_dists = getindex.(knn(KDTree(pos_data), pos_data, nn_id + 1, true)[2], nn_id + 1)

    df_spatial[!,:confidence] = interpolate_linear.(mean_dists, quantile(mean_dists, border_quantiles)...);
end

function encode_genes(gene_list)
    gene_names = unique(gene_list);
    gene_ids = Dict(zip(gene_names, 1:length(gene_names)))
    return [gene_ids[g] for g in gene_list], gene_names
end

assign_cells_to_centers(spatial_df::DataFrame, centers::DataFrame)::Array{Int, 1} =
    [v[1] for v in knn(KDTree(position_data(centers)), position_data(spatial_df), 1)[1]]


function position_data_by_assignment(spatial_df::DataFrame, assignment::Array{Int, 1})
    filt_mask = (assignment .> 0)
    spatial_df, assignment = spatial_df[filt_mask, :], assignment[filt_mask]

    pos_data = position_data(spatial_df);
    return [pos_data[:, ids] for ids in split(1:length(assignment), assignment)]
end

center_data_from_assignment(spatial_df::DataFrame, assignment_col::Symbol; kwargs...)  =
    center_data_from_assignment(spatial_df, Array(spatial_df[!, assignment_col]); kwargs...)

function center_data_from_assignment(spatial_df::DataFrame, assignment::Array{Int, 1}; cov_mult::Float64=0.5)
    cluster_centers = hcat([vec(mean(pos_data, dims=2)) for pos_data in position_data_by_assignment(spatial_df, assignment)]...)
    covs = covs_from_assignment(spatial_df, assignment)
    return CenterData(DataFrame(cluster_centers', [:x, :y]), [cov_mult .* m for m in covs], estimate_scale_from_centers(cluster_centers)...)
end

function covs_from_assignment(spatial_df::DataFrame, assignment::Array{Int, 1})
    pos_data_by_assignment = position_data_by_assignment(spatial_df, assignment)
    stds = [vec(std(pos_data, dims=2)) for pos_data in pos_data_by_assignment]
    mean_stds = vec(median(hcat(stds[size.(pos_data_by_assignment, 2) .> 1]...), dims=2))

    for i in findall(size.(pos_data_by_assignment, 2) .<= 1)
        stds[i] = deepcopy(mean_stds)
    end

    return [diagm(0 => s .^ 2) for s in stds]
end

function cell_centers_with_clustering(spatial_df::DataFrame, n_clusters::Int; scale::Union{Real, Nothing})
    n_clusters = min(n_clusters, size(spatial_df, 1))

    pos_data = position_data(spatial_df);
    cluster_centers = kshiftmedoids(pos_data, n_clusters)[1];
    cluster_labels = kshiftlabels(pos_data, cluster_centers);

    covs = (scale === nothing) ? covs_from_assignment(spatial_df, cluster_labels) : scale ^ 2
    return InitialParams(copy(cluster_centers'), covs, cluster_labels)
end

"""
    main function for initialization of bm_data
"""
function initial_distribution_arr(df_spatial::DataFrame, args...; n_frames::Int, n_frames_return::Int=0, n_cells_init::Union{Int, Nothing}=nothing,
                                  confidence_nn_id::Union{Int, Nothing}=nothing, confidence_border_quantiles::Tuple{Float64, Float64}=(0.3, 0.975), 
                                  min_molecules_per_cell::Union{Int, Nothing}=nothing, kwargs...)::Array{BmmData, 1}
    df_spatial = deepcopy(df_spatial)

    confidence_nn_id = something(confidence_nn_id, default_param_value(:confidence_nn_id, min_molecules_per_cell))

    if !(:confidence in names(df_spatial)) && (confidence_nn_id > 0)
        @info "Estimate confidence per molecule"
        append_confidence!(df_spatial; nn_id=confidence_nn_id, border_quantiles=confidence_border_quantiles)
        @info "Done"
    end

    dfs_spatial = n_frames > 1 ? split_spatial_data(df_spatial, n_frames) : [df_spatial]
    @info "#frames: $(length(dfs_spatial)); mean number of molecules per frame: $(median(size.(dfs_spatial, 1)))."

    if (n_frames_return > 0) && (n_frames_return < n_frames)
        dfs_spatial = dfs_spatial[1:n_frames_return]
    end

    n_cells_init = something(n_cells_init, default_param_value(:n_cells_init, min_molecules_per_cell, n_molecules=size(df_spatial, 1)))

    return initial_distribution_arr(dfs_spatial, args...; n_cells_init=n_cells_init, min_molecules_per_cell=min_molecules_per_cell, kwargs...)
end

function initial_distribution_arr(dfs_spatial::Array{DataFrame, 1}, centers::CenterData; n_cells_init::Int, 
                                  scale::Union{Number, Nothing}=nothing, scale_std::Union{Float64, String, Nothing}=nothing, 
                                  new_component_weight::Number=0.2, center_component_weight::Number=1.0, n_degrees_of_freedom_center::Union{Int, Nothing}=nothing,
                                  update_priors::Symbol=:no, min_molecules_per_cell::Union{Int, Nothing}=nothing, kwargs...)::Array{BmmData, 1}
    scale = something(scale, centers.scale_estimate)
    scale_std = parse_scale_std(something(scale_std, centers.scale_std_estimate), scale)

    n_degrees_of_freedom_center = something(n_degrees_of_freedom_center, default_param_value(:n_celln_degrees_of_freedom_centers_init, min_molecules_per_cell))
    centers.center_covs = something(centers.center_covs, [diagm(0 => [scale / 2, scale / 2] .^ 2) for i in 1:size(centers.centers, 1)])

    centers_per_frame = subset_by_coords.(Ref(centers), dfs_spatial);

    size_prior = ShapePrior(Float64[scale, scale], Float64[scale_std, scale_std], min_molecules_per_cell);

    if any([size(c.centers, 1) == 0 for c in centers_per_frame])
        error("Some frames don't contain cell centers. Try to reduce number of frames or provide better segmentation.")
    end

    @info "Initializing algorithm. Scale: $scale, scale std: $scale_std, initial #clusters: $n_cells_init."

    n_cells_init = max(div(n_cells_init, length(dfs_spatial)), 2)
    return initialize_bmm_data.(dfs_spatial, centers_per_frame; size_prior=size_prior, new_component_weight=new_component_weight,
            prior_component_weight=center_component_weight, default_std=scale, n_degrees_of_freedom_center=n_degrees_of_freedom_center,
            update_priors=update_priors, n_cells_init=n_cells_init, kwargs...);
end

function initial_distribution_arr(dfs_spatial::Array{DataFrame, 1}; scale::Number, scale_std::Union{Float64, String, Nothing}=nothing, n_cells_init::Int, 
                                  new_component_weight::Number=0.2, min_molecules_per_cell::Union{Int, Nothing}=nothing, kwargs...)::Array{BmmData, 1}
    scale_std = parse_scale_std(scale_std, scale)

    @info "Initializing algorithm. Scale: $scale, scale std: $scale_std, initial #clusters: $n_cells_init."
    size_prior = ShapePrior(Float64[scale, scale], Float64[scale_std, scale_std], min_molecules_per_cell);
    initial_params_per_frame = cell_centers_with_clustering.(dfs_spatial, n_cells_init; scale=scale)
    return initialize_bmm_data.(dfs_spatial, initial_params_per_frame; size_prior=size_prior, new_component_weight=new_component_weight, kwargs...)
end

function initialize_bmm_data(df_spatial::DataFrame, args...; update_priors::Symbol=:no, kwargs...)::BmmData
    adjacent_points, adjacent_dists = adjacency_list(df_spatial)
    real_edge_length = quantile(vcat(adjacent_dists...), 0.3)

    # Distance from point to itself equal to real_edge_length
    for i in 1:length(adjacent_points)
        for j in 1:length(adjacent_points[i])
            if adjacent_points[i][j] == i
                adjacent_dists[i][j] = real_edge_length
            end
        end
    end

    adjacent_weights = [1 ./ d for d in adjacent_dists]
    max_weight = quantile(vcat(adjacent_weights...), 0.975)
    adjacent_weights = [min.(w, max_weight) for w in adjacent_weights]

    components, sampler, assignment = initial_distributions(df_spatial, args...; kwargs...)
    return BmmData(components, df_spatial, adjacent_points, adjacent_weights, 1 ./ real_edge_length, sampler, assignment, update_priors=update_priors)
end

function initial_distribution_data(df_spatial::DataFrame, prior_centers::CenterData; n_degrees_of_freedom_center::Int, default_std::Union{Real, Nothing}=nothing,
           gene_num::Int=maximum(df_spatial[!,:gene]), noise_confidence_threshold::Float64=0.1, n_cells_init::Int=0)

    mtx_centers = position_data(prior_centers.centers)
    pos_data = position_data(df_spatial)
    can_be_dropped = falses(size(mtx_centers, 2))
    n_added_clusters = 0

    if n_cells_init > size(mtx_centers, 2)
        cluster_centers = kshiftmedoids(pos_data, n_cells_init)[1];
        dists_to_prior = getindex.(knn(KDTree(mtx_centers), cluster_centers, 1)[2], 1);

        additional_centers = cluster_centers[:, dists_to_prior .> sort(dists_to_prior)[size(mtx_centers, 2)]];
        n_added_clusters = size(additional_centers, 2)

        mtx_centers = hcat(mtx_centers, additional_centers);
        can_be_dropped = vcat(can_be_dropped, trues(n_added_clusters));
    end

    assignment = kshiftlabels(pos_data, mtx_centers);
    if :confidence in names(df_spatial)
        assignment[df_spatial[!,:confidence] .< noise_confidence_threshold] .= 0
    end

    covs = (default_std === nothing) ?
        covs_from_assignment(df_spatial, assignment) :
        [Float64[default_std 0; 0 default_std] .^ 2 for i in 1:size(mtx_centers, 2)]

    pos_distributions = [MvNormal(vec(mtx_centers[:, i]), cov) for (i, cov) in enumerate(covs)]

    center_priors = [CellCenter(pd.μ, deepcopy(cov), n_degrees_of_freedom_center)
        for (pd,cov) in zip(pos_distributions[1:length(prior_centers.center_covs)], prior_centers.center_covs)]
    center_priors = vcat(center_priors, [nothing for x in 1:n_added_clusters])

    return assignment, pos_distributions, center_priors, can_be_dropped
end

"""
    Creates `distributions with `prior_centers` centers and `default_std` standard deviation

    # Arguments
    - `df_spatial::DataFrame`: DataFrame with columns `:x`, `:y` and `:gene`
    - `prior_centers::CenterData`: info on centers, extracted from DAPIs. Must have `center_covs` filled.
    - `size_prior::ShapePrior`: shape prior for position_params
    - `new_component_weight::Float64`:
    - `prior_component_weight::Float64`:
    - `n_degrees_of_freedom_center::Int`:
    - `default_std::Union{Real, Nothing}=nothing`: initial std for position_params in components. Estimated from assignment if `nothing` is passed.
    - `gene_num::Int=maximum(df_spatial[:gene])`: total number of genes in the dataset
    - `noise_confidence_threshold::Float64=0.1`: all molecules with confidence < `noise_confidence_threshold` are assigned to noise class
    - `n_cells_init::Int=0`: number of clusters for initialization. Ignored if less than number of centers. Otherwise corresponding clusters are added with KShift clustering
"""
function initial_distributions(df_spatial::DataFrame, prior_centers::CenterData; size_prior::ShapePrior, new_component_weight::Float64,
                               prior_component_weight::Float64, n_degrees_of_freedom_center::Int, default_std::Union{Real, Nothing}=nothing,
                               gene_num::Int=maximum(df_spatial[!,:gene]), noise_confidence_threshold::Float64=0.1,
                               n_cells_init::Int=0)
    if new_component_weight < 1e-10
        n_cells_init = 0
    end

    assignment, pos_distributions, center_priors, can_be_dropped = initial_distribution_data(df_spatial, prior_centers;
        n_degrees_of_freedom_center=n_degrees_of_freedom_center, default_std=default_std, gene_num=gene_num,
        noise_confidence_threshold=noise_confidence_threshold, n_cells_init=n_cells_init)

    gene_prior = SingleTrialMultinomial(ones(Int, gene_num));
    n_mols_per_center = count_array(assignment .+ 1, max_value=length(pos_distributions) + 1)[2:end];
    components = [Component(pd, deepcopy(gene_prior), shape_prior=deepcopy(size_prior), center_prior=cp,
                            n_samples=n, prior_weight=prior_component_weight, can_be_dropped=dr)
                    for (pd, cp, n, dr) in zip(pos_distributions, center_priors, n_mols_per_center, can_be_dropped)];

    ids_per_comp = split(collect(1:length(assignment)), assignment .+ 1)[2:end]
    for (ids, comp) in zip(ids_per_comp, components)
        if comp.n_samples > 0
            comp.position_params = maximize(comp.position_params, position_data(df_spatial[ids,:]))
            comp.composition_params = maximize(comp.composition_params, composition_data(df_spatial[ids,:]))
        end
    end

    gene_sampler = SingleTrialMultinomial(ones(Int, gene_num))
    sampler = Component(MvNormal(zeros(2)), gene_sampler, shape_prior=deepcopy(size_prior), prior_weight=new_component_weight, can_be_dropped=false) # position_params are never used

    return components, sampler, assignment
end

function initial_distributions(df_spatial::DataFrame, initial_params::InitialParams; size_prior::ShapePrior, new_component_weight::Float64,
                               gene_smooth::Real=1.0, gene_num::Int=maximum(df_spatial[!,:gene]))
    gene_distributions = [SingleTrialMultinomial(ones(Int, gene_num), smooth=Float64(gene_smooth)) for i in 1:initial_params.n_comps]

    position_distrubutions = [MvNormal(initial_params.centers[i,:], initial_params.covs[i]) for i in 1:initial_params.n_comps]
    components = [Component(pd, gd, shape_prior=deepcopy(size_prior), prior_weight=new_component_weight, can_be_dropped=true)
                    for (pd, gd) in zip(position_distrubutions, gene_distributions)]

    gene_sampler = SingleTrialMultinomial(ones(Int, gene_num))
    sampler = Component(MvNormal(zeros(2)), gene_sampler, shape_prior=deepcopy(size_prior), prior_weight=new_component_weight, can_be_dropped=false) # position_params are never used

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

function load_df(data_path; x_col::Symbol=:x, y_col::Symbol=:y, gene_col::Symbol=:gene, min_molecules_per_gene::Int=0)
    df_spatial = read_spatial_df(data_path, x_col=x_col, y_col=y_col, gene_col=gene_col)

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
