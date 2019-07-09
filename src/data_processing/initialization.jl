using NearestNeighbors
using StatsBase: countmap, denserank

import CSV
import Distances

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
    center_data_from_assignment(spatial_df, Array(spatial_df[assignment_col]); kwargs...)

function center_data_from_assignment(spatial_df::DataFrame, assignment::Array{Int, 1}; cov_mult::Float64=0.5, scale_mult::Float64=0.75)
    cluster_centers = hcat([vec(mean(pos_data, dims=2)) for pos_data in position_data_by_assignment(spatial_df, assignment)]...)
    covs = covs_from_assignment(spatial_df, assignment)
    return CenterData(DataFrame(cluster_centers', [:x, :y]), [cov_mult .* m for m in covs],
        estimate_scale_from_centers(cluster_centers, scale_mult=scale_mult))
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
    - `shape_deg_freedom::Int`: number of degrees of freedom for `size_prior`. Ignored if `size_prior !== nothing`.
    - `noise_confidence_threshold::Float64=0.1`: all molecules with confidence < `noise_confidence_threshold` are assigned to noise class
"""
function initial_distributions(df_spatial::DataFrame, prior_centers::CenterData; size_prior::ShapePrior, new_component_weight::Float64,
                               prior_component_weight::Float64, n_degrees_of_freedom_center::Int, default_std::Union{Real, Nothing}=nothing,
                               gene_num::Int=maximum(df_spatial[:gene]), shape_deg_freedom::Int, noise_confidence_threshold::Float64=0.1)
    assignment = assign_cells_to_centers(df_spatial, prior_centers.centers);
    if :confidence in names(df_spatial)
        assignment[df_spatial[:confidence] .< noise_confidence_threshold] .= 0
    end

    mtx_centers = position_data(prior_centers.centers);

    covs = (default_std === nothing) ? covs_from_assignment(df_spatial, assignment) : [Float64[default_std 0; 0 default_std] .^ 2 for i in 1:size(mtx_centers, 2)]

    prior_distributions = [MvNormal(vec(mtx_centers[:, i]), cov) for (i, cov) in enumerate(covs)];
    gene_prior = SingleTrialMultinomial(ones(Int, gene_num));

    n_mols_per_center = count_array(assignment .+ 1, max_value=length(prior_distributions) + 1)[2:end];
    center_priors = [CellCenter(pd.μ, deepcopy(cov), n_degrees_of_freedom_center) for (pd,cov) in zip(prior_distributions, prior_centers.center_covs)]
    components = [Component(pd, deepcopy(gene_prior), shape_prior=deepcopy(size_prior), center_prior=cp,
                            n_samples=n, prior_weight=prior_component_weight, can_be_dropped=false)
                    for (pd, cp, n) in zip(prior_distributions, center_priors, n_mols_per_center)];

    ids_per_comp = split(collect(1:length(assignment)), assignment .+ 1)[2:end]
    for (ids, comp) in zip(ids_per_comp, components)
        if comp.n_samples > 0
            comp.position_params = maximize(comp.position_params, position_data(df_spatial)[:,ids])
            comp.composition_params = maximize(comp.composition_params, composition_data(df_spatial)[ids])
        end
    end

    gene_sampler = SingleTrialMultinomial(ones(Int, gene_num))
    sampler = Component(MvNormal(zeros(2)), gene_sampler, shape_prior=deepcopy(size_prior), prior_weight=new_component_weight, can_be_dropped=false) # position_params are never used

    return components, sampler, assignment
end

function initial_distributions(df_spatial::DataFrame, initial_params::InitialParams; size_prior::ShapePrior, new_component_weight::Float64,
                               gene_smooth::Real=1.0, gene_num::Int=maximum(df_spatial[:gene]))
    gene_distributions = [SingleTrialMultinomial(ones(Int, gene_num), smooth=Float64(gene_smooth)) for i in 1:initial_params.n_comps]

    position_distrubutions = [MvNormal(initial_params.centers[i,:], initial_params.covs[i]) for i in 1:initial_params.n_comps]
    components = [Component(pd, gd, shape_prior=deepcopy(size_prior), prior_weight=new_component_weight, can_be_dropped=true)
                    for (pd, gd) in zip(position_distrubutions, gene_distributions)]

    gene_sampler = SingleTrialMultinomial(ones(Int, gene_num))
    sampler = Component(MvNormal(zeros(2)), gene_sampler, shape_prior=deepcopy(size_prior), prior_weight=new_component_weight, can_be_dropped=false) # position_params are never used

    return components, sampler, initial_params.assignment
end

function initial_distribution_arr(df_spatial::DataFrame, args...; n_frames::Int, n_frames_return::Int=0, kwargs...)::Array{BmmData, 1}
    dfs_spatial = n_frames > 1 ? split_spatial_data(df_spatial, n_frames) : [df_spatial]
    @info "Mean number of molecules per frame: $(median(size.(dfs_spatial, 1)))"
    @info "Done."

    if (n_frames_return > 0) && (n_frames_return < n_frames)
        dfs_spatial = dfs_spatial[1:n_frames_return]
    end

    return initial_distribution_arr(dfs_spatial, args...; kwargs...)
end

function initial_distribution_arr(dfs_spatial::Array{DataFrame, 1}, centers::CenterData; shape_deg_freedom::Int, scale::Union{Number, Nothing}=nothing,
                                  new_component_weight::Number=0.2, center_component_weight::Number=1.0, n_degrees_of_freedom_center::Int=1000,
                                  update_priors::Union{Symbol, Nothing}=nothing, kwargs...)::Array{BmmData, 1}
    if scale === nothing
        scale = centers.scale_estimate
        if update_priors === nothing
            update_priors = :centers
        end
    end
    if update_priors === nothing
        update_priors = :no
    end

    if centers.center_covs === nothing
        centers.center_covs = [diagm(0 => [scale / 2, scale / 2] .^ 2) for i in 1:size(centers.centers, 1)]
    end
    centers_per_frame = subset_by_coords.(Ref(centers), dfs_spatial);

    size_prior = ShapePrior(shape_deg_freedom, [scale, scale].^2);

    if any([size(c.centers, 1) == 0 for c in centers_per_frame])
        error("Some frames don't contain cell centers. Try to reduce number of frames or provide better segmentation.")
    end

    return initialize_bmm_data.(dfs_spatial, centers_per_frame; size_prior=size_prior, new_component_weight=new_component_weight,
                prior_component_weight=center_component_weight, default_std=scale, n_degrees_of_freedom_center=n_degrees_of_freedom_center,
                shape_deg_freedom=shape_deg_freedom, update_priors=update_priors, kwargs...);
end

function initial_distribution_arr(dfs_spatial::Array{DataFrame, 1}; shape_deg_freedom::Int, scale::Number,
                                  n_cells_init::Int=1000, new_component_weight::Number=0.2, kwargs...)::Array{BmmData, 1}
    size_prior = ShapePrior(shape_deg_freedom, [scale, scale].^2)
    initial_params_per_frame = cell_centers_with_clustering.(dfs_spatial, max(div(n_cells_init, length(dfs_spatial)), 2); scale=scale)
    return initialize_bmm_data.(dfs_spatial, initial_params_per_frame; size_prior=size_prior, new_component_weight=new_component_weight, kwargs...)
end

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

    gene_counts = StatsBase.countmap(df_spatial[:gene]);
    large_genes = Set{String}(collect(keys(gene_counts))[collect(values(gene_counts)) .> min_molecules_per_gene]);
    df_spatial = df_spatial[in.(df_spatial[:gene], Ref(large_genes)),:];

    df_spatial[:x] = Array{Float64, 1}(df_spatial[:x])
    df_spatial[:y] = Array{Float64, 1}(df_spatial[:y])
    df_spatial[:gene], gene_names = encode_genes(df_spatial[:gene]);
    return df_spatial, gene_names;
    # return filter_background(df_spatial), gene_names;
end

function split_spatial_data(df::DataFrame, n::Int, key::Symbol)::Array{DataFrame, 1}
    factor = vec(sum(hcat([df[key] .<= quantile(df[key], q) for q in range(1 / n, stop=1.0, length=n)]...), dims=2))
    return split(df, factor)
end

split_spatial_data(df::DataFrame, n_hor::Int, n_ver::Int) = vcat(split_spatial_data.(split_spatial_data(df, n_ver, :y), n_hor, :x)...)
split_spatial_data(df::DataFrame, n::Int) = split_spatial_data(df, floor(Int, sqrt(n)), ceil(Int, sqrt(n))) # TODO: very approximate separation. Example: n=3.
split_spatial_data(df::DataFrame; mean_mols_per_frame::Int) = split_spatial_data(df, round(Int, size(df, 1) / mean_mols_per_frame))
