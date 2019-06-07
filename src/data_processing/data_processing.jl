using VoronoiDelaunay
using StatsBase: countmap, denserank

import CSV
import Distances
import GeometricalPredicates
import LightGraphs

import GeometricalPredicates.getx, GeometricalPredicates.gety

function encode_genes(gene_list)
    gene_names = unique(gene_list);
    gene_ids = Dict(zip(gene_names, 1:length(gene_names)))
    return [gene_ids[g] for g in gene_list], gene_names
end

function assign_cells_to_centers(spatial_df::DataFrame, centers::DataFrame)::Array{Int, 1}
    return [v[1] for v in knn(KDTree(position_data(centers)), position_data(spatial_df), 1)[1]]
end

function cell_centers_with_clustering(spatial_df::DataFrame, n_clusters::Int; scale::Real)
    n_clusters = min(n_clusters, size(spatial_df, 1))

    pos_data = position_data(spatial_df);
    cluster_centers = kshiftmedoids(pos_data, n_clusters)[1];
    cluster_labels = kshiftlabels(pos_data, cluster_centers);

    return InitialParams(copy(cluster_centers'), scale^2, cluster_labels)
end

function initial_distributions(df_spatial::DataFrame, prior_centers::DataFrame, center_std::Real; size_prior, new_component_weight::Float64, prior_component_weight::Float64,
                               n_degrees_of_freedom_center::Int, default_cov::Array{Float64, 2}=[1.0 0.0; 0.0 1.0], gene_num::Int=maximum(df_spatial[:gene]))
    adjacent_points = adjacency_list(df_spatial)
    assignment = assign_cells_to_centers(df_spatial, prior_centers);

    center_cov = Float64[center_std 0; 0 center_std] .^ 2;
    mtx_centers = Matrix{Float64}(prior_centers);
    prior_distributions = [MvNormal(mtx_centers[i,:], copy(default_cov)) for i in 1:size(mtx_centers, 1)];
    gene_prior = SingleTrialMultinomial(ones(Int, gene_num));

    n_mols_per_center = count_array(assignment, max_value=length(prior_distributions));
    center_priors = [CellCenter(pd.μ, deepcopy(center_cov), n_degrees_of_freedom_center) for pd in prior_distributions]
    components = [Component(pd, deepcopy(gene_prior), shape_prior=deepcopy(size_prior), center_prior=cp,
                            n_samples=n, prior_weight=prior_component_weight, can_be_dropped=false)
                    for (pd, cp, n) in zip(prior_distributions, center_priors, n_mols_per_center)];

    @assert all([c.n_samples for c in components] .== count_array(assignment, max_value=length(components)));

    ids_per_comp = split(collect(1:length(assignment)), assignment)
    for (ids, comp) in zip(ids_per_comp, components)
        if comp.n_samples > 0
            comp.position_params = maximize(comp.position_params, position_data(df_spatial)[:,ids])
            comp.composition_params = maximize(comp.composition_params, composition_data(df_spatial)[ids])
        end
    end

    shape_sampler = MvNormal(zeros(size(prior_centers, 2)), copy(default_cov))
    gene_sampler = SingleTrialMultinomial(ones(Int, gene_num))
    sampler = Component(shape_sampler, gene_sampler, shape_prior=deepcopy(size_prior), prior_weight=new_component_weight, can_be_dropped=false)

    return BmmData(components, df_spatial, adjacent_points, sampler, assignment)
end

function initial_distributions(df_spatial::DataFrame, initial_params::InitialParams; size_prior::ShapePrior, new_component_weight::Float64, 
                               gene_smooth::Real=1.0, gene_num::Int=maximum(df_spatial[:gene]))
    adjacent_points = adjacency_list(df_spatial)

    gene_distributions = [SingleTrialMultinomial(ones(Int, gene_num), smooth=Float64(gene_smooth)) for i in 1:initial_params.n_comps]

    position_distrubutions = [MvNormal(initial_params.centers[i,:], initial_params.covs[i]) for i in 1:initial_params.n_comps]
    params = Component[]
    params = [Component(pd, gd, shape_prior=deepcopy(size_prior), prior_weight=new_component_weight, can_be_dropped=true)
                    for (pd, gd) in zip(position_distrubutions, gene_distributions)]

    mean_std = reshape([mean([std[i] for std in initial_params.covs]) for i in 1:length(initial_params.covs[1])], size(initial_params.covs[1]));
    shape_sampler = MvNormal(zeros(size(initial_params.centers, 2)), mean_std)
    gene_sampler = SingleTrialMultinomial(ones(Int, gene_num))
    sampler = Component(shape_sampler, gene_sampler, shape_prior=deepcopy(size_prior), prior_weight=new_component_weight, can_be_dropped=false)

    return BmmData(params, df_spatial, adjacent_points, sampler, initial_params.assignment)
end

function initial_distribution_arr(df_spatial::DataFrame; n_frames::Int, shape_deg_freedom::Int, scale::Number, 
        n_cells_init::Int=1000, new_component_weight::Number=0.2, df_centers::Union{DataFrame, Nothing}=nothing,
        center_std::Union{Number, Nothing}=nothing, center_component_weight::Number=1.0, 
        n_degrees_of_freedom_center::Int=1000, min_molecules_per_cell::Int=1)::Array{BmmData, 1}
    dfs_spatial = n_frames > 1 ? split_spatial_data(df_spatial, n_frames) : [df_spatial]

    @info "Mean number of molecules per frame: $(median(size.(dfs_spatial, 1)))"

    @info "Done."

    size_prior = ShapePrior(shape_deg_freedom, [scale, scale].^2);

    bm_data_arr = nothing
    if df_centers !== nothing
        if center_std === nothing
            center_std = scale;
        end

        dfs_centers = subset_df_by_coords.(Ref(df_centers), dfs_spatial);

        if any(size.(dfs_centers, 1) .== 0)
            error("Some frames don't contain cell centers. Try to reduce number of frames or provide better segmentation.")
        end

        bm_data_arr = initial_distributions.(dfs_spatial, dfs_centers, center_std; size_prior=size_prior, new_component_weight=new_component_weight,
                    prior_component_weight=center_component_weight, default_cov=[scale 0.0; 0.0 scale].^2,
                    n_degrees_of_freedom_center=n_degrees_of_freedom_center);
    else
        initial_params_per_frame = cell_centers_with_clustering.(dfs_spatial, max(div(n_cells_init, length(dfs_spatial)), 2); scale=scale);
        bm_data_arr = initial_distributions.(dfs_spatial, initial_params_per_frame, size_prior=size_prior, new_component_weight=new_component_weight);
    end

    return bm_data_arr
end

## Triangulation

struct IndexedPoint2D <: GeometricalPredicates.AbstractPoint2D
    _x::Float64
    _y::Float64
    _index::Int
    IndexedPoint2D(x::Float64,y::Float64, index::Int) = new(x, y, index)
end

IndexedPoint2D() = IndexedPoint2D(0., 0., 0)
IndexedPoint2D(x::Float64, y::Float64) = IndexedPoint2D(x, y, 0)

getx(p::IndexedPoint2D) = p._x
gety(p::IndexedPoint2D) = p._y
geti(p::IndexedPoint2D) = p._index

function adjacency_list(points::AbstractArray{T, 2} where T <: Real; filter::Bool=true, n_mads::Real=2)::Array{Array{Int64,1},1}
    @assert size(points, 1) == 2

    points = deepcopy(points)
    points .-= minimum(points)
    points ./= maximum(points) * 1.1
    points .+= 1.01

    hashes = vec(mapslices(row -> "$(row[1]) $(row[2])", round.(points, digits=3), dims=1));
    is_duplicated = get.(Ref(countmap(hashes)), hashes, 0) .> 1;
    points[:, is_duplicated] .+= (rand(Float64, (2, sum(is_duplicated))) .- 0.5) .* 2e-3;

    points_g = [IndexedPoint2D(points[:,i]..., i) for i in 1:size(points, 2)];

    tess = DelaunayTessellation2D(length(points_g), IndexedPoint2D());
    push!(tess, points_g);

    edge_list = hcat([geti.([geta(v), getb(v)]) for v in delaunayedges(tess)]...);

    if filter
        adj_dists = log10.(vec(sum((points[:, edge_list[1,:]] .- points[:, edge_list[2,:]]) .^ 2, dims=1) .^ 0.5))
        d_threshold = median(adj_dists) + n_mads * mad(adj_dists, normalize=true)
        edge_list = edge_list[:, adj_dists .< d_threshold]
    end

    res = [vcat(v...) for v in zip(split(edge_list[2,:], edge_list[1,:], max_factor=size(points, 2)),
                                   split(edge_list[1,:], edge_list[2,:], max_factor=size(points, 2)))];

    for i in 1:length(res) # point is adjacent to itself
        push!(res[i], i)
    end
    return res
end

adjacency_list(spatial_df::DataFrame) = adjacency_list(position_data(spatial_df))

function connected_components(adjacent_points::Array{Array{Int, 1}, 1})
    g = LightGraphs.SimpleGraph(length(adjacent_points));
    for (v1, vs) in enumerate(adjacent_points)
        for v2 in vs
            LightGraphs.add_edge!(g, v1, v2)
        end
    end

    return LightGraphs.connected_components(g);
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

function read_spatial_df(data_path; x_col::Symbol=:x, y_col::Symbol=:y, gene_col::Union{Symbol, Nothing}=:gene)
    df_spatial = CSV.read(data_path);

    if gene_col === nothing
        df_spatial = df_spatial[[x_col, y_col]]
        DataFrames.rename!(df_spatial, x_col => :x, y_col => :y);
    else
        df_spatial = df_spatial[[x_col, y_col, gene_col]]
        DataFrames.rename!(df_spatial, x_col => :x, y_col => :y, gene_col => :gene);
    end

    return df_spatial
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

# Splitting

function split(df::DataFrame, factor::Array{Int, 1})
    res = Array{DataFrame, 1}(undef, maximum(factor))
    for i in unique(factor)
        res[i] = df[factor .== i, :]
    end

    return res
end

function split_spatial_data(df::DataFrame, n::Int, key::Symbol)::Array{DataFrame, 1}
    factor = vec(sum(hcat([df[key] .<= quantile(df[key], q) for q in range(1 / n, stop=1.0, length=n)]...), dims=2))
    return split(df, factor)
end

split_spatial_data(df::DataFrame, n_hor::Int, n_ver::Int) = vcat(split_spatial_data.(split_spatial_data(df, n_ver, :y), n_hor, :x)...)
split_spatial_data(df::DataFrame, n::Int) = split_spatial_data(df, floor(Int, sqrt(n)), ceil(Int, sqrt(n))) # TODO: very approximate separation. Example: n=3.
split_spatial_data(df::DataFrame; mean_mols_per_frame::Int) = split_spatial_data(df, round(Int, size(df, 1) / mean_mols_per_frame))

function subset_df_by_coords(subsetting_df::DataFrame, coord_df::DataFrame)
    pos_subs = position_data(subsetting_df)
    pos_coords = position_data(coord_df)
    ids = vec(all((pos_subs .>= minimum(pos_coords, dims=2)) .& (pos_subs .<= maximum(pos_coords, dims=2)), dims=1));

    return subsetting_df[ids,:]
end
