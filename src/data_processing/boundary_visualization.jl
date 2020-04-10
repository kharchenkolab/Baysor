import Distributions
import LightGraphs
import NearestNeighbors
import StatsBase
import FFTW

using DataFrames
using SimpleWeightedGraphs
using SparseArrays

using LightGraphs: src, dst

function estimate_density_norm(train_points::Array{Float64,2}, eval_points::Array{Float64,2})::Array{Float64, 1}
    dist = Distributions.MvNormal(vec(mean(train_points, dims=1)), adjust_cov_matrix!(cov(train_points)));
    return Distributions.pdf(dist, eval_points)
end

function grid_point_coords(min_x::Vector{T}, max_x::Vector{T}, step::T) where T <: Real
    grid_points = collect(Iterators.product([s:step:e for (s, e) in zip(min_x, max_x)]...));
    return map(collect, grid_points);
end

function grid_borders_per_label(grid_labels::Array{Int64,2})
    d_cols = [-1 0 1 0]
    d_rows = [0 1 0 -1]
    borders_per_label = [Array{Int, 1}[] for i in 1:maximum(grid_labels)]

    for row in 1:size(grid_labels, 1)
        for col in 1:size(grid_labels, 2)
            cur_label = grid_labels[row, col];
            if cur_label == 0
                continue
            end

            for (d_row, d_col) in zip(d_rows, d_cols)
                n_row = row + d_row
                n_col = col + d_col

                if (n_row < 1) || (n_col < 1) || (n_row > size(grid_labels, 1)) || (n_col > size(grid_labels, 2)) ||
                        (grid_labels[n_row, n_col] != cur_label)
                    push!(borders_per_label[cur_label], [row, col])
                    break
                end
            end
        end
    end

    return borders_per_label
end

function border_mst(border_points::Array{Array{T,1},1} where T<:Real; min_edge_length::Float64=0.1, max_edge_length::Float64=1.5)
    leafsize = min(size(border_points, 1), 10) # Some performance feature?
    tree = KDTree(Array{Float64, 2}(hcat(border_points...)), leafsize=leafsize);
    src_verts, dst_verts, weights = Int[], Int[], Float64[]
    edges = Set{Tuple{Int, Int}}()

    n_neighbs = min(length(border_points), 5)
    for (i, p) in enumerate(border_points)
        inds, dists = knn(tree, p, n_neighbs, true) # 5 to get border for sure
        filt_mask = (dists .> min_edge_length) .& (dists .< max_edge_length)

        for (ind, dist) in zip(inds[filt_mask], dists[filt_mask])
            edge = sort([i, ind])
            edge = (edge...,)
            if !(edge in edges)
                push!(src_verts, edge[1])
                push!(dst_verts, edge[2])
                push!(weights, dist)
                push!(edges, edge)
            end
        end
    end

    adj_mtx = sparse(src_verts, dst_verts, weights, size(border_points, 1), size(border_points, 1));
    return LightGraphs.kruskal_mst(SimpleWeightedGraph(adj_mtx + adj_mtx'));
end

function longest_paths(edges::Array{T, 1})::Array{Array{Int, 1}, 1} where T <: LightGraphs.AbstractEdge
    if length(edges) == 0
        return []
    end

    vert_counts = StatsBase.countmap(vcat(map(x -> [src(x), dst(x)], edges)...));
    n_verts = maximum(keys(vert_counts));
    leafs = collect(keys(vert_counts))[collect(values(vert_counts)) .== 1];

    adj_mtx = sparse(src.(edges), dst.(edges), weight.(edges), n_verts, n_verts)
    g_inc = SimpleWeightedGraph(adj_mtx + adj_mtx')

    conn_comps = LightGraphs.connected_components(g_inc);
    longest_paths = Array{Int, 1}[]
    for vert_inds in conn_comps
        if size(vert_inds, 1) < 3
            continue
        end

        comp_leafs = intersect(vert_inds, leafs)
        max_dist, max_source, max_target = -1, -1, -1
        paths = []
        for src_id in comp_leafs
            paths = LightGraphs.dijkstra_shortest_paths(g_inc, src_id);
            cur_target = comp_leafs[findmax(paths.dists[comp_leafs])[2]];
            if isinf(paths.dists[cur_target])
                @warn "Infinite distance for verteces $src_id and $cur_target"
                continue
            end

            if paths.dists[cur_target] > max_dist
                max_dist = paths.dists[cur_target]
                max_target = cur_target
                max_source = src_id
            end
        end

        if max_source < 0 | max_target < 0
            @warn "Source/target vertex id is 0"
            continue
        end

        cur_target = max_target
        path = Int[cur_target]
        paths = LightGraphs.dijkstra_shortest_paths(g_inc, max_source);
        if isinf(paths.dists[max_target])
            @warn "Infinite distance for verteces $src_id and $cur_target"
            continue
        end
        while cur_target != max_source
            cur_target = paths.parents[cur_target]
            push!(path, cur_target)
            if size(path, 1) > 2 * paths.dists[max_target]
                @warn "Can't build path"
                break
            end
        end
        push!(longest_paths, path)
    end

    return longest_paths
end

function order_points_to_polygon(vert_inds::Vector{Int}, border_coords::Matrix{T} where T <: Real; max_dev::T2 where T2 <: Real=10.0)
    polygon_inds = Vector{Int}[]
    while !isempty(vert_inds)
        ci = vert_inds[1]
        cur_ids = [ci]
        kdt = KDTree(Float64.(border_coords[:, vert_inds]));
        while true
            nn_ids = setdiff(vert_inds[inrange(kdt, Float64.(border_coords[:, ci]), max_dev)], cur_ids)
            if length(nn_ids) == 0
                break
            end

            nn_ids = nn_ids[sortperm(vec(sum((border_coords[:, ci] .- border_coords[:, nn_ids]) .^ 2, dims=1)))]

            ci = nn_ids[1]
            push!(cur_ids, ci)
        end

        push!(polygon_inds, cur_ids)
        vert_inds = setdiff(vert_inds, cur_ids)
    end

    return polygon_inds
end

boundary_polygons(bm_data::BmmData; kwargs...) = boundary_polygons(bm_data.x, bm_data.assignment; kwargs...)
boundary_polygons(spatial_df::DataFrame, cell_labels::Array{Int64,1}; kwargs...) =
    boundary_polygons(position_data(spatial_df), cell_labels; kwargs...)

function find_grid_point_labels_dens(pos_data::Matrix{T} where T <: Real, cell_labels::Array{Int64,1}, grid_points_mat::Matrix{T} where T <: Real;
        dens_threshold::Float64=1e-5, min_molecules_per_cell::Int=3, verbose::Bool=false, use_kde::Bool=true, bandwidth::Float64=1.0)::Vector{Int}
    coords_per_label = [pos_data[:, cell_labels .== i] for i in 1:maximum(cell_labels)];
    coords_per_label = coords_per_label[size.(coords_per_label, 2) .>= min_molecules_per_cell];

    if length(coords_per_label) == 0
        return Int[]
    end

    FFTW.set_num_threads(1); # see https://github.com/JuliaStats/KernelDensity.jl/issues/80

    dens_per_label = [Float64[] for i in 1:length(coords_per_label)]
    if use_kde
        p = verbose ? Progress(length(coords_per_label)) : nothing
        # @threads
        for i in 1:length(coords_per_label)
            dens_per_label[i] = estimate_density_kde(coords_per_label[i], grid_points_mat, bandwidth)
            if verbose
                next!(p)
            end
        end
    else
        dens_per_label = estimate_density_norm.(coords_per_label, Ref(grid_points_mat))
    end

    densities = hcat(dens_per_label...);

    grid_labels = mapslices(x -> findmax(x)[2], densities, dims=2)
    grid_labels[sum(densities, dims=2) .< dens_threshold] .= 0;
    return vec(grid_labels)
end

function find_grid_point_labels_knn(pos_data::Matrix{T} where T <: Real, cell_labels::Array{Int64,1}, grid_points_mat::Matrix{T} where T <: Real;
        grid_step::Float64, k::Int=10, verbose::Bool=false)::Vector{Int}
    real_ids = findall(cell_labels .> 0)
    knn_ids, knn_dists = knn(KDTree(grid_points_mat), pos_data[:, real_ids], k);
    grid_labels = zeros(Int, size(grid_points_mat, 2))
    grid_dits = zeros(Int, size(grid_points_mat, 2)) .+ Inf
    p = verbose ? Progress(length(knn_ids)) : nothing
    for ci in 1:length(knn_ids)
        c_ids = knn_ids[ci]
        c_dists = knn_dists[ci]
        for gi in 1:length(c_ids)
            if (c_dists[gi] .< grid_step) && (c_dists[gi] .< grid_dits[c_ids[gi]])
                grid_labels[c_ids[gi]] = cell_labels[real_ids[ci]]
                grid_dits[c_ids[gi]] = c_dists[gi]
            end
        end

        if verbose
            next!(p)
        end
    end

    return grid_labels
end

function boundary_polygons(pos_data::Matrix{T} where T <: Real, cell_labels::Array{Int64,1}; min_x::Union{Array, Nothing}=nothing, max_x::Union{Array, Nothing}=nothing,
                           grid_step::Float64=5.0, min_border_length::Int=3, method::Symbol=:knn, shape_method::Symbol=:path, max_dev::TD where TD <: Real=10.0,
                           bandwidth::T2 where T2 <: Real =grid_step / 2, exclude_labels::Vector{Int}=Int[], kwargs...)::Array{Array{Float64, 2}, 1}
    min_x = something(min_x, vec(mapslices(minimum, pos_data, dims=2)))
    max_x = something(max_x, vec(mapslices(maximum, pos_data, dims=2)))

    grid_points = grid_point_coords(min_x, max_x, grid_step);
    grid_points_mat = hcat(grid_points...)

    grid_labels = Int[]
    if method == :knn
        grid_labels = find_grid_point_labels_knn(pos_data, cell_labels, grid_points_mat; grid_step=grid_step, kwargs...)
    elseif method in (:norm, :kde)
        grid_labels = find_grid_point_labels_dens(pos_data, cell_labels, grid_points_mat; use_kde=(method == :kde), bandwidth=bandwidth, kwargs...)
    else
        error("Unknown method: $method")
    end

    if length(grid_labels) == 0
        return Matrix{Float64}[]
    end

    grid_labels_plane = reshape(grid_labels, [length(unique(grid_points_mat[i,:])) for i in 1:2]...);

    borders_per_label = grid_borders_per_label(grid_labels_plane);
    if !isempty(exclude_labels)
        borders_per_label = borders_per_label[setdiff(1:length(borders_per_label), exclude_labels)]
    end

    borders_per_label = borders_per_label[length.(borders_per_label) .>= min_border_length]
    grid_points_plane = reshape(grid_points, size(grid_labels_plane));

    paths = nothing
    edges_per_label = border_mst.(borders_per_label)
    if shape_method == :path
        paths = longest_paths.(edges_per_label);
    elseif shape_method == :order
        conn_comps = [LightGraphs.connected_components(LightGraphs.SimpleGraphFromIterator(LightGraphs.SimpleGraphs.SimpleEdge.(src.(edges), dst.(edges)))) for edges in edges_per_label]
        paths = [vcat(order_points_to_polygon.(cc, Ref(hcat(borders...)), max_dev=max_dev)...) for (borders, cc) in zip(borders_per_label, conn_comps)];
    else
        error("Unknown shape_method: $shape_method")
    end

    paths = [ps[length.(ps) .> min_border_length] for ps in paths]
    polygons = vcat([[borders[p] for p in cur_paths] for (borders, cur_paths) in zip(borders_per_label, paths)]...);
    polygons = [Array(hcat([grid_points_plane[c...] for c in coords]...)') for coords in polygons]
    return [vcat(cp, cp[1,:]') for cp in polygons]
end

## Futures

_is_inner(ps::Array{Float64, 2}, grid_step::Float64; k=min(5, size(ps, 2))) =
    [all(ds .< grid_step * 1.001) for ds in knn(KDTree(ps), ps, k)[2]]

# For future optimizations. Should replace boundary_polygons
function _boundary_polygons_opt(pos_data::Matrix{T} where T <: Real, cell_labels::Array{Int64,1}; min_x::Union{Array, Nothing}=nothing, max_x::Union{Array, Nothing}=nothing,
        grid_step::Float64=5.0, min_border_length::Int=3)
    real_ids = findall(cell_labels .> 0)
    pos_tree = KDTree(pos_data[:, real_ids])
    grid_labels = Int[]
    grid_points = Vector{Float64}[]
    for (x,y) in Iterators.product([s:grid_step:e for (s, e) in zip(min_x, max_x)]...)
        id,dist = knn(pos_tree, [x,y], 1)
        if dist[1] > grid_step * 2.5 # TODO: 2.5 is to account for holes in molecules inside a cell. Won't work in general case. Also makes polygons larger than it should be.
            continue
        end
        push!(grid_points, [x,y])
        push!(grid_labels, cur_df.cell[real_ids[id[1]]])
    end

    # t_pd = hcat(grid_points...);
    # Plots.scatter(t_pd[1,:], t_pd[2,:], color=Baysor.distinguishable_colors(grid_labels)[:colors], format=:png, legend=:none,
    #     ms=0.1, markerstrokewidth=0)

    grid_points_per_label = split(grid_points, grid_labels .+ 1)[2:end];
    # TODO: it's alternative to grid_borders_per_label, which picks only border molecules. As cells can have holes, it won't work in general case
    grid_points_per_label = [ps[.!_is_inner(hcat(ps...), grid_step)] for ps in grid_points_per_label if length(ps) >= min_border_length];
    grid_points_per_label = grid_points_per_label[length.(grid_points_per_label) .>= min_border_length];

    # t_pd = hcat([hcat(x...) for x in grid_points_per_label]...);
    # Plots.scatter(t_pd[1,:], t_pd[2,:], format=:png, legend=:none, ms=0.1, markerstrokewidth=0, size=(750, 750))
    paths = longest_paths.(border_mst.(grid_points_per_label, min_edge_length=grid_step * 0.1, max_edge_length=grid_step * 1.5));
    polygons = vcat([[borders[p] for p in cur_paths] for (borders, cur_paths) in zip(grid_points_per_label, paths)]...);
    return [copy(hcat(x...)') for x in polygons];
end