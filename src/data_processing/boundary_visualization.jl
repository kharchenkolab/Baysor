import Distributions
import LightGraphs
import NearestNeighbors
import StatsBase

using DataFrames
using SimpleWeightedGraphs
using SparseArrays

using LightGraphs: src, dst

function estimate_density_norm(train_points::Array{Float64,2}, eval_points::Array{Float64,2})::Array{Float64, 1}
    dist = Distributions.MvNormal(vec(mean(train_points, dims=1)), adjust_cov_matrix(cov(train_points)));
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
                if n_row < 1 || n_col < 1 || n_row > size(grid_labels, 1) || n_col > size(grid_labels, 2)
                    continue
                end

                if grid_labels[n_row, n_col] != cur_label
                    push!(borders_per_label[cur_label], [row, col])
                    break
                end
            end
        end
    end

    return borders_per_label
end

function border_mst(border_points::Array{Array{Int64,1},1}, min_edge_length::Float64=0.1, max_edge_length::Float64=1.5)
    leafsize = min(size(border_points, 1), 10) # Some performance feature?
    tree = NearestNeighbors.KDTree(Array{Float64, 2}(hcat(border_points...)), leafsize=leafsize);
    src_verts, dst_verts, weights = Int[], Int[], Float64[]
    edges = Set{Tuple{Int, Int}}()

    n_neighbs = min(length(border_points), 5)
    for (i, p) in enumerate(border_points)
        inds, dists = NearestNeighbors.knn(tree, p, n_neighbs, true) # 5 to get border for sure
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

    # conn_components = LightGraphs.connected_components.(LightGraphs.SimpleGraph.(adj_per_cell));
    adj_mtx = sparse(src.(edges), dst.(edges), weight.(edges), n_verts, n_verts)
    g_inc = SimpleWeightedGraph(adj_mtx + adj_mtx')
    # g_inc = SimpleWeightedGraph(length(vert_counts))
    # for edge in edges
    #     LightGraphs.add_edge!(g_inc, edge)
    # end

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

boundary_polygons(bm_data::BmmData; kwargs...) = boundary_polygons(bm_data.x, bm_data.assignment; kwargs...)
boundary_polygons(spatial_df::DataFrame, cell_labels::Array{Int64,1}; kwargs...) = 
    boundary_polygons(position_data(spatial_df), cell_labels; kwargs...)
function boundary_polygons(pos_data::Matrix{T} where T <: Real, cell_labels::Array{Int64,1}; min_x::Union{Array, Nothing}=nothing, max_x::Union{Array, Nothing}=nothing,
                           grid_step::Float64=5.0, dens_threshold::Float64=1e-5, min_border_length::Int=3, min_molecules_per_cell::Int=3, use_kde::Bool=true,
                           bandwidth::Float64=grid_step / 2)::Array{Array{Float64, 2}, 1}
    min_x = something(min_x, vec(mapslices(minimum, pos_data, dims=2)))
    max_x = something(max_x, vec(mapslices(maximum, pos_data, dims=2)))

    grid_points = grid_point_coords(min_x, max_x, grid_step);
    grid_points_mat = hcat(grid_points...)

    coords_per_label = [pos_data[:, cell_labels .== i] for i in 1:maximum(cell_labels)];
    coords_per_label = coords_per_label[size.(coords_per_label, 2) .>= min_molecules_per_cell];

    if length(coords_per_label) == 0
        return Array{Float64, 2}[]
    end

    dens_per_label = use_kde ? 
        estimate_density_kde.(coords_per_label, Ref(grid_points_mat), bandwidth) : 
        estimate_density_norm.(coords_per_label, Ref(grid_points_mat))

    densities = hcat(dens_per_label...);

    grid_labels = mapslices(x -> findmax(x)[2], densities, dims=2)
    grid_labels[sum(densities, dims=2) .< dens_threshold] .= 0;
    grid_labels = vec(grid_labels)
    grid_labels_plane = reshape(grid_labels, [length(unique(grid_points_mat[i,:])) for i in 1:2]...);

    borders_per_label = grid_borders_per_label(grid_labels_plane);
    borders_per_label = borders_per_label[map(length, borders_per_label) .>= min_border_length]
    grid_points_plane = reshape(grid_points, size(grid_labels_plane));
    paths = longest_paths.(border_mst.(borders_per_label));

    polygons = vcat([[borders[p] for p in cur_paths] for (borders, cur_paths) in zip(borders_per_label, paths)]...);
    return [Array(hcat([grid_points_plane[c...] for c in coords]...)') for coords in polygons]
end
