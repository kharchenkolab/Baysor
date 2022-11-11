@lazy import FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"

import JSON
import Graphs

using SimpleWeightedGraphs
using SparseArrays
using StatsBase: countmap

using Graphs: src, dst

function grid_borders_per_label(grid_labels::Matrix{<:Unsigned})
    d_cols = [-1 0 1 0]
    d_rows = [0 1 0 -1]
    borders_per_label = [Array{UInt32, 1}[] for i in 1:maximum(grid_labels)]

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
                    push!(borders_per_label[cur_label], [col, row])
                    break
                end
            end
        end
    end

    return borders_per_label
end

function border_mst(border_points::Array{<:Vector{<:Real},1}; min_edge_length::Float64=0.1, max_edge_length::Float64=1.5)
    leafsize = min(size(border_points, 1), 10) # Some performance feature?
    tree = KDTree(Matrix{Float64}(hcat(border_points...)), leafsize=leafsize);
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
    return Graphs.kruskal_mst(SimpleWeightedGraph(adj_mtx + adj_mtx'));
end

function find_longest_paths(edges::Array{<:Graphs.AbstractEdge, 1})::Array{Vector{Int}, 1}
    if length(edges) == 0
        return []
    end

    vert_counts = StatsBase.countmap(vcat(map(x -> [src(x), dst(x)], edges)...));
    n_verts = maximum(keys(vert_counts));
    leafs = collect(keys(vert_counts))[collect(values(vert_counts)) .== 1];

    adj_mtx = sparse(src.(edges), dst.(edges), weight.(edges), n_verts, n_verts)
    g_inc = SimpleWeightedGraph(adj_mtx + adj_mtx')

    conn_comps = Graphs.connected_components(g_inc);
    longest_paths = Array{Int, 1}[]
    for vert_inds in conn_comps
        if size(vert_inds, 1) < 3
            continue
        end

        comp_leafs = intersect(vert_inds, leafs)
        max_dist, max_source, max_target = -1, -1, -1
        paths = []
        for src_id in comp_leafs
            paths = Graphs.dijkstra_shortest_paths(g_inc, src_id);
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
        paths = Graphs.dijkstra_shortest_paths(g_inc, max_source);
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

function find_grid_point_labels_kde(pos_data::Matrix{Float64}, cell_labels::Vector{<:Integer}, min_x::Union{Vector{Float64}, Tuple{Float64, Float64}},
        max_x::Union{Vector{Float64}, Tuple{Float64, Float64}}; grid_step::Float64, bandwidth::Float64, dens_threshold::Float64=1e-5, min_molecules_per_cell::Int=3,
        verbose::Bool=false)::Matrix{<:Unsigned}
    coords_per_label = [pos_data[:, ids] for ids in split_ids(cell_labels .+ 1)];
    filt_mask = (size.(coords_per_label, 2) .>= min_molecules_per_cell)
    coords_per_label = coords_per_label[filt_mask];
    cell_labels = (0:maximum(cell_labels))[filt_mask]

    if length(coords_per_label) == 0
        return Matrix{Int}(undef,0,0)
    end

    FFTW.set_num_threads(1); # see https://github.com/JuliaStats/KernelDensity.jl/issues/80

    xw, yw = ceil.(Int, (max_x .- min_x) ./ grid_step)
    label_mat = nothing
    if maximum(cell_labels) < (2 ^ 16) - 1
        label_mat = zeros(UInt16, yw, xw)
    else
        label_mat = zeros(UInt32, yw, xw)
    end
    dens_mat = zeros(Float16, yw, xw)

    p = verbose ? Progress(length(coords_per_label)) : nothing
    for (ci, lab) in zip(1:length(coords_per_label), cell_labels)
        c_bandwidth = (ci == 1) ? bandwidth / 1.5 : bandwidth
        coords = coords_per_label[ci]
        sxi = max(floor(Int, (minimum(coords[1,:]) - min_x[1] - 3 * c_bandwidth) / grid_step), 1);
        ex = maximum(coords[1,:]) + 3 * c_bandwidth;

        syi = max(floor(Int, (minimum(coords[2,:]) - min_x[2] - 3 * c_bandwidth) / grid_step), 1);
        ey = maximum(coords[2,:]) + 3 * c_bandwidth;

        cxs = (sxi * grid_step + min_x[1]):grid_step:ex;
        cys = (syi * grid_step + min_x[2]):grid_step:ey;

        dens = KDE.kde((coords[2,:], coords[1,:]), (cys, cxs), bandwidth=(c_bandwidth, c_bandwidth)).density

        for dyi in 1:length(cys)
            yi = dyi + syi
            if yi > size(dens_mat, 1)
                break
            end

            for dxi in 1:length(cxs)
                xi = dxi + sxi
                if xi > size(dens_mat, 2)
                    break
                end
                d = dens[dyi, dxi] * size(coords, 2)
                if (d > dens_threshold) && (d > dens_mat[yi, xi])
                    dens_mat[yi, xi] = d
                    label_mat[yi, xi] = lab
                end
            end
        end

        if verbose
            next!(p)
        end
    end

    return ImageMorphology.closing(label_mat)
end

function extract_polygons_from_label_grid(grid_labels::Matrix{<:Unsigned}; min_border_length::Int=3, shape_method::Symbol=:path, max_dev::Float64=10.0,
        exclude_labels::Vector{Int}=Int[], offset::Vector{Float64}=[0.0, 0.0], grid_step::Float64=1.0)::Array{Matrix{Float64}, 1}
    borders_per_label = grid_borders_per_label(grid_labels);
    if !isempty(exclude_labels)
        borders_per_label = borders_per_label[setdiff(1:length(borders_per_label), exclude_labels)]
    end

    borders_per_label = borders_per_label[length.(borders_per_label) .>= min_border_length]

    paths = nothing
    edges_per_label = border_mst.(borders_per_label)
    if shape_method == :path
        paths = find_longest_paths.(edges_per_label);
    elseif shape_method == :order
        conn_comps = [Graphs.connected_components(Graphs.SimpleGraphFromIterator(Graphs.SimpleGraphs.SimpleEdge.(src.(edges), dst.(edges)))) for edges in edges_per_label]
        paths = [vcat(order_points_to_polygon.(cc, Ref(hcat(borders...)), max_dev=max_dev)...) for (borders, cc) in zip(borders_per_label, conn_comps)];
    else
        error("Unknown shape_method: $shape_method")
    end

    paths = [ps[length.(ps) .> min_border_length] for ps in paths]
    polygons = vcat([[borders[p] for p in cur_paths] for (borders, cur_paths) in zip(borders_per_label, paths)]...);
    polygons = [copy((hcat(coords...) .* grid_step .+ offset)') for coords in polygons]
    return [vcat(cp, cp[1,:]') for cp in polygons]
end

boundary_polygons(bm_data::BmmData; kwargs...) = boundary_polygons(bm_data.x, bm_data.assignment; kwargs...)
boundary_polygons(spatial_df::DataFrame, args...; kwargs...) =
    boundary_polygons(position_data(spatial_df), args...; kwargs...)

function boundary_polygons(pos_data::Matrix{Float64}, cell_labels::Vector{<:Integer}; min_x::Union{Array{Float64}, Nothing}=nothing, max_x::Union{Array{Float64}, Nothing}=nothing,
                           grid_step::Float64=5.0, min_border_length::Int=3, shape_method::Symbol=:path, max_dev::Float64=10.0,
                           bandwidth::Float64=(grid_step / 2), exclude_labels::Vector{Int}=Int[], kwargs...)::Array{Matrix{Float64}, 1}
    pos_data = pos_data[1:2,:]
    if min_x === nothing
        min_x = vec(mapslices(minimum, pos_data, dims=2))
    end

    if max_x === nothing
        max_x = vec(mapslices(maximum, pos_data, dims=2))
    end

    grid_labels_plane = find_grid_point_labels_kde(pos_data, cell_labels, min_x, max_x; grid_step=grid_step, bandwidth=bandwidth, kwargs...)
    GC.gc()

    if length(grid_labels_plane) == 0
        return Matrix{Float64}[]
    end

    return extract_polygons_from_label_grid(grid_labels_plane; min_border_length=min_border_length, shape_method=shape_method, max_dev=max_dev,
        exclude_labels=exclude_labels, offset=min_x, grid_step=grid_step)
end

function polygons_to_geojson(polygons::Vector{Matrix{T}}) where T <:Real
    geoms = [Dict("type" => "Polygon", "coordinates" => [collect.(eachrow(p))]) for p in polygons]
    return Dict("type" => "GeometryCollection", "geometries" => geoms);
end

polygons_to_geojson(polygons::Dict{String, Vector{Matrix{T}}}) where T <:Real =
    [merge!(polygons_to_geojson(poly), Dict("z" => k)) for (k,poly) in polygons]

function save_polygons_to_geojson(polygons::Union{Vector{Matrix{T}}, Dict{String, Vector{Matrix{T}}}}, path::String) where T <: Real
    open(path, "w") do f
        print(f, JSON.json(polygons_to_geojson(polygons)))
    end
end
