import GeometricalPredicates
import LightGraphs

using DataFrames
using Statistics
using VoronoiDelaunay

import GeometricalPredicates.getx, GeometricalPredicates.gety

using StatsBase: countmap

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

function adjacency_list(points::AbstractArray{T, 2} where T <: Real; filter::Bool=true, n_mads::Real=2)
    @assert size(points, 1) == 2

    points = deepcopy(points)
    points .-= minimum(points)
    points ./= maximum(points) * 1.1
    points .+= 1.01

    hashes = vec(mapslices(row -> "$(row[1]) $(row[2])", round.(points, digits=10), dims=1));
    is_duplicated = get.(Ref(countmap(hashes)), hashes, 0) .> 1;
    points[:, is_duplicated] .+= (rand(Float64, (2, sum(is_duplicated))) .- 0.5) .* 2e-10;

    points_g = [IndexedPoint2D(points[:,i]..., i) for i in 1:size(points, 2)];

    tess = DelaunayTessellation2D(length(points_g), IndexedPoint2D());
    push!(tess, points_g);

    edge_list = hcat([geti.([geta(v), getb(v)]) for v in delaunayedges(tess)]...);
    adj_dists = vec(sum((points[:, edge_list[1,:]] .- points[:, edge_list[2,:]]) .^ 2, dims=1)) .^ 0.5

    if filter
        adj_dists_log = log10.(adj_dists)
        d_threshold = median(adj_dists_log) + n_mads * mad(adj_dists_log, normalize=true)

        filt_mask = (adj_dists_log .< d_threshold)
        edge_list = edge_list[:, filt_mask]
        adj_dists = adj_dists[filt_mask]
    end

    return edge_list, adj_dists
end

adjacency_list(spatial_df::DataFrame; kwargs...) = adjacency_list(position_data(spatial_df); kwargs...)

function connected_components(adjacent_points::Array{Array{Int, 1}, 1})
    g = LightGraphs.SimpleGraph(length(adjacent_points));
    for (v1, vs) in enumerate(adjacent_points)
        for v2 in vs
            LightGraphs.add_edge!(g, v1, v2)
        end
    end

    return LightGraphs.connected_components(g);
end