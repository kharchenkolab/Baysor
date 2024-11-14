using StatsBase: countmap
using StaticArrays
using CategoricalArrays: cut, levels, levelcode
using OrderedCollections

"""
Check if an array of `points` are inside a triangle defined by its vertex coordinates `tri`.

Returns the number of points that are inside the triangle.
"""
function get_n_points_in_triangle(points::AbstractMatrix{Float64}, tri::SMatrix{2, 3, Float64})
    size(points, 2) > 0 || return 0
    a = tri[:,1]
    v0 = tri[:,3] .- a
    v1 = tri[:,2] .- a

    dot00 = dot(v0, v0)
    dot01 = dot(v0, v1)
    dot11 = dot(v1, v1)

    n_inner_points = 0

    # Calculate barycentric coordinates for each point
    for i in 1:size(points, 2)
        v2 = SVector{2}(points[1, i] - a[1], points[2, i] - a[2])

        dot02 = dot(v0, v2)
        dot12 = dot(v1, v2)

        # Compute barycentric coordinates
        inv_denom = 1.0 / (dot00 * dot11 - dot01 * dot01)
        u = (dot11 * dot02 - dot01 * dot12) * inv_denom
        v = (dot00 * dot12 - dot01 * dot02) * inv_denom

        # Check if point is in the triangle
        n_inner_points += (u >= 0) && (v >= 0) && (u + v < 1)
    end

    return n_inner_points
end

function grid_borders_per_label(grid_labels::Matrix{<:Unsigned})
    d_cols = [-1 0 1 0]
    d_rows = [0 1 0 -1]
    borders_per_label = [Vector{UInt32}[] for _ in 1:maximum(grid_labels)]

    for row in 1:size(grid_labels, 1)
        for col in 1:size(grid_labels, 2)
            cur_label = grid_labels[row, col];
            cur_label > 0 || continue

            for (d_row, d_col) in zip(d_rows, d_cols)
                n_row = row + d_row
                n_col = col + d_col

                if (
                        (n_row < 1) || (n_col < 1) ||
                        (n_row > size(grid_labels, 1)) || (n_col > size(grid_labels, 2)) ||
                        (grid_labels[n_row, n_col] != cur_label)
                    )
                    push!(borders_per_label[cur_label], [col, row])
                    break
                end
            end
        end
    end

    return borders_per_label
end


function extract_triangle_verts(tess::VD.DelaunayTessellation2D; closed::Bool=true)
    triangles = Vector{Int}[]
    vert_fs = closed ? [VD.geta, VD.getb, VD.getc, VD.geta] : [VD.geta, VD.getb, VD.getc]
    for tr in tess
        push!(triangles, [geti(f(tr)) for f in vert_fs])
    end
    return triangles
end

function extract_border_edges(triangle_point_ids::Vector{Vector{Int}})
    !isempty(triangle_point_ids) || return Vector{Int}[]

    edges::Vector{Vector{Int}} = vcat([[ps[ids] for ids in [[1, 2], [2, 3], [3, 1]]] for ps in triangle_point_ids]...)
    edge_counts = countmap(sort.(edges))
    border_edges = [e for (e,n) in edge_counts if n == 1];
    border_edges = [e for e in edges if (e in border_edges) || ([e[2], e[1]] in border_edges)];
    return border_edges
end

function find_border_without_admixture(
        triangles::Vector{Vector{Int}}, pos_data::Matrix{Float64}, non_cell_pos::Matrix{Float64};
        max_iters::Int=100
    )

    edges_per_tri = [[fsort(ps[i1], ps[i2]) for (i1, i2) in [(1, 2), (2, 3), (3, 1)]] for ps in triangles];
    edge_counts = countmap(vcat(edges_per_tri...));

    n_borders_per_node = Dict{Int, Int}() # Used to avoid loops after triangle filtering
    for ((s,e), n) in edge_counts
        n == 1 || continue
        n_borders_per_node[s] = get(n_borders_per_node, s, 0) + 1
        n_borders_per_node[e] = get(n_borders_per_node, e, 0) + 1
    end

    border_triangles = Int[]
    is_excluded = falses(length(edges_per_tri))

    for i in 1:max_iters
        converged = true
        empty!(border_triangles)
        for (i,(tri,es)) in enumerate(zip(triangles, edges_per_tri))
            is_excluded[i] && continue

            n_borders = 0
            for e in es
                n_borders += (edge_counts[e] == 1) # border triangle
            end

            (n_borders > 0) || continue
            if (n_borders > 1)
                push!(border_triangles, i)
                continue
            end

            @views n_admix_points = get_n_points_in_triangle(non_cell_pos, SMatrix{2, 3, Float64}(pos_data[:,tri]))
            (n_admix_points > 0) || continue
            skip = false
            for e in es
                edge_counts[e] == 2 || continue
                if sum((get(n_borders_per_node, k, 0) == 2) for k in e) == 2
                    # The condition shows that both nodes of an internal edge already have other border edges
                    # It means that if we exclude the triangle, we'd have a loop in our polygon
                    push!(border_triangles, i)
                    skip = true
                end
            end
            skip && continue

            converged = false
            is_excluded[i] = true
            for e in es
                edge_counts[e] -= 1
                ne = edge_counts[e]

                for k in e
                    if ne == 0
                        n_borders_per_node[k] -= 1
                    elseif ne == 1
                        n_borders_per_node[k] = get(n_borders_per_node, k, 0) + 1
                    end
                end
            end
        end

        converged && break
        if i == max_iters
            @warn "Polygon filtering didn't converge for cell $cell_id"
        end
    end

    border_edges = [e for (e,n) in edge_counts if n == 1];
    return border_edges
end

function extract_cell_border_ids(
        points::Vector{IndexedPoint2D}, pos_data::Matrix{Float64}, cell_labels::Vector{Int}, cell_id::Int;
        max_iters::Int=100, id_map::Union{Dict{Int, Int}, Nothing}=nothing
    )
    cell_mask = (cell_labels .== cell_id)
    cell_points = points[cell_mask]
    non_cell_pos = pos_data[:, .!cell_mask]

    c_tess = VD.DelaunayTessellation2D(length(cell_points), IndexedPoint2D());
    push!(c_tess, cell_points);

    triangles = extract_triangle_verts(c_tess; closed=false)
    if id_map !== nothing
        triangles = [[id_map[i] for i in t] for t in triangles]
    end

    border_edges = find_border_without_admixture(triangles, pos_data, non_cell_pos; max_iters)
    
    !isempty(border_edges) || return Int[]
    return border_edges_to_poly(border_edges)
end

function border_edges_to_poly(border_edges::Vector{<:Union{Vector{Int}, Tuple{Int, Int}}}; max_border_len::Int=10000)
    (length(border_edges) > 2) || return Int[]

    bes = border_edges
    border_edges = Dict{Int, Tuple{Int, Int}}()
    for (a,b) in bes
        if a in keys(border_edges)
            border_edges[a] = (border_edges[a][1], b)
        else
            border_edges[a] = (b, 0)
        end

        if b in keys(border_edges)
            border_edges[b] = (border_edges[b][1], a)
        else
            border_edges[b] = (a, 0)
        end
    end

    ce = start = first(keys(border_edges))
    pe = 0
    border_poly = [start]

    for _ in 1:max_border_len
        a,b = border_edges[ce]
        ce, pe = ((a == pe) ? b : a), ce

        push!(border_poly, ce)
        (ce == start) && return border_poly
    end

    @show bes
    @warn "Can't build border for a polygon of size $(length(border_edges))" # should never happen
    return Int[]
end

function get_boundary_box_per_cell(pos_data::Matrix{Float64}, assignment::Vector{Int})
    ids_per_cell = split_ids(assignment .+ 1, drop_zero=true)[2:end];
    @views bbox_per_cell = [
        (
            isempty(ids) ?
            Tuple{Float64, Float64}[] :
            [val_range(pos_data[ri, ids]) for ri in 1:2]
        ) for ids in ids_per_cell
    ];

    return bbox_per_cell
end

function extract_ids_per_bbox(xv::Vector{Float64}, yv::Vector{Float64}, bbox_per_cell::Vector{Vector{Tuple{Float64, Float64}}})
    xord = sortperm(xv)
    xv = xv[xord]
    yv = yv[xord]

    ids_per_cell = [Int[] for _ in 1:length(bbox_per_cell)]
    @threads for ci in 1:length(bbox_per_cell)
        bb = bbox_per_cell[ci]
        !isempty(bb) || continue

        (xs, xe), (ys, ye) = bb
        x_range_start = searchsortedfirst(xv, xs)
        for mi in x_range_start:length(xv)
            x = xv[mi]
            y = yv[mi]

            (x <= xe) || break

            if (y >= ys) && (y <= ye)
                push!(ids_per_cell[ci], xord[mi])
            end
        end
    end

    return ids_per_cell
end

function boundary_polygons_from_grid(grid_labels::Matrix{<:Unsigned}; grid_step::Int=1)
    borders_per_label = grid_borders_per_label(grid_labels);
    polys = Matrix{Float64}[]
    for bords in borders_per_label
        if size(bords, 1) < 2
            push!(polys, Float64[;;])
            continue
        end
        c_coords =Float64.(hcat(bords...));
        tess = adjacency_list(c_coords, adjacency_type=:triangulation, filter=false, return_tesselation=true)[1];
        poly_ids = extract_triangle_verts(tess) |> extract_border_edges |> border_edges_to_poly;
        push!(polys, Matrix((c_coords[:,poly_ids] .* grid_step)'))
    end

    return polys
end

boundary_polygons(bm_data::BmmData) = boundary_polygons(position_data(bm_data), bm_data.assignment)

function boundary_polygons(
        pos_data::Matrix{Float64}, cell_labels::Vector{<:Integer}; 
        cell_names::Union{Vector{String}, Nothing}=nothing, offset_rel::Float64=0.01
    )
    if size(pos_data, 1) == 3
        pos_data = pos_data[1:2,:]
    elseif size(pos_data, 1) != 2
        @error "Only 2D and 3D data is supported"
    end

    if size(pos_data, 2) < 3
        return Dict{Int, Matrix{Float64}}()
    end

    points = normalize_points(pos_data)
    points = [IndexedPoint2D(p[1], p[2], i) for (i,p) in enumerate(eachcol(points))];

    bbox_per_cell = get_boundary_box_per_cell(pos_data, cell_labels);
    ids_per_bbox = extract_ids_per_bbox(pos_data[1,:], pos_data[2,:], bbox_per_cell);

    mids_per_cell = split_ids(cell_labels, drop_zero=true);
    dists_to_nn = @p KDTree(pos_data) |>
        knn_parallel(_, pos_data, 2, sorted=true)[2] |>
        getindex.(_, 2);
    offset = mean(dists_to_nn) * offset_rel;

    cell_borders = [Matrix{Float64}([;;]) for _ in 1:length(ids_per_bbox)]
    @threads for cid in 1:length(ids_per_bbox)
        if length(mids_per_cell[cid]) == 1
            cpd = pos_data[:, mids_per_cell[cid]]
            poly = hcat(cpd .+ [offset, 0], cpd .- [offset, 0], cpd .+ [0, offset], cpd .- [0, offset])
        elseif length(mids_per_cell[cid]) == 2
            cpd = pos_data[:, mids_per_cell[cid]]
            center = mean(cpd, dims=2)
            poly = hcat(cpd, center .+ [offset, 0], center .- [offset, 0])
        else
            mids = ids_per_bbox[cid]
            !isempty(mids) || continue

            id_map = Dict(si => ti for (ti,si) in enumerate(mids))
            cpd = pos_data[:, mids]
            cbord = extract_cell_border_ids(points[mids], cpd, cell_labels[mids], cid; id_map)
            @views poly = cpd[:,cbord]
            !isempty(cbord) || continue
        end
        cell_borders[cid] = Matrix(poly')
    end

    if cell_names == nothing
        cell_names = 1:length(cell_borders)
    end

    return Dict(cid => cb for (cid,cb) in zip(cell_names, cell_borders) if !isempty(cb))
end

function boundary_polygons_auto(
        pos_data::Matrix{Float64}, assignment::Vector{<:Integer};
        estimate_per_z::Bool, cell_names::Union{Vector{String}, Nothing}=nothing, verbose::Bool=true,
        max_z_slices::Int=10
    )
    verbose && @info "Estimating boundary polygons"

    poly_joined = boundary_polygons(pos_data, assignment; cell_names);

    if !estimate_per_z || size(pos_data, 1) == 2
        return poly_joined, poly_joined
    end

    z_coords = @view pos_data[3,:]
    z_vals = sort(unique(z_coords))
    if length(z_vals) > max_z_slices
        @warn "To many values of z ($(length(z_vals))). Binning z-stack into $max_z_slices layers for polygon estimation."
        clip = min(1 / max_z_slices / 4, 0.025);
        breaks = @p quantile(z_coords, [clip, 1 - clip]) |> range(_..., length=(max_z_slices + 1));
        z_cat = cut(z_coords, breaks[2:end-1], extend=true)
        z_vals = levels(z_cat)
        z_coords = levelcode.(z_cat)
    end

    mask_per_z = [(z_coords .≈ z) for z in sort(unique(z_coords))] # can't use z_vals here because of binning

    pmap = verbose ? progress_map : map
    poly_per_z = pmap(
        mask -> boundary_polygons(pos_data[:, mask], assignment[mask]; cell_names),
        mask_per_z
    );
    poly_per_z = OrderedDict("$k" => p for (k,p) in zip(z_vals, poly_per_z))
    poly_per_z["2d"] = poly_joined

    return poly_joined, poly_per_z
end

