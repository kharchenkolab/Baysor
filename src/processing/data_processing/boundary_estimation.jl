using StatsBase: countmap

function grid_borders_per_label(grid_labels::Matrix{<:Unsigned})
    d_cols = [-1 0 1 0]
    d_rows = [0 1 0 -1]
    borders_per_label = [Array{UInt32, 1}[] for _ in 1:maximum(grid_labels)]

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


function extract_triangle_verts(tess::VD.DelaunayTessellation2D)
    triangles = Vector{Int}[]
    for tr in tess
        push!(triangles, [geti(f(tr)) for f in [VD.geta, VD.getb, VD.getc, VD.geta]])
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

function border_edges_to_poly(border_edges::Vector{Vector{Int}})
    if !isempty(border_edges)
        border_edges = Dict(a => b for (a,b) in border_edges)
        ce = start = first(keys(border_edges))
        border_poly = [start]
        for _ in 1:10000
            ce = border_edges[ce]
            push!(border_poly, ce)
            (ce != start) || return border_poly
        end
    end

    @warn "Can't build border for a polygon of size $(length(border_edges))" # should never happen
    return Int[]
end

function boundary_polygons_from_grid(grid_labels::Matrix{<:Unsigned}, grid_step::Int=1)
    borders_per_label = grid_borders_per_label(grid_labels);
    polys = Matrix{Float64}[]
    for bords in borders_per_label
        c_coords =Float64.(hcat(bords...));
        tess = adjacency_list(c_coords, adjacency_type=:triangulation, filter=false, return_tesselation=true)[1];
        poly_ids = extract_triangle_verts(tess) |> extract_border_edges |> border_edges_to_poly;
        push!(polys, Matrix((c_coords[:,poly_ids] .* grid_step)'))
    end

    return polys
end

boundary_polygons(bm_data::BmmData) = boundary_polygons(position_data(bm_data), bm_data.assignment)

function boundary_polygons(pos_data::Matrix{Float64}, cell_labels::Vector{<:Integer})
    pos_data = pos_data[1:2,:]
    tess, points_g = adjacency_list(pos_data, adjacency_type=:triangulation, filter=false, return_tesselation=true);
    points_g = points_g[sortperm(geti.(points_g))];

    triangles_per_cell = Dict(k => Vector{Int}[] for k in unique(cell_labels) if k != 0)
    for pids in extract_triangle_verts(tess)
        cids = cell_labels[pids]
        (cids[1] != 0) && (length(unique(cids)) == 1) || continue
        push!(triangles_per_cell[cids[1]], pids)
    end

    borders_per_cell = Dict(k => extract_border_edges(trs) for (k,trs) in triangles_per_cell);

    # fall back to independent triangulagion for cells with loops in the border
    for (cid,bids) in borders_per_cell
        if isempty(bids)
            vids = findall(cell_labels .== cid)
            if length(vids) <= 2
                borders_per_cell[cid] = (length(vids) == 1) ? [[vids; vids]] : [vids, [vids[2], vids[1]]]
                continue
            end
            # for :knn field, there is a chance that a cell will not have a single inner triangle
        else
            any(length(unique(getindex.(bids, d))) != length(bids) for d in 1:2) || continue
        end

        cpoints = points_g[cell_labels .== cid]

        if length(cpoints) == 3
            # Triangulation sometimes fails for this case
            cpoints = geti.(cpoints)
            borders_per_cell[cid] = [[cpoints[1], cpoints[2]], [cpoints[2], cpoints[3]], [cpoints[3], cpoints[1]]]
            continue
        end

        c_tess = VD.DelaunayTessellation2D(length(cpoints), IndexedPoint2D());
        push!(c_tess, cpoints);

        borders_per_cell[cid] = extract_triangle_verts(c_tess) |> extract_border_edges
        if isempty(borders_per_cell[cid])
            @warn "Can't find border for cell $cid" # should never happen
            continue
        end
    end

    return Dict(k => Matrix(pos_data[:, border_edges_to_poly(bes)]') for (k,bes) in borders_per_cell);
end

function boundary_polygons_auto(pos_data::Matrix{Float64}, assignment::Vector{<:Integer}; estimate_per_z::Bool)
    @info "Estimating boundary polygons"

    poly_joint = boundary_polygons(pos_data, assignment);

    if !estimate_per_z || size(pos_data, 1) == 2
        return poly_joint, poly_joint
    end

    z_coords = @view pos_data[3,:]
    z_vals = sort(unique(z_coords))
    if length(z_vals) > (size(pos_data, 1) / 100)
        @warn "To many values of z. Using 2D polygons"
        return poly_joint, poly_joint
    end

    mask_per_z = [(z_coords .≈ z) for z in z_vals]
    poly_per_z = progress_map(
        mask -> boundary_polygons(pos_data[mask,:], assignment[mask]),
        mask_per_z
    );
    poly_per_z = Dict("$k" => p for (k,p) in zip(z_vals, poly_per_z))
    poly_per_z["joint"] = poly_joint

    return poly_joint, poly_per_z
end

