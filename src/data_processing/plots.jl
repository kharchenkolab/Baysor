using Random
using Colors
using Measures
using DataFrames
using DataFramesMeta
using NearestNeighbors
using Statistics

import MultivariateStats
import Plots

plot_cell_borders_density(bmm_data::BmmData; kwargs...) = plot_cell_borders_density(bmm_data.x, bmm_data.assignment; kwargs...)
function plot_cell_borders_density(df_spatial::DataFrame, cell_labels::Array{Int, 1}; min_n_molecules::Int=3, kwargs...)
    polygons = boundary_polygons(df_spatial, cell_labels, min_molecules_per_cell=min_n_molecules);
    return plot_cell_borders_polygons(df_spatial, polygons; kwargs...)
end

plot_cell_borders_polygons!(args...; kwargs...) =
    plot_cell_borders_polygons(args...; append=true, kwargs...)

plot_cell_borders_polygons(df_spatial::DataFrame, df_centers::DataFrame; kwargs...) =
    plot_cell_borders_polygons(df_spatial, Array{Float64, 2}[], df_centers; kwargs...)

function plot_cell_borders_polygons(df_spatial::DataFrame, polygons::Array{Array{Float64, 2}, 1}=Array{Float64, 2}[], df_centers=nothing; point_size=2,
                                    color::Union{Vector, Symbol}=:gene, center_size::Real=3.0, polygon_line_width=1, polygon_line_color="black", polygon_alpha::Float64=1.0,
                                    size=(800, 800), xlims=nothing, ylims=nothing, append::Bool=false, alpha=0.5, offset=(0, 0), subplot::Int=1,
                                    is_noise::Union{Vector, BitArray, Symbol, Nothing}=nothing, annotation::Union{TA, Nothing} where TA <: AbstractVector = nothing,
                                    ann_colors::Union{Nothing, Dict} = nothing, legend=(annotation !== nothing), legend_bg_alpha::Float64=0.85, fontsize=8,
                                    noise_ann = nothing, format::Symbol=:png, noise_kwargs::Union{Dict, Nothing}=nothing, shuffle_colors::Bool=false, kwargs...)
    noise_args_default = Dict(:markershape => :xcross, :alpha => alpha, :markersize => point_size / 2, :legend => legend, :markerstrokewidth => 0, :color => "black");
    if noise_kwargs === nothing
        noise_kwargs = noise_args_default
    else
        noise_kwargs = Dict{Symbol, Any}(noise_kwargs)
        for k in keys(noise_args_default)
            if !(k in keys(noise_kwargs))
                noise_kwargs[k] = noise_args_default[k]
            end
        end

        for k in keys(kwargs)
            if !(k in keys(noise_kwargs))
                noise_kwargs[k] = kwargs[k]
            end
        end
    end

    if typeof(color) === Symbol
        color = df_spatial[!,color]
    end

    if xlims === nothing
        xlims = (minimum(df_spatial.x), maximum(df_spatial.x))
    end

    if ylims === nothing
        ylims = (minimum(df_spatial.y), maximum(df_spatial.y))
    end

    fig = append ? Plots.plot!(format=format) : Plots.plot(format=format, size=size, xtickfontsize=fontsize, ytickfontsize=fontsize)

    df_noise = nothing

    if typeof(is_noise) === Symbol
        is_noise = df_spatial[!,is_noise]
    end

    if is_noise !== nothing
        df_noise = df_spatial[is_noise,:]
        df_spatial = df_spatial[.!is_noise,:]
        color = color[.!is_noise]
    end

    if is_noise !== nothing
        Plots.scatter!(df_noise.x .+ offset[1], df_noise.y .+ offset[2]; subplot=subplot, noise_kwargs...)
    end

    if annotation === nothing
        fig = Plots.scatter!(df_spatial.x .+ offset[1], df_spatial.y .+ offset[2]; color=color, markerstrokewidth=0, markersize=point_size,
                             alpha=alpha, legend=false, subplot=subplot, kwargs...)
    else
        ann_vals = annotation[annotation .!= noise_ann] |> unique |> sort
        c_map = Colors.distinguishable_colors(length(ann_vals), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
        # c_map = Colors.distinguishable_colors(length(ann_vals), Colors.colorant"#007a10")
        if shuffle_colors
            Random.shuffle!(c_map)
        end
        for (color, ann) in zip(c_map, ann_vals)
            style_dict = (ann_colors === nothing) ? Dict() : Dict(:color => ann_colors[ann])
            fig = Plots.scatter!(df_spatial.x[annotation .== ann] .+ offset[1], df_spatial.y[annotation .== ann] .+ offset[2];
                                 markerstrokewidth=0, markersize=point_size, alpha=alpha, label=ann, legend=legend, color=color,
                                 bg_legend=Colors.RGBA(1.0, 1.0, 1.0, legend_bg_alpha), subplot=subplot, style_dict..., kwargs...)
        end

        if noise_ann in annotation
            fig = Plots.scatter!(df_spatial.x[annotation .== noise_ann] .+ offset[1], df_spatial.y[annotation .== noise_ann] .+ offset[2];
                                 label=noise_ann, subplot=subplot, noise_kwargs...)
        end
    end

    shapes = [Plots.Shape(pg[:,1] .+ offset[1], pg[:,2] .+ offset[2]) for pg in polygons]
    Plots.plot!(shapes, fill=(0, 0.0), linewidth=polygon_line_width, linecolor=polygon_line_color, alpha=polygon_alpha, subplot=subplot, label="")

    if df_centers !== nothing
        Plots.scatter!(df_centers[!,:x] .+ offset[1], df_centers[!,:y] .+ offset[2], color=colorant"#cc1300", markerstrokewidth=1, markersize=center_size,
            subplot=subplot, label="")
    end

    Plots.xlims!(xlims .+ offset[1], subplot=subplot)
    Plots.ylims!(ylims .+ offset[2], subplot=subplot)

    return fig
end

### Tracing

function plot_num_of_cells_per_iterarion(tracer::Dict{Symbol, Any}; kwargs...)
    if !(:n_components in keys(tracer)) || length(tracer[:n_components]) == 0
        error("No data about #components per iteration was stored")
    end

    n_components_per_iter = hcat(collect.(values.(tracer[:n_components]))...);
    labels = collect(keys(tracer[:n_components][1]));

    n_components_per_iter = n_components_per_iter[sortperm(labels),:]
    labels = sort(labels)

    p = Plots.plot(; legendtitle="Min #molecules", title="Convergence", xlabel="Iteration", ylabel="#Cells",
                   background_color_legend=Colors.RGBA(1.0, 1.0, 1.0, 0.5), legend=:topleft, kwargs...)
    for i in 1:size(n_components_per_iter, 1)
        p = Plots.plot!(n_components_per_iter[i,:], label="$(labels[i])")
    end

    p = Plots.plot!(xlabel="Iteration", ylabel="#Cells", inset = (1, Plots.bbox(0.05,0.05,0.35,0.35,:top,:right)), subplot=2,
        bg_inside=nothing, legend=nothing, kwargs...)

    subplot_start = ceil(Int, 0.75 * size(n_components_per_iter, 2))
    iter_vals = subplot_start:size(n_components_per_iter, 2)
    for i in 1:size(n_components_per_iter, 1)
        p = Plots.plot!(iter_vals, n_components_per_iter[i, iter_vals], subplot=2, ylims=[0, maximum(n_components_per_iter[:, iter_vals]) * 1.05])
    end

    return p
end

function plot_prior_shape_per_iteration(tracer::Dict{Symbol, Any})
    Plots.plot(get.(tracer[:prior_shape], 1, 0) .^ 0.5, label="(eigenvalue 1)^0.5",
        xlabel="Iteration", ylabel="Eigenvalue", title="Shape prior")
    Plots.plot!(get.(tracer[:prior_shape], 2, 0) .^ 0.5, label="(eigenvalue 2)^0.5")
end


### Gene composition visualization

neighborhood_count_matrix(bmm_data::BmmData, k::Int) = neighborhood_count_matrix(bmm_data.x, k, maximum(bmm_data.x.gene))
neighborhood_count_matrix(df::DataFrame, k::Int, args...; kwargs...) =
    neighborhood_count_matrix(position_data(df), df.gene, k, args...; kwargs...)

function neighborhood_count_matrix(pos_data::Matrix{T} where T <: Real, genes::Vector{Int}, k::Int, n_genes::Int=maximum(genes); normalize_by_dist::Bool=true, normalize::Bool=true)
    if k < 3
        @warn "Too small value of k: $k. Setting it to 3."
        k = 3
    end

    if !normalize
        normalize_by_dist = false
    end

    k = min(k, size(pos_data, 2))

    neighbors, dists = knn(KDTree(pos_data), pos_data, k, true);

    n_cm = zeros((normalize ? Float64 : Int), n_genes, size(pos_data, 2));

    if !normalize_by_dist
        @threads for (i,ids) in collect(enumerate(neighbors))
            if !normalize
                count_array!(view(n_cm, :, i), genes[ids])
            else
                prob_array!(view(n_cm, :, i), genes[ids])
            end
        end

        return n_cm
    end

    # normalize_by_dist doesn't affect result much, but neither it slow the work down. In theory, should work better for border cases.
    med_closest_dist = median(getindex.(dists, 2))
    @threads for (i,(ids, dists)) in collect(enumerate(zip(neighbors, dists)))
        c_probs = view(n_cm, :, i)
        for (gene, dist) in zip(genes[ids], dists)
            c_probs[gene] += 1 / max(dist, med_closest_dist)
        end
    end

    return n_cm ./ sum(n_cm, dims=1);
end

function gene_composition_transformation(count_matrix::Array{Float64, 2}, confidence::Vector{Float64}=ones(size(count_matrix, 2));
        sample_size::Int=10000, seed::Int=42, method::Symbol=:umap, kwargs...)
    sample_size = min(sample_size, size(count_matrix, 2))
    Random.seed!(seed)

    pc2 = transform(fit(MultivariateStats.PCA, count_matrix, maxoutdim=2), count_matrix);
    sample_ids = select_ids_uniformly(pc2[1,:], pc2[2,:], confidence, sample_size)

    count_matrix_sample = count_matrix[:,sample_ids]

    if method == :umap
        return fit(UmapFit, count_matrix_sample, n_components=3; kwargs...);
    end

    if method != :pca
        error("Unknown method: '$method'")
    end

    return fit(MultivariateStats.PCA, count_matrix_sample, maxoutdim=3; kwargs...);
end

function gene_composition_colors(count_matrix::Array{Float64, 2}, transformation::UmapFit; color_range::T where T<:Real =2000.0)
    mtx_colors = MultivariateStats.transform(transformation, count_matrix);
    min_vals = minimum(transformation.embedding, dims=2)
    mtx_colors .-= min_vals
    mtx_colors ./= (maximum(transformation.embedding, dims=2) - min_vals);

    mtx_colors[1,:] .*= 100
    mtx_colors[2:3,:] .-= 0.5
    mtx_colors[2:3,:] .*= color_range

    return vec(mapslices(col -> Colors.Lab(col...), mtx_colors, dims=1))
end

function shuffle_labels(labels::Array{Int})
    new_labs = deepcopy(labels)
    mask = (new_labs .!= 0)
    new_labs[mask] = shuffle(1:maximum(labels))[new_labs[mask]]
    return new_labs
end

function plot_expression_vectors(vecs...; gene_names::Vector{String}, min_expr_frac::Float64=0.05, alpha::Float64=0.5, fontsize::Int=5, text_offset::Float64=0.005,
        labels::Vector{String}=["y$i" for i in 1:length(vecs)], kwargs...)
    p = Plots.plot(;kwargs...)
    for (v,l) in zip(vecs, labels)
        p = Plots.bar!(v, alpha=alpha, label=l)
    end

    y_vals = maximum(hcat(vecs...), dims=2) |> vec
    scale = sum(y_vals)
    ann_genes = findall(y_vals .>= min_expr_frac * scale)
    p = Plots.annotate!(collect(zip(ann_genes, y_vals[ann_genes] .+ text_offset * scale, Plots.text.(gene_names[ann_genes], fontsize))))
    return p
end

### Summary plots

subset_df(df_spatial::DataFrame, x_start::Real, y_start::Real, frame_size::Real) =
    @where(df_spatial, :x .>= x_start, :y .>= y_start, :x .< (x_start + frame_size), :y .< (y_start + frame_size));

function plot_cell_boundary_polygons_all(df_res::DataFrame, assignment::Array{Int, 1}, df_centers::Union{DataFrame, Nothing};
                                         gene_composition_neigborhood::Int, frame_size::Int, grid_size::Int=500, return_raw::Bool=false,
                                         min_molecules_per_cell::Int, plot_width::Int=800, margin=5*Plots.mm)
    df_res = @transform(df_res, cell=assignment)

    frame_size = min(frame_size, max(maximum(df_res.x) - minimum(df_res.x), maximum(df_res.y) - minimum(df_res.y)))
    neighb_cm = neighborhood_count_matrix(df_res, gene_composition_neigborhood);
    transformation = gene_composition_transformation(neighb_cm, df_res.confidence)

    borders = [(minimum(df_res[!, s]), maximum(df_res[!, s])) for s in [:x, :y]];
    borders = [collect(range(b[1], b[1] + floor((b[2] - b[1]) / frame_size) * frame_size, step=frame_size)) for b in borders]
    borders = hcat(collect.(Iterators.product(borders...))...);

    df_subsets = subset_df.(Ref(df_res), borders[1,:], borders[2,:], frame_size);
    filt_mask = size.(df_subsets, 1) .> max(gene_composition_neigborhood, min_molecules_per_cell)

    borders = borders[:, filt_mask]
    df_subsets = df_subsets[filt_mask];

    assignments = [df.cell for df in df_subsets];
    genes_per_frame = [df.gene for df in df_subsets];
    grid_step = frame_size / grid_size

    plot_info = @showprogress "Extracting plot info..." pmap(zip(df_subsets, genes_per_frame, assignments)) do (cdf, g, a)
        pd = position_data(cdf)
        pol = boundary_polygons(pd, a; grid_step=grid_step)
        col = gene_composition_colors(neighborhood_count_matrix(pd, g, gene_composition_neigborhood, maximum(df_res.gene)), transformation)
        pol, col
    end;

    df_centers = (df_centers === nothing) ? fill(nothing, length(df_subsets)) : subset_by_coords.(Ref(df_centers), df_subsets);

    if return_raw
        return df_subsets, plot_info, df_centers, borders
    end

    @info "Plotting..."

    plots_col = [plot_cell_borders_polygons(dfs, p, dfc; color=col, xlims=(xs, xs + frame_size), ylims=(ys, ys + frame_size), size=(plot_width, plot_width), margin=margin)
        for (dfs, p, dfc, col, xs, ys) in zip(df_subsets, getindex.(plot_info, 1), df_centers, getindex.(plot_info, 2), borders[1,:], borders[2,:])]

    if !in(:cluster, names(df_res))
        return plots_col, nothing
    end

    plots_clust = [plot_cell_borders_polygons(dfs, p, dfc; annotation=dfs.cluster, xlims=(xs, xs + frame_size), ylims=(ys, ys + frame_size), size=(plot_width, plot_width), margin=margin)
        for (dfs, p, dfc, col, xs, ys) in zip(df_subsets, getindex.(plot_info, 1), df_centers, getindex.(plot_info, 2), borders[1,:], borders[2,:])]

    return plots_col, plots_clust
end

### Colormaps

function map_to_colors(vals::Array{T, 1} where T; lims=nothing, palette=Colors.sequential_palette(0, 11))
    offset = (lims === nothing) ? minimum(vals) : lims[1]
    scale = (lims === nothing) ? maximum(vals) - offset : (lims[2] - lims[1])

    if lims !== nothing
        vals = min.(max.(vals, lims[1]), lims[2])
    end

    color_ticks = collect(range(0.0, 1.0, length=length(palette))) .* scale .+ offset
    colors = palette[floor.(Int, ((vals .- offset) ./ scale) .* (length(palette) - 1) .+ 1)]

    return Dict(:colors=>colors, :ticks=>color_ticks, :palette=>palette)
end

function distinguishable_colors(vals::Array{T, 1} where T, args...; kwargs...)
    id_per_type = Dict(c => i for (i,c) in enumerate(unique(vals)));
    colors_uniq = Colors.distinguishable_colors(length(id_per_type), args...; kwargs...);
    colors = colors_uniq[get.(Ref(id_per_type), vals, 0)];

    return Dict(:colors=>colors, :ticks=>unique(vals), :palette=>colors_uniq)
end

plot_colorbar(colors; kwargs...) = plot_colorbar(colors[:ticks], colors[:palette]; kwargs...)

function plot_colorbar(color_ticks, palette; size=(500, 60), rotation=0)
    p = Plots.bar(color_ticks, ones(length(palette)), color=palette, size=size, legend=false, yticks=false)
    p.subplots[1][:xaxis][:rotation] = rotation;

    return p
end