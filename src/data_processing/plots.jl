using Colors
using DataFrames
using DataFramesMeta
using NearestNeighbors
using Random
using Statistics
using VegaLite

import AbstractPlotting
import Base.show
import Base64
import CairoMakie as MK
import ColorSchemes
import MultivariateStats
import Plots

plot_molecules!(args...; kwargs...) = plot_molecules(args...; append=true, kwargs...)

function plot_molecules(df_spatial::DataFrame, polygons::Array{Matrix{Float64}, 1}=Matrix{Float64}[]; markersize=2,
        color::Union{Vector, Symbol, String}=:gene, size=(800, 800), poly_strokewidth=1, xlims=nothing, ylims=nothing, offset=(0, 0),
        is_noise::Union{Vector, BitArray, Symbol, Nothing}=nothing, annotation::Union{<:AbstractVector, Symbol, Nothing} = nothing,
        ann_colors::Union{Nothing, Dict} = nothing, legend=(annotation !== nothing), fontsize=8,
        noise_ann = nothing, shuffle_colors::Bool=false, append::Bool=false, 
        polygon_kwargs::KWArgT=nothing, axis_kwargs::KWArgT=nothing, noise_kwargs::KWArgT=nothing, legend_kwargs::KWArgT=nothing, kwargs...)

    noise_args_default = (marker=:xcross, markersize=(0.75 * markersize), strokewidth=0, color="black");
    axis_args_default = (xticklabelsize=12, yticklabelsize=12, xticksvisible=false, xticklabelsvisible=false, yticksvisible=false, yticklabelsvisible=false);
    legend_args_default = (bgcolor=Colors.RGBA(1, 1, 1, 0.85),);
    polygon_args_default = (strokecolor="black", color="transparent", strokewidth=poly_strokewidth, label="")

    noise_kwargs = update_args(update_args(noise_args_default, Dict(kwargs...)), noise_kwargs)
    axis_kwargs = update_args(axis_args_default, axis_kwargs)
    legend_kwargs = update_args(legend_args_default, legend_kwargs)
    polygon_kwargs = update_args(polygon_args_default, polygon_kwargs)

    if typeof(color) === Symbol
        color = df_spatial[!,color]
    end

    if annotation !== nothing
        if typeof(annotation) === Symbol
            annotation = df_spatial[!,annotation]
        end

        annotation = ["$a" for a in annotation]
        if noise_ann !== nothing
            noise_ann = "$noise_ann"
        end
    end

    xlims = something(xlims, val_range(df_spatial.x))
    ylims = something(ylims, val_range(df_spatial.y))

    if !append
        fig = MK.Figure(resolution=size)
        fig[1, 1] = MK.Axis(fig; axis_kwargs...);
    end

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
        MK.scatter!(df_noise.x .+ offset[1], df_noise.y .+ offset[2]; noise_kwargs...)
    end

    if annotation === nothing
        MK.scatter!(df_spatial.x .+ offset[1], df_spatial.y .+ offset[2]; color=color, 
            strokewidth=0, markersize=markersize, kwargs...)
    else
        ann_vals = annotation[annotation .!= noise_ann] |> unique |> sort
        c_map = Colors.distinguishable_colors(length(ann_vals), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
        if shuffle_colors
            Random.shuffle!(c_map)
        end

        for (color, ann) in zip(c_map, ann_vals)
            style_dict = (ann_colors === nothing) ? Dict() : Dict(:color => ann_colors[ann])
            MK.scatter!(df_spatial.x[annotation .== ann] .+ offset[1], df_spatial.y[annotation .== ann] .+ offset[2];
                strokewidth=0, markersize=markersize, label=ann, color=color, style_dict..., kwargs...)
        end

        if noise_ann in annotation
            MK.scatter!(df_spatial.x[annotation .== noise_ann] .+ offset[1], df_spatial.y[annotation .== noise_ann] .+ offset[2];
                label=noise_ann, noise_kwargs...)
        end
        
        if legend
            MK.axislegend(;legend_kwargs...)
        end
    end

    if length(polygons) > 0
        MK.poly!([MK.Point2.(eachrow(p .+ [offset[1] offset[2]])) for p in polygons]; polygon_kwargs...)
    end

    MK.xlims!(MK.current_axis(), xlims .+ offset[1])
    MK.ylims!(MK.current_axis(), ylims .+ offset[2])

    return MK.current_figure()
end

function shuffle_labels(labels::Array{Int})
    new_labs = deepcopy(labels)
    mask = (new_labs .!= 0)
    new_labs[mask] = shuffle(1:maximum(labels))[new_labs[mask]]
    return new_labs
end

function shuffle_colors(colors::Vector)
    uniq_cols = unique(colors);
    col_ord = Dict(Pair.(uniq_cols, shuffle(1:length(uniq_cols))));
    return [uniq_cols[col_ord[tc]] for tc in colors]
end

function plot_expression_vectors(vecs...; gene_names::Vector{String}, min_expr_frac::Float64=0.05, alpha::Float64=0.5, fontsize::Int=5, text_offset::Float64=0.005,
        labels::Vector{String}=["y$i" for i in 1:length(vecs)], xrotation=60, xticks::Bool=false, kwargs...)
    y_vals = maximum(hcat(vecs...), dims=2) |> vec
    scale = sum(y_vals)

    p = Plots.plot(;widen=false, ylims=(0, maximum(y_vals) + text_offset * scale), kwargs...)
    for (v,l) in zip(vecs, labels)
        p = Plots.bar!(v, alpha=alpha, label=l)
    end

    ann_genes = findall(y_vals .>= min_expr_frac * scale)
    p = Plots.annotate!(collect(zip(ann_genes, y_vals[ann_genes] .+ text_offset * scale, Plots.text.(gene_names[ann_genes], fontsize))))
    if xticks
        Plots.xticks!(1:length(gene_names), gene_names, xrotation=xrotation, xgrid=false)
    end
    return p
end

function clustermap(mtx::T where T <: AbstractMatrix{Float64}, gene_names::Vector{String}; gene_ord::Union{<:AbstractVector{<:Integer}, Bool}=true,
        cell_ord::Union{<:AbstractVector{<:Integer},Bool}=true, diag_genes::Bool=false, kwargs...)
    if typeof(cell_ord) === Bool
        if cell_ord
            cell_dists = 1 .- cor(mtx);
            cell_ord = Clustering.hclust(cell_dists, linkage=:ward).order;
        else
            cell_ord = 1:size(mtx, 2)
        end
    end

    if typeof(gene_ord) === Bool
        if gene_ord
            if diag_genes
                gene_ord = sortperm(vec(mapslices(x -> findmax(x)[2], mtx[:, cell_ord], dims=2)), rev=true);
            else
                gene_dists = 1 .- cor(mtx');
                gene_ord = Clustering.hclust(gene_dists, linkage=:ward).order;
            end
        else
            gene_ord = 1:size(mtx, 1)
        end
    end

    Plots.heatmap(mtx[gene_ord, cell_ord]; yticks=(1:length(gene_names), gene_names[gene_ord]), kwargs...), cell_ord, gene_ord
end

### Tracing

function plot_num_of_cells_per_iterarion(tracer::Dict{Symbol, Any})
    if !(:n_components in keys(tracer)) || length(tracer[:n_components]) == 0
        error("No data about #components per iteration was stored")
    end

    n_components_per_iter = hcat(collect.(values.(tracer[:n_components]))...);
    labels = collect(keys(tracer[:n_components][1]));

    n_components_per_iter = n_components_per_iter[sortperm(labels),:]
    labels = sort(labels)

    fig = MK.Figure(resolution=(600, 400));
    axis_args = (xticksvisible=false, yticksvisible=false, xticklabelsize=14, yticklabelsize=14, xlabel="Iteration", ylabel="Num. cells", 
        xlabelpadding=0, ylabelpadding=0, font="Helvetica")
    fig[1, 1] = MK.Axis(fig; title="Convergence", axis_args...);

    for i in 1:size(n_components_per_iter, 1)
        p = MK.lines!(n_components_per_iter[i,:], label="$(labels[i])", color=ColorSchemes.Dark2_8[mod(i - 1, 8) + 1])
    end
    MK.axislegend("Min #molecules"; position=(0.01, 1), nbanks = 2, titlesize=16, labelsize=14)
    MK.xlims!(MK.current_axis(), 0, size(n_components_per_iter, 2))
    MK.ylims!(MK.current_axis(), 0, maximum(n_components_per_iter))

    ax2 = MK.Axis(fig; bbox=MK.BBox(370, 550, 180, 330), axis_args...)
    subplot_start = ceil(Int, 0.75 * size(n_components_per_iter, 2))
    iter_vals = subplot_start:size(n_components_per_iter, 2)
    for i in 1:size(n_components_per_iter, 1)
        p = MK.lines!(iter_vals, n_components_per_iter[i, iter_vals], color=ColorSchemes.Dark2_8[mod(i - 1, 8) + 1])
    end

    MK.ylims!(MK.current_axis(), 0, maximum(n_components_per_iter[:, iter_vals]))

    return fig
end

### Diagnostics

function plot_clustering_convergence(clust_res::NamedTuple, title::String="")::VegaLite.VLSpec
    plt = DataFrame(:x => 1:length(clust_res.diffs), :diff => 100 .* clust_res.diffs, :change_frac => 100 .* clust_res.change_fracs) |>
        @vlplot(x={:x, title="Iteration"}, title=title) +
        @vlplot(:line, y={:diff, title="Change, %"}, color={datum="Max prob. difference"}) +
        @vlplot(:line, y={:change_frac}, color={datum="Molecules changed"})
end

### Colormaps

function map_to_colors(vals::Array{T, 1} where T; lims=nothing, palette::Vector=Colors.sequential_palette(0, 11))
    offset = (lims === nothing) ? minimum(vals) : lims[1]
    scale = (lims === nothing) ? maximum(vals) - offset : (lims[2] - lims[1])

    if lims !== nothing
        vals = min.(max.(vals, lims[1]), lims[2])
    end

    color_ticks = collect(range(0.0, 1.0, length=length(palette))) .* scale .+ offset
    colors = palette[floor.(Int, ((vals .- offset) ./ scale) .* (length(palette) - 1) .+ 1)]

    return Dict(:colors=>colors, :ticks=>color_ticks, :palette=>palette)
end

function distinguishable_colors(vals::T where T <: AbstractVector, args...; kwargs...)
    id_per_type = Dict(c => i for (i,c) in enumerate(unique(vals)));
    colors_uniq = Colors.distinguishable_colors(length(id_per_type), args...; kwargs...);
    colors = colors_uniq[get.(Ref(id_per_type), vals, 0)];

    return Dict(:colors=>colors, :ticks=>unique(vals), :palette=>colors_uniq)
end

plot_colorbar(colors; kwargs...) = plot_colorbar(colors[:ticks], colors[:palette]; kwargs...)

function plot_colorbar(color_ticks, palette; size=(500, 110), xticklabelsize=12, 
        xlabelsize=14, xlabelpadding=0, kwargs...)
    fig = MK.Figure(resolution=size)
    fig[1, 1] = MK.Axis(fig; yticksvisible=false, yticklabelsvisible=false, xticksvisible=false, xticklabelsize=12, 
        xlabelsize=14, xlabelpadding=0, kwargs...)
    MK.barplot!(color_ticks, ones(length(color_ticks)), color=palette)
    MK.xlims!(MK.current_axis(), val_range(color_ticks))

    return fig
end

### Utils

function Base.show(io::IO, m::MIME"text/html", fig::AbstractPlotting.Scene)
    io64 = Base64.IOBuffer();
    iob64 = Base64.Base64EncodePipe(io64)
    show(iob64, MIME("image/png"), fig);
    close(iob64)
    str = String(take!(io64))
    close(io64)

    print(io, "<img src=\"data:image/png;base64,$str\" />")
    # show(io, MIME("text/plain"), "<img src=\"data:image/png;base64,$str\" />")
end
