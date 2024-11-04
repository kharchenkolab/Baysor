using NearestNeighbors
using Random

import Base.show
import Base64

@lazy import ColorSchemes = "35d6a980-a343-548e-a6ea-1d62b119f2f4"

plot_molecules!(args...; kwargs...) = plot_molecules(args...; append=true, kwargs...)

plot_molecules(df_spatial::DataFrame, polygons::AbstractDict{<:Union{Int, String}, Matrix{Float64}}; kwargs...) =
    plot_molecules(df_spatial, collect(values(polygons)); kwargs...)

"""
```julia
plot_molecules(df_spatial, polygons=Matrix{Float64}[]; markersize=2, color=:gene, size=(800, 800), poly_strokewidth=1, xlims=nothing, ylims=nothing, offset=(0, 0), is_noise=nothing, annotation=nothing, ann_colors=nothing, legend=(annotation !== nothing), fontsize=8, ticksvisible=false, noise_ann=nothing, shuffle_colors=false, append=false, polygon_kwargs=nothing, axis_kwargs=nothing, noise_kwargs=nothing, legend_kwargs=nothing, kwargs...)
```

Plot molecules from spatial transcriptomics data with optional polygon overlays.

### Arguments
- `df_spatial`: DataFrame containing spatial information of molecules.
- `polygons=Matrix{Float64}[]`: Array of matrices representing polygons to be plotted.

### Keyword Arguments
- `markersize=2`: Size of markers representing molecules.
- `color=:gene`: Color of markers. Can be a vector of colors, a symbol referring to a column in `df_spatial`, or a string.
- `size=(800, 800)`: Size of the plot.
- `poly_strokewidth=1`: Stroke width of polygon outlines.
- `xlims=nothing`: x-axis limits. If `nothing`, automatically determined.
- `ylims=nothing`: y-axis limits. If `nothing`, automatically determined.
- `offset=(0, 0)`: Offset for the plot.
- `is_noise=nothing`: Indicator for noise points. Can be a vector, BitArray, symbol referring to a column in `df_spatial`, or `nothing`.
- `annotation=nothing`: Annotations for points. Can be a vector, symbol referring to a column in `df_spatial`, or `nothing`.
- `ann_colors=nothing`: Dictionary of colors for annotations.
- `legend=(annotation !== nothing)`: Whether to show a legend.
- `fontsize=8`: Font size for annotations.
- `ticksvisible=false`: Whether to show axis ticks.
- `noise_ann=nothing`: Annotation for noise points.
- `shuffle_colors=false`: Whether to shuffle colors of annotations.
- `append=false`: Whether to append to an existing plot.
- `polygon_kwargs=nothing`: Additional keyword arguments for polygon plotting.
- `axis_kwargs=nothing`: Additional keyword arguments for axis (`MK.Axis`).
- `noise_kwargs=nothing`: Additional keyword arguments for noise points (`MK.scatter!`).
- `legend_kwargs=nothing`: Additional keyword arguments for legend (`MK.axislegend`).
- `kwargs...`: Additional keyword arguments for `MK.scatter!` for molecule plotting.

### Returns
A plot of molecules with optional polygon overlays.
"""
function plot_molecules(
        df_spatial::DataFrame, polygons::Array{Matrix{Float64}, 1}=Matrix{Float64}[]; markersize::Union{<:Real, Vector{<:Real}}=2,
        color::Union{Vector, Symbol, String}=:gene, size=(800, 800), poly_strokewidth=1, xlims=nothing, ylims=nothing, offset=(0, 0),
        is_noise::Union{Vector, BitArray, Symbol, Nothing}=nothing, annotation::Union{<:AbstractVector, Symbol, Nothing} = nothing,
        ann_colors::Union{Nothing, Dict} = nothing, legend=(annotation !== nothing), fontsize=8, ticksvisible::Bool=false,
        noise_ann = nothing, shuffle_colors::Bool=false, append::Bool=false,
        polygon_kwargs::KWArgT=nothing, axis_kwargs::KWArgT=nothing, noise_kwargs::KWArgT=nothing, legend_kwargs::KWArgT=nothing, kwargs...
    )

    noise_args_default = (marker=:xcross, markersize=(0.75 * markersize), strokewidth=0, color="black");
    axis_args_default = (xticklabelsize=12, yticklabelsize=12, xticksvisible=ticksvisible, xticklabelsvisible=ticksvisible, yticksvisible=ticksvisible, yticklabelsvisible=ticksvisible);
    legend_args_default = (backgroundcolor=Colors.RGBA(1, 1, 1, 0.85),);
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
        fig = MK.Figure(;size)
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
        MK.scatter!(
            df_spatial.x .+ offset[1], df_spatial.y .+ offset[2];
            color=color, strokewidth=0, markersize=markersize, kwargs...
        )
    else
        ann_vals = annotation[annotation .!= noise_ann] |> unique |> sort
        c_map = Colors.distinguishable_colors(length(ann_vals), Colors.colorant"#007a10", lchoices=range(20, stop=70, length=15))
        if shuffle_colors
            Random.shuffle!(c_map)
        end

        for (color, ann) in zip(c_map, ann_vals)
            style_dict = (ann_colors === nothing) ? Dict() : Dict(:color => ann_colors[ann])
            mask = (annotation .== ann)
            ms = isa(markersize, Vector) ? markersize[mask] : markersize
            MK.scatter!(
                df_spatial.x[mask] .+ offset[1], df_spatial.y[mask] .+ offset[2];
                strokewidth=0, markersize=ms, label=ann, color=color,
                style_dict..., kwargs...
            )
        end

        if noise_ann in annotation
            mask = (annotation .== noise_ann)
            ms = isa(markersize, Vector) ? markersize[mask] : markersize
            MK.scatter!(
                df_spatial.x[mask] .+ offset[1], df_spatial.y[mask] .+ offset[2];
                strokewidth=0, markersize=ms, label=noise_ann, noise_kwargs...
            )
        end

        if legend
            MK.axislegend(;legend_kwargs...)
        end
    end

    if length(polygons) > 0
        MK.poly!([MK.Point2.(eachrow(p .+ [offset[1] offset[2]])) for p in polygons if Base.size(p, 1) > 1]; polygon_kwargs...)
        # We can also do `for p in polygons MK.lines!(p[:,1], p[:,2], color="black") end`, but this is 10+ times slower
    end

    MK.xlims!(MK.current_axis(), xlims .+ offset[1])
    MK.ylims!(MK.current_axis(), ylims .+ offset[2])

    return MK.current_figure()
end

function plot_expression_vectors(vecs...; gene_names::Vector{String}, min_expr_frac::Float64=0.025, alpha::Float64=0.5, text_offset::Float64=0.005,
        labels::Vector{String}=["y$i" for i in 1:length(vecs)], xrotation=-60, xticks::Bool=false, ylabel::String="Expression")
    y_vals = maximum(hcat(vecs...), dims=2) |> vec
    scale = sum(y_vals)

    p_df = DataFrame(:gene => gene_names)
    layers = []

    for (l,cnt) in zip(labels, vecs)
        p_df[!,l] = cnt
        layer = Mark(:bar, opacity=alpha) *
            Encoding(x="gene", y=(;field=l, type="quantitative", title=ylabel), color=(;datum=l))
        push!(layers, layer)
    end
    p_df[!, :_yf] = y_vals
    p_df[!, :_yt] = y_vals .+ text_offset * scale

    push!(layers, Mark(:text) * Encoding(x="gene", y="_yt:q", text="gene") * transform_filter("datum._yf > $(min_expr_frac * scale)"))

    return Data(p_df) * sum(layers) * config(view=(width=700, height=300))
end

### Tracing

function plot_num_of_cells_per_iterarion(
        n_comps_dict::Vector{Dict{Int64, Int64}}; legend_position::Symbol=:lt,
        figsize=(500, 350), inset_bbox::Tuple{Int, Int, Int, Int}=(300, 480, 150, 310),
        legend_titlesize::Int=16, legend_labelsize::Int=14
    )
    (length(n_comps_dict) > 0) || error("No data about #components per iteration was stored")

    n_components_per_iter = hcat(collect.(values.(n_comps_dict))...);
    labels = collect(keys(n_comps_dict[1]));

    n_components_per_iter = n_components_per_iter[sortperm(labels),:]
    labels = sort(labels)

    fig = MK.Figure(size=figsize);
    axis_args = (
        xticksvisible=false, yticksvisible=false, xticklabelsize=14, yticklabelsize=14,
        xlabel="Iteration", ylabel="Num. cells", xlabelpadding=0, ylabelpadding=0
    )
    fig[1, 1] = MK.Axis(fig; title="Convergence", axis_args...);

    for i in 1:size(n_components_per_iter, 1)
        MK.lines!(n_components_per_iter[i,:], label="$(labels[i])", color=ColorSchemes.Dark2_8[mod(i - 1, 8) + 1])
    end
    MK.axislegend("Min #molecules"; position=legend_position, nbanks = 2, titlesize=legend_titlesize, labelsize=legend_labelsize)
    MK.xlims!(MK.current_axis(), 0, size(n_components_per_iter, 2))
    MK.ylims!(MK.current_axis(), 0, maximum(n_components_per_iter))

    ax2 = MK.Axis(fig; bbox=MK.BBox(inset_bbox...), axis_args...)
    subplot_start = ceil(Int, 0.75 * size(n_components_per_iter, 2))
    iter_vals = subplot_start:size(n_components_per_iter, 2)
    for i in 1:size(n_components_per_iter, 1)
        MK.lines!(iter_vals, n_components_per_iter[i, iter_vals], color=ColorSchemes.Dark2_8[mod(i - 1, 8) + 1])
    end

    MK.ylims!(MK.current_axis(), 0, maximum(n_components_per_iter[:, iter_vals]))

    return fig
end

### Diagnostics

function plot_clustering_convergence(clust_res::NamedTuple, title::String=""; digits::Int=2)
    p_df = DataFrame(
        :x => 1:length(clust_res.diffs),
        :diff => round.(100 .* clust_res.diffs; digits=digits),
        :change_frac => round.(100 .* clust_res.change_fracs; digits=digits)
    )

    return Data(p_df) * Mark(:line) * Encoding(x=(;field="x", type="quantitative", title="Iteration")) *(
        Encoding(y=(;field="diff", type="quantitative", title="Change, %"), color=(;datum="Max prob. difference")) +
        Encoding(y=(;field="change_frac", type="quantitative", title="Change, %"), color=(;datum="Molecules changed"))
    ) * Deneb.title(title) * config(view=(width=300, height=250))
end

### Colormaps

function map_to_colors(vals::Vector{T} where T; lims=nothing, palette::Vector=Colors.sequential_palette(0, 11))
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
    fig = MK.Figure(;size)
    fig[1, 1] = MK.Axis(fig; yticksvisible=false, yticklabelsvisible=false, xticksvisible=false, xticklabelsize=12,
        xlabelsize=14, xlabelpadding=0, kwargs...)
    MK.barplot!(color_ticks, ones(length(color_ticks)), color=palette)
    MK.xlims!(MK.current_axis(), val_range(color_ticks))

    return fig
end

### Utils

function makie_to_base64(fig::Union{MK.Scene, MK.Figure})
    io64 = Base64.IOBuffer();
    iob64 = Base64.Base64EncodePipe(io64)
    show(iob64, MIME("image/png"), fig);
    close(iob64)
    str = String(take!(io64))
    close(io64)

    return "<img src=\"data:image/png;base64,$str\" />"
end
