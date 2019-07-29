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

function plot_cell_borders_polygons(df_spatial::DataFrame, polygons::Array{Array{Float64, 2}, 1}=Array{Float64, 2}[], df_centers=nothing; point_size=2, 
                                    color::Union{Vector, Symbol}=:gene, center_size::Real=3.0, polygon_line_width=1, size=(800, 800), xlims=nothing, ylims=nothing, 
                                    append::Bool=false, alpha=0.5, offset=(0, 0), is_noise::Union{Vector, BitArray, Symbol, Nothing}=nothing, kwargs...)
    if typeof(color) === Symbol
        color = df_spatial[!,color]
    end

    if xlims === nothing
        xlims = (minimum(df_spatial.x), maximum(df_spatial.x))
    end

    if ylims === nothing
        ylims = (minimum(df_spatial.y), maximum(df_spatial.y))
    end

    scat_func = append ? (a...; kw...) -> Plots.scatter!(a...; kw...) : (a...; kw...) -> Plots.scatter(a...; kw..., size=size)
    df_noise = nothing

    if typeof(is_noise) === Symbol
        is_noise = df_spatial[!,is_noise]
    end

    if is_noise !== nothing
        df_noise = df_spatial[is_noise,:]
        df_spatial = df_spatial[.!is_noise,:]
        color = color[.!is_noise]
    end

    fig = scat_func(df_spatial.x .+ offset[1], df_spatial.y .+ offset[2]; color=color, markerstrokewidth=0, markersize=point_size,
                    alpha=alpha, legend=false, format=:png, kwargs...)
    
    if is_noise !== nothing
        Plots.scatter!(df_noise.x .+ offset[1], df_noise.y .+ offset[2]; color="black", 
            markerstrokewidth=0, markersize=point_size, alpha=alpha, legend=false, kwargs...)
    end

    for pg in polygons
        Plots.plot!(Plots.Shape(pg[:,1] .+ offset[1], pg[:,2] .+ offset[2]), fill=(0, 0.0), linewidth=polygon_line_width)
    end

    if df_centers !== nothing
        Plots.scatter!(df_centers[!,:x] .+ offset[1], df_centers[!,:y] .+ offset[2], color=colorant"#cc1300", markerstrokewidth=1, markersize=center_size, legend=false)
    end

    Plots.xlims!(xlims)
    Plots.ylims!(ylims)

    return fig
end

### Tracing

function plot_num_of_cells_per_iterarion(tracer::Dict{String, Any})
    if !("n_components" in keys(tracer)) || length(tracer["n_components"]) == 0
        error("No data about #components per iteration was stored")
    end

    n_components_per_iter = hcat(collect.(values.(tracer["n_components"]))...);
    labels = collect(keys(tracer["n_components"][1]));

    n_components_per_iter = n_components_per_iter[sortperm(labels),:]
    labels = sort(labels)

    p = Plots.plot(legendtitle="Min #molecules", title="Convergence", xlabel="Iteration", ylabel="#Cells",
                   background_color_legend=Colors.RGBA(1.0, 1.0, 1.0, 0.5), legend=:topleft)
    for i in 1:size(n_components_per_iter, 1)
        p = Plots.plot!(n_components_per_iter[i,:], label="$(labels[i])")
    end

    return p
end

function plot_prior_shape_per_iteration(tracer::Dict{String, Any})
    Plots.plot(get.(tracer["prior_shape"], 1, 0) .^ 0.5, label="(eigenvalue 1)^0.5",
        xlabel="Iteration", ylabel="Eigenvalue", title="Shape prior")
    Plots.plot!(get.(tracer["prior_shape"], 2, 0) .^ 0.5, label="(eigenvalue 2)^0.5")
end


### Gene composition visualization

neighborhood_count_matrix(bmm_data::BmmData, k::Int) = neighborhood_count_matrix(bmm_data.x, k, maximum(bmm_data.x[!, :gene]))
function neighborhood_count_matrix(df::DataFrame, k::Int, n_genes::Int=maximum(df[!, :gene]); normalize_by_dist::Bool=true)
    if k < 3
        @warn "Too small value of k: $k. Setting it to 3."
        k = 3
    end

    points = position_data(df);
    neighbors, dists = knn(KDTree(points), points, k, true);

    n_cm = zeros(maximum(df.gene), size(df, 1));

    if !normalize_by_dist
        for (i,ids) in enumerate(neighbors)
            prob_array!(view(n_cm, :, i), df[ids, :gene])
        end

        return n_cm
    end

    # normalize_by_dist doesn't affect result much, but neither it slow the work down. In theory, should work better for border cases.
    med_closest_dist = median(getindex.(dists, 2))
    for (i,(ids, dists)) in enumerate(zip(neighbors, dists))
        c_probs = view(n_cm, :, i)
        for (gene, dist) in zip(df[ids, :gene], dists)
            c_probs[gene] += 1 / max(dist, med_closest_dist)
        end
    end

    return n_cm ./ sum(n_cm, dims=1);
end

function gene_composition_transformation(count_matrix::Array{Float64, 2}; sample_size::Int=10000, seed::Int=42, method::Symbol=:umap, kwargs...)
    sample_size = min(sample_size, size(count_matrix, 2))
    Random.seed!(seed)
    count_matrix_sample = count_matrix[:,randperm(size(count_matrix, 2))[1:sample_size]]

    if method == :umap
        return fit(UmapFit, count_matrix_sample, n_components=3, kwargs...);
    end

    if method != :pca
        error("Unknown method: '$method'")
    end

    return fit(MultivariateStats.PCA, count_matrix_sample, maxoutdim=3, kwargs...);
end

function gene_composition_colors(count_matrix::Array{Float64, 2}, transformation)
    mtx_trans = MultivariateStats.transform(transformation, count_matrix);

    mtx_colors = mtx_trans .- minimum(mtx_trans, dims=2)
    mtx_colors ./= maximum(mtx_colors, dims=2);
    mtx_colors .*= 100;

    return vec(mapslices(col -> Colors.Lab(col...), mtx_colors, dims=1))
end

function shuffle_labels(labels::Array{Int})
    new_labs = deepcopy(labels)
    mask = (new_labs .!= 0)
    new_labs[mask] = shuffle(1:maximum(labels))[new_labs[mask]]
    return new_labs
end

### Summary plots

function extract_plot_information(df_spatial::DataFrame, df_centers::DataFrame, x_start::Real, y_start::Real; color_transformation, frame_size::Int=5000,
                                  min_molecules_per_cell::Int=5, grid_step::Union{Float64, Nothing}=nothing, dens_threshold::Float64=1e-10, k::Int=20, 
                                  plot=true, center_size::Real=3.0)::Dict{Symbol,Any}
    x_end, y_end = [x_start, y_start] .+ frame_size
    cur_df = @where(df_spatial, :x .>= x_start, :y .>= y_start, :x .< x_end, :y .< y_end);
    if size(cur_df, 1) < k + 1
        return Dict{Symbol,Any}()
    end

    if grid_step === nothing
        grid_step = frame_size / 300
    end

    cur_df_centers = subset_df_by_coords(df_centers, cur_df);
    polygons = boundary_polygons(cur_df, cur_df[!,:cell], min_molecules_per_cell=min_molecules_per_cell, grid_step=grid_step, dens_threshold=dens_threshold);
    gene_colors = gene_composition_colors(neighborhood_count_matrix(cur_df, k, maximum(df_spatial[!,:gene])), color_transformation);

    res = Dict(:df => cur_df, :polygons => polygons, :gene_colors => gene_colors, :centers => cur_df_centers)
    if plot
        res[:plot] = plot_cell_borders_polygons(res[:df], res[:polygons], res[:centers]; color=res[:gene_colors], center_size=center_size,
                                                xlims=(x_start, x_end), ylims=(y_start, y_end))
    end

    return res
end

### Colormaps

function map_to_colors(vals::Array{T, 1} where T; h=0, n::Int=11, lims=nothing)
    palette = Colors.sequential_palette(h, n)
    offset = (lims === nothing) ? minimum(vals) : lims[1]
    scale = (lims === nothing) ? maximum(vals) - offset : (lims[2] - lims[1])

    color_ticks = collect(range(0.0, 1.0, length=n)) .* scale .+ offset
    colors = palette[floor.(Int, ((vals .- offset) ./ scale) .* (n - 1) .+ 1)]

    return Dict(:colors=>colors, :ticks=>color_ticks, :palette=>palette)
end

function distinguishable_colors(vals::Array{T, 1} where T)
    id_per_type = Dict(c => i for (i,c) in enumerate(unique(vals)));
    colors_uniq = Colors.distinguishable_colors(length(id_per_type));
    colors = colors_uniq[get.(Ref(id_per_type), vals, 0)];

    return Dict(:colors=>colors, :ticks=>unique(vals), :palette=>colors_uniq)
end

plot_colorbar(colors; kwargs...) = plot_colorbar(colors[:ticks], colors[:palette]; kwargs...)

function plot_colorbar(color_ticks, palette; size=(500, 60), rotation=0)
    p = Plots.bar(color_ticks, ones(length(palette)), color=palette, size=size, legend=false, yticks=false)
    p.subplots[1][:xaxis][:rotation] = rotation;

    return p
end