using Random
using Colors
using Measures
using DataFramesMeta

import MultivariateStats
import Plots

plot_cell_borders_density(bmm_data::BmmData; kwargs...) = plot_cell_borders_density(bmm_data.x, bmm_data.assignment; kwargs...)
function plot_cell_borders_density(df_spatial::DataFrame, cell_labels::Array{Int, 1}; min_n_molecules::Int=3, kwargs...)
    polygons = boundary_polygons(df_spatial, cell_labels, min_molecules_per_cell=min_n_molecules);
    return plot_cell_borders_polygons(df_spatial, polygons; kwargs...)
end

function plot_cell_borders_polygons(df_spatial::DataFrame, polygons::Array{Array{Float64, 2}, 1}=Array{Float64, 2}[], df_centers=nothing; point_size=2, color=:gene,
                                    center_size::Real=50.0, polygon_line_width=1, size=(800, 600), xlims=nothing, ylims=nothing, kwargs...)
    if typeof(color) === Symbol
        color = df_spatial[color]
    end

    fig = Plots.scatter(df_spatial[:x], df_spatial[:y]; color=color, markerstrokewidth=0, markersize=point_size,
                        alpha=0.5, legend=false, size=size, format=:png, kwargs...)

    for pg in polygons
        Plots.plot!(Plots.Shape(pg[:,1], pg[:,2]), fill=(0, 0.0), linewidth=polygon_line_width)
    end

    if df_centers !== nothing
        Plots.scatter!(df_centers[:x], df_centers[:y], color=colorant"#cc1300", markerstrokewidth=1, markersize=center_size, legend=false)
    end

    if xlims !== nothing
        Plots.xlims!(xlims)
    end

    if ylims !== nothing
        Plots.ylims!(ylims)
    end

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

neighborhood_count_matrix(bmm_data::BmmData, k::Int) = neighborhood_count_matrix(bmm_data.x, k, maximum(bmm_data.x[:gene]))
function neighborhood_count_matrix(df::DataFrame, k::Int, n_genes::Int=maximum(df[:gene]))
    points = position_data(df);
    neighbors = knn(KDTree(points), points, k)[1];

    return hcat([prob_vec(df[:gene][ids], n_genes) for ids in neighbors]...)
end

function prob_vec(gene_ids::Array{Int, 1}, max_gene::Int, smooth::Int=0)
    probs = Array{Int}(zeros(max_gene))
    for id in gene_ids
        probs[id] += 1
    end

    probs .+= smooth
    return probs ./ sum(probs)
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

    mtx_colors = mtx_trans .- mapslices(minimum, mtx_trans, dims=2)
    mtx_colors ./= mapslices(maximum, mtx_colors, dims=2);
    mtx_colors .*= 100;

    # return vec(mapslices(col -> "#" * hex(Colors.Lab(col...)), mtx_colors, dims=1))
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
                                  min_molecules_per_cell::Int=5, grid_step::Float64=10.0, dens_threshold::Float64=1e-6, k::Int=20, plot=true, center_size::Real=5.0)::Dict{Symbol,Any}
    x_end, y_end = [x_start, y_start] .+ frame_size
    cur_df = @where(df_spatial, :x .>= x_start, :y .>= y_start, :x .< x_end, :y .< y_end);
    if size(cur_df, 1) < k + 1
        return Dict{Symbol,Any}()
    end

    cur_df_centers = subset_df_by_coords(df_centers, cur_df);
    polygons = boundary_polygons(cur_df, cur_df[:assignment], min_molecules_per_cell=min_molecules_per_cell, grid_step=grid_step, dens_threshold=dens_threshold);
    gene_colors = gene_composition_colors(neighborhood_count_matrix(cur_df, k, maximum(df_spatial[:gene])), color_transformation);

    res = Dict(:df => cur_df, :polygons => polygons, :gene_colors => gene_colors, :centers => cur_df_centers)
    if plot
        res[:plot] = plot_cell_borders_polygons(res[:df], res[:polygons], res[:centers]; color=res[:gene_colors], center_size=center_size,
                                                xlims=(x_start, x_end), ylims=(y_start, y_end))
    end

    return res
end