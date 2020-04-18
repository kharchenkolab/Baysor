import UMAP

function plot_noise_estimation_diagnostics(edge_length::Vector{Float64}, confidences::Vector{Float64}, d1::T, d2::T;
        confidence_nn_id::Union{Int, String}="k", linewidth::Float64=2.0, bins::Int=200, kwargs...) where T <: Distributions.UnivariateDistribution
    x_max = quantile(edge_length, 0.99)
    Plots.histogram(edge_length[edge_length .< x_max]; bins=bins, normalize=true, label="Observed distribution",
        xlabel="Distance to $(confidence_nn_id)'th nearest neighbor", ylabel="Density", title="Noise estimation",
        linewidth=linewidth, kwargs...)
    x_vals = range(0, x_max, length=1000)
    n1 = sum(confidences)
    n2 = length(confidences) - n1
    Plots.plot!(x_vals, n1 / (n1 + n2) .* pdf.(d1, x_vals), label="'Real' distributon", linewidth=linewidth)
    Plots.plot!(x_vals, n2 / (n1 + n2) .* pdf.(d2, x_vals), label="Noise distribution", linewidth=linewidth)
end

function plot_num_transcript_overview(genes::Vector{Int}, confidences::Vector{Float64}, gene_names::Vector; size=(800, 300), kwargs...)
    order = sortperm(gene_names)
    plot_expression_vectors(
        count_array(genes[confidences .>= 0.5], max_value=length(gene_names))[order],
        count_array(genes[confidences .< 0.5], max_value=length(gene_names))[order],
        gene_names=gene_names[order]; labels=["Real", "Noise"], xlabel="Gene ID", ylabel="#Transcripts",
        min_expr_frac=0.01, legend=:topleft, alpha=0.3, xlims=(0, length(gene_names)), size=size, kwargs...)
end

function plot_gene_structure(df_spatial::DataFrame, gene_names::Vector, confidence::Vector{Float64}=df_spatial.confidence; kwargs...)
    adjacent_points, adjacent_weights = build_molecule_graph(df_spatial, filter=false);
    cor_mat = pairwise_gene_spatial_cor(df_spatial.gene, df_spatial.confidence, adjacent_points, adjacent_weights);
    cor_vals = vec(cor_mat)
    min_cor = quantile(cor_vals[cor_vals .> 0], 0.001) ./ 5;
    max_cor = quantile(cor_vals, 0.99);
    p_dists = 1 .- max.(min.(cor_mat, max_cor), min_cor) ./ max_cor;
    p_dists[diagind(p_dists)] .= 0.0;

    embedding = UMAP.umap(p_dists, 2; metric=:precomputed, spread=1.0, min_dist=0.1, n_epochs=5000);
    marker_sizes = log.(count_array(df_spatial.gene));
    marker_sizes ./= median(marker_sizes) ./ 2;

    Plots.scatter(embedding[1,:], embedding[2,:]; ms=marker_sizes, legend=false, alpha=0.1, size=(800, 800),
        ticks=false, xlabel="UMAP-1", ylabel="UMAP-2", title="Gene local structure", kwargs...)

    Plots.annotate!(collect(zip(embedding[1,:], embedding[2,:], Plots.text.(gene_names, ceil.(Int, marker_sizes .* 1.5)))), color=Colors.RGBA(0.0, 0.0, 0.0, 0.75))
end

function estimate_panel_plot_size(df_spatial::DataFrame, min_molecules_per_cell::Int, min_pixels_per_cell::Int=7)
    n_cells_per_side = sqrt(size(df_spatial, 1) / min_molecules_per_cell)
    plot_size = min_pixels_per_cell * n_cells_per_side
    x_rng = val_range(df_spatial.x)
    y_rng = val_range(df_spatial.y)
    y_ratio = (y_rng[2] - y_rng[1]) / (x_rng[2] - x_rng[1])

    y_ratio = y_ratio^0.5;
    x_ratio = 1 / y_ratio^0.5;
    plot_size = (x_ratio * plot_size, y_ratio * plot_size)
    return plot_size, n_cells_per_side
end

function plot_dataset_colors(df_spatial::DataFrame, colors::Vector; min_molecules_per_cell::Int, min_pixels_per_cell::Int=7,
        ms::Float64=0.1, alpha::Union{Float64, Vector{Float64}}=0.1, polygons::Array{Array{Float64, 2}, 1}=Array{Float64, 2}[], kwargs...)
    plot_size = estimate_panel_plot_size(df_spatial, min_molecules_per_cell, min_pixels_per_cell)[1]

    return plot_cell_borders_polygons(df_spatial, polygons; color=colors, ms=ms, alpha=alpha, size=plot_size,
        minorgrid=true, fgaxis=Colors.RGBA(0.0, 0.0, 0.0, 0.0), kwargs...)
end