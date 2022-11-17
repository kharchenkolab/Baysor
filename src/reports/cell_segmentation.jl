using ProgressMeter: progress_map

function plot_diagnostics_panel(df_res::DataFrame, assignment::Array{Int, 1}, tracer::Dict;
        file::String, clust_res::Union{NamedTuple, Nothing}=nothing, comp_segs::Union{NamedTuple, Nothing}=nothing)
    @info "Plot diagnostics"
    open(file, "w") do io
        vega_plots = Dict{String, VL.VLSpec}()
        println(io, "<html>")
        println(io, vega_header("Report"))
        println(io, vega_style())
        println(io, "<body>")
        # Molecule clustering convergence
        if clust_res !== nothing
            println(io, "<div id='vg_clust_conv'></div>")
            vega_plots["vg_clust_conv"] = plot_clustering_convergence(clust_res, "Molecule clustering convergence")
        end

        if comp_segs !== nothing
            println(io, "<div id='vg_compart_conv'></div>")
            vega_plots["vg_compart_conv"] = plot_clustering_convergence(comp_segs, "Compartment segmentation convergence")
        end

        println(io, "<br><br>")
        # Main algorithm convergence
        if (:n_components in keys(tracer)) && length(tracer[:n_components]) != 0
            p_conv = plot_num_of_cells_per_iterarion(tracer);
            show(io, MIME("text/html"), p_conv)
        end

        println(io, "<br><br>")

        # Confidence per molecule
        if :confidence in propertynames(df_res)
            println(io, "<div id='vg_mol_confidence'></div>")
            vega_plots["vg_mol_confidence"] = plot_confidence_distribution(df_res.confidence, assignment, size=(500, 250))
        end

        # Assignment confidence
        if :assignment_confidence in propertynames(df_res)
            println(io, "<div id='vg_assignment_confidence'></div>")
            vega_plots["vg_assignment_confidence"] = plot_assignment_confidence_distribution(df_res.assignment_confidence[assignment .> 0])
        end

        println(io, "<br><br>")

        # Num. of molecules per cell
        n_mols_per_cell = BPR.count_array(assignment, drop_zero=true)
        n_mols_per_cell = n_mols_per_cell[(n_mols_per_cell .> 1) .& (n_mols_per_cell .< quantile(n_mols_per_cell, 0.99) / 0.99)]

        println(io, "<div id='vg_num_transcripts'></div>")
        vega_plots["vg_num_transcripts"] = plot_n_molecules_per_cell(n_mols_per_cell)

        println(io, "</body>")
        if length(vega_plots) > 0
            println(io, vega_plot_html(vega_plots))
        end
        println(io, "</html>")
    end
end

function plot_transcript_assignment_panel(df_res::DataFrame, assignment::Vector{Int}, args::Dict;
        clusters::Vector{Int}, prior_polygons::Array{Matrix{Float64}, 1}, file::String,
        gene_colors::Union{Vector{<:Colors.ColorTypes.Color}, Nothing}=nothing)
    if gene_colors === nothing
        @info "Estimating local colors"
        gene_colors = BPR.gene_composition_colors(df_res, args["gene-composition-neigborhood"])
    end

    @info "Estimating boundary polygons" # For some parameters, this step takes a lot of time and memory
    grid_step = args["scale"] / args["min-pixels-per-cell"];
    polygons = poly_joint = BPR.boundary_polygons(df_res, assignment; grid_step=grid_step, bandwidth=args["scale"]/10);

    if (args["save-polygons"] !== nothing) && (:z in propertynames(df_res))
        if length(unique(df_res.z)) > (size(df_res, 1) / args["min-molecules-per-cell"])
            @warn "To many values of z. Using 2D polygons"
        else
            z_vals = sort(unique(df_res.z))
            mask_per_z = [(df_res.z .≈ z) for z in z_vals]
            poly_per_z = progress_map(
                mask -> boundary_polygons(df_res[mask,:], assignment[mask]; grid_step=grid_step, bandwidth=args["scale"]/10),
                mask_per_z
            );
            poly_per_z = Dict("$k" => p for (k,p) in zip(z_vals, poly_per_z))
            poly_per_z["joint"] = polygons
            polygons = poly_per_z
        end
    end

    if gene_colors !== nothing
        gene_colors = Colors.alphacolor.(gene_colors, 0.5)
    end

    @info "Plot transcript assignment"
    gc_plot = plot_dataset_colors(
        df_res, gene_colors; polygons=poly_joint, prior_polygons=prior_polygons,
        min_molecules_per_cell=args["min-molecules-per-cell"], min_pixels_per_cell=args["min-pixels-per-cell"],
        title="Local expression similarity"
    )

    clust_plot = nothing
    if !isempty(clusters)
        clust_plot = plot_dataset_colors(
            df_res, gene_colors; polygons=poly_joint, prior_polygons=prior_polygons, annotation=clusters,
            min_molecules_per_cell=args["min-molecules-per-cell"], min_pixels_per_cell=args["min-pixels-per-cell"],
            title="Molecule clustering"
        )
    end

    open(file, "w") do io
        show(io, MIME("text/html"), gc_plot)

        if clust_plot !== nothing
            show(io, MIME("text/html"), clust_plot)
        end
    end

    return polygons
end
