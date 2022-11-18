function plot_diagnostics_panel(df_res::DataFrame, assignment::Vector{Int}, tracer::Dict{Symbol};
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
        n_mols_per_cell = count_array(assignment, drop_zero=true)
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

function plot_transcript_assignment_panel(df_res::DataFrame; clusters::Vector{Int}, gene_colors::Vector{<:Colors.ColorTypes.Color}, file::String, kwargs...)
    @info "Plot transcript assignment"

    gene_colors = Colors.alphacolor.(gene_colors, 0.5)

    gc_plot = plot_dataset_colors(df_res, gene_colors; title="Local expression similarity", kwargs...)

    clust_plot = nothing
    if !isempty(clusters)
        clust_plot = plot_dataset_colors(df_res, gene_colors; annotation=clusters, title="Molecule clustering", kwargs...)
    end

    open(file, "w") do io
        show(io, MIME("text/html"), gc_plot)

        if clust_plot !== nothing
            show(io, MIME("text/html"), clust_plot)
        end
    end
end

function plot_segmentation_report(
        segmented_df::DataFrame; tracer::Dict{Symbol},
        clust_res::Union{NamedTuple, Nothing}=nothing, comp_segs::Union{NamedTuple, Nothing}=nothing,
        plot_transcripts::Bool, diagnostic_file::String, molecule_file::String, kwargs...
    )
    plot_diagnostics_panel(
        segmented_df, segmented_df.cell, tracer;
        clust_res=clust_res, comp_segs=comp_segs, file=diagnostic_file
    )

    if plot_transcripts
        clusters = (:cluster in propertynames(segmented_df)) ? segmented_df.cluster : Int[]
        plot_transcript_assignment_panel(segmented_df; clusters=clusters, file=molecule_file, kwargs...)
    end
end