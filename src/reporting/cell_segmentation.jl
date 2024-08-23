using StatsBase: countmap

get_diagnostic_report_template(main_conv_plot::String, vega_plot_js::String) = """
<!DOCTYPE html>
<html>
$(vega_header("Report"))
$(vega_style())
<style>
figure {
    display: flex;
    flex-flow: column;
    padding: 5px;
    max-width: 720px;
    margin: auto;
}

figcaption {
    font: italic smaller sans-serif;
    padding: 3px;
    text-align: center;
}
</style>
<body>
    <section>
        <h1>Algorithm convergence</h1>
        <div id='vg_clust_conv'></div>
        <div id='vg_compart_conv'></div>
        <div id='main_conv'>
            $(main_conv_plot)
        </div>
    </section>
    <section>
        <h1>Molecule confidence</h1>
        <div id='vg_mol_confidence'></div>
        <div id='vg_assignment_confidence'></div>
    </section>
    <section>
        <h1>Number of molecules per cell</h1>
        <div id='vg_num_transcripts'></div>
    </section>
</body>
$(vega_plot_js)
</html>
"""

function plot_diagnostics_panel(
        df_res::DataFrame, assignment::Vector{<:Union{Int, String}}, tracer::Dict{Symbol};
        file::String, clust_res::Union{NamedTuple, Nothing}=nothing, comp_segs::Union{NamedTuple, Nothing}=nothing
    )
    @info "Plot diagnostics"
    vega_plots = Dict{String, Deneb.AbstractSpec}()

    if clust_res !== nothing
        vega_plots["vg_clust_conv"] = plot_clustering_convergence(clust_res, "Molecule clustering convergence")
    end

    if comp_segs !== nothing
        vega_plots["vg_compart_conv"] = plot_clustering_convergence(comp_segs, "Compartment segmentation convergence")
    end

    main_conf_plot = ""
    if (:n_components in keys(tracer)) && length(tracer[:n_components]) != 0
        # Main algorithm convergence
        p_conv = plot_num_of_cells_per_iterarion(tracer[:n_components]);
        io = IOBuffer()
        print(io, makie_to_base64(p_conv))
        main_conf_plot = String(take!(io))
    end

    is_noise = Utils.isnoise.(assignment)
    if :confidence in propertynames(df_res)
        # Confidence per molecule
        vega_plots["vg_mol_confidence"] = plot_confidence_distribution(df_res.confidence, is_noise, size=(500, 250))
    end

    if :assignment_confidence in propertynames(df_res)
        # Assignment confidence
        vega_plots["vg_assignment_confidence"] = plot_assignment_confidence_distribution(df_res.assignment_confidence[.!is_noise])
    end

    # Num. of molecules per cell
    n_mols_per_cell = countmap(assignment[.!is_noise]) |> values |> collect
    n_mols_per_cell = n_mols_per_cell[(n_mols_per_cell .> 1) .& (n_mols_per_cell .< quantile(n_mols_per_cell, 0.99) / 0.99)]
    vega_plots["vg_num_transcripts"] = plot_n_molecules_per_cell(n_mols_per_cell)

    vega_str = (length(vega_plots) > 0) ? vega_plot_html(vega_plots) : ""

    # Write to file
    report_str = get_diagnostic_report_template(main_conf_plot, vega_str)
    open(file, "w") do io
        println(io, report_str)
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
        print(io, makie_to_base64(gc_plot))

        if clust_plot !== nothing
            print(io, makie_to_base64(clust_plot))
        end
    end
end

function plot_segmentation_report(
        segmented_df::DataFrame; tracer::Dict{Symbol}, gene_colors::Symbol=:ncv_color,
        clust_res::Union{NamedTuple, Nothing}=nothing, comp_segs::Union{NamedTuple, Nothing}=nothing,
        plot_transcripts::Bool=true, diagnostic_file::String, molecule_file::String, kwargs...
    )
    plot_diagnostics_panel(
        segmented_df, segmented_df.cell, tracer;
        clust_res=clust_res, comp_segs=comp_segs, file=diagnostic_file
    )

    if plot_transcripts
        clusters = (:cluster in propertynames(segmented_df)) ? segmented_df.cluster : Int[]
        plot_transcript_assignment_panel(
            segmented_df; gene_colors=parse.(Colors.Colorant, segmented_df[!, gene_colors]),
            clusters=clusters, file=molecule_file, kwargs...
        )
    end
end