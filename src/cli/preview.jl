using Statistics

"""
Plot an html with the dataset preview.

# Args

- `coordinates`:            CSV or Parquet file with coordinates of molecules and gene type

# Options

- `-c, --config=<config.toml>`:         TOML file with a config. The function uses `[data]` and `[plotting]` sections.
- `-x, --x-column=<x>`:                 Name of x column. Overrides the config value.
- `-y, --y-column=<y>`:                 Name of y column. Overrides the config value.
- `-z, --z-column=<z>`:                 Name of z column. Overrides the config value.
- `-g, --gene-column=<gene>`:           Name of gene column. Overrides the config value.
- `-m, --min-molecules-per-cell=<m>`:   Minimal number of molecules for a cell to be considered as real.
                                        It's an important parameter, as it's used to infer several other parameters.
                                        Overrides the config value.
- `-o, --output=<path>`:                Name of the output file or path to the output directory (default: "preview.html")
"""
@cast function preview(
        coordinates::String;
        config::RunOptions=RunOptions(),
        x_column::Symbol=config.data.x, y_column::Symbol=config.data.y, z_column::Symbol=config.data.z,
        gene_column::Symbol=config.data.gene, min_molecules_per_cell::Int=config.data.min_molecules_per_cell,
        output::String="preview.html"
    )
    # Parse options

    opts = config;
    opts.data = from_dict(DataOptions,
        merge(to_dict(opts.data), Dict(
            "x" => x_column, "y" => y_column, "z" => z_column, "gene" => gene_column,
            "min_molecules_per_cell" => min_molecules_per_cell
        ))
    )

    fill_and_check_options!(opts.data)
    if isdir(output) || isdirpath(output)
        output = joinpath(output, "preview.html")
    end

    Random.seed!(1)
    log_file = setup_logger(output, "preview_log.log")

    @info "# CLI params: `$(join(ARGS[2:end], " "))`"
    @info get_baysor_run_str()

    # Estimate preview

    @info "Loading data..."
    df_spatial, gene_names = DAT.load_df(coordinates, opts.data)

    @info "Estimating noise level"
    # TODO: use segmentation mask if available here
    edge_lengths, confidences, d1, d2 = BPR.append_confidence!(df_spatial, nn_id=opts.data.confidence_nn_id)
    @info "Done"

    @info "Estimating local colors"

    fill_and_check_options!(opts.plotting, opts.data.min_molecules_per_cell, length(gene_names))
    gene_colors = BPR.gene_composition_colors(
        df_spatial, opts.plotting.gene_composition_neigborhood; method=Symbol(opts.plotting.ncv_method)
    )

    # Prepare plots

    @info "Building transcript plots"
    gc_plot = REP.plot_dataset_colors(
        df_spatial, gene_colors; title="Local expression similarity",
        opts.data.min_molecules_per_cell, opts.plotting.min_pixels_per_cell, opts.plotting.max_plot_size
    )

    conf_colors = REP.map_to_colors(confidences, lims=(0.0, 1.0), palette=Colors.diverging_palette(10, 250, s=0.75, w=1.0));
    cc_plot = REP.plot_dataset_colors(
        df_spatial, conf_colors[:colors]; title="Transcript confidence",
        opts.data.min_molecules_per_cell, opts.plotting.min_pixels_per_cell, opts.plotting.max_plot_size
    )

    @info "Building gene structure plot"
    vega_plots = Dict{String, Any}()

    gene_emb = BPR.estimate_gene_structure_embedding(df_spatial, gene_names)
    vega_plots["vg_gene_structure"] = REP.plot_gene_structure(gene_emb)
    vega_plots["vg_num_trans"] = REP.plot_num_transcript_overview(df_spatial, confidences, gene_names)
    vega_plots["vg_noise_dist"] = REP.plot_noise_estimation_diagnostics(
        edge_lengths, confidences, d1, d2, confidence_nn_id=opts.data.confidence_nn_id
    )

    # Save plots

    @info "Plotting"

    open(output, "w") do io
        print(io, """
            <!DOCTYPE html>
            <html>
            $(REP.vega_header("Report", "<style> body {background-color: #ffffff;}</style>"))
            $(REP.vega_style())

            <body>
            <div>
            <h2>Content</h2>
            <ul>
                <li><a href="#Transcript_Plots">Transcript plots</a></li>
                <li><a href="#Noise_Level">Noise level</a></li>
                <li><a href="#Gene_Structure">Gene structure</a></li>
            </ul>
            </div>
            """)

        println(io, "<h1 id=\"Transcript_Plots\">Transcript plots</h1><br>")
        println(io, REP.makie_to_base64(gc_plot))
        println(io, "<br>")

        println(io, "<hr>\n<h1 id=\"Noise_Level\">Noise level</h1><br>")
        println(io, REP.makie_to_base64(cc_plot))

        println(io, "<br>")
        println(io, "<div id='vg_noise_dist'></div>")

        println(io, "<br>")
        println(io, "Minimal noise level=$(round(100 * mean(confidences .< 0.01), sigdigits=3))%. ",
            "Expected noise level=$(round(100 * mean(1 .- confidences), sigdigits=2))%.")

        println(io, "<br>")
        println(io, "<hr>\n<h1 id=\"Gene_Structure\">Gene structure</h1><br>")
        println(io, "<div id='vg_num_trans'></div>")

        println(io, "<br>")
        println(io, "<div id='vg_gene_structure'></div>")

        println(io, "</body>")

        println(io, REP.vega_plot_html(vega_plots))
        println(io, "</html>")
    end

    @info "All done!"

    close(log_file)

    return 0
end
