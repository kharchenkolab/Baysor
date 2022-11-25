using ArgParse
using Statistics

function parse_preview_commandline(args::Union{Nothing, Array{String, 1}}=nothing)
    s = ArgParseSettings(prog="baysor preview")
    @add_arg_table! s begin
        "--config", "-c"
            help = "TOML file with config"
        "--x-column", "-x"
            help = "Name of x column. Overrides the config value."
        "--y-column", "-y"
            help = "Name of gene column. Overrides the config value."
        "--z-column", "-z"
            help = "Name of gene column. Overrides the config value."
        "--gene-column", "-g"
            help = "Name of gene column. Overrides the config value."
        "--min-molecules-per-cell", "-m"
            help = "Minimal number of molecules for a cell to be considered as real. It's an important parameter, as it's used to infer several other parameters. Overrides the config value."
            arg_type = Int
        "--min-pixels-per-cell"
            help = "Minimal number of pixels per cell. Used to estimate size of the dataset plot."
            arg_type = Int
            default = 15
        "--gene-composition-neigborhood"
            help = "Number of neighbors (i.e. 'k' in k-NN), which is used for gene composition visualization. Larger numbers leads to more global patterns. Default: estimate from min-molecules-per-cell"
            arg_type = Int
        "--exclude-genes"
            help = "Comma-separated list of genes to ignore during segmentation"
        "--output", "-o"
            help = "Name of the output file or path to the output directory"
            default = "preview.html"
        "coordinates"
            help = "CSV file with coordinates of transcripts and gene type"
            required = true
        # TODO: add verbosity level
    end

    return (args === nothing) ? parse_args(s) : parse_args(args, s)
end

function parse_preview_configs(args::Union{Nothing, Array{String, 1}}=nothing)
    r = parse_preview_commandline(args)
    if r === nothing
        return nothing
    end

    if r["config"] !== nothing
        extend_params_with_config!(r, parse_toml_config(r["config"]))
    else
        extend_params_with_config!(r, get_default_config())
    end

    for k in ["gene-column", "x-column", "y-column", "z-column"]
        r[k] = Symbol(r[k])
    end

    if isdir(r["output"]) || isdirpath(r["output"])
        r["output"] = joinpath(r["output"], "preview.html")
    end

    return r
end

function run_cli_preview(args::Union{Nothing, Array{String, 1}}=nothing)
    Random.seed!(1)
    args_str = join(args, " ")
    args = parse_preview_configs(args)
    (args !== nothing) || return 1

    log_file = setup_logger(args["output"], "preview_log.log")

    @info "# CLI params: `$args_str`"
    @info get_baysor_run_str()

    # Run preview

    @info "Loading data..."
    args["min-molecules-per-gene"] = 0
    df_spatial, gene_names = DAT.load_df(args)

    @info "Estimating noise level"
    confidence_nn_id = default_param_value(:confidence_nn_id, args["min-molecules-per-cell"])
    edge_lengths, confidences, d1, d2 = BPR.append_confidence!(df_spatial, nn_id=confidence_nn_id) # TODO: use segmentation mask if available here
    @info "Done"

    @info "Estimating local colors"

    if args["gene-composition-neigborhood"] === nothing
        args["gene-composition-neigborhood"] = default_param_value(
            :composition_neighborhood, args["min-molecules-per-cell"], n_genes=length(gene_names)
        )
    end

    gene_colors = BPR.gene_composition_colors(df_spatial, args["gene-composition-neigborhood"])

    ## Plot

    @info "Building transcript plots"
    gc_plot = REP.plot_dataset_colors(
        df_spatial, gene_colors; min_molecules_per_cell=args["min-molecules-per-cell"],
        min_pixels_per_cell=args["min-pixels-per-cell"], title="Local expression similarity"
    )

    conf_colors = REP.map_to_colors(confidences, lims=(0.0, 1.0), palette=Colors.diverging_palette(10, 250, s=0.75, w=1.0));
    cc_plot = REP.plot_dataset_colors(
        df_spatial, conf_colors[:colors]; min_molecules_per_cell=args["min-molecules-per-cell"],
        min_pixels_per_cell=args["min-pixels-per-cell"], title="Transcript confidence"
    )

    @info "Building gene structure plot"
    vega_plots = Dict{String, Any}()

    gene_emb = BPR.estimate_gene_structure_embedding(df_spatial, gene_names)
    vega_plots["vg_gene_structure"] = REP.plot_gene_structure(gene_emb)
    vega_plots["vg_num_trans"] = REP.plot_num_transcript_overview(df_spatial, confidences, gene_names)
    vega_plots["vg_noise_dist"] = REP.plot_noise_estimation_diagnostics(edge_lengths, confidences, d1, d2, confidence_nn_id=confidence_nn_id)

    @info "Plotting"

    open(args["output"], "w") do io
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