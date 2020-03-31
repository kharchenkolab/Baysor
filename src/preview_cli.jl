using ArgParse
using DataFrames
using DataFramesMeta
using Distributed
using ProgressMeter
using Statistics

import CSV
import Pkg.TOML
import Plots

function parse_preview_commandline(args::Union{Nothing, Array{String, 1}}=nothing) # TODO: add verbosity level
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--config", "-c"
            help = "TOML file with config"
        "--x-column", "-x"
            help = "Name of x column. Overrides the config value."
        "--y-column", "-y"
            help = "Name of gene column. Overrides the config value."
        "--gene-column"
            help = "Name of gene column. Overrides the config value."
        "--min-molecules-per-cell", "-m"
            help = "Minimal number of molecules for a cell to be considered as real. It's an important parameter, as it's used to infer several other parameters. Overrides the config value."
            arg_type = Int
        "--min-pixels-per-cell"
            help = "Minimal number of pixels per cell. Used to estimate size of the dataset plot."
            arg_type = Int
            default = 7
        "--output", "-o"
            help = "Name of the output file or path to the output directory"
            default = "preview.html"
        "coordinates"
            help = "CSV file with coordinates of transcripts and gene type"
            required = true
        # TODO: add gene-composition-neighborhood
    end

    return (args === nothing) ? parse_args(s) : parse_args(args, s)
end

function parse_preview_configs(args::Union{Nothing, Array{String, 1}}=nothing)
    r = parse_preview_commandline(args)
    if r["config"] !== nothing
        extend_params_with_config!(r, parse_toml_config(r["config"]))
    else
        extend_params_with_config!(r, get_default_config())
    end

    for k in ["gene-column", "x-column", "y-column"]
        r[k] = Symbol(r[k])
    end

    if isdir(r["output"]) || isdirpath(r["output"])
        r["output"] = joinpath(r["output"], "preview.html")
    end

    return r
end

function run_preview_cli(args::Union{Nothing, Array{String, 1}, String}=nothing)
    if args == "build" # need to call this function during build without any actual work
        return 0
    end

    args = parse_preview_configs(args)

    # Set up logger

    log_file = open(append_suffix(args["output"], "preview_log.log"), "w")
    Base.CoreLogging.global_logger(DoubleLogger(log_file, stdout; force_flush=true))

    # Set up plotting

    ENV["GKSwstype"] = "100"; # Disable output device
    Plots.gr()

    # Run preview

    @info "Run"
    @info "Load data..."
    args["min-molecules-per-gene"] = 0
    df_spatial, gene_names = load_df(args)

    @info "Estimating noise level"
    confidence_nn_id = default_param_value(:confidence_nn_id, args["min-molecules-per-cell"])
    edge_lengths, confidences, d1, d2 = append_confidence!(df_spatial, nn_id=confidence_nn_id) # TODO: use segmentation mask if available here
    @info "Done"

    @info "Estimating local neighborhoods"

    composition_neighborhood = default_param_value(:composition_neighborhood, args["min-molecules-per-cell"], n_genes=length(gene_names))
    neighb_cm = neighborhood_count_matrix(df_spatial, composition_neighborhood);
    color_transformation = gene_composition_transformation(neighb_cm[:, confidences .> 0.95]);

    @info "Estimating local colors"
    gene_colors = Baysor.gene_composition_colors(neighb_cm, color_transformation)

    @info "Building transcript plots"
    gc_plot, cc_plot = plot_dataset_colors(df_spatial, gene_colors; min_molecules_per_cell=args["min-molecules-per-cell"],
        min_pixels_per_cell=args["min-pixels-per-cell"], confidence=confidences)

    @info "Building gene structure plot"
    gene_structure_plot = plot_gene_structure(df_spatial, gene_names)
    ## Plots

    n_tr_plot = plot_num_transcript_overview(df_spatial.gene, confidences, gene_names)
    noise_dist_plot = plot_noise_estimation_diagnostics(edge_lengths, confidences, d1, d2, confidence_nn_id=confidence_nn_id)

    @info "Plotting"

    open(args["output"], "w") do io
        print(io, """
            <!DOCTYPE html>
            <html>
            <head><style> body {background-color: #ffffff;}</style></head>

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
        show(io, MIME("text/html"), gc_plot)
        println(io, "<br>")

        println(io, "<hr>\n<h1 id=\"Noise_Level\">Noise level</h1><br>")
        show(io, MIME("text/html"), cc_plot)
        println(io, "<br>")
        show(io, MIME("text/html"), noise_dist_plot)
        println(io, "<br>")
        println(io, "Minimal noise level=$(round(100 * mean(confidences .< 0.01), sigdigits=3))%. ",
            "Expected noise level=$(round(100 * mean(1 .- confidences), sigdigits=2))%.")
        println(io, "<br>")

        println(io, "<hr>\n<h1 id=\"Gene_Structure\">Gene structure</h1><br>")
        show(io, MIME("text/html"), n_tr_plot)
        println(io, "<br>")
        show(io, MIME("text/html"), gene_structure_plot)
        println(io, "</body>")
        println(io, "</html>")
    end

    @info "All done!"

    close(log_file)

    return 0
end