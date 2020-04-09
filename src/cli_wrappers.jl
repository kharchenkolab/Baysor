using ArgParse
using DataFrames
using DataFramesMeta
using Distributed
using ProgressMeter
using Statistics

import CSV
import Pkg.TOML
import Plots

## Common

parse_toml_config(config::T where T <: AbstractString) =
    parse_toml_config(TOML.parsefile(config))

function get_default_config()
    return deepcopy(Dict{String, Any}(
        "Data" => Dict{String, Any}(
            "x-column" => "x",
            "y-column" => "y",
            "gene-column" => "gene",
            "min-molecules-per-gene" => 1,
            "min-molecules-per-cell" => 3,
            "estimate-scale-from-centers" => true,
            "scale" => nothing,
            "scale-std" => "25%",
            "min-center-size" => 200
        ),
        "Sampling" => Dict{String, Any}(
            "new-component-weight" => 0.2,
            "new-component-fraction" => 0.3,
            "center-component-weight" => 1.0,
            "n-degrees-of-freedom-center" => nothing
        ),
        "Plotting" => Dict{String, Any}(
            "gene-composition-neigborhood" => 20, # TODO: should it be nothing to use default estimation?
            "plot-frame-size" => 5000
        )
    ))
end

function parse_toml_config(config::Dict{AS, Any}) where AS <: AbstractString
    res_config = get_default_config()
    for (k,v) in config
        if !(k in keys(res_config))
            error("Unexpected value in the config: '$k'")
        end

        cur_def = res_config[k]

        for (k2,v2) in v
            if !(k2 in keys(cur_def))
                error("Unexpected value in the config: '$k' -> '$k2'")
            end

            cur_def[k2] = v2
        end
    end

    return res_config
end

function extend_params_with_config!(params::Dict, config::Dict)
    for sub_cfg in values(config)
        for (k, v) in sub_cfg
            if !(k in keys(params)) || params[k] === nothing
                params[k] = v
            end
        end
    end
end

load_df(args::Dict) = load_df(args["coordinates"]; x_col=args["x-column"], y_col=args["y-column"], gene_col=args["gene-column"], min_molecules_per_gene=args["min-molecules-per-gene"])

append_suffix(output::String, suffix) = "$(splitext(output)[1])_$suffix"

## Main run

function parse_commandline(args::Union{Nothing, Array{String, 1}}=nothing) # TODO: add verbosity level
    s = ArgParseSettings(prog="baysor run")
    @add_arg_table! s begin
        "--config", "-c"
            help = "TOML file with config"
        "--x-column", "-x"
            help = "Name of x column. Overrides the config value."
        "--y-column", "-y"
            help = "Name of gene column. Overrides the config value."
        "--gene-column"
            help = "Name of gene column. Overrides the config value."

        "--iters", "-i"
            help = "Number of iterations"
            arg_type = Int
            default = 500
        "--min-molecules-per-cell", "-m"
            help = "Minimal number of molecules for a cell to be considered as real. It's an important parameter, as it's used to infer several other parameters. Overrides the config value."
            arg_type = Int
        "--n-clusters"
            help = "Number of molecule clusters, i.e. major cell types. Depends on protocol resolution, but should not be too high. In most cases something between 3 and 15 should work well."
            arg_type = Int
            default=4
        "--n-frames", "-n"
            help = "Number of frames, which is the same as number of processes. Algorithm data is splitted by frames to allow parallel run over frames."
            arg_type = Int
            default=1
        "--num-cells-init"
            help = "Initial number of cells."
            arg_type = Int
        "--output", "-o"
            help = "Name of the output file or path to the output directory"
            default = "segmentation.csv"
        "--plot", "-p"
            help = "Save pdf with plot of the segmentation"
            action = :store_true

        "--scale", "-s"
            help = "Scale parameter, which suggest approximate cell radius for the algorithm. Overrides the config value. Set 'estimate-scale-from-centers' to false."
            arg_type = Float64

        "coordinates"
            help = "CSV file with coordinates of transcripts and gene type"
            required = true
        "centers"
            help = "Either CSV file with coordinates of cell centers, extracted from DAPI staining or image with segmentation mask (either boolean or component indexing)"
    end

    return (args === nothing) ? parse_args(s) : parse_args(args, s)
end

function parse_configs(args::Union{Nothing, Array{String, 1}}=nothing)
    r = parse_commandline(args)
    if r["scale"] !== nothing
        r["estimate-scale-from-centers"] = false
    end

    if r["config"] !== nothing
        extend_params_with_config!(r, parse_toml_config(r["config"]))
    else
        @warn "No config file provided. Back-up to default parameters."
        extend_params_with_config!(r, get_default_config())
    end

    for k in ["gene-column", "x-column", "y-column"]
        if !(k in keys(r)) || (r[k] === nothing) # should never be the case as we have defaults
            error("$k must be specified")
        end
        r[k] = Symbol(r[k])
    end

    if r["centers"] === nothing && r["scale"] === nothing
        println("Either `centers` or `scale` must be provided.")
        exit(1)
    end

    if r["n-degrees-of-freedom-center"] === nothing
        r["n-degrees-of-freedom-center"] = default_param_value(:n_degrees_of_freedom_center, r["min-molecules-per-cell"])
    end

    if isdir(r["output"]) || isdirpath(r["output"])
        r["output"] = joinpath(r["output"], "segmentation.csv")
    end

    return r
end

function plot_diagnostics_panel(df_res::DataFrame, assignment::Array{Int, 1}, tracer::Dict, args::Dict; margin=5*Plots.mm,
        max_diffs::Union{Vector{Float64}, Nothing}=nothing, change_fracs::Union{Vector{Float64}, Nothing}=nothing)
    @info "Plot diagnostics"
    open(append_suffix(args["output"], "diagnostics.html"), "w") do io
        # Molecule clustering convergence
        if max_diffs !== nothing
            p_mol_conv = Plots.plot(max_diffs, xlabel="Iteration", ylabel="Maximal probability change", title="Molcecule clustering convergence", label="Maximal change");
            if change_fracs !== nothing
                p_mol_conv = Plots.plot!(change_fracs, label="Fraction of molecules changed")
            end
            show(io, MIME("text/html"), p_mol_conv)
        end

        # Main algorithm convergence
        if (:n_components in keys(tracer)) && length(tracer[:n_components]) != 0
            p_conv = plot_num_of_cells_per_iterarion(tracer);
            show(io, MIME("text/html"), p_conv)
        end

        println(io, "<br><br>")

        # Confidence per molecule
        if :confidence in names(df_res)
            bins = 0.0:0.025:1.0
            p_conf = Plots.histogram(df_res.confidence[assignment .!= 0], bins=bins, label="Assigned molecules",
                xlabel="Confidence", ylabel="#Molecules", title="Confidence per molecule", margin=margin)
            p_conf = Plots.histogram!(df_res.confidence[assignment .== 0], alpha=0.5, bins=bins, label="Noise molecules")
            show(io, MIME("text/html"), p_conf)
        end

        # Assignment confidence
        if :assignment_confidence in names(df_res)
            p_conf = Plots.histogram(df_res.assignment_confidence[assignment .> 0], bins=50, legend=false)
            p_conf = Plots.vline!([0.95], xlabel="Assignment confidence", ylabel="#Molecules", xlims=(-0.01, 1.03),
                title="Assignment confidence per real molecules")
            show(io, MIME("text/html"), p_conf)
        end

        println(io, "<br><br>")

        # Num. of molecules per cell
        n_mols_per_cell = count_array(assignment .+ 1)[2:end]
        p_n_mols = Plots.histogram(n_mols_per_cell[(n_mols_per_cell .> 1) .& (n_mols_per_cell .< quantile(n_mols_per_cell, 0.99) / 0.99)],
            title="Num. molecules per cell", xlabel="Num. molecules per cell", ylabel="Num. cells")
        show(io, MIME("text/html"), p_n_mols)
    end
end

function plot_transcript_assignment_panel(df_res::DataFrame, assignment::Array{Int, 1}, df_centers::Union{DataFrame, Nothing}, args::Dict;
                                          plot_width::Int=800, margin=5*Plots.mm)
    @info "Plot transcript assignment"
    plots_col, plots_clust = plot_cell_boundary_polygons_all(df_res, assignment, df_centers; gene_composition_neigborhood=args["gene-composition-neigborhood"],
        frame_size=args["plot-frame-size"], min_molecules_per_cell=args["min-molecules-per-cell"], plot_width=plot_width, margin=margin)

    open(append_suffix(args["output"], "borders.html"), "w") do io
        for p in plots_col
            show(io, MIME("text/html"), p)
        end

        if plots_clust !== nothing
            println(io, "<br>")
            for p in plots_clust
                show(io, MIME("text/html"), p)
            end
        end
    end
end

function run_cli_main(args::Union{Nothing, Array{String, 1}}=nothing)
    args_str = join(args, " ")
    args = parse_configs(args)

    # Dump parameters

    if args["config"] !== nothing
        dump_dst = abspath(append_suffix(args["output"], "config.toml"))
        if abspath(args["config"]) != dump_dst
            cp(args["config"], dump_dst, force=true)
        end
    end

    open(append_suffix(args["output"], "params.dump"), "w") do f
        println(f, "# CLI params: `$args_str`")
        TOML.print(f, Dict(k => (v !== nothing) ? ((typeof(v) === Symbol) ? String(v) : v) : "" for (k,v) in args))
    end

    # Set up logger

    log_file = open(append_suffix(args["output"], "log.log"), "w")
    Base.CoreLogging.global_logger(DoubleLogger(log_file, stdout; force_flush=true))

    # Set up plotting

    ENV["GKSwstype"] = "100"; # Disable output device
    Plots.gr()

    # Run algorithm

    @info "Run"
    @info "Load data..."
    df_spatial, gene_names = load_df(args)

    df_centers = nothing
    bm_data_arr = BmmData[]
    confidence_nn_id = default_param_value(:confidence_nn_id, args["min-molecules-per-cell"])

    @info "Estimating noise level"
    append_confidence!(df_spatial, nn_id=confidence_nn_id) # TODO: use segmentation mask if available here
    @info "Done"

    max_diffs, change_fracs = nothing, nothing
    if args["n-clusters"] > 1
        @info "Clustering molecules..."
        adjacent_points, adjacent_weights = build_molecule_graph(df_spatial, filter=false);
        for i in 1:length(adjacent_weights)
            cur_points = adjacent_points[i]
            cur_weights = adjacent_weights[i]
            for j in 1:length(cur_weights)
                cur_weights[j] *= df_spatial.confidence[cur_points[j]]
            end
        end

        mol_clusts = cluster_molecules_on_mrf(df_spatial, adjacent_points, adjacent_weights; n_clusters=args["n-clusters"], weights_pre_adjusted=true)

        df_spatial[!, :cluster] = mol_clusts.assignment;
        max_diffs, change_fracs = mol_clusts.diffs, mol_clusts.change_fracs

        # TODO: store mol_clusts.exprs at bm_data.misc?
        @info "Done"
    end

    if args["centers"] !== nothing
        centers = load_centers(args["centers"], x_col=args["x-column"], y_col=args["y-column"], min_segment_size=args["min-center-size"])
        df_centers = centers.centers

        scale = args["estimate-scale-from-centers"] ? nothing : args["scale"]
        scale_std = args["estimate-scale-from-centers"] ? nothing : args["scale-std"]
        bm_data_arr = initial_distribution_arr(df_spatial, centers; n_frames=args["n-frames"], scale=scale, scale_std=scale_std,
            n_cells_init=args["num-cells-init"], new_component_weight=args["new-component-weight"],
            center_component_weight=args["center-component-weight"], n_degrees_of_freedom_center=args["n-degrees-of-freedom-center"],
            min_molecules_per_cell=args["min-molecules-per-cell"], confidence_nn_id=confidence_nn_id);
    else
        bm_data_arr = initial_distribution_arr(df_spatial; n_frames=args["n-frames"], scale=args["scale"], scale_std=args["scale-std"],
            n_cells_init=args["num-cells-init"], new_component_weight=args["new-component-weight"],
            min_molecules_per_cell=args["min-molecules-per-cell"], confidence_nn_id=confidence_nn_id);
    end

    history_depth = round(Int, args["iters"] * (args["iters"] >= 1000 ? 0.05 : 0.01))
    bm_data = run_bmm_parallel!(bm_data_arr, args["iters"], new_component_frac=args["new-component-fraction"],
                                min_molecules_per_cell=args["min-molecules-per-cell"], assignment_history_depth=history_depth);

    @info "Processing complete."

    # Save results

    segmentated_df = get_segmentation_df(bm_data, gene_names)
    cell_stat_df = get_cell_stat_df(bm_data; add_qc=true)

    @info "Save data to $(args["output"])"
    CSV.write(args["output"], segmentated_df);
    CSV.write(append_suffix(args["output"], "cell_stats.csv"), cell_stat_df);

    if args["plot"]
        plot_diagnostics_panel(segmentated_df, bm_data.assignment, bm_data.tracer, args; max_diffs=max_diffs, change_fracs=change_fracs)
        plot_transcript_assignment_panel(bm_data.x, bm_data.assignment, df_centers, args)
    end

    @info "All done!"

    close(log_file)

    return 0
end

## Pre-view


function parse_preview_commandline(args::Union{Nothing, Array{String, 1}}=nothing)
    s = ArgParseSettings(prog="baysor preview")
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
        # TODO: add verbosity level
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

function run_cli_preview(args::Union{Nothing, Array{String, 1}}=nothing)
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
    color_transformation = gene_composition_transformation(neighb_cm, confidences);

    @info "Estimating local colors"
    gene_colors = gene_composition_colors(neighb_cm, color_transformation)

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

## All

function run_cli(args::Vector{String}=ARGS)::Cint
    help_message = "Usage: baysor <command> [options]\n\nCommands:\n\trun\t\trun segmentation of the dataset\n\tpreview\t\tgenerate preview diagnostics of the dataset\n"

    debug = false
    if "--debug" in args
        args = args[args .!= "--debug"]
        debug = true
    end

    try
        if (length(args) == 0) || (length(args) == 1) && (args[1] == "-h" || args[1] == "--help")
            println(help_message)
            return 0
        end

        if args[1] == "run"
            return run_cli_main(args[2:end])
        end

        if args[1] == "preview"
            return run_cli_preview(args[2:end])
        end

        @error "Can't parse argument $(args[1])"
        println(help_message)
        return 1
    catch err
        if debug
            rethrow()
        else
            @error("$err\n\n" * join(["$s" for s in stacktrace(catch_backtrace())], "\n"))
        end
    end

    return 2
end

julia_main()::Cint = run_cli()