using ArgParse
using DataFrames
using DataFramesMeta
using Distributed
using ProgressMeter
using Statistics

import CSV
import Pkg.TOML
import Plots

parse_toml_config(config::T where T <: AbstractString) =
    parse_toml_config(TOML.parsefile(config))

function parse_toml_config(config::Dict{AS, Any}) where AS <: AbstractString
    res_config = Dict{AS, Any}(
        "Data" => Dict{AS, Any}(
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
        "Sampling" => Dict{AS, Any}(
            "new-component-weight" => 0.2,
            "new-component-fraction" => 0.3,
            "center-component-weight" => 1.0,
            "n-degrees-of-freedom-center" => nothing
        ),
        "Plotting" => Dict{AS, Any}(
            "gene-composition-neigborhood" => 20,
            "plot-frame-size" => 5000
        )
    )

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

function parse_commandline(args::Union{Nothing, Array{String, 1}}=nothing) # TODO: add verbosity level
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--config", "-c"
            help = "TOML file with config"
        "--x-column", "-x"
            help = "Name of x column. Overrides the config value."
            default = "x"
        "--y-column", "-y"
            help = "Name of gene column. Overrides the config value."
            default = "y"
        "--gene-column"
            help = "Name of gene column. Overrides the config value."
            default = "gene"

        "--iters", "-i"
            help = "Number of iterations"
            arg_type = Int
            default = 500
        "--min-molecules-per-cell"
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
        cfg = parse_toml_config(r["config"])
        for sub_cfg in values(cfg)
            for (k, v) in sub_cfg
                if !(k in keys(r)) || r[k] === nothing
                    r[k] = v
                end
            end
        end
    else
        @warn "No config file provided. Back-up to default parameters."
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

load_df(args::Dict) = load_df(args["coordinates"]; x_col=args["x-column"], y_col=args["y-column"], gene_col=args["gene-column"], min_molecules_per_gene=args["min-molecules-per-gene"])

append_suffix(output::String, suffix) = "$(splitext(output)[1])_$suffix"

function plot_diagnostics_panel(df_res::DataFrame, assignment::Array{Int, 1}, tracer::Dict, args::Dict; margin=5*Plots.mm)
    @info "Plot diagnostics"
    open(append_suffix(args["output"], "diagnostics.html"), "w") do io
        # Convergence
        if ("n_components" in keys(tracer)) && length(tracer["n_components"]) != 0
            p_cov = plot_num_of_cells_per_iterarion(tracer, margin=margin);
            show(io, MIME("text/html"), p_cov)
        end

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
            for p in plots_clust
                show(io, MIME("text/html"), p)
            end
        end
    end
end

function run_cli(args::Union{Nothing, Array{String, 1}, String}=nothing)
    if args == "build" # need to call this function during build without any actual work
        return 0
    end

    args = parse_configs(args)

    # Dump parameters

    if args["config"] !== nothing
        dump_dst = abspath(append_suffix(args["output"], "config.toml"))
        if abspath(args["config"]) != dump_dst
            cp(args["config"], dump_dst, force=true)
        end
    end

    open(append_suffix(args["output"], "params.dump"), "w") do f
        println(f, "# CLI params: `$(join(ARGS, " "))`")
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

    if args["n-clusters"] > 1
        @info "Clustering molecules..."
        adjacent_points, adjacent_weights = build_molecule_graph(df_spatial, filter=false);
        mol_cluster_centers, cluster_per_molecule = cluster_molecules_on_mrf(df_spatial, adjacent_points, adjacent_weights; n_clusters=args["n-clusters"],
            nn_num=confidence_nn_id)[1:2] # TODO: use graph from confidence estimation

        df_spatial[!, :cluster] = cluster_per_molecule;

        # TODO: store mol_cluster_centers at bm_data.misc?
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
        plot_diagnostics_panel(segmentated_df, bm_data.assignment, bm_data.tracer, args)
        plot_transcript_assignment_panel(bm_data.x, bm_data.assignment, df_centers, args)
    end

    @info "All done!"

    close(log_file)

    return 0
end