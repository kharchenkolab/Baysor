using ArgParse
using DataFrames
using DataFramesMeta
using Distributed
using ProgressMeter
using Statistics

import CSV
import TOML
import Plots

parse_toml_config(config::T where T <: AbstractString) =
    parse_toml_config(TOML.parsefile(config))

function parse_toml_config(config::Dict{AbstractString, Any})
    res_config = Dict{AbstractString, Any}(
        "Data" => Dict{AbstractString, Any}(
            "x-column" => "x",
            "y-column" => "y",
            "gene-column" => "gene",
            "min-molecules-per-gene" => 1,
            "min-molecules-per-cell" => 3,
            "estimate-scale-from-centers" => true,
            "scale" => nothing,
            "min-center-size" => 200
        ),
        "Sampling" => Dict{AbstractString, Any}(
            "new-component-weight" => 0.2,
            "new-component-fraction" => 0.3,
            "center-component-weight" => 1.0,
            "n-degrees-of-freedom-center" => nothing,
            "shape-deg-freedom" => nothing       
        ),
        "Plotting" => Dict{AbstractString, Any}(
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

function parse_commandline(args::Union{Nothing, Array{String, 1}}=nothing)
    s = ArgParseSettings()
    @add_arg_table s begin
        "--config", "-c"
            help = "TOML file with config"
        "--x-column", "-x"
            help = "Name of x column. Overrides JSON value."
        "--y-column", "-y"
            help = "Name of gene column. Overrides JSON value."
        "--gene-column"
            help = "Name of gene column. Overrides JSON value."

        "--iters", "-i"
            help = "Number of iterations"
            arg_type = Int
            default = 500
        "--min-molecules-per-gene"
            help = "Minimal number of molecules per gene. Overrides JSON value."
            arg_type = Int
        "--n-frames", "-n"
            help = "Number of frames, which is the same as number of processes. Algorithm data is splitted by frames to allow parallel run over frames."
            arg_type = Int
            default=1
        "--num-cells-init"
            help = "Initial number of cells. Ignored if CSV with centers is provided. Overrides JSON value."
            arg_type = Int
        "--output", "-o"
            help = "Name of the output file or path to the output directory"
            default = "segmentation.csv"
        "--plot", "-p"
            help = "Save pdf with plot of the segmentation"
            action = :store_true
        "--refinement-iters"
            help = "Number of iterations for refinement of results. In most cases, default is enough."
            arg_type = Int
            default = 50

        "--scale", "-s"
            help = "Scale parameter, which suggest approximate cell radius for the algorithm. Overrides JSON value."
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
        if !(k in keys(r)) || (r[k] === nothing)
            error("$k must be specified")
        end
        r[k] = Symbol(r[k])
    end

    if r["shape-deg-freedom"] === nothing
        r["shape-deg-freedom"] = default_param_value(:shape_deg_freedom, r["min-molecules-per-cell"])
    end

    if r["n-degrees-of-freedom-center"] === nothing
        r["n-degrees-of-freedom-center"] = default_param_value(:n_degrees_of_freedom_center, r["min-molecules-per-cell"])
    end

    if r["centers"] === nothing && r["scale"] === nothing
        print("Either `centers` or `scale` must be provided.\n" * usage_string(s) * "\n")
        exit(1)
    end

    if isdir(r["output"]) || isdirpath(r["output"])
        r["output"] = joinpath(r["output"], "segmentation.csv")
    end

    return r
end

load_df(args::Dict) = load_df(args["coordinates"]; x_col=args["x-column"], y_col=args["y-column"], gene_col=args["gene-column"], min_molecules_per_gene=args["min-molecules-per-gene"])

append_suffix(output::String, suffix) = "$(splitext(output)[1])_$suffix"

function plot_results(df_res::DataFrame, assignment::Array{Int, 1}, df_centers::Union{DataFrame, Nothing}, tracer::Dict, args::Dict{String,T} where T)
    # Convergence
    p_cov = plot_num_of_cells_per_iterarion(tracer);
    Plots.savefig(append_suffix(args["output"], "convergence.pdf"))

    # Transcripts
    plots = plot_cell_boundary_polygons_all(df_res, assignment, df_centers; gene_composition_neigborhood=args["gene-composition-neigborhood"], 
        frame_size=args["plot-frame-size"], min_molecules_per_cell=args["min-molecules-per-cell"])

    open(append_suffix(args["output"], "borders.html"), "w") do io
        for (i,p) in enumerate(plots)
            show(io, MIME("text/html"), p)
        end
    end

    @info "Done!"
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
        TOML.print(f, Dict(k => (v === nothing) ? ((typeof(v) === Symbol) ? String(v) : "") : v for (k,v) in args))
    end

    # Run algorithm

    @info "Run"
    @info "Load data..."
    df_spatial, gene_names = load_df(args)

    df_centers = nothing
    bm_data_arr = BmmData[]
    confidence_nn_id = max(div(args["min-molecules-per-cell"], 2) + 1, 3)

    if args["centers"] !== nothing
        centers = load_centers(args["centers"], x_col=args["x-column"], y_col=args["y-column"], min_segment_size=args["min-center-size"])
        df_centers = centers.centers

        bm_data_arr = initial_distribution_arr(df_spatial, centers; n_frames=args["n-frames"],
            shape_deg_freedom=args["shape-deg-freedom"], scale=(args["estimate-scale-from-centers"] ? nothing : args["scale"]), n_cells_init=args["num-cells-init"],
            new_component_weight=args["new-component-weight"], center_component_weight=args["center-component-weight"], 
            n_degrees_of_freedom_center=args["n-degrees-of-freedom-center"], min_molecules_per_cell=args["min-molecules-per-cell"], confidence_nn_id=confidence_nn_id);
    else
        bm_data_arr = initial_distribution_arr(df_spatial; n_frames=args["n-frames"],
            shape_deg_freedom=args["shape-deg-freedom"], scale=args["scale"], n_cells_init=args["num-cells-init"],
            new_component_weight=args["new-component-weight"], min_molecules_per_cell=args["min-molecules-per-cell"], confidence_nn_id=confidence_nn_id);
    end

    if length(bm_data_arr) > 1
        addprocs(length(bm_data_arr) - nprocs())
        eval(:(@everywhere using Baysor))
    end

    bm_data = run_bmm_parallel(bm_data_arr, args["iters"], new_component_frac=args["new-component-fraction"],
                               min_molecules_per_cell=args["min-molecules-per-cell"], n_refinement_iters=args["refinement-iters"]);

    @info "Processing complete."

    # Save results

    segmentated_df = get_segmentation_df(bm_data, gene_names)
    cell_stat_df = get_cell_stat_df(bm_data; add_qc=true)

    @info "Save data to $(args["output"])"
    CSV.write(args["output"], segmentated_df);
    CSV.write(append_suffix(args["output"], "cell_stats.csv"), cell_stat_df);

    if args["plot"]
        @info "Plot results"
        plot_results(bm_data.x, bm_data.assignment, df_centers, bm_data.tracer, args)
    end

    @info "All done!"

    return 0
end