using ArgParse
using Statistics

import Pkg.TOML

## Configs and parameters

function parse_commandline(args::Union{Nothing, Array{String, 1}}=nothing) # TODO: add verbosity level
    s = ArgParseSettings(prog="baysor run")
    @add_arg_table! s begin
        "--config", "-c"
            help = "TOML file with config"
        "--x-column", "-x"
            help = "Name of x column. Overrides the config value."
        "--y-column", "-y"
            help = "Name of y column. Overrides the config value."
        "--z-column", "-z"
            help = "Name of z column. Overrides the config value."
        "--gene-column", "-g"
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
        "--num-cells-init"
            help = "Initial number of cells."
            arg_type = Int
        "--output", "-o"
            help = "Name of the output file or path to the output directory"
            default = "segmentation.csv"
        "--plot", "-p"
            help = "Save pdf with plot of the segmentation"
            action = :store_true
        "--save-polygons"
            help = "Save estimated cell boundary polygons to a file with a specified FORMAT. Only 'GeoJSON' format is currently supported. The option requires setting '-p' to work."
            range_tester = (x -> in(lowercase(x), ["geojson"]))
            metavar = "FORMAT"
        "--force-2d"
            help = "Ignores z-column in the data if it is provided"
            action = :store_true
        "--exclude-genes"
            help = "Comma-separated list of genes or regular expressions to ignore during segmentation. Example: --exclude-genes='Blank*,MALAT1'"
        "--nuclei-genes"
            help = "Comma-separated list of nuclei-specific genes. If provided, `cyto-genes` has to be set, as well."
        "--cyto-genes"
            help = "Comma-separated list of cytoplasm-specific genes. If provided, `nuclei-genes` has to be set, as well."
        "--scale-std"
            help = "Standard deviation of scale across cells. Can be either number, which means absolute value of the std, or string ended with '%' to set it relative to scale. Default: 25%"
            default = "25%"

        "--scale", "-s"
            help = "Scale parameter, which suggest approximate cell radius for the algorithm. Must be in the same units as 'x' and 'y' molecule coordinates. Overrides the config value. Sets 'estimate-scale-from-centers' to false."
            arg_type = Float64
        "--prior-segmentation-confidence"
            help = "Confidence of the `prior_segmentation` results. Value in [0; 1]. If you want the final segmentation not contradicting to prior_segmentation, set it to 1. Otherwise, if you assume errors in prior_segmentation, values in [0.2-0.7] allow flexibility for the algorithm."
            arg_type = Float64
            default = 0.2

        "--no-ncv-estimation"
            dest_name = "estimate-ncvs"
            help = "Turns off neighborhood composition vectors estimation"
            action = :store_false

        "coordinates"
            help = "CSV file with coordinates of molecules and gene type"
            required = true
        "prior_segmentation"
            help = "Image or a MAT file with segmentation mask (either boolean or component indexing) or CSV column with integer segmentation labels. If it's the column name, it should be preceded ':' symbol (e.g. :cell)"
    end

    return (args === nothing) ? parse_args(s) : parse_args(args, s)
end

function parse_configs(args::Union{Nothing, Array{String, 1}}=nothing)
    r = parse_commandline(args)
    if r === nothing
        return nothing
    end

    if r["scale"] !== nothing
        r["estimate-scale-from-centers"] = false
    end

    if r["config"] !== nothing
        extend_params_with_config!(r, parse_toml_config(r["config"]))
    else
        @warn "No config file provided. Using default parameters."
        extend_params_with_config!(r, get_default_config())
    end

    for k in ["gene-column", "x-column", "y-column"]
        if !(k in keys(r)) || (r[k] === nothing) # should never be the case as we have defaults
            @warn "$k must be specified"
            return nothing
        end
        r[k] = Symbol(r[k])
    end

    r["z-column"] = Symbol(r["z-column"])

    if r["prior_segmentation"] === nothing && r["scale"] === nothing
        @warn "Either `prior_segmentation` or `scale` must be provided."
        return nothing
    end

    if r["min-molecules-per-segment"] === nothing
        r["min-molecules-per-segment"] = default_param_value(:min_molecules_per_segment, r["min-molecules-per-cell"])
    end

    if isdir(r["output"]) || isdirpath(r["output"])
        r["output"] = joinpath(r["output"], "segmentation.csv")
    end

    if (r["save-polygons"] !== nothing) && (!r["plot"])
        @warn "--plot option is required for saving polygons (--save-polygons). The polygons will not be saved."
    end

    if xor(r["nuclei-genes"] === nothing, r["cyto-genes"] === nothing)
        @warn "Only one of `nuclei-genes` and `cyto-genes` is provided. It has to be either both or none."
        return nothing
    end

    if (r["nuclei-genes"] !== nothing) && (r["n-clusters"] > 1)
        @warn "Setting n-clusters > 1 is not recommended with compartment-specific expression patterns (nuclei- and cyto-genes parameters)."
    end

    return r
end

function dump_parameters(args::Dict{String, Any}, args_str::String)
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
end

function load_prior_segmentation!(df_spatial, args::Dict{String, Any})
    prior_polygons = Matrix{Float64}[]
    if args["prior_segmentation"] !== nothing
        prior_seg_labels, scale, scale_std = DAT.load_prior_segmentation!(
            args["prior_segmentation"], df_spatial, BPR.position_data(df_spatial);
            min_mols_per_cell=args["min-molecules-per-cell"], min_molecules_per_segment=args["min-molecules-per-segment"]
        )

        if args["estimate-scale-from-centers"]
            args["scale"], args["scale-std"] = scale, scale_std
        end

        if (prior_seg_labels !== nothing) && args["plot"]
            @info "Estimating prior segmentation polygons..."
            prior_polygons = BPR.extract_polygons_from_label_grid(Matrix{UInt32}(prior_seg_labels[1:5:end, 1:5:end]); grid_step=5.0) # subset to save memory and time
            @info "Done"
        end
    end
    GC.gc()

    return prior_polygons
end

## CLI

function run_cli_main(args::Union{Nothing, Array{String, 1}}=nothing)
    Random.seed!(1)
    args_str = join(args, " ")
    args = parse_configs(args)
    (args !== nothing) || return 1

    dump_parameters(args, args_str)

    log_file = setup_logger(args["output"], "log.log")

    @info get_baysor_run_str()

    @info "Loading data..."
    df_spatial, gene_names = DAT.load_df(args, filter_cols=false)

    if args["gene-composition-neigborhood"] === nothing
        args["gene-composition-neigborhood"] = default_param_value(
            :composition_neighborhood, args["min-molecules-per-cell"], n_genes=length(gene_names)
        )
    end

    prior_polygons = load_prior_segmentation!(df_spatial, args)
    BPR.append_confidence!(df_spatial, args)

    (segmented_df, tracer, mol_clusts, comp_segs, poly_joint, cell_stat_df, cm, polygons) = BPR.run_segmentation(
        df_spatial, gene_names, args
    )

    @info "Saving results to $(args["output"])"

    out_paths = get_output_paths(args["output"])
    DAT.save_segmentation_results(
        segmented_df, cell_stat_df, cm, polygons, out_paths; poly_format=args["save-polygons"]
    )

    if args["plot"]
        @info "Plotting results"
        REP.plot_segmentation_report(
            segmented_df; tracer=tracer, clust_res=mol_clusts, comp_segs=comp_segs,
            prior_polygons=prior_polygons, polygons=poly_joint,
            diagnostic_file=out_paths.diagnostic_report, molecule_file=out_paths.molecule_plot,
            plot_transcripts=args["estimate-ncvs"],
            gene_colors=:ncv_color,
            min_molecules_per_cell=args["min-molecules-per-cell"], min_pixels_per_cell=args["min-pixels-per-cell"]
        )
    end

    @info "All done!"

    close(log_file)

    return 0
end
