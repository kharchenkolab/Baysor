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
        r["min-molecules-per-segment"] = BPR.default_param_value(:min_molecules_per_segment, r["min-molecules-per-cell"])
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

## Utils

struct OutputPaths
    # The class is here for clarity of what files are saved
    segmented_df::String;
    cell_stats::String;
    counts::String;
    diagnostic_report::String;
    molecule_plot::String;
    polygons::String;
end
OutputPaths(;segmented_df::String, cell_stats::String, counts::String, diagnostic_report::String, molecule_plot::String, polygons::String) =
    OutputPaths(segmented_df, cell_stats, counts, diagnostic_report, molecule_plot, polygons)

function load_and_preprocess_data!(args::Dict{String, Any})
    @info "Loading data..."
    df_spatial, gene_names = load_df(args, filter_cols=false)
    df_spatial[!, :molecule_id] = 1:size(df_spatial, 1)

    @info "Loaded $(size(df_spatial, 1)) transcripts"
    if size(df_spatial, 1) != size(unique(df_spatial), 1)
        @warn "$(size(df_spatial, 1) - size(unique(df_spatial), 1)) records are duplicates. You may need to filter them beforehand."
    end

    if args["gene-composition-neigborhood"] === nothing
        args["gene-composition-neigborhood"] = BPR.default_param_value(:composition_neighborhood, args["min-molecules-per-cell"], n_genes=length(gene_names))
    end

    prior_polygons = Matrix{Float64}[]
    if args["prior_segmentation"] !== nothing
        df_spatial[!, :prior_segmentation], prior_polygons, args["scale"], args["scale-std"] = DAT.load_prior_segmentation(df_spatial, args)
    end
    GC.gc()

    confidence_nn_id = BPR.default_param_value(:confidence_nn_id, args["min-molecules-per-cell"])

    @info "Estimating noise level"
    prior_seg = (args["prior_segmentation"]===nothing) ? nothing : df_spatial.prior_segmentation
    BPR.append_confidence!(df_spatial, prior_seg, nn_id=confidence_nn_id, prior_confidence=args["prior-segmentation-confidence"])
    @info "Done"

    return df_spatial, gene_names, prior_polygons
end

function get_output_paths(segmented_df_path::String)
    return OutputPaths(;
        segmented_df=segmented_df_path,
        Dict(k => append_suffix(segmented_df_path, v) for (k,v) in [
            :cell_stats => "cell_stats.csv",
            :counts => "counts.tsv",
            :diagnostic_report => "diagnostics.html",
            :molecule_plot => "borders.html",
            :polygons => "polygons.json"
        ])...
    )
end

## CLI

function run_cli_main(args::Union{Nothing, Array{String, 1}}=nothing)
    Random.seed!(1)
    args_str = join(args, " ")
    args = parse_configs(args)
    (args !== nothing) || return 1

    dump_parameters(args, args_str)

    log_file = setup_logger(args["output"], "log.log")

    # run_id = get_run_id()
    # @info "Run $run_id"
    # TODO: add run_id to cell ids

    @info get_baysor_run_str()

    # Load data
    df_spatial, gene_names, prior_polygons = load_and_preprocess_data!(args)

    # Region-based segmentations

    comp_segs, comp_genes = nothing, Vector{Int}[]
    adjacent_points, adjacent_weights = BPR.build_molecule_graph(df_spatial; use_local_gene_similarities=false, adjacency_type=:triangulation)[1:2]
    if args["nuclei-genes"] !== nothing
        comp_segs, comp_genes, df_spatial[!, :compartment] = BPR.estimate_molecule_compartments(
            df_spatial, gene_names; nuclei_genes=args["nuclei-genes"], cyto_genes=args["cyto-genes"], scale=args["scale"]
        )
        df_spatial[!, :nuclei_probs] = 1 .- comp_segs.assignment_probs[2,:];

        adjacent_points, adjacent_weights = BPR.adjust_mrf_with_compartments(
            adjacent_points, adjacent_weights,
            comp_segs.assignment_probs[1,:], comp_segs.assignment_probs[2,:]
        );
    end

    mol_clusts = nothing
    if args["n-clusters"] > 1
        mol_clusts = BPR.estimate_molecule_clusters(df_spatial, args["n-clusters"])
        df_spatial[!, :cluster] = mol_clusts.assignment;
    end

    # Cell segmentation

    bm_data = BPR.initialize_bmm_data(
        df_spatial; scale=args["scale"], scale_std=args["scale-std"], n_cells_init=args["num-cells-init"],
        prior_seg_confidence=args["prior-segmentation-confidence"], min_molecules_per_cell=args["min-molecules-per-cell"],
        confidence_nn_id=0, adjacent_points=adjacent_points, adjacent_weights=adjacent_weights,
        na_genes=Vector{Int}(vcat(comp_genes...))
    );

    @info "Using $(size(BPR.position_data(bm_data), 1))D coordinates"

    history_depth = round(Int, args["iters"] * 0.1)
    bm_data = BPR.bmm!(
        bm_data; n_iters=args["iters"],
        new_component_frac=args["new-component-fraction"], new_component_weight=args["new-component-weight"],
        min_molecules_per_cell=args["min-molecules-per-cell"], assignment_history_depth=history_depth
    );

    @info "Processing complete."

    # Save results

    ## Extract results

    segmented_df, cell_stat_df, cm = BPR.get_segmentation_results(bm_data, gene_names)
    gene_colors = nothing
    if args["estimate-ncvs"]
        @info "Estimating local colors"
        gene_colors = BPR.gene_composition_colors(bm_data.x, args["gene-composition-neigborhood"])
        segmented_df[!, :ncv_color] = "#" .* Colors.hex.(gene_colors)
    end

    ## Save results

    @info "Saving results to $(args["output"])"
    out_paths = get_output_paths(args["output"])

    DAT.save_segmented_df(segmented_df, out_paths.segmented_df);
    DAT.save_cell_stat_df(cell_stat_df, out_paths.cell_stats);
    DAT.save_molecule_counts(cm, out_paths.counts)

    poly_joint, polygons = nothing, nothing
    if args["save-polygons"] !== nothing || args["plot"] && args["estimate-ncvs"]
        poly_joint, polygons = BPR.boundary_polygons_auto(
            BPR.position_data(bm_data), bm_data.assignment; scale=args["scale"], min_pixels_per_cell=args["min-pixels-per-cell"],
            estimate_per_z=(args["save-polygons"] !== nothing)
        )
    end

    if args["save-polygons"] !== nothing
        DAT.save_polygons(polygons; format=args["save-polygons"], file=out_paths.polygons)
    end

    # Plot report

    if args["plot"]
        REP.plot_segmentation_report(
            segmented_df; tracer=bm_data.tracer, plot_transcripts=args["estimate-ncvs"],
            gene_colors=gene_colors, prior_polygons=prior_polygons, polygons=poly_joint,
            min_molecules_per_cell=args["min-molecules-per-cell"], min_pixels_per_cell=args["min-pixels-per-cell"],
            diagnostic_file=out_paths.diagnostic_report, molecule_file=out_paths.molecule_plot
        )
    end

    # Finish!

    @info "All done!"

    close(log_file)

    return 0
end
