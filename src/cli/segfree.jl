using ArgParse
using Statistics

function parse_segfree_commandline(args::Union{Nothing, Array{String, 1}}=nothing)
    s = ArgParseSettings(prog="baysor preview")
    @add_arg_table! s begin
        "--k", "-k"
            help = "Number of neighbors for segmentation-free pseudo-cells. It's suggested to set it to the expected minimal number of molecules per cell."
        "-n"
            help = "Number of NCVs to save. Default: 0 (all of them). If set, selects NCVs uniformly across the first 2 PCs"
            default = 0
            arg_type = Int
        "--min-molecules-per-cell", "-m"
            help = "Minimal number of molecules for a cell to be considered as real. Only used to estimate other parameters."
            arg_type = Int
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
        "--exclude-genes"
            help = "Comma-separated list of genes to ignore during segmentation"
        "--output", "-o"
            help = "Name of the output file or path to the output directory"
            default = "ncvs.loom"

        "coordinates"
            help = "CSV file with coordinates of transcripts and gene type"
            required = true
    end

    return (args === nothing) ? parse_args(s) : parse_args(args, s)
end

function parse_segfree_configs(args::Union{Nothing, Array{String, 1}}=nothing)
    r = parse_segfree_commandline(args)
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
        r["output"] = joinpath(r["output"], "ncvs.loom")
    end

    return r
end

function run_cli_segfree(args::Union{Nothing, Array{String, 1}}=nothing)
    Random.seed!(1)
    args_str = join(args, " ")
    args = parse_segfree_configs(args)
    (args !== nothing) || return 1

    log_file = setup_logger(args["output"], "segfree_log.log")

    run_id = get_run_id()
    @info "Run $run_id"
    @info "# CLI params: `$args_str`"
    @info get_baysor_run_str()

    # Run NCV estimates

    @info "Loading data..."
    args["min-molecules-per-gene"] = 0
    df_spatial, gene_names = load_df(args)

    @info "Loaded $(size(df_spatial, 1)) transcripts"
    if size(df_spatial, 1) != size(unique(df_spatial), 1)
        @warn "$(size(df_spatial, 1) - size(unique(df_spatial), 1)) records are duplicates. You may need to filter them beforehand."
    end

    if args["k"] === nothing
        if args["min-molecules-per-cell"] === nothing
            error("Either `min-molecules-per-cell` or `k` must be provided")
        end
        args["k"] = default_param_value(:composition_neighborhood, args["min-molecules-per-cell"], n_genes=length(gene_names))
    end

    @info "Estimating neighborhoods..."
    neighb_cm = BPR.neighborhood_count_matrix(df_spatial, args["k"]);

    @info "Estimating molecule confidences..."
    confidence_nn_id = default_param_value(:confidence_nn_id, args["min-molecules-per-cell"])
    BPR.append_confidence!(df_spatial, nn_id=confidence_nn_id)

    @info "Estimating gene colors..."
    transformation = BPR.gene_composition_transformation(neighb_cm, df_spatial.confidence)
    gene_colors = BPR.gene_composition_colors(neighb_cm, transformation)
    gene_colors = "#" .* Colors.hex.(gene_colors)

    @info "Saving results..."
    DAT.save_matrix_to_loom(
        neighb_cm; gene_names=gene_names, cell_names=["$run_id-$i" for i in 1:size(neighb_cm, 1)],
        col_attrs=Dict("ncv_color" => gene_colors), file_path=args["output"]
    )

    @info "Done!"
    close(log_file)

    return 0
end