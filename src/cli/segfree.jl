using ArgParse
using DataFrames
using Statistics

import Colors

function parse_segfree_commandline(args::Union{Nothing, Array{String, 1}}=nothing)
    s = ArgParseSettings(prog="baysor preview")
    @add_arg_table! s begin
        "--k", "-k"
            help = "Number of neighbors for segmentation-free pseudo-cells. It's suggested to set it to the expected minimal number of molecules per cell."
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
            default = "preview.html"

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

    return r
end

function run_cli_segfree(args::Union{Nothing, Array{String, 1}}=nothing)
    args = parse_segfree_configs(args)

    # Set up logger

    log_file = open(append_suffix(args["output"], "segfree_log.log"), "w")
    Base.CoreLogging.global_logger(DoubleLogger(log_file, stdout; force_flush=true))

    # Run preview

    @info get_baysor_run_str()
    @info "Loading data..."
    args["min-molecules-per-gene"] = 0
    df_spatial, gene_names = load_df(args)

    @info "Loaded $(size(df_spatial, 1)) transcripts"
    if size(df_spatial, 1) != size(unique(df_spatial), 1)
        @warn "$(size(df_spatial, 1) - size(unique(df_spatial), 1)) records are duplicates. You may need to filter them beforehand."
    end

    @error "NOT IMPLEMENTED" # TODO

    close(log_file)

    return 0
end