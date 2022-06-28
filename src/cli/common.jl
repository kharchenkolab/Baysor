using DataFrames
using Statistics

import HDF5
import Dates
import LibGit2
import Pkg
import Pkg.TOML
import UUIDs

parse_toml_config(config::T where T <: AbstractString) =
    parse_toml_config(TOML.parsefile(config))

function get_default_config()
    return deepcopy(Dict{String, Any}(
        "Data" => Dict{String, Any}(
            "x-column" => "x",
            "y-column" => "y",
            "z-column" => "z",
            "gene-column" => "gene",
            "min-molecules-per-gene" => 1,
            "min-molecules-per-cell" => 3,
            "estimate-scale-from-centers" => true,
            "scale" => nothing,
            "scale-std" => "25%",
            "min-molecules-per-segment" => nothing
        ),
        "Sampling" => Dict{String, Any}(
            "new-component-weight" => 0.2,
            "new-component-fraction" => 0.3
        ),
        "Plotting" => Dict{String, Any}(
            "gene-composition-neigborhood" => nothing,
            "min-pixels-per-cell" => 15
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

get_run_id() = "$(UUIDs.uuid1())"[25:end] # must be independent of Random.seed

function setup_logger(prefix::String, file_name::String)
    log_file = open(append_suffix(prefix, file_name), "w")
    Base.CoreLogging.global_logger(DoubleLogger(log_file, stdout; force_flush=true))
    return log_file
end

function load_df(args::Dict; kwargs...)
    exc_genes = (args["exclude-genes"] === nothing) ? String[] : String.(strip.(Base.split(args["exclude-genes"], ",")))

    df_spatial, gene_names = load_df(args["coordinates"]; x_col=args["x-column"], y_col=args["y-column"], z_col=args["z-column"], gene_col=args["gene-column"],
        min_molecules_per_gene=args["min-molecules-per-gene"], exclude_genes=exc_genes, kwargs...)

    if :z in propertynames(df_spatial)
        if ("force-2d" in keys(args)) && args["force-2d"] || (length(unique(df_spatial.z)) < 2)
            select!(df_spatial, Not(:z))
        end
    elseif args["z-column"] != :z
        error("z-column $(args["z-column"]) not found in the data")
    end

    return df_spatial, gene_names
end

append_suffix(output::String, suffix) = "$(splitext(output)[1])_$suffix"

function get_baysor_run_str()::String
    pkg_str = "v$pkg_version"
    try
        pkg_info = Pkg.dependencies()[Base.UUID("cc9f9468-1fbe-11e9-0acf-e9460511877c")]
        repo = LibGit2.GitRepo(pkg_info.source) # In case it's a development version with a commit instead of a release
        hash_str = "$(LibGit2.GitShortHash(LibGit2.peel(LibGit2.GitCommit, LibGit2.head(repo))))"
        pkg_str = "$pkg_str [$hash_str]"
    catch err
    end

    return "($(Dates.Date(Dates.now()))) Run Baysor $pkg_str"
end

function save_matrix_to_loom(matrix; gene_names::Vector{String}, cell_names::Vector{String}, file_path::String,
        row_attrs::Union{Dict{String, T1}, Nothing} where T1=nothing, col_attrs::Union{Dict{String, T2}, Nothing} where T2=nothing)
    # Specification: https://linnarssonlab.org/loompy/format/index.html
    HDF5.h5open(file_path, "w") do fid
        fid["matrix", chunk=(64,64), compress=3] = matrix
        HDF5.create_group(fid, "row_attrs")
        HDF5.create_group(fid, "col_attrs")
        HDF5.create_group(fid, "attrs")
        fid["row_attrs"]["Name"] = gene_names
        fid["col_attrs"]["CellID"] = cell_names
        if row_attrs !== nothing
            for (k,v) in row_attrs
                fid["row_attrs"][k] = v
            end
        end

        if col_attrs !== nothing
            for (k,v) in col_attrs
                fid["col_attrs"][k] = v
            end
        end
    end;
end

run_cli(args::String) = run_cli(String.(Base.split(args)))

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

        if args[1] == "segfree"
            return run_cli_segfree(args[2:end])
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