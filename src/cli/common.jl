using Statistics

import Dates
import LibGit2
import Pkg
import UUIDs

function default_param_value(param::Symbol, min_molecules_per_cell::Union{Int, Nothing};
                             n_molecules::Union{Int, Nothing}=nothing, n_genes::Union{Int, Nothing}=nothing)
    if min_molecules_per_cell === nothing
        error("Either `$param` or `min_molecules_per_cell` must be provided")
    end

    min_molecules_per_cell = max(min_molecules_per_cell, 3)

    if param == :min_molecules_per_segment
        return max(round(Int, min_molecules_per_cell / 4), 2)
    end

    if param == :confidence_nn_id
        return max(div(min_molecules_per_cell, 2) + 1, 5)
    end

    if param == :composition_neighborhood
        if n_genes === nothing
            return max(min_molecules_per_cell, 3)
        end

        return max(div(n_genes, 10), min_molecules_per_cell, 3)
    end

    if param == :n_gene_pcs
        (n_genes !== nothing) || error("Either `$param` or `n_genes` must be provided")
        return min(max(div(n_genes, 3), 30), 100, n_genes)
    end

    if param == :n_cells_init
        (n_molecules !== nothing) || error("Either `$param` or `n_molecules` must be provided")
        return div(n_molecules, min_molecules_per_cell) * 2
    end
end

get_run_id() = "R" * "$(UUIDs.uuid1())"[25:33] # must be independent of Random.seed

function setup_logger(prefix::String, file_name::String)
    log_file = open(append_suffix(prefix, file_name), "w")
    Base.CoreLogging.global_logger(DoubleLogger(log_file, stdout; force_flush=true))
    return log_file
end

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

function fill_and_check_options!(opts::DataOptions)
    opts.min_molecules_per_cell > 0 || cmd_error("`min_molecules_per_cell` must be positive")

    if opts.min_molecules_per_segment == 0
        opts.min_molecules_per_segment = default_param_value(:min_molecules_per_segment, opts.min_molecules_per_cell)
    end

    if opts.confidence_nn_id == 0
        opts.confidence_nn_id = default_param_value(:confidence_nn_id, opts.min_molecules_per_cell)
    end

    return opts
end

function fill_and_check_options!(opts::PlottingOptions, min_molecules_per_cell::Int, n_genes::Union{Int, Nothing}=nothing)
    if opts.gene_composition_neigborhood <= 0
        opts.gene_composition_neigborhood = default_param_value(
            :composition_neighborhood, min_molecules_per_cell, n_genes=n_genes
        )
    end

    (opts.ncv_method in ["ri", "dense", "sparse"]) || error("`ncv_method` must be one of \"ri\", \"dense\", or \"sparse\"")

    return opts
end

Comonicon.from_dict(::Type{DataOptions}, ::Type{Symbol}, s) = Symbol(s)
Comonicon.to_dict(::Type, s::Symbol) = "$s"
Base.tryparse(::Type{Symbol}, s::AbstractString) = Symbol(s)