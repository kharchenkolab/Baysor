module Baysor

include("LazySubmodules.jl")
using .LazySubmodules

# export
#     BmmData, bmm!,
#     build_molecule_graph, load_df, initialize_bmm_data,
#     boundary_polygons,
#     neighborhood_count_matrix, gene_composition_transformation, gene_composition_colors,
#     append_confidence!,
#     plot_molecules, plot_molecules!, plot_expression_vectors,
#     cluster_molecules_on_mrf

export load_module

LazySubmodules.__init__() # Somehow without it __init__ is only called after all @lazy_submodule macroses

# Utils: Minimal functions with zero compilation time shared across submodules
include("utils/Utils.jl")

# DataLoading: work with different input and output formats. Isolated as it requires additional file-type specific libraries
@lazy_submodule DAT = "data_loading/DataLoading.jl"

# Reporting: all plotting code isolated, as plotting libraries are the slowest in Julia
@lazy_submodule REP = "reporting/Reporting.jl"

# Processing: all the actual logic
@lazy_submodule BPR = "processing/Processing.jl"

# CLI: minimal CLI code required for fast responsive CLI
@lazy_submodule CLI = "cli/CLI.jl"

command_main(args...; kwargs...) = CLI.command_main(args...; kwargs...)
julia_main(args...; kwargs...) = CLI.julia_main(args...; kwargs...)

end # module
