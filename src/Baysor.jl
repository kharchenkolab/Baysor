module Baysor

# export
#     BmmData, bmm!,
#     build_molecule_graph, load_df, initialize_bmm_data,
#     boundary_polygons,
#     neighborhood_count_matrix, gene_composition_transformation, gene_composition_colors,
#     append_confidence!,
#     plot_molecules, plot_molecules!, plot_expression_vectors,
#     cluster_molecules_on_mrf

# Utils: Minimal functions with zero compilation time shared across submodules
include("utils/Utils.jl")

DAT = include("data_loading/DataLoading.jl")
REP = include("reporting/Reporting.jl")
BPR = include("processing/Processing.jl")
CLI = include("cli/CLI.jl")

using Comonicon: @cast, cmd_error
import Comonicon: @main as @app
# using .CLI: command_main, julia_main

using .CLI: CASTED_COMMANDS
@app

end # module
