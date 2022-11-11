module Baysor

# export
#     BmmData, bmm!,
#     build_molecule_graph, load_df, initialize_bmm_data,
#     boundary_polygons,
#     neighborhood_count_matrix, gene_composition_transformation, gene_composition_colors,
#     append_confidence!,
#     plot_molecules, plot_molecules!, plot_expression_vectors,
#     cluster_molecules_on_mrf


include("LazySubmodules.jl")
using .LazySubmodules

@lazy_submodule BPR = "Processing.jl"
@lazy_submodule CLI = "CLI.jl"

end # module
