module Baysor

using Distributed
using Distributions
using LinearAlgebra
using ProgressMeter
using Statistics

import Distributions.pdf
import Distributions.logpdf
import Statistics.rand

export
    BmmData, bmm!,
    build_molecule_graph, load_df, initial_distribution_arr,
    boundary_polygons,
    neighborhood_count_matrix, gene_composition_transformation, gene_composition_colors,
    append_confidence!,
    plot_cell_borders_polygons, plot_cell_borders_polygons!, plot_expression_vectors,
    cluster_molecules_on_mrf

include("utils/utils.jl")
include("utils/kmeans.jl")
include("utils/logging.jl")
include("utils/convex_hull.jl")
include("utils/spatial.jl")

include("distributions/MvNormal.jl")
include("distributions/CategoricalSmoothed.jl")

include("models/InitialParams.jl")
include("models/Component.jl")
include("models/BmmData.jl")

include("data_processing/prior_segmentation.jl")
include("data_processing/watershed.jl")
include("data_processing/triangulation.jl")
include("data_processing/umap_wrappers.jl")
include("data_processing/expression_local_vectors.jl")
include("data_processing/noise_estimation.jl")
include("data_processing/initialization.jl")
include("data_processing/plots.jl")
include("data_processing/boundary_estimation.jl")
include("data_processing/preview_diagnostics.jl")

include("bmm_algorithm/molecule_clustering.jl")
include("bmm_algorithm/tracing.jl")
include("bmm_algorithm/distribution_samplers.jl")
include("bmm_algorithm/history_analysis.jl")
include("bmm_algorithm/doublet_splitting.jl")
include("bmm_algorithm/bmm_algorithm.jl")
include("bmm_algorithm/smoothing.jl")

include("data_processing/validation.jl")

include("cli_wrappers.jl")

end # module
