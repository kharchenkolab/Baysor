module Processing

using LazyModules

using DataFrames
@lazy import CairoMakie as MK = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
@lazy import CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
@lazy import VegaLite = "112f6efa-9a02-5b7d-90c0-432ed331239a"
@lazy import KernelDensity as KDE = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
@lazy import VegaLite as VL = "112f6efa-9a02-5b7d-90c0-432ed331239a"
@lazy import Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
@lazy import ImageMorphology = "787d08f9-d448-5407-9aad-5290dd7ab264"

using Distributions
using LinearAlgebra
using ProgressMeter
using Statistics

import Distributions.pdf
import Distributions.logpdf
import Statistics.rand

include("utils/utils.jl")
include("utils/convex_hull.jl")
include("utils/spatial.jl")
include("utils/vega_wrappers.jl")

include("distributions/MvNormal.jl")
include("distributions/CategoricalSmoothed.jl")

include("models/InitialParams.jl")
include("models/Component.jl")
include("models/BmmData.jl")

include("data_processing/prior_segmentation.jl")
include("data_processing/triangulation.jl")
include("data_processing/umap_wrappers.jl")
include("data_processing/neighborhood_composition.jl")
include("data_processing/noise_estimation.jl")
include("data_processing/initialization.jl")
include("data_processing/boundary_estimation.jl")
include("data_processing/plots.jl")
include("data_processing/diagnostic_plots.jl")

include("bmm_algorithm/molecule_clustering.jl")
include("bmm_algorithm/compartment_segmentation.jl")
include("bmm_algorithm/tracing.jl")
include("bmm_algorithm/distribution_samplers.jl")
include("bmm_algorithm/history_analysis.jl")
include("bmm_algorithm/bmm_algorithm.jl")

include("data_loading/cli_wrappers.jl")
include("reports/cell_segmentation.jl")

end