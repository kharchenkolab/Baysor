module Processing

using LazyModules
using Compat

using DataFrames

@lazy import Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
@lazy import ImageMorphology = "787d08f9-d448-5407-9aad-5290dd7ab264"

using Distributions
using LinearAlgebra
using ProgressMeter
using Statistics
using Pipe: @pipe as @p

import Distributions.pdf
import Distributions.logpdf
import Statistics.rand

using ..Baysor.Utils
using ..Baysor.Utils: SegmentationOptions, PlottingOptions, split_ids, fsort
import ..Baysor.Utils.split

include("utils/utils.jl")
include("utils/convex_hull.jl")
include("utils/cli_wrappers.jl")

include("distributions/MvNormal.jl")
include("distributions/CategoricalSmoothed.jl")

include("models/AdjList.jl")
include("models/InitialParams.jl")
include("models/Component.jl")
include("models/BmmData.jl")

include("data_processing/triangulation.jl")
include("data_processing/umap_wrappers.jl")
include("data_processing/neighborhood_composition.jl")
include("data_processing/noise_estimation.jl")
include("data_processing/initialization.jl")
include("data_processing/boundary_estimation.jl")

include("bmm_algorithm/molecule_clustering.jl")
include("bmm_algorithm/compartment_segmentation.jl")
include("bmm_algorithm/tracing.jl")
include("bmm_algorithm/history_analysis.jl")
include("bmm_algorithm/bmm_algorithm.jl")

end