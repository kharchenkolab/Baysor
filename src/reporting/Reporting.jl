module Reporting

using LazyModules
using DataFrames

@lazy import Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
@lazy import CairoMakie as MK = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"

using Deneb

using ..Baysor.Utils

include("utils.jl")
include("vega_wrappers.jl")
include("plots.jl")
include("diagnostic_plots.jl")
include("cell_segmentation.jl")

end
