module Reporting

using LazyModules
using DataFrames

import ..Baysor: BPR

@lazy import Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"

import CairoMakie as MK # lazy loading doesn't help, as it's used as an argument type in Base.show
import VegaLite as VL # lazy loading doesn't help, as macroses trigger compilation

include("reports/utils.jl")
include("reports/vega_wrappers.jl")
include("reports/plots.jl")
include("reports/diagnostic_plots.jl")
include("reports/cell_segmentation.jl")

end
