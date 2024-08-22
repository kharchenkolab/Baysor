module CommandLine

import Pkg
import Random
using LazyModules

@lazy import Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"

import ..Baysor: BPR, DAT, REP
using ..Utils
using Comonicon
using Configurations

pkg_version = Pkg.TOML.parsefile(
    joinpath(dirname(dirname(dirname(@__FILE__))), "Project.toml")
)["version"] # Pre-estimated for binary builds

include("logging.jl")

include("common.jl")
include("main.jl")
include("preview.jl")
include("segfree.jl")


end