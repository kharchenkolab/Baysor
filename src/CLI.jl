module CommandLine

import Pkg
import Random
using LazyModules

@lazy import Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"

import ..Baysor: BPR, DAT, REP

pkg_version = Pkg.TOML.parsefile(joinpath(dirname(dirname(@__FILE__)), "Project.toml"))["version"] # Pre-estimated for binary builds

include("cli/logging.jl")

include("cli/common.jl")
include("cli/main.jl")
include("cli/preview.jl")
include("cli/segfree.jl")

end