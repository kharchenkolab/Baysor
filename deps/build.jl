ENV["LazyModules_lazyload"] = "false"

using Baysor
using Comonicon: Builder
using Pkg

Pkg.activate(pkgdir(Baysor))
Pkg.instantiate() # Comonicon doesn't pull non-registered dependencies, such as VegaLite.jl

# delete!(ENV, "SHELL") # Disable autocompletion installation, as it doesn't work with sub-modules
Builder.command_main(Baysor)