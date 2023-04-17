ENV["LazyModules_lazyload"] = "false"

using Baysor
using Comonicon: Builder
using Pkg

# Just to make sure they're initialized. Failsafe for `LazyModules_lazyload`.
Baysor.load_module(Baysor.DAT)
Baysor.load_module(Baysor.REP)
Baysor.load_module(Baysor.BPR)
Baysor.load_module(Baysor.CLI)

Pkg.activate(pkgdir(Baysor))
Pkg.instantiate() # Comonicon doesn't pull non-registered dependencies, such as VegaLite.jl

delete!(ENV, "SHELL") # Disable autocompletion installation, as it doesn't work with sub-modules. v1.0.1 doesn't allow to disable it in TOML
Builder.command_main(Baysor)