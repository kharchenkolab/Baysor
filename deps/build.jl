ENV["LazyModules_lazyload"] = "false"

using Baysor
using Comonicon: Builder
using Pkg

Pkg.activate(pkgdir(Baysor))
Pkg.instantiate() # Comonicon doesn't pull non-registered dependencies if there are any

# delete!(ENV, "SHELL") # Disable autocompletion installation, as it doesn't work with sub-modules. v1.0.1 doesn't allow to disable it in TOML
Builder.command_main(Baysor)