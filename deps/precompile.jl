ENV["LazyModules_lazyload"] = "false"

import Baysor, Pkg

Baysor.command_main(["-h"]);
for cmd in keys(Baysor.CLI.CASTED_COMMANDS)
    if cmd != "main"
        Baysor.command_main([cmd, "-h"]);
    end
end

Pkg.activate(pkgdir(Baysor))
include(joinpath(pkgdir(Baysor), "test", "runtests.jl"))
