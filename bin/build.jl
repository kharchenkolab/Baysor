import Baysor as B

B.run_cli(["--help"]);
B.run_cli(["run", "--help"]);
B.run_cli(["preview", "--help"]);

# To build, run create_sysimage(:Baysor; precompile_execution_file="$(dirname(pathof(Baysor)))/../bin/build.jl", replace_default=true)
# If you need to build app, comment @require part in __init__ function of ImageCore (~/.julia/packages/ImageCore/BTvmx/src/ImageCore.jl:167)
# and use create_app instead of create_sysimage
