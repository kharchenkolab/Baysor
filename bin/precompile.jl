import Baysor as B

include("../test/runtests.jl")

B.run_cli(["--help"]);
B.run_cli(["run", "--help"]);