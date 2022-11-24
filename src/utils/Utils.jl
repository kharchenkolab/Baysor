module Utils

export
    append_suffix, get_output_paths, OutputPaths,
    default_param_value,
    count_array, count_array!, fmax, fmin

include("cli.jl")
include("utils.jl")

end