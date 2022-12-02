module Utils

export
    append_suffix, get_output_paths, OutputPaths,
    default_param_value, split_string_list,
    count_array, count_array!, fmax, fmin,
    DataOptions, SegmentationOptions, PlottingOptions, RunOptions

include("options.jl")
include("cli.jl")
include("utils.jl")

end