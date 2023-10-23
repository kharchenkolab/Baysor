module Utils

export
    append_suffix, get_output_paths, OutputPaths,
    default_param_value, split_string_list, get_cell_name, isnoise,
    count_array, count_array!, fmax, fmin, fsort, val_range,
    DataOptions, SegmentationOptions, PlottingOptions, RunOptions

include("options.jl")
include("cli.jl")
include("general.jl")

end
