append_suffix(output::String, suffix) = "$(splitext(output)[1])_$suffix"

Base.@kwdef struct OutputPaths
    # The class is here for clarity of what files are saved
    segmented_df::String;
    cell_stats::String;
    counts::String;
    diagnostic_report::String;
    molecule_plot::String;
    polygons::String;
end

get_output_paths(segmented_df_path::String) = OutputPaths(;
    segmented_df=segmented_df_path,
    Dict(k => append_suffix(segmented_df_path, v) for (k,v) in [
        :cell_stats => "cell_stats.csv",
        :counts => "counts.tsv",
        :diagnostic_report => "diagnostics.html",
        :molecule_plot => "borders.html",
        :polygons => "polygons.json"
    ])...
)