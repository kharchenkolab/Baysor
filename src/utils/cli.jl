append_suffix(output::String, suffix) = "$(splitext(output)[1])_$suffix"

Base.@kwdef struct OutputPaths
    # The class is here for clarity of what files are saved
    segmented_df::String;
    cell_stats::String;
    counts::String;
    diagnostic_report::String;
    molecule_plot::String;
    polygons_2d::String;
    polygons_3d::String;
end

function get_output_paths(segmented_df_path::String; count_matrix_format::String)
    if isdir(segmented_df_path)
        segmented_df_path = joinpath(segmented_df_path, "segmentation.csv")
    end

    segmented_df = splitext(segmented_df_path)[2] == "" ? "$(segmented_df_path)_segmentation.csv" : segmented_df_path
    return OutputPaths(;
        segmented_df=segmented_df,
        Dict(k => append_suffix(segmented_df_path, v) for (k,v) in [
            :cell_stats => "cell_stats.csv",
            :counts => "counts." * count_matrix_format,
            :diagnostic_report => "diagnostics.html",
            :molecule_plot => "borders.html",
            :polygons_2d => "polygons_2d.json",
            :polygons_3d => "polygons_3d.json"
        ])...
    )
end