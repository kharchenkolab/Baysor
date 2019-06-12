using DataFrames
using NearestNeighbors

struct CenterData
    centers::DataFrame
    scale_estimate::Float64
end

estimate_scale_from_centers(centers::Array{Float64, 2}) = median(maximum.(knn(KDTree(centers), centers, 2)[2]))

function load_centers(path::String; kwargs...)::CenterData
    file_ext = splitext(path)[2]
    if file_ext == ".csv"
        df_centers = read_spatial_df(path, gene_col=nothing, kwargs...) |> unique;
        return CenterData(df_centers, position_data(df_centers))
    end

    error("Unsupported file extension: '$file_ext'")
end