using DataFrames
using NearestNeighbors

mutable struct CenterData
    centers::DataFrame
    center_covs::Union{Array{Array{Float64,2},1}, Nothing}
    scale_estimate::Float64

    CenterData(centers::DataFrame, scale_estimate::Float64) = new(centers, nothing, scale_estimate)
    CenterData(centers::DataFrame, center_covs::Array{Array{Float64,2},1}, scale_estimate::Float64) = new(centers, center_covs, scale_estimate)
end

estimate_scale_from_centers(centers::Array{Float64, 2}) = median(maximum.(knn(KDTree(centers), centers, 2)[2])) / 2

function load_centers(path::String; kwargs...)::CenterData
    file_ext = splitext(path)[2]
    if file_ext == ".csv"
        df_centers = read_spatial_df(path; gene_col=nothing, kwargs...) |> unique;
        return CenterData(df_centers, estimate_scale_from_centers(position_data(df_centers)))
    end

    error("Unsupported file extension: '$file_ext'")
end

# function subset_by_coords(subsetting_df::DataFrame, coord_df::DataFrame)
#     pos_subs = position_data(subsetting_df)
#     pos_coords = position_data(coord_df)
#     ids = vec(all((pos_subs .>= minimum(pos_coords, dims=2)) .& (pos_subs .<= maximum(pos_coords, dims=2)), dims=1));

#     return subsetting_df[ids,:]
# end

function subset_by_coords(centers::CenterData, coord_df::DataFrame)
    pos_subs = position_data(centers.centers)
    pos_coords = position_data(coord_df)
    ids = vec(all((pos_subs .>= minimum(pos_coords, dims=2)) .& (pos_subs .<= maximum(pos_coords, dims=2)), dims=1));

    center_covs = centers.center_covs === nothing ? nothing : deepcopy(centers.center_covs[ids])

    return CenterData(deepcopy(centers.centers[ids,:]), center_covs, centers.scale_estimate)
end