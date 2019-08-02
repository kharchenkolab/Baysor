using DataFrames
using NearestNeighbors

import Images

mutable struct CenterData
    centers::DataFrame
    center_covs::Union{Array{Array{Float64,2},1}, Nothing}
    scale_estimate::Float64

    CenterData(centers::DataFrame, scale_estimate::Float64) =
        new(centers, nothing, scale_estimate)
    CenterData(centers::DataFrame, center_covs::Union{Array{Array{Float64,2},1}, Nothing}, scale_estimate::Float64) =
        new(centers, center_covs, scale_estimate)
end

"""
Estimates scale as a 0.5 * median distance between two nearest centers multiplied by `scale_mult`
"""
estimate_scale_from_centers(centers::Array{Float64, 2}; scale_mult::Float64=0.75) =
    scale_mult * median(maximum.(knn(KDTree(centers), centers, 2)[2])) / 2

function load_centers(path::String; min_segment_size::Int=0, scale_mult::Float64=0.75, kwargs...)::CenterData
    file_ext = splitext(path)[2]
    if file_ext == ".csv"
        df_centers = read_spatial_df(path; gene_col=nothing, kwargs...) |> unique;
        scale = estimate_scale_from_centers(position_data(df_centers), scale_mult=scale_mult)
        return CenterData(df_centers, scale)
    end

    if file_ext in [".png", ".jpg", ".tiff"]
        segmentation_labels = Images.load(path) |> Images.label_components |> Array{Int}
        coords_per_label = Baysor.coords_per_segmentation_label(segmentation_labels);
        coords_per_label = coords_per_label[size.(coords_per_label, 1) .>= min_segment_size]

        centers = hcat(vec.(mean.(coords_per_label, dims=1))...);
        center_covs = cov.(coords_per_label);

        for i in findall([any(isnan.(c)) for c in center_covs])
            center_covs[i] = copy(Float64[1. 0.; 0. 1.])
        end

        center_covs = adjust_cov_matrix.(center_covs)

        is_pos_def = isposdef.(center_covs)
        centers = centers[:, is_pos_def]
        center_covs = center_covs[is_pos_def]

        scale = estimate_scale_from_centers(centers, scale_mult=scale_mult)
        return CenterData(DataFrame(centers', [:y, :x])[:,[:x, :y]], center_covs, scale)
    end

    error("Unsupported file extension: '$file_ext'")
end

function subset_by_coords(subsetting_df::DataFrame, coord_df::DataFrame)
    pos_subs = position_data(subsetting_df)
    pos_coords = position_data(coord_df)
    ids = vec(all((pos_subs .>= minimum(pos_coords, dims=2)) .& (pos_subs .<= maximum(pos_coords, dims=2)), dims=1));

    return subsetting_df[ids,:]
end

function subset_by_coords(centers::CenterData, coord_df::DataFrame)
    pos_subs = position_data(centers.centers)
    pos_coords = position_data(coord_df)
    ids = vec(all((pos_subs .>= minimum(pos_coords, dims=2)) .& (pos_subs .<= maximum(pos_coords, dims=2)), dims=1));

    center_covs = centers.center_covs === nothing ? nothing : deepcopy(centers.center_covs[ids])

    return CenterData(deepcopy(centers.centers[ids,:]), center_covs, centers.scale_estimate)
end

function coords_per_segmentation_label(segmentation_mask::Array{Int, 2})
    coords = [Array{Float64, 1}[] for v in 1:maximum(segmentation_mask)];

    @inbounds for r in 1:size(segmentation_mask, 1)
        for c in 1:size(segmentation_mask, 2)
            label = segmentation_mask[r, c]
            if label > 0
                push!(coords[label], Float64[r, c])
            end
        end
    end

    return [hcat(c...)' for c in coords]
end
