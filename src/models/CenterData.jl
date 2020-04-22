using DataFrames
using NearestNeighbors
using StatsBase

import Images

mutable struct CenterData
    centers::DataFrame
    center_covs::Union{Array{CovMat,1}, Nothing}
    scale_estimate::Float64
    scale_std_estimate::Float64
end

CenterData(centers::DataFrame, scale_estimate::Float64, scale_std_estimate::Float64) =
    CenterData(centers, nothing, scale_estimate, scale_std_estimate)

estimate_scale_from_centers(center_scales::Vector{Float64}) =
    (median(center_scales), mad(center_scales; normalize=true))

"""
Estimates scale as a 0.5 * median distance between two nearest centers
"""
estimate_scale_from_centers(centers::Matrix{Float64}) =
    estimate_scale_from_centers(maximum.(knn(KDTree(centers), centers, 2)[2]) ./ 2)

function estimate_scale_from_centers(seg_mask::Matrix{Int}; min_segment_size::Int=0)
    comp_lengths = Images.component_lengths(seg_mask)[2:end]
    return estimate_scale_from_centers(sqrt.(comp_lengths[comp_lengths .>= min_segment_size] / π))
end

extract_centers_from_mask(segmentation::BitMatrix; kwargs...) =
    extract_centers_from_mask(Images.label_components(segmentation); kwargs...)

function extract_centers_from_mask(segmentation::Matrix{Float64}; kwargs...)
    segmentation = round.(Int, Float64.(segmentation) .* 2^16)
    uniq_labels = sort(unique(segmentation))
    if length(uniq_labels) == 2
        return extract_centers_from_mask(Images.label_components(segmentation); kwargs...)
    end

    rank_per_label = Dict(Pair.(0:(length(uniq_labels)-1), uniq_labels))
    for i in 1:length(segmentation)
        segmentation[i] = rank_per_label[segmentation[i]]
    end

    return extract_centers_from_mask(segmentation; kwargs...)
end

function extract_centers_from_mask(segmentation_labels::Matrix{Int}; min_segment_size::Int=5)
    coords_per_label = coords_per_segmentation_label(segmentation_labels);
    coords_per_label = coords_per_label[size.(coords_per_label, 1) .>= min_segment_size]

    centers = hcat(vec.(mean.(coords_per_label, dims=1))...);
    center_covs = CovMat.(cov.(coords_per_label));

    for i in findall([any(isnan.(c)) for c in center_covs])
        center_covs[i] = CovMat(Float64[1. 0.; 0. 1.])
    end

    adjust_cov_matrix!.(center_covs);

    is_pos_def = isposdef.(center_covs)
    centers = centers[:, is_pos_def]
    center_covs = center_covs[is_pos_def]

    scale, scale_std = estimate_scale_from_centers(segmentation_labels; min_segment_size=min_segment_size)
    return CenterData(DataFrame(centers', [:y, :x])[:,[:x, :y]], center_covs, scale, scale_std)
end

function load_centers(path::String; min_segment_size::Int=5, kwargs...)::CenterData
    file_ext = splitext(path)[2]
    if file_ext == ".csv"
        df_centers = read_spatial_df(path; gene_col=nothing, kwargs...) |> unique;
        scale, scale_std = estimate_scale_from_centers(position_data(df_centers))
        return CenterData(df_centers, scale, scale_std)
    end

    if file_ext in [".png", ".jpg", ".tiff", ".tif"]
        return extract_centers_from_mask(Images.load(path); min_segment_size=min_segment_size)
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

    return CenterData(deepcopy(centers.centers[ids,:]), center_covs, centers.scale_estimate, centers.scale_std_estimate)
end

coords_per_segmentation_label(labels::Array{Int, 2}) =
    [hcat(mod.(ids .- 1, size(labels, 1)) .+ 1, div.(ids .- 1, size(labels, 1)) .+ 1) for ids in Images.component_indices(labels)[2:end]]
