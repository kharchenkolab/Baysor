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


load_segmentation_mask(path::String) =
    load_segmentation_mask(Float64.(Images.load(path)))

function load_segmentation_mask(mask::Matrix{Float64})::Matrix{Int}
    segmentation = round.(Int, mask .* 2^16)
    if length(unique(segmentation)) == 2
        return Images.label_components(segmentation)
    end

    return segmentation
end

"""
Estimates scale as a 0.5 * median distance between two nearest centers
"""
estimate_scale_from_centers(centers::Matrix{Float64}) =
    estimate_scale_from_centers(maximum.(knn(KDTree(centers), centers, 2)[2]) ./ 2)

estimate_scale_from_centers(seg_mask::Matrix{Int}) =
    estimate_scale_from_centers(sqrt.(Images.component_lengths(seg_mask)[2:end] ./ π))

extract_centers_from_mask(segmentation::BitMatrix, args...; kwargs...) =
    extract_centers_from_mask!(Images.label_components(segmentation), args...; kwargs...)

extract_centers_from_mask(segmentation::Matrix{Float64}, args...; kwargs...) =
    extract_centers_from_mask!(load_segmentation_mask(segmentation), args...; kwargs...)

function extract_centers_from_mask!(segmentation_labels::Matrix{Int}, df_spatial::Union{DataFrame, Nothing}=nothing; min_transcripts_per_center::Int=2)
    if df_spatial !== nothing
        n_mols_per_label = count_array(staining_value_per_transcript(df_spatial, segmentation_labels), max_value=maximum(segmentation_labels), drop_zero=true)

        for i in 1:length(segmentation_labels)
            lab = segmentation_labels[i]
            if (lab > 0) && (n_mols_per_label[lab] < min_transcripts_per_center)
                segmentation_labels[i] = 0
            end
        end
    elseif min_transcripts_per_center > 0
        @warn "df_spatial not provided, but min_transcripts_per_center > 0. Centers weren't filtered by min_transcripts_per_center."
    end

    uniq_labels = sort(unique(segmentation_labels))
    rank_per_label = Dict(Pair.(uniq_labels, 0:(length(uniq_labels)-1)))
    for i in 1:length(segmentation_labels)
        segmentation_labels[i] = rank_per_label[segmentation_labels[i]]
    end

    coords_per_label = coords_per_segmentation_label(segmentation_labels);

    centers = hcat(vec.(mean.(coords_per_label, dims=1))...);
    center_covs = CovMat.(cov.(coords_per_label)); # TODO: shoule it be MLE cov?

    for i in findall([any(isnan.(c)) for c in center_covs])
        center_covs[i] = CovMat(Float64[1. 0.; 0. 1.])
    end

    adjust_cov_matrix!.(center_covs);

    is_pos_def = isposdef.(center_covs)
    centers = centers[:, is_pos_def]
    center_covs = center_covs[is_pos_def]

    scale, scale_std = estimate_scale_from_centers(segmentation_labels)
    # return centers
    return CenterData(DataFrame(centers', [:y, :x])[:,[:x, :y]], center_covs, scale, scale_std)
end

function load_centers(path::String, df_spatial::Union{DataFrame, Nothing}=nothing; min_transcripts_per_center::Int=2, kwargs...)#::CenterData
    file_ext = splitext(path)[2]
    if file_ext == ".csv"
        df_centers = read_spatial_df(path; gene_col=nothing, kwargs...) |> unique;
        scale, scale_std = estimate_scale_from_centers(position_data(df_centers))
        return CenterData(df_centers, scale, scale_std)
    end

    if file_ext in [".png", ".jpg", ".jpeg", ".tiff", ".tif"]
        return extract_centers_from_mask(Float64.(Images.load(path)), df_spatial; min_transcripts_per_center=min_transcripts_per_center)
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
