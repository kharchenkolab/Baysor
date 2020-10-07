using DataFrames
using NearestNeighbors
using SparseArrays
using StatsBase

import Images
import MAT

function load_segmentation_mask(path::String)::SparseMatrixCSC
    if lowercase(splitext(path)[end]) == ".mat"
        labels = first(MAT.matread(path))[2]
        if !(typeof(labels) <: Integer)
            for v in labels
                if (v < -1e-5) || (abs(round(Int, v) - v) > 1e-5)
                    error(".mat file must have non-negative integer array with labels in its first slot, but it contains value '$v'")
                end
            end
        end

        return SparseMatrixCSC{Int, Int}(labels)
    end

    labels = Images.load(path) |> Images.channelview |> Images.rawview |> sparse |> dropzeros!
    if length(unique(nonzeros(labels))) == 1
        return BitMatrix(labels) |> Images.label_components |> sparse
    end

    return labels
end

estimate_scale_from_centers(center_scales::Vector{Float64}) =
    (median(center_scales), mad(center_scales; normalize=true))

"""
Estimates scale as a 0.5 * median distance between two nearest centers
"""
estimate_scale_from_centers(centers::Matrix{Float64}) =
    estimate_scale_from_centers(maximum.(knn(KDTree(centers), centers, 2)[2]) ./ 2)

estimate_scale_from_centers(seg_labels::Matrix{<:Integer}) =
    estimate_scale_from_centers(sqrt.(Images.component_lengths(seg_labels)[2:end] ./ π))

function estimate_scale_from_centers(seg_labels::SparseMatrixCSC{<:Integer})
    nz_vals = nonzeros(seg_labels);
    nz_vals = nz_vals[nz_vals .> 0]
    if isempty(nz_vals)
        error("No transcripts detected inside the segmented regions. Please, check that transcript coordinates match those in the segmentation mask.")
    end
    return estimate_scale_from_centers(sqrt.(filter(x -> x > 0, count_array(nz_vals)) ./ π))
end

filter_segmentation_labels!(segmentation_labels::MT where MT <: AbstractMatrix{<:Integer}, df_spatial::DataFrame; quiet::Bool=false, kwargs...) =
    filter_segmentation_labels!(segmentation_labels, staining_value_per_transcript(df_spatial, segmentation_labels, quiet=quiet); kwargs...)

function filter_segmentation_labels!(segmentation_labels::SparseMatrixCSC{<:Integer}, segment_per_transcript::Vector{<:Integer}; kwargs...)
    filter_segmentation_labels!(segmentation_labels.nzval, segment_per_transcript; kwargs...)
    dropzeros!(segmentation_labels)
    return segmentation_labels
end

function filter_segmentation_labels!(segmentation_labels::MT where MT <: AbstractArray{<:Integer}, segment_per_transcript::Vector{<:Integer}; min_transcripts_per_segment::Int)
    n_mols_per_label = count_array(segment_per_transcript, max_value=maximum(segmentation_labels), drop_zero=true)

    for i in 1:length(segmentation_labels)
        lab = segmentation_labels[i]
        if (lab > 0) && (n_mols_per_label[lab] < min_transcripts_per_segment)
            segmentation_labels[i] = 0
        end
    end
    return segmentation_labels
end
