using NearestNeighbors
using SparseArrays
using StatsBase

using FileIO: load # reqiored for `load` to work properly
import ImageIO

@lazy import ImageCore = "a09fc81d-aa75-5fe9-8630-4744c3626534"
@lazy import ImageMorphology = "787d08f9-d448-5407-9aad-5290dd7ab264"
@lazy import MAT = "23992714-dd62-5051-b70f-ba57cb901cac"


staining_value_per_transcript(df_spatial::DataFrame, args...; kwargs...) =
    staining_value_per_transcript(df_spatial.x, df_spatial.y, args...; kwargs...)

staining_value_per_transcript(pos_data::Matrix{<:Real}, args...; kwargs...) =
    staining_value_per_transcript(pos_data[1,:], pos_data[2,:], args...; kwargs...)

function staining_value_per_transcript(
        x_vals::Vector{T2}, y_vals::Vector{T2}, staining::MT where MT <: AbstractMatrix{T}; quiet::Bool=false
    )::Vector{T} where T <: Real where T2 <: Real

    x_vals = round.(Int, x_vals)
    y_vals = round.(Int, y_vals)
    if !quiet && ((maximum(x_vals) > size(staining, 2)) || (maximum(y_vals) > size(staining, 1)))
        @warn "Maximum transcript coordinates are $((maximum(y_vals), maximum(x_vals))), which is larger than the prior image size: $(size(staining)). Filling it with 0."
    end

    if !quiet && ((minimum(x_vals) < 1) || (minimum(y_vals) < 1))
        @warn "Minimum transcript coordinates are < 1: $((minimum(y_vals), minimum(x_vals))). Filling it with 0."
    end

    if !quiet && ((maximum(x_vals) < 0.5 * size(staining, 2)) || (maximum(y_vals) < 0.5 * size(staining, 1)))
        @warn "Maximum transcript coordinates are $((maximum(y_vals), maximum(x_vals))), which is much smaller than the prior image size: $(size(staining)). May be result of an error."
    end

    inds = findall((x_vals .> 0) .& (x_vals .<= size(staining, 2)) .& (y_vals .> 0) .& (y_vals .<= size(staining, 1)))
    staining_vals = zeros(T, length(x_vals))
    for i in inds
        staining_vals[i] = staining[y_vals[i], x_vals[i]];
    end

    return staining_vals
end

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

        (length(size(labels)) == 2) || error("Segmentation mask must be a 2D image, but it has $(length(size(labels))) dimensions.")
        return SparseMatrixCSC{Int, Int}(labels)
    end

    labels = load(path) |> ImageCore.channelview |> ImageCore.rawview;
    (length(size(labels)) == 2) || error("Segmentation mask must be a 2D image, but it has $(length(size(labels))) dimensions.")

    labels = labels |> sparse |> dropzeros!
    if length(unique(nonzeros(labels))) == 1
        return BitMatrix(labels .> 0) |> ImageMorphology.label_components |> sparse
    end

    return labels
end

estimate_scale_from_centers(radius_per_segment::Vector{Float64}) =
    (median(radius_per_segment), mad(radius_per_segment; normalize=true))

"""
Estimates scale as a 0.5 * median distance between two nearest centers
"""
estimate_scale_from_centers(centers::Matrix{Float64}) =
    estimate_scale_from_centers(maximum.(knn(KDTree(centers), centers, 2)[2]) ./ 2)

estimate_scale_from_centers(seg_labels::Matrix{<:Integer}) =
    estimate_scale_from_centers(sqrt.(ImageMorphology.component_lengths(seg_labels)[2:end] ./ π))

function estimate_scale_from_centers(seg_labels::SparseMatrixCSC{<:Integer})
    nz_vals = nonzeros(seg_labels);
    nz_vals = nz_vals[nz_vals .> 0]
    if isempty(nz_vals)
        error("No transcripts detected inside the segmented regions. Please, check that transcript coordinates match those in the segmentation mask.")
    end
    return estimate_scale_from_centers(sqrt.(filter(x -> x > 0, count_array(nz_vals)) ./ π))
end

function estimate_scale_from_assignment(pos_data::Matrix{Float64}, assignment::Vector{Int}; min_molecules_per_cell::Int)
    pd_per_cell = [pos_data[:,ids] for ids in Utils.split_ids(assignment, drop_zero=true) if length(ids) >= min_molecules_per_cell];
    if length(pd_per_cell) < 3
        error("Not enough prior cells pass the min_molecules_per_cell=$(min_molecules_per_cell) threshold. Please, specify scale manually.")
    end
    return estimate_scale_from_centers(hcat(mean.(pd_per_cell, dims=2)...))
end

filter_segmentation_labels!(segment_per_transcript::Vector{<:Integer}; kwargs...) =
    filter_segmentation_labels!(segment_per_transcript, segment_per_transcript; kwargs...)[1]

filter_segmentation_labels!(segmentation_labels::MT where MT <: AbstractMatrix{<:Integer}, df_spatial::DataFrame; quiet::Bool=false, kwargs...) =
    filter_segmentation_labels!(segmentation_labels, staining_value_per_transcript(df_spatial, segmentation_labels, quiet=quiet); kwargs...)

function filter_segmentation_labels!(segmentation_labels::SparseMatrixCSC{<:Integer}, segment_per_transcript::Vector{<:Integer}; kwargs...)
    filter_segmentation_labels!(segmentation_labels.nzval, segment_per_transcript; kwargs...)
    dropzeros!(segmentation_labels)
    return segmentation_labels
end

function filter_segmentation_labels!(segmentation_labels::MT where MT <: AbstractArray{<:Integer}, segment_per_transcript::Vector{<:Integer}; min_molecules_per_segment::Int)
    n_mols_per_label = count_array(segment_per_transcript, max_value=maximum(segmentation_labels), drop_zero=true)

    for labs in (segmentation_labels, segment_per_transcript)
        for i in 1:length(labs)
            lab = labs[i]
            if (lab > 0) && (n_mols_per_label[lab] < min_molecules_per_segment)
                labs[i] = 0
            end
        end
    end

    return (segmentation_labels, segment_per_transcript)
end
