using DataStructures
using KernelDensity

import ImageSegmentation
import Images

# modification of the function from the ImageSegmentation package
function watershed(img::Array{T,2}, markers::Array{S,2}, img_bin::Union{Array{Bool,2}, BitArray{2}, Nothing}=nothing) where {T<:Images.NumberLike, S<:Integer}

    if axes(img) != axes(markers)
        error("image size doesn't match marker image size")
    end

    if img_bin !== nothing
        markers = markers .* img_bin
    end
    segments = copy(markers)
    pq = PriorityQueue{CartesianIndex{2}, ImageSegmentation.PixelKey{T}}()
    time_step = 0

    R = CartesianIndices(axes(img))
    Istart, Iend = first(R), last(R)
    for i in R
        if markers[i] > 0
            for j in CartesianIndices(ImageSegmentation._colon(max(Istart,i-one(i)), min(i+one(i),Iend)))
                if (img_bin !== nothing) && (!img_bin[j])
                    continue
                end

                if segments[j] == 0
                    segments[j] = markers[i]
                    enqueue!(pq, j, ImageSegmentation.PixelKey(img[i], time_step))
                    time_step += 1
                end
            end
        end
    end

    while !isempty(pq)
        current = dequeue!(pq)
        segments_current = segments[current]
        img_current = img[current]
        for j in CartesianIndices(ImageSegmentation._colon(max(Istart,current-one(current)), min(current+one(current),Iend)))
            if (img_bin !== nothing) && (!img_bin[j])
                continue
            end

            if segments[j] == 0
                segments[j] = segments_current
                enqueue!(pq, j, ImageSegmentation.PixelKey(img_current, time_step))
                time_step += 1
            end
        end
    end

    TM = ImageSegmentation.meantype(T)
    num_segments        = Int(maximum(segments))
    labels              = Array(1:num_segments)
    region_means        = Dict{Int, TM}()
    region_pix_count    = Dict{Int, Int}()

    for i in R
        region_pix_count[segments[i]] = get(region_pix_count, segments[i], 0) + 1
        region_means[segments[i]] = get(region_means, segments[i], zero(TM)) + (img[i] - get(region_means, segments[i], zero(TM)))/region_pix_count[segments[i]]
    end
    return ImageSegmentation.SegmentedImage(segments, labels, region_means, region_pix_count)
end

function estimate_kde_grid(coords::Array{Int, 2}; bandwidth::Real)
    coords = Array{Float64, 2}(coords)
    min_vals = minimum(coords, dims=2)
    max_vals = maximum(coords, dims=2)

    midpoints = Tuple(KernelDensity.kde_range((s,e), Int(e - s)) for (s,e) in zip(min_vals, max_vals));

    return kde((coords[1,:], coords[2,:]), midpoints, bandwidth=Tuple(bandwidth for i in 1:2));
end

"""
    segment image, starting from local minima, using Otsu thresholding
"""
function segment_watershed(image::Array{Float64, 2}; verbose::Bool = false)
    if verbose
        @info "Initialize watershed"
    end

    local_minima = Images.findlocalminima(image);

    img_bin = (image .< Images.otsu_threshold(image));

    markers = zeros(Int, size(image));
    markers[local_minima] .= 1:length(local_minima);

    if verbose
        @info "Run watershed"
    end

    segments = watershed(image, markers, img_bin);

    labels = Images.label_components(segments.image_indexmap);

    if verbose
        @info "Done!"
    end

    return labels
end

function assign_molecules_to_kde_segmentation(labels::Array{Int, 2}, df_spatial::DataFrame, kde_res::KernelDensity.BivariateKDE)::Array{Int, 1}
    nn_ids = [vcat(knn(KDTree(Array(collect(rng)')), Array(df_spatial[!,s]'), 1)[1]...) for (s, rng) in zip([:x, :y], [kde_res.x, kde_res.y])];
    return labels[CartesianIndex.(nn_ids...)];
end

function segment_molecules_watershed(df_spatial::DataFrame; bandwidth::Real, verbose=false)
    if verbose
        @info "Estimate KDE density"
    end

    c_kde = estimate_kde_grid(round.(Int, position_data(df_spatial)), bandwidth=bandwidth);
    img_labels = segment_watershed(.-c_kde.density; verbose=verbose);
    return assign_molecules_to_kde_segmentation(img_labels, df_spatial, c_kde);
end