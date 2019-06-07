import Images

function extract_cell_centers(stain::Array{Float64,2}; min_x::Array{Int64,1}=nothing, brightness_threshold::Float64=0.15, pixel_threshold::Int=50)::Array{Float64,2}
    max_brightness = maximum(stain);
    if max_brightness < 1e-5
        return Array{Float64,2}(undef, 0, 0)
    end

    labels = Images.label_components(stain ./ max_brightness .> brightness_threshold);
    pixels_per_label = count_array(vec(labels .+ 1))[2:end];

    coord_array = [[] for l in unique(labels)[2:end]]
    for y in 1:size(labels, 1)
        for x in 1:size(labels, 2)
            l = labels[y, x]
            if l == 0
                continue
            end
            push!(coord_array[l], [x, y])
        end
    end

    coord_array = [hcat(a...) for a in coord_array[length.(coord_array) .> pixel_threshold]];
    if length(coord_array) == 0
        return Array{Float64,2}(undef, 0, 0)
    end

    centers = copy(hcat(mean.(coord_array, dims=2)...)');
    if min_x !== nothing
        centers .+= min_x'
    end

    return centers
end

function coords_per_label(segmentation_mask::Array{Int, 2})
    coords = [Array{Float64, 1}[] for v in 1:maximum(segmentation_mask)];

    @inbounds for r in 1:size(segmentation_mask, 1)
        for c in 1:size(segmentation_mask, 2)
            label = segmentation_mask[r, c]
            if label > 0
                push!(coords[label], Float64[r, c])
            end
        end
    end

    return coords
end

function parse_segmentation_mask(segmentation_mask::Array{Int, 2})
    coords = coords_per_label(segmentation_mask)
    centers = DataFrame(hcat([vec(mean(hcat(c...), dims=2)) for c in coords]...)', [:x, :y])
    polygons = [copy(p') for p in convex_hull.(coords)]

    return centers, polygons
end