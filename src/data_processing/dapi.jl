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
