function select_ids_uniformly(xs::Vector{T}, ys::Vector{T}, confidence::Vector{Float64}, n::Int; confidence_threshold::Float64=0.95)::Vector{Int} where T <: Real
    dense_ids = findall(confidence .>= confidence_threshold)
    if length(dense_ids) < n
        return dense_ids
    end

    xs = xs[dense_ids] .- minimum(xs[dense_ids])
    ys = ys[dense_ids] .- minimum(ys[dense_ids])

    n_y_bins = round(Int, sqrt(n) * maximum(ys) / maximum(xs))
    n_x_bins = div(n, n_y_bins)

    index_vals = round.(Int, xs .* (n_x_bins / maximum(xs))) .* n_y_bins .+ round.(Int, ys .* (n_y_bins / maximum(ys)));
    return dense_ids[sortperm(index_vals)[unique(round.(Int, range(1, length(dense_ids), length=n)))]];
end

neighborhood_count_matrix(data::Union{BmmData, T} where T <: AbstractDataFrame, k::Int; kwargs...) =
    neighborhood_count_matrix(position_data(data), composition_data(data), k; kwargs...)

function neighborhood_count_matrix(pos_data::Matrix{T} where T <: Real, genes::Vector{Int}, k::Int;
    n_genes::Int=maximum(genes), normalize_by_dist::Bool=true, normalize::Bool=true)
    if k < 3
        @warn "Too small value of k: $k. Setting it to 3."
        k = 3
    end

    if !normalize
        normalize_by_dist = false
    end

    k = min(k, size(pos_data, 2))

    neighbors, dists = knn(KDTree(pos_data), pos_data, k, true);

    n_cm = zeros((normalize ? Float64 : Int), n_genes, size(pos_data, 2));

    if !normalize_by_dist
        @threads for (i,ids) in collect(enumerate(neighbors))
            if !normalize
                count_array!(view(n_cm, :, i), genes[ids])
            else
                prob_array!(view(n_cm, :, i), genes[ids])
            end
        end

        return n_cm
    end

    # normalize_by_dist doesn't affect result much, but neither it slow the work down. In theory, should work better for border cases.
    med_closest_dist = median(getindex.(dists, 2))
    @threads for (i,(ids, dists)) in collect(enumerate(zip(neighbors, dists)))
        c_probs = view(n_cm, :, i)
        for (gene, dist) in zip(genes[ids], dists)
            c_probs[gene] += 1 / max(dist, med_closest_dist)
        end
    end

    return n_cm ./ sum(n_cm, dims=1);
end

function gene_composition_transformation(count_matrix::Array{Float64, 2}, confidence::Vector{Float64}=ones(size(count_matrix, 2));
        sample_size::Int=10000, seed::Int=42, method::Symbol=:umap, kwargs...)
    sample_size = min(sample_size, size(count_matrix, 2))
    Random.seed!(seed)

    pc2 = transform(fit(MultivariateStats.PCA, count_matrix, maxoutdim=2), count_matrix);
    sample_ids = select_ids_uniformly(pc2[1,:], pc2[2,:], confidence, n=sample_size)

    count_matrix_sample = count_matrix[:,sample_ids]

    if method == :umap
        return fit(UmapFit, count_matrix_sample, n_components=3; kwargs...);
    end

    if method != :pca
        error("Unknown method: '$method'")
    end

    return fit(MultivariateStats.PCA, count_matrix_sample, maxoutdim=3; kwargs...);
end

function gene_composition_colors(count_matrix::Array{Float64, 2}, transformation::UmapFit; color_range::T where T<:Real =2000.0)
    mtx_colors = MultivariateStats.transform(transformation, count_matrix);
    min_vals = minimum(transformation.embedding, dims=2)
    mtx_colors .-= min_vals
    mtx_colors ./= (maximum(transformation.embedding, dims=2) - min_vals);

    mtx_colors[1,:] .*= 100
    mtx_colors[2:3,:] .-= 0.5
    mtx_colors[2:3,:] .*= color_range

    return vec(mapslices(col -> Colors.Lab(col...), mtx_colors, dims=1))
end