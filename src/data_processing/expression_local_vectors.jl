using Base.Threads

neighborhood_count_matrix(data::Union{BmmData, T} where T <: AbstractDataFrame, k::Int; kwargs...) =
    neighborhood_count_matrix(position_data(data), composition_data(data), k; kwargs...)

function neighborhood_count_matrix(pos_data::Matrix{<:Real}, genes::Vector{Int}, k::Int;
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

    # normalize_by_dist doesn't affect result much, but neither it slows the work down. In theory, should work better for border cases.
    med_closest_dist = median([d[d .> 1e-15][1] for d in dists if any(d .> 1e-15)]) # account for problems with points with duplicating coordinates
    @threads for (i,(ids, dists)) in collect(enumerate(zip(neighbors, dists)))
        c_probs = view(n_cm, :, i)
        for (gene, dist) in zip(genes[ids], dists)
            c_probs[gene] += 1 / max(dist, med_closest_dist)
        end
    end

    return n_cm ./ sum(n_cm, dims=1);
end

truncate_pca(pca::MultivariateStats.PCA, outdim::Int) =
    MultivariateStats.PCA(pca.mean, pca.proj[:, 1:outdim], pca.prinvars[1:outdim], pca.tvar)

function gene_composition_transformation(count_matrix::Matrix{Float64}, confidence::Vector{Float64}=ones(size(count_matrix, 2));
        sample_size::Int=10000, seed::Int=42, method::Symbol=:umap, return_all::Bool=false, n_pcs::Int=15, kwargs...)
    sample_size = min(sample_size, size(count_matrix, 2))
    Random.seed!(seed)

    if n_pcs <= 3
        n_pcs = 3
        method = :pca
    end

    pca = fit(MultivariateStats.PCA, count_matrix, maxoutdim=n_pcs)
    pca_res = transform(truncate_pca(pca, 2), count_matrix);
    sample_ids = select_ids_uniformly(pca_res', confidence, n=sample_size)

    count_matrix_sample = count_matrix[:,sample_ids]

    emb = nothing
    if method == :umap
        emb = fit(UmapFit, count_matrix_sample, pca; n_components=3, kwargs...);
    elseif method == :pca
        emb = truncate_pca(pca, 3);
    else
        error("Unknown method: '$method'")
    end

    if return_all
        return (emb, sample_ids, pca)
    end

    return emb
end

gene_composition_colors(count_matrix::Matrix{Float64}, transformation::UmapFit; kwargs...) =
    gene_composition_colors!(MultivariateStats.transform(transformation, count_matrix); kwargs...)

gene_composition_colors(embedding::Matrix{Float64}; kwargs...) =
    gene_composition_colors!(deepcopy(embedding); kwargs...)

function gene_composition_colors!(embedding::Matrix{Float64}; lrange::Tuple{<:Real, <:Real}=(10, 90), log_colors::Bool=false, trim_frac::Float64=0.0125)
    @assert size(embedding, 1) == 3 "Color embedding must have exactly 3 rows"

    embedding .-= mapslices(x -> quantile(x, trim_frac), embedding, dims=2)
    embedding .= max.(embedding, 0.0)
    max_val = quantile(vec(embedding), 1.0 - trim_frac)
    embedding ./= max_val;
    embedding .= min.(embedding, 1.0)

    if log_colors
        embedding .= log10.(embedding .+ max(quantile(vec(embedding), 0.05), 1e-3))
        embedding .-= minimum(embedding, dims=2)
        embedding ./= maximum(embedding, dims=2)
    end

    embedding[1,:] .*= lrange[2] - lrange[1]
    embedding[1,:] .+= lrange[1]

    embedding[2:3,:] .-= 0.5
    embedding[2:3,:] .*= 200

    return vec(mapslices(col -> Colors.Lab(col...), embedding, dims=1))
end

function gene_composition_colors(df_spatial::T where T <: AbstractDataFrame, k::Int; kwargs...)
    neighb_cm = neighborhood_count_matrix(df_spatial, k);
    confidence = (:confidence in propertynames(df_spatial)) ? df_spatial.confidence : ones(size(df_spatial, 1))
    transformation = gene_composition_transformation(neighb_cm, confidence)
    return gene_composition_colors(neighb_cm, transformation; kwargs...)
end