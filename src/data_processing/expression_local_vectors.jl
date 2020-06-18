using Base.Threads

# This function is a determenistic analogue of sampling. It picks points in a manner that preserves the distributions across x and y.
function select_ids_uniformly(vals::Union{Vector{<:Real}, <:AbstractMatrix{<:Real}}, confidence::Union{Vector{Float64}, Nothing}=nothing; n::Int, confidence_threshold::Float64=0.95)::Vector{Int}
    if n <= 1
        error("n must be > 1")
    end

    high_conf_ids = (confidence===nothing) ? (1:size(vals, 1)) : findall(confidence .>= confidence_threshold)
    if length(high_conf_ids) < n
        @warn "n=$n, which is > length(high_conf_ids) ($(length(high_conf_ids)))"
        return high_conf_ids
    end

    vals = sum(vals, dims=2)[high_conf_ids]
    return high_conf_ids[sortperm(vals)[unique(round.(Int, range(1, length(high_conf_ids), length=n)))]];
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

function gene_composition_transformation(count_matrix::Array{Float64, 2}, confidence::Vector{Float64}=ones(size(count_matrix, 2));
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

function extract_filtered_local_vectors(df_spatial::DataFrame, adjacent_points::Array{Vector{Int}, 1}, k::Int;
        confidence::Union{Vector{Float64}, Nothing}=df_spatial.confidence, clust_per_mol::Union{Vector{Int}, Nothing}=nothing,
        n_vectors::Int=0, confidence_threshold::Float64=0.95, kwargs...)
    if clust_per_mol === nothing
        clust_per_mol = ones(Int, size(df_spatial, 1))
    end

    if confidence !== nothing
        clust_per_mol[confidence .< confidence_threshold] .= 0
    end

    conn_comps, rci, mol_ids_per_cell = get_connected_components_per_label(clust_per_mol, adjacent_points, 1)
    mol_ids_per_comp = vcat([[mol_ids_per_cell[i][ids] for ids in conn_comps[i]] for i in 1:length(conn_comps)]...)
    mol_ids_per_comp = mol_ids_per_comp[length.(mol_ids_per_comp) .>= k];

    pos_data = position_data(df_spatial)
    n_genes = maximum(df_spatial.gene)
    neighb_mat_filt = hcat([neighborhood_count_matrix(pos_data[:, ids], df_spatial.gene[ids], k; n_genes=n_genes, kwargs...) for ids in mol_ids_per_comp]...)

    mol_id_per_vec = vcat(mol_ids_per_comp...)
    if n_vectors > 0
        ids = select_ids_uniformly(Matrix(df_spatial[[:x, :y], mol_id_per_vec]), confidence[mol_id_per_vec]; n=n_vectors)
        mol_id_per_vec = mol_id_per_vec[ids]
        neighb_mat_filt = neighb_mat_filt[:, ids]
    end

    order = sortperm(mol_id_per_vec)
    return neighb_mat_filt[:, order], mol_id_per_vec[order]
end