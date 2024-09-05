using Base.Threads
import FastRandPCA

neighborhood_count_matrix(data::Union{BmmData, DataFrame}, k::Int; kwargs...) =
    neighborhood_count_matrix(position_data(data), composition_data(data), k; kwargs...)

function neighborhood_count_matrix(
        pos_data::Matrix{Float64}, genes::Vector{<:Union{Int, Missing}}, k::Int;
        confidences::Union{Vector{Float64}, Nothing}=nothing,
        n_genes::Int=maximum(skipmissing(genes)), normalize_by_dist::Bool=true, normalize::Bool=true
    )
    if k < 3
        @warn "Too small value of k: $k. Setting it to 3."
        k = 3
    end

    if !normalize
        normalize_by_dist = false
    end

    if (confidences !== nothing) && normalize_by_dist
        @warn "`confidences` are not supported with `normalize_by_dist=true`. Ignoring them."
        confidences = nothing
    end

    k = min(k, size(pos_data, 2))

    neighbors, dists = knn_parallel(KDTree(pos_data), pos_data, k; sorted=true);

    if normalize_by_dist
        # account for problems with points with duplicating coordinates
        med_closest_dist = median(d[findfirst(d .> 1e-15)] for d in dists if any(d .> 1e-15));

        return sphstack([ # Not sure if making it parallel will have large effect, as we have a lot of allocations here
            count_array_sparse(Float32, genes[nns], 1 ./ max.(ds, med_closest_dist); total=n_genes, normalize=normalize)
            for (nns,ds) in zip(neighbors, dists)
        ]);
    end

    s_vecs = Vector{SparseArrays.SparseVector{Float32, Int64}}(undef, length(neighbors))
    if !normalize || (confidences === nothing)
        @threads for i in eachindex(neighbors)
            s_vecs[i] = count_array_sparse(Float32, view(genes, neighbors[i]); total=n_genes, normalize=normalize)
        end
    else
        @threads for i in eachindex(neighbors)
            s_vecs[i] = count_array_sparse(Float32, view(genes, neighbors[i]), view(confidences, neighbors[i]); total=n_genes, normalize=normalize)
        end
    end

    return sphstack(s_vecs);
end

function estimate_gene_vectors(
        count_matrix::AbstractMatrix{<:Real}, gene_ids::Vector{Int};
        n_components::Int, method::Symbol=:ri, per_molecule::Bool=false
    )::Matrix{<:Real}
    if method in (:dense, :sparse)
        mol_vecs, gene_vecs = gene_pca(count_matrix, n_components; method=method)
        return per_molecule ? mol_vecs : gene_vecs
    elseif method == :ri
        gene_vecs = generate_randomized_gene_vectors(count_matrix, gene_ids; n_components=n_components)
        return per_molecule ? Matrix((count_matrix' * gene_vecs)') : Matrix(gene_vecs')
    end

    error("Unknown method: $method. Only :dense, :sparse or :ri are supported")
end

function gene_pca(count_matrix::AbstractMatrix{<:Real}, n_pcs::Int; method::Symbol=:auto)::Tuple{Matrix{<:Real}, Matrix{<:Real}}
    if method == :auto
        method = (prod(size(count_matrix)) < 1e9) ? :dense : :sparse
    end

    if method == :dense
        count_matrix = Matrix(count_matrix)
        pca_trans = MultivariateStats.fit(MultivariateStats.PCA, count_matrix, maxoutdim=n_pcs)
        pca = MultivariateStats.predict(pca_trans, count_matrix)
        return pca, Matrix(pca_trans.proj')
    end

    if method == :sparse
        pca = FastRandPCA.pca(count_matrix, n_pcs)
        proj, _, pca = FastRandPCA.pca(count_matrix, n_pcs);
        return Matrix(pca'), Matrix(proj)
    end

    error("Unknown method: $method. Only :dense, :sparse or :auto are supported")
end

function generate_randomized_gene_vectors(
        neighb_cm::SparseArrays.SparseMatrixCSC{<:Real, Int64}, gene_ids::Vector{Int};
        n_components::Int=50, seed::Int=42, var_clip::Float64=0.05
    )::Matrix{Float64}
    Random.seed!(seed)
    random_vectors_init = randn(size(neighb_cm, 1), n_components);

    coexpr_mat = neighb_cm * neighb_cm'

    if var_clip > 0
        # Genes that mostly co-variate with themselves don't get updated by other genes
        # Som we clip their covariance, which greatly improves convergence
        diag_vals = diag(coexpr_mat)
        total_var = sum(coexpr_mat, dims=1)[:]
        diag_frac = diag_vals ./ total_var

        coexpr_mat[diagind(coexpr_mat)] .= min.(
            (quantile(diag_frac, 1 - var_clip) .* total_var),
            diag_vals
        )
    end

    gene_emb = (coexpr_mat * random_vectors_init) ./ sum(coexpr_mat, dims=2)
    return gene_emb
end

normalize_embedding_to_lab_range(embedding::AbstractMatrix{<:Real}; kwargs...) =
    normalize_embedding_to_lab_range!(deepcopy(embedding); kwargs...)

function normalize_embedding_to_lab_range!(
        embedding::AbstractMatrix{<:Real}; lrange::Tuple{<:Real, <:Real}=(10, 90), log_colors::Bool=false, trim_frac::Float64=0.0125
    )
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

    return embedding
end

function gene_composition_color_embedding(
        pca::AbstractMatrix{<:Real}, confidence::Vector{<:Real}=ones(size(pca, 2));
        normalize::Bool=true, sample_size::Int=10000, seed::Int=42, kwargs...
    )

    sample_size = min(sample_size, size(pca, 2))
    Random.seed!(seed)

    (size(pca, 1) > 3) || error("pca must have at least 3 components")

    sample_ids = select_ids_uniformly(pca[1:2,:]', confidence, n=sample_size)

    ump = fit(UmapFit, pca[:,sample_ids]; n_components=3, kwargs...);
    emb = MultivariateStats.transform(ump, pca)

    if normalize
        # It's useful if we want to average embedding across clusters or take a subsample
        normalize_embedding_to_lab_range!(emb)
    end

    return emb
end

embedding_to_lab(embedding::AbstractMatrix{<:Real}) =
    vec(mapslices(col -> Colors.Lab(col...), embedding, dims=1))

embedding_to_hex(embedding::AbstractMatrix{<:Real}) =
    "#" .* Colors.hex.(embedding_to_lab(embedding))

# It's used in segfree.jl:73
function gene_composition_colors(count_matrix::Matrix{Float64}, transformation::UmapFit; kwargs...)
    gene_composition_colors!(MultivariateStats.transform(transformation, count_matrix); kwargs...)
end

# This is an important function for notebook usage
gene_composition_colors(embedding::AbstractMatrix{<:Real}; kwargs...) =
    gene_composition_colors!(deepcopy(embedding); kwargs...)

function gene_composition_colors!(embedding::AbstractMatrix{<:Real}; kwargs...)
    normalize_embedding_to_lab_range!(embedding; kwargs...)
    return vec(mapslices(col -> Colors.Lab(col...), embedding, dims=1))
end

# It's used in preview.jl:62 and cli_wrappers.jl:76
function gene_composition_colors(df_spatial::DataFrame, k::Int; method::Symbol=:ri, n_components::Int=20, log_counts::Bool=true, kwargs...)
    neighb_cm = neighborhood_count_matrix(df_spatial, k);

    if log_counts
        neighb_cm.nzval .= log.(neighb_cm.nzval .* 10000 .+ 1e-5);
    end

    # TODO: embed gene vectors and then estimate embedded molecule vectors
    mol_vecs = estimate_gene_vectors(neighb_cm, df_spatial.gene; n_components, method, per_molecule=true);

    confidence = (:confidence in propertynames(df_spatial)) ? df_spatial.confidence : ones(size(df_spatial, 1))
    col_emb = gene_composition_color_embedding(mol_vecs, confidence; kwargs...);

    return embedding_to_lab(col_emb)
end
