import CSV
using DataFrames

function read_spatial_df(data_path; x_col::Symbol=:x, y_col::Symbol=:y, gene_col::Union{Symbol, Nothing}=:gene, filter_cols::Bool=false)
    df_spatial = CSV.read(data_path) |> DataFrame;

    for (cn, co) in zip((:x, :y, :gene), (x_col, y_col, gene_col))
        if co == nothing
            continue
        end

        if (cn in names(df_spatial)) & (cn != co)
            cr = Symbol(String(cn) * "_reserved")
            if cr in names(df_spatial)
                DataFrames.select!(df_spatial, DataFrames.Not(cr))
            end
            DataFrames.rename!(df_spatial, cn => cr);
        end

        DataFrames.rename!(df_spatial, co => cn);
    end

    if gene_col === nothing
        if filter_cols
            df_spatial = df_spatial[:, [:x, :y]]
        end
    else
        if filter_cols
            df_spatial = df_spatial[:, [:x, :y, :gene]]
        end
        df_spatial[!, :gene] = ["$g" for g in df_spatial.gene]
    end

    return df_spatial
end


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