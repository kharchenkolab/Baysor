import CSV
using DataFrames

function read_spatial_df(data_path; x_col::Symbol=:x, y_col::Symbol=:y, gene_col::Union{Symbol, Nothing}=:gene, filter_cols::Bool=false)
    df_spatial = DataFrame!(CSV.File(data_path));

    for (cn, co) in zip((:x, :y, :gene), (x_col, y_col, gene_col))
        if co == nothing
            continue
        end

        if (cn in propertynames(df_spatial)) & (cn != co)
            cr = Symbol(String(cn) * "_reserved")
            if cr in propertynames(df_spatial)
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
