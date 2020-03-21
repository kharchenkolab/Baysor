import CSV
using DataFrames

function read_spatial_df(data_path; x_col::Symbol=:x, y_col::Symbol=:y, gene_col::Union{Symbol, Nothing}=:gene, filter_cols::Bool=false)
    df_spatial = CSV.read(data_path) |> DataFrame;

    if (:x in names(df_spatial)) & (x_col != :x)
        if :x_reserved in names(df_spatial)
            DataFrames.select!(df_spatial, DataFrames.Not(:x_reserved))
        end
        DataFrames.rename!(df_spatial, :x => :x_reserved);
    end

    if (:y in names(df_spatial)) & (y_col != :y)
        if :y_reserved in names(df_spatial)
            DataFrames.select!(df_spatial, DataFrames.Not(:y_reserved))
        end
        DataFrames.rename!(df_spatial, :y => :y_reserved);
    end

    if gene_col === nothing
        if filter_cols
            df_spatial = df_spatial[:, [x_col, y_col]]
        end
        DataFrames.rename!(df_spatial, x_col => :x, y_col => :y);
    else
        if filter_cols
            df_spatial = df_spatial[:, [x_col, y_col, gene_col]]
        end
        DataFrames.rename!(df_spatial, x_col => :x, y_col => :y, gene_col => :gene);
        df_spatial[!, :gene] = ["$g" for g in df_spatial.gene]
    end

    return df_spatial
end