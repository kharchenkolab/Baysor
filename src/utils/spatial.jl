import CSV
using DataFrames

function read_spatial_df(data_path; x_col::Symbol=:x, y_col::Symbol=:y, gene_col::Union{Symbol, Nothing}=:gene, filter_cols::Bool=false)
    df_spatial = CSV.read(data_path);

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
    end

    return df_spatial
end