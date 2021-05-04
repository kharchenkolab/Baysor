import CSV
using DataFrames

function read_spatial_df(data_path::String; x_col::Symbol=:x, y_col::Symbol=:y, z_col::Union{Symbol, Nothing}=:z, gene_col::Union{Symbol, Nothing}=:gene, filter_cols::Bool=false)
    df_spatial = DataFrame(CSV.File(data_path), copycols=false);

    if (z_col === :z) && !(:z in propertynames(df_spatial))
        z_col = nothing
    end

    for (cn, co) in zip((:x, :y, :z, :gene), (x_col, y_col, z_col, gene_col))
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

function convert_segmentation_to_counts(genes::Vector{Int}, cell_assignment::Vector{Int}; drop_empty_labels::Bool=false, gene_names::Union{Vector{String}, Nothing}=nothing)
    if drop_empty_labels
        if minimum(cell_assignment) == 0
            cell_assignment = denserank(cell_assignment) .- 1
        else
            cell_assignment = denserank(cell_assignment)
        end
    end

    cm = zeros(Int, maximum(genes), maximum(cell_assignment))
    for i in 1:length(genes)
        if cell_assignment[i] == 0
            continue
        end
        cm[genes[i], cell_assignment[i]] += 1
    end

    if gene_names !== nothing
        cm = DataFrame(cm, [Symbol("$c") for c in 1:size(cm, 2)])
        cm[!, :gene] = gene_names
        cm = cm[:, vcat(end, 1:end-1)]
    end

    return cm
end
