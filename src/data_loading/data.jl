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

function read_spatial_df(data_path::String; x_col::Symbol=:x, y_col::Symbol=:y, z_col::Union{Symbol, Nothing}=:z, gene_col::Union{Symbol, Nothing}=:gene, filter_cols::Bool=false)
    df_spatial = DataFrame(CSV.File(data_path), copycols=false);

    if (z_col === :z) && !(:z in propertynames(df_spatial))
        z_col = nothing
    end

    for (cn, co) in zip((:x, :y, :z, :gene), (x_col, y_col, z_col, gene_col))
        if co === nothing
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

function load_df(data_path::String; min_molecules_per_gene::Int=0, exclude_genes::Vector{String}=String[], kwargs...)
    df_spatial = read_spatial_df(data_path; kwargs...)

    gene_counts = StatsBase.countmap(df_spatial[!, :gene]);
    large_genes = Set{String}(collect(keys(gene_counts))[collect(values(gene_counts)) .>= min_molecules_per_gene]);
    df_spatial = df_spatial[in.(df_spatial.gene, Ref(large_genes)),:];

    if length(exclude_genes) > 0
        exclude_genes = match_gene_names(exclude_genes, unique(df_spatial.gene))
        df_spatial = df_spatial[.!in.(df_spatial.gene, Ref(exclude_genes)),:];
        @info "Excluding genes: " * join(sort(collect(exclude_genes)), ", ")
    end

    df_spatial[!, :x] = Array{Float64, 1}(df_spatial[!, :x])
    df_spatial[!, :y] = Array{Float64, 1}(df_spatial[!, :y])
    df_spatial[!, :gene], gene_names = encode_genes(df_spatial[!, :gene]);
    return df_spatial, gene_names;
end