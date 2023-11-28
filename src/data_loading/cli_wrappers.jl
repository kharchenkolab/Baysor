import HDF5
import JSON
import LinearAlgebra: Adjoint
using ProgressMeter
using StatsBase: denserank

PolygonCollection = Dict{TI, Matrix{TV}} where {TV <: Real, TI <: Union{String, Integer}}
PolygonStack = Dict{String, Dict{TI, Matrix{TV}}} where {TV <: Real, TI <: Union{String, Integer}}
Polygons = Union{PolygonCollection, PolygonStack}

function parse_prior_assignment(pos_data::Matrix{Float64}, prior_segmentation::Vector; col_name::Symbol, min_molecules_per_segment::Int, min_mols_per_cell::Int)
    try
        prior_segmentation = Int.(prior_segmentation);
    catch
        error("The prior segmentation column '$col_name' must be of integer type")
    end
    if minimum(prior_segmentation) < 0
        error("The prior segmentation column '$col_name' must not contain negative numbers")
    end

    prior_segmentation = denserank(prior_segmentation) .- 1
    filter_segmentation_labels!(prior_segmentation, min_molecules_per_segment=min_molecules_per_segment)
    scale, scale_std = estimate_scale_from_assignment(pos_data, prior_segmentation; min_mols_per_cell=min_mols_per_cell)

    return prior_segmentation, scale, scale_std
end

# TODO: move it to CLI
function load_prior_segmentation(file::String, pos_data::Matrix{Float64}; min_molecules_per_segment::Int)
    @info "Loading segmentation mask..."

    prior_seg_labels = load_segmentation_mask(file)
    prior_segmentation = Int.(staining_value_per_transcript(pos_data, prior_seg_labels));
    filter_segmentation_labels!(prior_seg_labels, prior_segmentation; min_molecules_per_segment=min_molecules_per_segment)
    GC.gc()

    scale, scale_std = estimate_scale_from_centers(prior_seg_labels)
    @info "Done"

    return prior_segmentation, prior_seg_labels, scale, scale_std
end

### Data

save_segmented_df(segmented_df::DataFrame, file::String) =
    CSV.write(file, segmented_df[sortperm(segmented_df.molecule_id), :])

save_cell_stat_df(cell_stat_df::DataFrame, file::String) =
    CSV.write(file, cell_stat_df)

save_molecule_counts(counts::DataFrame, file::Union{String, Nothing}) =
    save_molecule_counts(counts, names(counts), file)

function save_molecule_counts(counts::Union{AbstractDataFrame, AbstractMatrix}, gene_names::Vector{String}, file::Union{String, Nothing})::String
    st = isnothing(file) ? IOBuffer() : open(file, "w")

    println(st, join(gene_names, "\t"))
    for r in eachrow(counts)
        println(st, join(["$v" for v in r], '\t'))
    end

    (file === nothing) && return String(take!(st))

    close(st)
    return file
end


function polygons_to_geojson(polygons::PolygonCollection)
    geoms = [Dict("type" => "Polygon", "coordinates" => [collect.(eachrow(p))], "cell" => c) for (c,p) in polygons]
    return Dict("type" => "GeometryCollection", "geometries" => geoms);
end

polygons_to_geojson(polygons::PolygonStack) =
    [merge!(polygons_to_geojson(poly), Dict("z" => k)) for (k,poly) in polygons]

function save_polygons_to_geojson(polygons::Polygons, file::String)
    open(file, "w") do f
        print(f, JSON.json(polygons_to_geojson(polygons)))
    end
end

function save_polygons(polygons::Polygons; format::String, file::String)
    if lowercase(format) == "geojson"
        save_polygons_to_geojson(polygons, file)
    else
        @warn "Unknown polygon format: $(format)"
    end
end

function save_matrix_to_loom!(matrix::AbstractMatrix{<:Real}, fid::HDF5.File; chunk=min.(size(matrix), 64), compress::Int=3)
    fid["matrix", chunk=chunk, compress=compress] = matrix
end

function save_matrix_to_loom!(
        matrix::Union{Adjoint{T, SparseMatrixCSC{T, T2}}, SparseMatrixCSC{T, T2}} where T<:Real where T2 <: Real, fid::HDF5.File;
        chunk=min.(size(matrix), 64), compress::Int=3
    )
    HDF5.create_dataset(
        fid,
        "matrix",
        eltype(matrix),
        size(matrix);
        chunk=chunk,
        filters=[HDF5.Filters.Shuffle(), HDF5.Filters.Deflate(compress)]
    )

    # @showprogress for (r,c,v) in collect(zip(findnz(SparseMatrixCSC(matrix))...))
    #     fid["matrix"][r, c] = v
    # end
    @showprogress for (i,c) in enumerate(eachcol(SparseMatrixCSC(matrix')))
        fid["matrix"][i,:] = c
    end
end

function save_matrix_to_loom(
        matrix::AbstractMatrix{<:Real}; gene_names::Vector{String}, cell_names::Vector{String}, file_path::String,
        row_attrs::Union{Dict{String, T1}, Nothing} where T1=nothing, col_attrs::Union{Dict{String, T2}, Nothing} where T2=nothing,
        kwargs...
    )
    # Specification: https://linnarssonlab.org/loompy/format/index.html

    if length(gene_names) != size(matrix, 2)
        error("Row attribute gene_names has wrong length ($(length(gene_names)) != $(size(matrix, 2)))")
    end

    if length(cell_names) != size(matrix, 1)
        error("Col attribute cell_names has wrong length ($(length(cell_names)) != $(size(matrix, 1)))")
    end

    HDF5.h5open(file_path, "w") do fid
        save_matrix_to_loom!(matrix, fid; kwargs...)
        HDF5.create_group(fid, "row_attrs")
        HDF5.create_group(fid, "col_attrs")
        HDF5.create_group(fid, "attrs")
        fid["attrs"]["LOOM_SPEC_VERSION"] = "3.0.0"
        fid["row_attrs"]["Name"] = gene_names
        fid["col_attrs"]["CellID"] = Float64.(1:length(cell_names))
        fid["col_attrs"]["Name"] = cell_names
        if row_attrs !== nothing
            for (k,v) in row_attrs
                if length(v) == size(matrix, 1)
                    error("Row attribute $(k) has wrong length ($(length(v)) != $(size(matrix, 1)))")
                end
                fid["row_attrs"][k] = v
            end
        end

        if col_attrs !== nothing
            for (k,v) in col_attrs
                if length(v) == size(matrix, 2)
                    error("Row attribute $(k) has wrong length ($(length(v)) != $(size(matrix, 2)))")
                end
                fid["col_attrs"][k] = v
            end
        end
    end;
end

save_matrix_to_loom(cm::AbstractMatrix{<:Real}, gene_names::Vector{String}, cell_stat_df::DataFrame, file_path::String) =
    save_matrix_to_loom(
        cm'; gene_names=gene_names, cell_names=cell_stat_df.cell,
        col_attrs=Dict(String(cn) => cell_stat_df[!,cn] for cn in propertynames(cell_stat_df) if cn != :cell),
        file_path=file_path
    )

using ..Utils: DataOptions

# function load_df(args::Dict; kwargs...)
function load_df(coordinates::String, data_opts::DataOptions; kwargs...)
    exclude_genes = split_string_list(data_opts.exclude_genes)

    # TODO: fix Symbol when Configurations are fixed. Also replace this whole call with `to_dict(data_opts)...`
    df_spatial, gene_names = load_df(
        coordinates; x_col=Symbol(data_opts.x), y_col=Symbol(data_opts.y), z_col=Symbol(data_opts.z), gene_col=Symbol(data_opts.gene),
        min_molecules_per_gene=data_opts.min_molecules_per_gene, exclude_genes=exclude_genes,
        drop_z=data_opts.force_2d, kwargs...
    )

    @info "Loaded $(size(df_spatial, 1)) transcripts"

    if size(df_spatial, 1) != size(unique(df_spatial), 1)
        @warn "$(size(df_spatial, 1) - size(unique(df_spatial), 1)) records are duplicates. You may need to filter them beforehand."
    end

    return df_spatial, gene_names, size(df_spatial, 1)
end

function save_segmentation_results(
        segmented_df::DataFrame, cell_stat_df::DataFrame, cm::SparseMatrixCSC{<:Real, Int}, polygons::Union{Polygons, Nothing},
        out_paths::OutputPaths; poly_format::String, gene_names::Vector{String}, matrix_format::Symbol=:loom
    )
    isempty(out_paths.segmented_df) || save_segmented_df(segmented_df, out_paths.segmented_df);
    isempty(out_paths.cell_stats) || save_cell_stat_df(cell_stat_df, out_paths.cell_stats);

    if matrix_format == :loom
        save_matrix_to_loom(cm, gene_names, cell_stat_df, out_paths.counts)
    elseif matrix_format == :tsv
        save_molecule_counts(cm, gene_names, out_paths.counts)
    else
        error("Unknown matrix format: $(matrix_format). Only :loom and :tsv are supported.")
    end

    if (poly_format !== "false") && (polygons !== nothing)
        save_polygons(polygons; format=poly_format, file=out_paths.polygons)
    end
end

function load_prior_segmentation!(
        path::String, df_spatial::DataFrame, pos_data::Matrix{Float64};
        min_molecules_per_segment::Int, min_mols_per_cell::Int
    )

    length(path) > 0 || error("Prior segmentation file path is empty")
    if path[1] == ':'
        prior_col = Symbol(path[2:end])
        prior_seg, scale, scale_std = parse_prior_assignment(
            pos_data, df_spatial[!, prior_col]; col_name=prior_col,
            min_molecules_per_segment=min_molecules_per_segment, min_mols_per_cell=min_mols_per_cell
        )
        prior_seg_labels = nothing
    else
        prior_seg, prior_seg_labels, scale, scale_std = load_prior_segmentation(
            path, pos_data; min_molecules_per_segment=min_molecules_per_segment
        )
    end

    df_spatial[!, :prior_segmentation] = prior_seg

    return prior_seg_labels, scale, scale_std
end
