import HDF5
import JSON

# TODO: move prior segmentation loading logic to DataProcessing

function parse_prior_assignment(df_spatial::DataFrame, args::Dict{String, Any})
    prior_col = Symbol(args["prior_segmentation"][2:end])
    prior_segmentation = df_spatial[!, prior_col]
    try
        prior_segmentation = Int.(prior_segmentation);
    catch
        error("The prior segmentation column '$prior_col' must be of integer type")
    end
    if minimum(prior_segmentation) < 0
        error("The prior segmentation column '$prior_col' must not contain negative numbers")
    end

    filter_segmentation_labels!(prior_segmentation, min_molecules_per_segment=args["min-molecules-per-segment"])
    if args["scale"] === nothing
        pos_data = position_data(df_spatial)
        scale, scale_std = estimate_scale_from_assignment(pos_data, prior_segmentation;
            min_mols_per_cell=args["min-molecules-per-cell"])
    else
        scale, scale_std = args["scale"], args["scale-std"]
    end

    return prior_segmentation, Matrix{Float64}[], scale, scale_std
end

# TODO: move it to CLI
function load_prior_segmentation(df_spatial::DataFrame, args::Dict{String, Any})
    if args["prior_segmentation"][1] == ':'
        return parse_prior_assignment(df_spatial, args)
    end

    @info "Loading segmentation mask..."

    prior_seg_labels = load_segmentation_mask(args["prior_segmentation"])
    prior_segmentation = Int.(staining_value_per_transcript(df_spatial, prior_seg_labels));
    filter_segmentation_labels!(prior_seg_labels, prior_segmentation; min_molecules_per_segment=args["min-molecules-per-segment"])
    GC.gc()

    scale, scale_std = args["scale"], args["scale-std"]
    if args["estimate-scale-from-centers"]
        scale, scale_std = estimate_scale_from_centers(prior_seg_labels)
    end
    @info "Done"

    prior_polygons = nothing
    if args["plot"]
        @info "Estimating prior segmentation polygons..."
        prior_polygons = extract_polygons_from_label_grid(Matrix{UInt32}(prior_seg_labels[1:5:end, 1:5:end]); grid_step=5.0) # subset to save memory and time
        @info "Done"
    end

    return prior_segmentation, prior_polygons, scale, scale_std
end

### Data

save_segmented_df(segmented_df::DataFrame, file::String) =
    CSV.write(file, segmented_df[sortperm(segmented_df.molecule_id), :])

save_cell_stat_df(cell_stat_df::DataFrame, file::String) =
    CSV.write(file, cell_stat_df)

function save_molecule_counts(counts::DataFrame, file::String)
    count_str = join(names(counts), "\t") * "\n" * join([join(["$v" for v in r], '\t') for r in eachrow(counts)], '\n');
    open(file, "w") do f; print(f, count_str) end
end


function polygons_to_geojson(polygons::Vector{Matrix{T}}) where T <:Real
    geoms = [Dict("type" => "Polygon", "coordinates" => [collect.(eachrow(p))]) for p in polygons]
    return Dict("type" => "GeometryCollection", "geometries" => geoms);
end

polygons_to_geojson(polygons::Dict{String, Vector{Matrix{T}}}) where T <:Real =
    [merge!(polygons_to_geojson(poly), Dict("z" => k)) for (k,poly) in polygons]

function save_polygons_to_geojson(polygons::Union{Vector{Matrix{T}}, Dict{String, Vector{Matrix{T}}}}, file::String) where T <: Real
    open(file, "w") do f
        print(f, JSON.json(polygons_to_geojson(polygons)))
    end
end

function save_polygons(polygons::Union{Vector{Matrix{T}}, Dict{String, Vector{Matrix{T}}}}; format::String, file::String) where T <: Real
    if lowercase(format) == "geojson"
        save_polygons_to_geojson(polygons, file)
    else
        @warn "Unknown polygon format: $(format)"
    end
end

function save_matrix_to_loom(matrix; gene_names::Vector{String}, cell_names::Vector{String}, file_path::String,
        row_attrs::Union{Dict{String, T1}, Nothing} where T1=nothing, col_attrs::Union{Dict{String, T2}, Nothing} where T2=nothing)
    # Specification: https://linnarssonlab.org/loompy/format/index.html
    HDF5.h5open(file_path, "w") do fid
        fid["matrix", chunk=(64,64), compress=3] = matrix
        HDF5.create_group(fid, "row_attrs")
        HDF5.create_group(fid, "col_attrs")
        HDF5.create_group(fid, "attrs")
        fid["row_attrs"]["Name"] = gene_names
        fid["col_attrs"]["CellID"] = cell_names
        if row_attrs !== nothing
            for (k,v) in row_attrs
                fid["row_attrs"][k] = v
            end
        end

        if col_attrs !== nothing
            for (k,v) in col_attrs
                fid["col_attrs"][k] = v
            end
        end
    end;
end
