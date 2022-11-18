import HDF5
import JSON

function parse_prior_assignment(pos_data::Matrix{Float64}, prior_segmentation::Vector; col_name::Symbol, min_molecules_per_segment::Int, min_mols_per_cell::Int)
    try
        prior_segmentation = Int.(prior_segmentation);
    catch
        error("The prior segmentation column '$col_name' must be of integer type")
    end
    if minimum(prior_segmentation) < 0
        error("The prior segmentation column '$col_name' must not contain negative numbers")
    end

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
