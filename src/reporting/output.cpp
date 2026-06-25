#include "baysor/reporting/output.h"
#include "baysor/data_loading/data.h"

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#include <hdf5.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>

namespace baysor {

namespace {

template<class T>
T arrow_unwrap(arrow::Result<T>&& result, const char* ctx) {
    if (!result.ok()) {
        throw std::runtime_error(std::string(ctx) + ": " + result.status().ToString());
    }
    return std::move(*result);
}

void arrow_check(const arrow::Status& status, const char* ctx) {
    if (!status.ok()) {
        throw std::runtime_error(std::string(ctx) + ": " + status.ToString());
    }
}

std::vector<std::array<double, 2>> polygon_vertices(const Eigen::MatrixXd& poly) {
    std::vector<std::array<double, 2>> vertices;
    if (poly.rows() == 2) {
        vertices.reserve(static_cast<size_t>(poly.cols()));
        for (int j = 0; j < poly.cols(); ++j) {
            vertices.push_back({poly(0, j), poly(1, j)});
        }
    } else if (poly.cols() == 2) {
        vertices.reserve(static_cast<size_t>(poly.rows()));
        for (int i = 0; i < poly.rows(); ++i) {
            vertices.push_back({poly(i, 0), poly(i, 1)});
        }
    }
    if (vertices.size() >= 3) {
        const auto& first = vertices.front();
        const auto& last = vertices.back();
        if (std::abs(first[0] - last[0]) > 1e-9 || std::abs(first[1] - last[1]) > 1e-9) {
            vertices.push_back(first);
        }
    }
    return vertices;
}

template<class T>
void append_le(std::string& out, T value) {
    const char* ptr = reinterpret_cast<const char*>(&value);
    out.append(ptr, ptr + sizeof(T));
}

std::string polygon_to_wkb(const Eigen::MatrixXd& poly) {
    const auto vertices = polygon_vertices(poly);
    if (vertices.size() < 4) {
        return {};
    }

    std::string out;
    out.reserve(1 + 4 + 4 + 4 + vertices.size() * 16);
    append_le<std::uint8_t>(out, 1);          // little endian
    append_le<std::uint32_t>(out, 3);         // Polygon
    append_le<std::uint32_t>(out, 1);         // one ring
    append_le<std::uint32_t>(out, static_cast<std::uint32_t>(vertices.size()));
    for (const auto& xy : vertices) {
        append_le<double>(out, xy[0]);
        append_le<double>(out, xy[1]);
    }
    return out;
}

std::shared_ptr<arrow::Table> make_table(
    const std::vector<std::shared_ptr<arrow::Field>>& fields,
    const std::vector<std::shared_ptr<arrow::Array>>& arrays,
    std::shared_ptr<arrow::KeyValueMetadata> metadata = nullptr
) {
    auto schema = std::make_shared<arrow::Schema>(fields, metadata);
    return arrow::Table::Make(schema, arrays);
}

void write_parquet_table(const std::shared_ptr<arrow::Table>& table, const std::string& path) {
    auto sink = arrow_unwrap(arrow::io::FileOutputStream::Open(path), "Open parquet output");
    auto writer = arrow_unwrap(
        parquet::arrow::FileWriter::Open(*table->schema(), arrow::default_memory_pool(), sink),
        "Open parquet writer");
    if (table->schema()->metadata()) {
        arrow_check(writer->AddKeyValueMetadata(table->schema()->metadata()),
                    "Attach parquet metadata");
    }
    arrow_check(writer->WriteTable(*table, std::max<int64_t>(1, std::min<int64_t>(table->num_rows(), 65536))),
                "Write parquet table");
    arrow_check(writer->Close(), "Close parquet writer");
    arrow_check(sink->Close(), "Close parquet output");
}

std::shared_ptr<arrow::Array> build_string_array(const std::vector<std::string>& values) {
    arrow::StringBuilder builder;
    arrow_check(builder.Reserve(static_cast<int64_t>(values.size())), "Reserve string builder");
    for (const auto& v : values) arrow_check(builder.Append(v), "Append string");
    return arrow_unwrap(builder.Finish(), "Finish string array");
}

std::shared_ptr<arrow::Array> build_double_array(const std::vector<double>& values) {
    arrow::DoubleBuilder builder;
    arrow_check(builder.Reserve(static_cast<int64_t>(values.size())), "Reserve double builder");
    for (double v : values) arrow_check(builder.Append(v), "Append double");
    return arrow_unwrap(builder.Finish(), "Finish double array");
}

std::shared_ptr<arrow::Array> build_int32_array(const std::vector<int>& values) {
    arrow::Int32Builder builder;
    arrow_check(builder.Reserve(static_cast<int64_t>(values.size())), "Reserve int32 builder");
    for (int v : values) arrow_check(builder.Append(v), "Append int32");
    return arrow_unwrap(builder.Finish(), "Finish int32 array");
}

std::shared_ptr<arrow::Array> build_bool_array(const std::vector<bool>& values) {
    arrow::BooleanBuilder builder;
    arrow_check(builder.Reserve(static_cast<int64_t>(values.size())), "Reserve bool builder");
    for (bool v : values) arrow_check(builder.Append(v), "Append bool");
    return arrow_unwrap(builder.Finish(), "Finish bool array");
}

std::shared_ptr<arrow::Array> build_binary_array(const std::vector<std::string>& values) {
    arrow::BinaryBuilder builder;
    arrow_check(builder.Reserve(static_cast<int64_t>(values.size())), "Reserve binary builder");
    for (const auto& v : values) {
        arrow_check(builder.Append(reinterpret_cast<const std::uint8_t*>(v.data()),
                                   static_cast<int32_t>(v.size())), "Append binary");
    }
    return arrow_unwrap(builder.Finish(), "Finish binary array");
}

std::shared_ptr<arrow::KeyValueMetadata> make_geoparquet_metadata(const std::string& geometry_name) {
    nlohmann::json geo = {
        {"version", "1.1.0"},
        {"primary_column", geometry_name},
        {"columns", {
            {geometry_name, {
                {"encoding", "WKB"},
                {"geometry_types", nlohmann::json::array({"Polygon"})},
                {"crs", nullptr}
            }}
        }}
    };
    return arrow::key_value_metadata({"geo"}, {geo.dump()});
}

} // namespace

OutputStyle parse_output_style(const std::string& style) {
    if (style == "legacy") return OutputStyle::Legacy;
    if (style == "parquet") return OutputStyle::Parquet;
    throw std::invalid_argument("Unknown output style: " + style);
}

std::string to_string(OutputStyle style) {
    switch (style) {
        case OutputStyle::Legacy: return "legacy";
        case OutputStyle::Parquet: return "parquet";
    }
    return "legacy";
}

static nlohmann::json polygons_to_geojson_json(
    const PolygonCollection& polygons,
    const std::string& format
) {
    nlohmann::json out;

    const bool is_feature = (format != "GeometryCollection");
    if (is_feature) {
        out["type"] = "FeatureCollection";
        out["features"] = nlohmann::json::array();
    } else {
        out["type"] = "GeometryCollection";
        out["geometries"] = nlohmann::json::array();
    }

    for (const auto& [cell_name, poly] : polygons) {
        const int nv = (poly.cols() > 0) ? static_cast<int>(poly.cols()) : 0;
        if (is_feature && nv < 3) {
            continue;
        }

        nlohmann::json ring = nlohmann::json::array();
        for (int j = 0; j < nv; ++j) {
            ring.push_back({poly(0, j), poly(1, j)});
        }
        if (nv > 0) {
            const bool already_closed =
                (std::abs(poly(0, nv - 1) - poly(0, 0)) < 1e-9) &&
                (std::abs(poly(1, nv - 1) - poly(1, 0)) < 1e-9);
            if (!already_closed) {
                ring.push_back({poly(0, 0), poly(1, 0)});
            }
        }

        nlohmann::json geom = {
            {"type", "Polygon"},
            {"coordinates", nlohmann::json::array({ring})}
        };

        if (is_feature) {
            out["features"].push_back({
                {"type", "Feature"},
                {"id", cell_name},
                {"geometry", geom},
                {"properties", {
                    {"cell", cell_name}
                }}
            });
        } else {
            out["geometries"].push_back({
                {"type", "Polygon"},
                {"coordinates", nlohmann::json::array({ring})},
                {"cell", cell_name}
            });
        }
    }

    return out;
}

// ============================================================================
// save_matrix_to_loom  —  C API implementation (avoids C++ ABI issues)
// ============================================================================
//
// Loom spec: https://linnarssonlab.org/loompy/format/index.html
//
// Layout produced (mirrors Julia's save_matrix_to_loom):
//   /matrix            float32  n_genes × n_cells  (chunked, shuffle+deflate)
//   /row_attrs/Name    variable-length UTF-8 strings  [n_genes]
//   /col_attrs/Name    variable-length UTF-8 strings  [n_cells]
//   /col_attrs/CellID  float64  [n_cells]   (1-based indices)
//   /col_attrs/<key>   float64[] or vlen-string[] per extra attribute
//   /attrs/LOOM_SPEC_VERSION  "3.0.0"

// RAII wrappers so we don't leak handles on exceptions.
struct HidGuard {
    hid_t id;
    explicit HidGuard(hid_t h) : id(h) {}
    ~HidGuard() { if (id >= 0) H5Idec_ref(id); }
    operator hid_t() const { return id; }
};

static void check(hid_t id, const char* ctx) {
    if (id < 0) throw std::runtime_error(std::string("HDF5 error: ") + ctx);
}
static void check(herr_t err, const char* ctx) {
    if (err < 0) throw std::runtime_error(std::string("HDF5 error: ") + ctx);
}

// Write a 1-D array of variable-length UTF-8 strings.
static void write_vlen_strings(hid_t grp, const char* name,
                                const std::vector<std::string>& strs) {
    hid_t stype = H5Tcopy(H5T_C_S1);
    H5Tset_size(stype, H5T_VARIABLE);
    H5Tset_cset(stype, H5T_CSET_UTF8);

    hsize_t dim = strs.size();
    hid_t space = H5Screate_simple(1, &dim, nullptr);
    hid_t ds = H5Dcreate2(grp, name, stype, space,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    check(ds, name);

    std::vector<const char*> ptrs(strs.size());
    for (size_t i = 0; i < strs.size(); ++i) ptrs[i] = strs[i].c_str();
    check(H5Dwrite(ds, stype, H5S_ALL, H5S_ALL, H5P_DEFAULT, ptrs.data()), name);

    H5Dclose(ds);
    H5Sclose(space);
    H5Tclose(stype);
}

// Write a 1-D array of float64.
static void write_doubles(hid_t grp, const char* name,
                           const std::vector<double>& vals) {
    hsize_t dim = vals.size();
    hid_t space = H5Screate_simple(1, &dim, nullptr);
    hid_t ds = H5Dcreate2(grp, name, H5T_NATIVE_DOUBLE, space,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    check(ds, name);
    check(H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals.data()), name);
    H5Dclose(ds);
    H5Sclose(space);
}

void save_matrix_to_loom(
    const Eigen::SparseMatrix<float, Eigen::RowMajor>& matrix,   // n_cells × n_genes
    const std::vector<std::string>& gene_names,
    const std::vector<std::string>& cell_names,
    const std::string& path,
    const LoomColAttrs& col_attrs
) {
    int n_cells = static_cast<int>(matrix.rows());
    int n_genes = static_cast<int>(matrix.cols());

    if (static_cast<int>(gene_names.size()) != n_genes)
        throw std::runtime_error("save_matrix_to_loom: gene_names length mismatch");
    if (static_cast<int>(cell_names.size()) != n_cells)
        throw std::runtime_error("save_matrix_to_loom: cell_names length mismatch");

    // Open file
    hid_t fid = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    check(fid, "H5Fcreate");

    // ---- /matrix (float32, chunked, shuffle + deflate-3) ----
    // We transpose it here to match the required Loom spec 
    Eigen::SparseMatrix<float, Eigen::RowMajor> gene_major = matrix.transpose();
    {
        hsize_t dims[2]  = { static_cast<hsize_t>(n_genes),
                             static_cast<hsize_t>(n_cells) };
        constexpr int row_block = 128;
        hsize_t chunk[2] = { static_cast<hsize_t>(std::min(n_genes, row_block)),
                             static_cast<hsize_t>(std::min(n_cells, 256)) };

        hid_t space  = H5Screate_simple(2, dims, nullptr);
        hid_t plist  = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist, 2, chunk);
        H5Pset_shuffle(plist);
        H5Pset_deflate(plist, 3);

        hid_t ds = H5Dcreate2(fid, "matrix", H5T_NATIVE_FLOAT, space,
                               H5P_DEFAULT, plist, H5P_DEFAULT);
        check(ds, "create /matrix");

        // Write in gene-row blocks to reduce HDF5 call overhead without
        // materializing the full dense matrix.
        std::vector<float> block_buf(static_cast<size_t>(row_block) * static_cast<size_t>(n_cells), 0.0f);
        hsize_t block_offset[2] = {0, 0};

        for (int block_start = 0; block_start < n_genes; block_start += row_block) {
            int n_block_rows = std::min(row_block, n_genes - block_start);
            std::fill(block_buf.begin(), block_buf.begin() + static_cast<size_t>(n_block_rows) * static_cast<size_t>(n_cells), 0.0f);
            for (int local_row = 0; local_row < n_block_rows; ++local_row) {
                int row = block_start + local_row;
                float* row_ptr = block_buf.data() + static_cast<size_t>(local_row) * static_cast<size_t>(n_cells);
                for (Eigen::SparseMatrix<float, Eigen::RowMajor>::InnerIterator it(gene_major, row); it; ++it) {
                    row_ptr[it.col()] = it.value();
                }
            }

            hsize_t block_dims[2] = {
                static_cast<hsize_t>(n_block_rows),
                static_cast<hsize_t>(n_cells)
            };
            block_offset[0] = static_cast<hsize_t>(block_start);

            hid_t file_space = H5Dget_space(ds);
            H5Sselect_hyperslab(file_space, H5S_SELECT_SET,
                                block_offset, nullptr, block_dims, nullptr);
            hid_t mem_space = H5Screate_simple(2, block_dims, nullptr);
            H5Dwrite(ds, H5T_NATIVE_FLOAT, mem_space, file_space, H5P_DEFAULT,
                     block_buf.data());
            H5Sclose(mem_space);
            H5Sclose(file_space);
        }

        H5Dclose(ds);
        H5Pclose(plist);
        H5Sclose(space);
    }

    // ---- /attrs ----
    {
        hid_t grp = H5Gcreate2(fid, "attrs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(grp, "create /attrs");
        write_vlen_strings(grp, "LOOM_SPEC_VERSION", {"3.0.0"});
        H5Gclose(grp);
    }

    // ---- /row_attrs ----
    {
        hid_t grp = H5Gcreate2(fid, "row_attrs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(grp, "create /row_attrs");
        write_vlen_strings(grp, "Name", gene_names);
        H5Gclose(grp);
    }

    // ---- /col_attrs ----
    {
        hid_t grp = H5Gcreate2(fid, "col_attrs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(grp, "create /col_attrs");

        write_vlen_strings(grp, "Name", cell_names);

        std::vector<double> cell_ids(n_cells);
        for (int i = 0; i < n_cells; ++i) cell_ids[i] = i + 1.0;
        write_doubles(grp, "CellID", cell_ids);

        for (const auto& [key, val] : col_attrs) {
            if (std::holds_alternative<std::vector<std::string>>(val)) {
                write_vlen_strings(grp, key.c_str(),
                                   std::get<std::vector<std::string>>(val));
            } else {
                write_doubles(grp, key.c_str(),
                              std::get<std::vector<double>>(val));
            }
        }

        H5Gclose(grp);
    }

    H5Fclose(fid);
}

void save_matrix_to_loom(
    const Eigen::SparseMatrix<float>& matrix,
    const std::vector<std::string>& gene_names,
    const std::vector<std::string>& cell_names,
    const std::string& path,
    const LoomColAttrs& col_attrs
) {
    Eigen::SparseMatrix<float, Eigen::RowMajor> row_major(matrix);
    save_matrix_to_loom(row_major, gene_names, cell_names, path, col_attrs);
}

void save_matrix_to_10x_h5(
    const Eigen::SparseMatrix<float>& matrix,   // n_cells x n_genes
    const std::vector<std::string>& gene_names,
    const std::vector<std::string>& cell_names,
    const std::string& path
) {
    const int n_cells = static_cast<int>(matrix.rows());
    const int n_genes = static_cast<int>(matrix.cols());

    Eigen::SparseMatrix<float> feature_by_cell = matrix.transpose(); // n_genes x n_cells
    feature_by_cell.makeCompressed();

    std::vector<int32_t> data(feature_by_cell.nonZeros());
    std::vector<int32_t> indices(feature_by_cell.nonZeros());
    std::vector<int64_t> indptr(n_cells + 1);
    for (int k = 0; k < feature_by_cell.nonZeros(); ++k) {
        data[k] = static_cast<int32_t>(std::lround(feature_by_cell.valuePtr()[k]));
        indices[k] = static_cast<int32_t>(feature_by_cell.innerIndexPtr()[k]);
    }
    for (int k = 0; k < n_cells + 1; ++k) {
        indptr[k] = static_cast<int64_t>(feature_by_cell.outerIndexPtr()[k]);
    }
    std::vector<int64_t> shape = {
        static_cast<int64_t>(n_genes), static_cast<int64_t>(n_cells)
    };

    hid_t fid = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    check(fid, "H5Fcreate 10x h5");

    hid_t matrix_grp = H5Gcreate2(fid, "matrix", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    check(matrix_grp, "create /matrix");
    hid_t features_grp = H5Gcreate2(matrix_grp, "features", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    check(features_grp, "create /matrix/features");

    write_vlen_strings(matrix_grp, "barcodes", cell_names);
    write_vlen_strings(features_grp, "id", gene_names);
    write_vlen_strings(features_grp, "name", gene_names);
    write_vlen_strings(features_grp, "feature_type", std::vector<std::string>(gene_names.size(), "Gene Expression"));
    write_vlen_strings(features_grp, "genome", std::vector<std::string>(gene_names.size(), ""));

    auto write_int32 = [](hid_t grp, const char* name, const std::vector<int32_t>& vals) {
        hsize_t dim = vals.size();
        hid_t space = H5Screate_simple(1, &dim, nullptr);
        hid_t ds = H5Dcreate2(grp, name, H5T_NATIVE_INT32, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(ds, name);
        check(H5Dwrite(ds, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals.data()), name);
        H5Dclose(ds);
        H5Sclose(space);
    };
    auto write_int64 = [](hid_t grp, const char* name, const std::vector<int64_t>& vals) {
        hsize_t dim = vals.size();
        hid_t space = H5Screate_simple(1, &dim, nullptr);
        hid_t ds = H5Dcreate2(grp, name, H5T_NATIVE_INT64, space,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        check(ds, name);
        check(H5Dwrite(ds, H5T_NATIVE_INT64, H5S_ALL, H5S_ALL, H5P_DEFAULT, vals.data()), name);
        H5Dclose(ds);
        H5Sclose(space);
    };

    write_int32(matrix_grp, "data", data);
    write_int32(matrix_grp, "indices", indices);
    write_int64(matrix_grp, "indptr", indptr);
    write_int64(matrix_grp, "shape", shape);

    H5Gclose(features_grp);
    H5Gclose(matrix_grp);
    H5Fclose(fid);
}

// ============================================================================
// get_output_paths
// ============================================================================

OutputPaths get_output_paths(const std::string& base_path,
                              OutputStyle style,
                              const std::string& count_matrix_format) {
    // Strip trailing slash if present
    std::string base = base_path;
    while (!base.empty() && (base.back() == '/' || base.back() == '\\')) base.pop_back();

    OutputPaths p;
    if (style == OutputStyle::Parquet) {
        p.segmented_df       = base + "/molecules.parquet";
        p.cell_stats         = base + "/cells.parquet";
        p.diagnostic_report  = base + "/diagnostic_report.html";
        p.molecule_plot      = base + "/segmentation_plot.html";
        p.polygons_2d        = base + "/cell_boundaries.parquet";
        p.polygons_3d        = base + "/cell_boundaries_3d.parquet";
        p.params_dump        = base + "/run_params.toml";
        p.log_file           = base + "/run.log";
        p.counts             = base + "/feature_matrix.h5";
        return p;
    }
    p.segmented_df       = base + "/segmentation.csv";
    p.cell_stats         = base + "/segmentation_cell_stats.csv";
    p.diagnostic_report  = base + "/diagnostic_report.html";
    p.molecule_plot      = base + "/segmentation_plot.html";
    p.polygons_2d        = base + "/segmentation_polygons_2d.json";
    p.polygons_3d        = base + "/segmentation_polygons_3d.json";
    p.params_dump        = base + "/segmentation_params.dump.toml";
    p.log_file           = base + "/segmentation_log.log";

    if (count_matrix_format == "tsv") {
        p.counts = base + "/segmentation_counts.tsv";
    } else {
        p.counts = base + "/segmentation_counts.loom";
    }
    return p;
}

// ============================================================================
// save_segmented_df — per-molecule CSV
// ============================================================================

void save_segmented_df(const MoleculeData& data,
                       const std::vector<int>& assignment,
                       const std::vector<std::string>& gene_names,
                       const std::string& path,
                       const std::vector<std::string>* ncv_color,
                       const std::vector<double>* assignment_confidence,
                       const std::vector<int>* cluster) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("save_segmented_df: cannot open " + path);

    bool has_z    = data.is_3d();
    bool has_conf = !data.confidence.empty();
    bool has_col  = ncv_color && !ncv_color->empty();
    bool has_ac   = assignment_confidence && !assignment_confidence->empty();
    bool has_cl   = cluster && !cluster->empty();
    bool has_tid  = data.source_transcript_id.size() == static_cast<size_t>(data.n_molecules());

    // Header
    if (has_tid) f << "transcript_id,";
    f << "cell,gene,x,y";
    if (has_z)    f << ",z";
    if (has_conf) f << ",confidence";
    if (has_cl)   f << ",cluster";
    if (has_col)  f << ",ncv_color";
    if (has_ac)   f << ",assignment_confidence";
    f << ",is_noise\n";

    int n = data.n_molecules();
    for (int i = 0; i < n; ++i) {
        int cell = assignment[i];
        std::string cell_name = (cell > 0) ? ("cell_" + std::to_string(cell)) : "0";
        int g = data.gene[i];
        std::string gene_name = (g > 0 && g <= static_cast<int>(gene_names.size()))
                                ? gene_names[g - 1] : std::to_string(g);

        if (has_tid) f << data.source_transcript_id[i] << ',';
        f << cell_name << ',' << gene_name << ',' << data.x[i] << ',' << data.y[i];
        if (has_z)    f << ',' << data.z[i];
        if (has_conf) f << ',' << data.confidence[i];
        if (has_cl)   f << ',' << (*cluster)[i];
        if (has_col)  f << ',' << (*ncv_color)[i];
        if (has_ac)   f << ',' << (*assignment_confidence)[i];
        if (has_tid) {
            f << ',' << (cell == 0 ? "true" : "false") << '\n';
        } else {
            f << ',' << (cell == 0 ? 1 : 0) << '\n';
        }
    }
}

void save_segmented_df_parquet(const MoleculeData& data,
                               const std::vector<int>& assignment,
                               const std::vector<std::string>& gene_names,
                               const std::string& path,
                               const std::vector<std::string>* ncv_color,
                               const std::vector<double>* assignment_confidence,
                               const std::vector<int>* cluster) {
    const bool has_z = data.is_3d();
    const bool has_conf = !data.confidence.empty();
    const bool has_col = ncv_color && !ncv_color->empty();
    const bool has_ac = assignment_confidence && !assignment_confidence->empty();
    const bool has_cl = cluster && !cluster->empty();
    const int n = data.n_molecules();

    std::vector<std::string> cells(n);
    std::vector<std::string> genes(n);
    std::vector<bool> is_noise(n);
    for (int i = 0; i < n; ++i) {
        int cell = assignment[i];
        cells[i] = (cell > 0) ? ("cell_" + std::to_string(cell)) : "0";
        int g = data.gene[i];
        genes[i] = (g > 0 && g <= static_cast<int>(gene_names.size()))
                 ? gene_names[g - 1] : std::to_string(g);
        is_noise[i] = (cell == 0);
    }

    std::vector<std::shared_ptr<arrow::Field>> fields = {
        arrow::field("cell", arrow::utf8()),
        arrow::field("gene", arrow::utf8()),
        arrow::field("x", arrow::float64()),
        arrow::field("y", arrow::float64())
    };
    std::vector<std::shared_ptr<arrow::Array>> arrays = {
        build_string_array(cells),
        build_string_array(genes),
        build_double_array(data.x),
        build_double_array(data.y)
    };
    if (has_z) {
        fields.push_back(arrow::field("z", arrow::float64()));
        arrays.push_back(build_double_array(data.z));
    }
    if (has_conf) {
        fields.push_back(arrow::field("confidence", arrow::float64()));
        arrays.push_back(build_double_array(data.confidence));
    }
    if (has_cl) {
        fields.push_back(arrow::field("cluster", arrow::int32()));
        arrays.push_back(build_int32_array(*cluster));
    }
    if (has_col) {
        fields.push_back(arrow::field("ncv_color", arrow::utf8()));
        arrays.push_back(build_string_array(*ncv_color));
    }
    if (has_ac) {
        fields.push_back(arrow::field("assignment_confidence", arrow::float64()));
        arrays.push_back(build_double_array(*assignment_confidence));
    }
    fields.push_back(arrow::field("is_noise", arrow::boolean()));
    arrays.push_back(build_bool_array(is_noise));

    write_parquet_table(make_table(fields, arrays), path);
}

// ============================================================================
// save_cell_stat_df — per-cell CSV
// ============================================================================

void save_cell_stat_df(const Eigen::MatrixXd& stats,
                       const std::vector<std::string>& cell_names,
                       const std::vector<std::string>& col_names,
                       const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("save_cell_stat_df: cannot open " + path);

    f << "cell";
    for (const auto& c : col_names) f << ',' << c;
    f << '\n';

    int nr = static_cast<int>(stats.rows());
    int nc = static_cast<int>(stats.cols());
    for (int r = 0; r < nr; ++r) {
        f << cell_names[r];
        for (int c = 0; c < nc; ++c) f << ',' << stats(r, c);
        f << '\n';
    }
}

void save_cell_stat_df_parquet(const Eigen::MatrixXd& stats,
                               const std::vector<std::string>& cell_names,
                               const std::vector<std::string>& col_names,
                               const std::string& path) {
    std::vector<std::shared_ptr<arrow::Field>> fields;
    std::vector<std::shared_ptr<arrow::Array>> arrays;

    fields.push_back(arrow::field("cell", arrow::utf8()));
    arrays.push_back(build_string_array(cell_names));

    for (int c = 0; c < stats.cols(); ++c) {
        std::vector<double> col(stats.rows());
        for (int r = 0; r < stats.rows(); ++r) col[r] = stats(r, c);
        fields.push_back(arrow::field(col_names[c], arrow::float64()));
        arrays.push_back(build_double_array(col));
    }

    write_parquet_table(make_table(fields, arrays), path);
}

// ============================================================================
// save_matrix_to_tsv
// ============================================================================

void save_matrix_to_tsv(const Eigen::SparseMatrix<double>& matrix,
                         const std::vector<std::string>& gene_names,
                         const std::vector<std::string>& cell_names,
                         const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("save_matrix_to_tsv: cannot open " + path);

    const int n_cells = static_cast<int>(matrix.rows());
    const int n_genes = static_cast<int>(matrix.cols());
    if (static_cast<int>(gene_names.size()) != n_genes)
        throw std::runtime_error("save_matrix_to_tsv: gene_names length mismatch");
    if (static_cast<int>(cell_names.size()) != n_cells)
        throw std::runtime_error("save_matrix_to_tsv: cell_names length mismatch");

    // Header: cell names
    f << "gene";
    for (const auto& cn : cell_names) f << '\t' << cn;
    f << '\n';

    // Write one row per gene. Internally counts are stored as n_cells x n_genes,
    // so transpose first to get a gene x cell sparse matrix with cheap row access.
    Eigen::SparseMatrix<double, Eigen::RowMajor> rm(matrix.transpose());
    for (int gi = 0; gi < n_genes; ++gi) {
        f << gene_names[gi];
        // Dense row output (sparse matrix — iterate nonzeros + fill zeros).
        std::vector<double> row(cell_names.size(), 0.0);
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(rm, gi); it; ++it) {
            row[it.col()] = it.value();
        }
        for (double v : row) f << '\t' << v;
        f << '\n';
    }
}

// ============================================================================
// save_polygons_geojson
// ============================================================================

void save_polygons_geojson(const PolygonCollection& polygons,
                            const std::string& path,
                            const std::string& format) {
    if (format == "none" || polygons.empty()) return;

    std::ofstream f(path);
    if (!f) throw std::runtime_error("save_polygons_geojson: cannot open " + path);
    f << polygons_to_geojson_json(polygons, format).dump();
}

void save_polygons_geoparquet(const PolygonCollection& polygons,
                              const std::string& path,
                              const std::string& geometry_name) {
    if (polygons.empty()) return;

    std::vector<std::string> cells;
    std::vector<int> n_vertices;
    std::vector<std::string> wkbs;
    cells.reserve(polygons.size());
    n_vertices.reserve(polygons.size());
    wkbs.reserve(polygons.size());

    for (const auto& [cell_name, poly] : polygons) {
        auto vertices = polygon_vertices(poly);
        if (vertices.size() < 4) continue;
        std::string wkb = polygon_to_wkb(poly);
        if (wkb.empty()) continue;
        cells.push_back(cell_name);
        n_vertices.push_back(static_cast<int>(vertices.size() - 1));
        wkbs.push_back(std::move(wkb));
    }
    if (cells.empty()) return;

    auto metadata = make_geoparquet_metadata(geometry_name);
    auto table = make_table(
        {
            arrow::field("cell", arrow::utf8()),
            arrow::field("n_vertices", arrow::int32()),
            arrow::field(geometry_name, arrow::binary())
        },
        {
            build_string_array(cells),
            build_int32_array(n_vertices),
            build_binary_array(wkbs)
        },
        metadata
    );
    write_parquet_table(table, path);
}

void save_polygon_stack_geojson(const PolygonStack& polygons,
                                const OutputPaths& out_paths,
                                const std::string& format) {
    if (format == "none" || polygons.empty()) return;

    nlohmann::json by_layer = nlohmann::json::object();
    bool has_3d = false;

    for (const auto& [layer_name, poly] : polygons) {
        if (layer_name == "2d") {
            save_polygons_geojson(poly, out_paths.polygons_2d, format);
            continue;
        }
        by_layer[layer_name] = polygons_to_geojson_json(poly, format);
        has_3d = true;
    }

    if (!has_3d) return;

    std::ofstream f(out_paths.polygons_3d);
    if (!f) throw std::runtime_error("save_polygon_stack_geojson: cannot open " + out_paths.polygons_3d);
    f << by_layer.dump();
}

void save_polygon_stack_geoparquet(const PolygonStack& polygons,
                                   const OutputPaths& out_paths,
                                   const std::string& geometry_name) {
    if (polygons.empty()) return;

    PolygonCollection combined_2d;
    std::vector<std::string> cells;
    std::vector<std::string> layers;
    std::vector<std::string> wkbs;
    std::vector<int> n_vertices;

    for (const auto& [layer_name, poly] : polygons) {
        if (layer_name == "2d") {
            combined_2d = poly;
            continue;
        }
        for (const auto& [cell_name, geom] : poly) {
            auto vertices = polygon_vertices(geom);
            if (vertices.size() < 4) continue;
            auto wkb = polygon_to_wkb(geom);
            if (wkb.empty()) continue;
            cells.push_back(cell_name);
            layers.push_back(layer_name);
            n_vertices.push_back(static_cast<int>(vertices.size() - 1));
            wkbs.push_back(std::move(wkb));
        }
    }

    if (!combined_2d.empty()) {
        save_polygons_geoparquet(combined_2d, out_paths.polygons_2d, geometry_name);
    }
    if (cells.empty()) return;

    auto metadata = make_geoparquet_metadata(geometry_name);
    auto table = make_table(
        {
            arrow::field("cell", arrow::utf8()),
            arrow::field("layer", arrow::utf8()),
            arrow::field("n_vertices", arrow::int32()),
            arrow::field(geometry_name, arrow::binary())
        },
        {
            build_string_array(cells),
            build_string_array(layers),
            build_int32_array(n_vertices),
            build_binary_array(wkbs)
        },
        metadata
    );
    write_parquet_table(table, out_paths.polygons_3d);
}

} // namespace baysor
