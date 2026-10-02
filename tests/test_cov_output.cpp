// Tests for the writers of src/reporting/output.cpp: molecule table and cell
// statistics (CSV + Parquet), count matrices (Loom, TSV), GeoJSON/GeoParquet
// polygons, output paths and error paths.

#include <gtest/gtest.h>

#include <nlohmann/json.hpp>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/reader.h>
#include <parquet/file_reader.h>
#include <hdf5.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <csignal>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <map>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#ifndef _WIN32
#include <sys/resource.h>
#endif

#include "baysor/data_loading/data.h"
#include "baysor/reporting/output.h"

#include "test_cov_helpers.h"

namespace {

namespace fs = std::filesystem;

using baysor_test::read_text_file;

class Cov4OutputFiles : public ::testing::Test {
protected:
    std::string path(const std::string& name) const { return tmp_.file(name); }

    baysor_test::TempDir tmp_{"cov4_output"};
};

#ifndef _WIN32
/// Temporarily cap the maximum file size of this process and ignore SIGXFSZ
/// so that writes past the cap fail with an error instead of killing it.
/// Restore happens in the destructor (RAII); ok() reports whether both the
/// rlimit query and the update succeeded.
class ScopedFileSizeLimit {
public:
    explicit ScopedFileSizeLimit(rlim_t bytes) {
        ok_ = (::getrlimit(RLIMIT_FSIZE, &old_) == 0);
        if (!ok_) return;
        keep_ = old_;
        keep_.rlim_cur = bytes;
        ok_ = (::setrlimit(RLIMIT_FSIZE, &keep_) == 0);
        if (!ok_) return;
        old_handler_ = std::signal(SIGXFSZ, SIG_IGN);
    }
    ~ScopedFileSizeLimit() {
        if (!ok_) return;
        std::signal(SIGXFSZ, old_handler_);
        ::setrlimit(RLIMIT_FSIZE, &old_);
    }
    ScopedFileSizeLimit(const ScopedFileSizeLimit&) = delete;
    ScopedFileSizeLimit& operator=(const ScopedFileSizeLimit&) = delete;
    bool ok() const { return ok_; }

private:
    struct rlimit old_ {};
    struct rlimit keep_ {};
    bool ok_ = false;
    void (*old_handler_)(int) = nullptr;
};
#endif  // !_WIN32

template <typename T>
T unwrap_arrow(arrow::Result<T>&& result) {
    EXPECT_TRUE(result.ok()) << result.status().ToString();
    return std::move(*result);
}

std::shared_ptr<arrow::Table> read_parquet(const std::string& path) {
    auto infile = unwrap_arrow(arrow::io::ReadableFile::Open(path));
    auto reader = unwrap_arrow(parquet::arrow::OpenFile(infile, arrow::default_memory_pool()));
    std::shared_ptr<arrow::Table> table;
    auto status = reader->ReadTable(&table);
    EXPECT_TRUE(status.ok()) << status.ToString();
    return table;
}

// Values of column `name` (first chunk).
template <class ArrayT>
auto column(const arrow::Table& t, const std::string& name) {
    using Value = decltype(std::declval<ArrayT>().Value(0));
    std::vector<std::conditional_t<std::is_same_v<Value, std::string_view>, std::string, Value>> out;
    const auto col = t.GetColumnByName(name);
    EXPECT_NE(col, nullptr) << "missing column " << name;
    if (!col) return out;
    const auto& arr = static_cast<const ArrayT&>(*col->chunk(0));
    for (int64_t i = 0; i < arr.length(); ++i) out.emplace_back(arr.Value(i));
    return out;
}

auto string_column(const arrow::Table& t, const std::string& name) { return column<arrow::StringArray>(t, name); }
auto double_column(const arrow::Table& t, const std::string& name) { return column<arrow::DoubleArray>(t, name); }
auto int32_column(const arrow::Table& t, const std::string& name) { return column<arrow::Int32Array>(t, name); }
auto bool_column(const arrow::Table& t, const std::string& name) { return column<arrow::BooleanArray>(t, name); }

std::vector<std::string> field_names(const arrow::Table& t) {
    std::vector<std::string> names;
    for (const auto& f : t.schema()->fields()) names.push_back(f->name());
    return names;
}

/// Read a key/value metadata entry from the parquet file footer.
std::string parquet_footer_kv(const std::string& path, const std::string& key) {
    auto infile = unwrap_arrow(arrow::io::ReadableFile::Open(path));
    auto pq_reader = parquet::ParquetFileReader::Open(infile);
    auto kv = pq_reader->metadata()->key_value_metadata();
    if (!kv) return {};
    for (int i = 0; i < kv->size(); ++i) {
        if (kv->key(i) == key) return kv->value(i);
    }
    return {};
}

/// RAII wrapper for HDF5 handles: closes with the matching H5*close function
/// on scope exit, including early returns.
struct H5Handle {
    hid_t id = -1;
    herr_t (*closer)(hid_t) = nullptr;

    H5Handle() = default;
    H5Handle(hid_t handle, herr_t (*close_fn)(hid_t)) : id(handle), closer(close_fn) {}
    ~H5Handle() {
        if (id >= 0 && closer != nullptr) closer(id);
    }
    H5Handle(const H5Handle&) = delete;
    H5Handle& operator=(const H5Handle&) = delete;
    operator hid_t() const { return id; }
};

std::vector<std::string> h5_read_vlen_strings(hid_t fid, const std::string& ds_path) {
    H5Handle ds(H5Dopen2(fid, ds_path.c_str(), H5P_DEFAULT), H5Dclose);
    EXPECT_GE(ds.id, 0) << ds_path;
    if (ds.id < 0) return {};
    H5Handle space(H5Dget_space(ds), H5Sclose);
    hsize_t dim = 0;
    H5Sget_simple_extent_dims(space, &dim, nullptr);
    H5Handle type(H5Dget_type(ds), H5Tclose);
    std::vector<char*> buf(static_cast<size_t>(dim), nullptr);
    EXPECT_GE(H5Dread(ds, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data()), 0) << ds_path;
    std::vector<std::string> out;
    for (hsize_t i = 0; i < dim; ++i) out.emplace_back(buf[i] ? buf[i] : "");
    // The strings were allocated by the HDF5 library: reclaim them (otherwise
    // AddressSanitizer reports the leak); the RAII handles close afterwards.
#if H5_VERSION_GE(1, 12, 0)
    H5Treclaim(type, space, H5P_DEFAULT, buf.data());
#else
    H5Dvlen_reclaim(type, space, H5P_DEFAULT, buf.data());
#endif
    return out;
}

std::vector<double> h5_read_doubles(hid_t fid, const std::string& ds_path) {
    H5Handle ds(H5Dopen2(fid, ds_path.c_str(), H5P_DEFAULT), H5Dclose);
    EXPECT_GE(ds.id, 0) << ds_path;
    if (ds.id < 0) return {};
    H5Handle space(H5Dget_space(ds), H5Sclose);
    hsize_t dim = 0;
    H5Sget_simple_extent_dims(space, &dim, nullptr);
    std::vector<double> out(static_cast<size_t>(dim), 0.0);
    EXPECT_GE(H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data()), 0)
        << ds_path;
    return out;
}

/// The loom /matrix dataset (genes x cells, row-major).
struct LoomMatrix {
    hsize_t rows = 0, cols = 0;
    std::vector<float> values;
};

LoomMatrix h5_read_matrix(hid_t fid) {
    LoomMatrix m;
    H5Handle ds(H5Dopen2(fid, "/matrix", H5P_DEFAULT), H5Dclose);
    EXPECT_GE(ds.id, 0);
    if (ds.id < 0) return m;
    H5Handle space(H5Dget_space(ds), H5Sclose);
    hsize_t dims[2] = {0, 0};
    EXPECT_EQ(H5Sget_simple_extent_dims(space, dims, nullptr), 2);
    m.rows = dims[0];
    m.cols = dims[1];
    m.values.assign(static_cast<size_t>(m.rows * m.cols), -1.0f);
    EXPECT_GE(H5Dread(ds, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, m.values.data()), 0);
    return m;
}

baysor::PolygonCollection make_triangle_collection(const std::string& cell_name) {
    baysor::PolygonCollection coll;
    Eigen::MatrixXd tri(2, 3);
    tri << 0.0, 1.0, 0.0,
           0.0, 0.0, 1.0;
    coll[cell_name] = tri;
    return coll;
}

} // namespace

// ============================================================================
// Output style / paths
// ============================================================================

TEST(Cov4OutputStyle, ToStringAndParseRoundTrip) {
    EXPECT_EQ(baysor::to_string(baysor::OutputStyle::Legacy), "legacy");
    EXPECT_EQ(baysor::to_string(baysor::OutputStyle::Parquet), "parquet");
    EXPECT_EQ(baysor::to_string(baysor::parse_output_style("legacy")),
              baysor::to_string(baysor::OutputStyle::Legacy));
    EXPECT_EQ(baysor::to_string(baysor::parse_output_style("parquet")),
              baysor::to_string(baysor::OutputStyle::Parquet));

    EXPECT_THROW_MSG(baysor::parse_output_style("csv"), std::invalid_argument, "Unknown output style: csv");
}

TEST_F(Cov4OutputFiles, GetOutputPathsSupportsLoomCountMatrixAndBackslashBase) {
    auto p = baysor::get_output_paths("out", baysor::OutputStyle::Legacy, "loom");
    EXPECT_EQ(p.counts, "out/segmentation_counts.loom");
    EXPECT_EQ(p.diagnostic_report, "out/diagnostic_report.html");
    EXPECT_EQ(p.molecule_plot, "out/segmentation_plot.html");

    // Trailing backslash is stripped as well as '/'.
    auto q = baysor::get_output_paths("out\\", baysor::OutputStyle::Legacy, "tsv");
    EXPECT_EQ(q.counts, "out/segmentation_counts.tsv");

    auto r = baysor::get_output_paths("out", baysor::OutputStyle::Parquet, "loom");
    EXPECT_EQ(r.counts, "out/feature_matrix.h5");
    EXPECT_EQ(r.diagnostic_report, "out/diagnostic_report.html");
    EXPECT_EQ(r.molecule_plot, "out/segmentation_plot.html");
}

// ============================================================================
// save_segmented_df (CSV)
// ============================================================================

TEST_F(Cov4OutputFiles, SaveSegmentedDf2DWithoutOptionalColumns) {
    baysor::MoleculeData data;
    data.x = {1.5, 2.5, 3.5};
    data.y = {4.5, 5.5, 6.5};
    data.gene = {1, 2, 5};
    data.gene_names = {"A", "B"};

    std::vector<int> assignment = {1, 0, 2};
    const std::string out = path("seg.csv");
    baysor::save_segmented_df(data, assignment, data.gene_names, out);

    const std::string content = read_text_file(out);
    EXPECT_EQ(content,
              "cell,gene,x,y,is_noise\n"
              "cell_1,A,1.5,4.5,0\n"
              "0,B,2.5,5.5,1\n"
              "cell_2,5,3.5,6.5,0\n");
}

TEST_F(Cov4OutputFiles, SaveSegmentedDf3DWithAllOptionalColumns) {
    baysor::MoleculeData data;
    data.x = {1.0, 2.0};
    data.y = {2.0, 3.0};
    data.z = {10.0, 20.0};
    data.gene = {1, 1};
    data.gene_names = {"A"};
    data.confidence = {0.95, 0.4};
    data.source_transcript_id = {101ULL, 102ULL};

    std::vector<int> assignment = {1, 0};
    std::vector<std::string> ncv_color = {"#ff0000", "#00ff00"};
    std::vector<double> assign_conf = {0.9, 0.75};
    std::vector<int> cluster = {3, 0};

    const std::string out = path("seg3d.csv");
    baysor::save_segmented_df(data, assignment, data.gene_names, out,
                              &ncv_color, &assign_conf, &cluster);

    const std::string content = read_text_file(out);
    EXPECT_EQ(content,
              "transcript_id,cell,gene,x,y,z,confidence,cluster,ncv_color,"
              "assignment_confidence,is_noise\n"
              "101,cell_1,A,1,2,10,0.95,3,#ff0000,0.9,false\n"
              "102,0,A,2,3,20,0.4,0,#00ff00,0.75,true\n");
}

// ============================================================================
// save_segmented_df_parquet
// ============================================================================

TEST_F(Cov4OutputFiles, SaveSegmentedDfParquetWritesAllColumnsIn3D) {
    baysor::MoleculeData data;
    data.x = {1.0, 2.0};
    data.y = {2.0, 3.0};
    data.z = {10.0, 20.0};
    data.gene = {1, 2};
    data.gene_names = {"A", "B"};
    data.confidence = {0.95, 0.4};

    std::vector<int> assignment = {1, 0};
    std::vector<std::string> ncv_color = {"#ff0000", "#00ff00"};
    std::vector<double> assign_conf = {0.9, 0.75};
    std::vector<int> cluster = {3, 0};

    const std::string out = path("molecules.parquet");
    baysor::save_segmented_df_parquet(data, assignment, data.gene_names, out,
                                      &ncv_color, &assign_conf, &cluster);

    auto table = read_parquet(out);
    ASSERT_NE(table, nullptr);
    EXPECT_EQ(table->num_rows(), 2);
    EXPECT_EQ(field_names(*table),
              (std::vector<std::string>{"cell", "gene", "x", "y", "z", "confidence",
                                        "cluster", "ncv_color", "assignment_confidence",
                                        "is_noise"}));
    EXPECT_EQ(string_column(*table, "cell"), (std::vector<std::string>{"cell_1", "0"}));
    EXPECT_EQ(string_column(*table, "gene"), (std::vector<std::string>{"A", "B"}));
    EXPECT_EQ(string_column(*table, "ncv_color"),
              (std::vector<std::string>{"#ff0000", "#00ff00"}));
    EXPECT_EQ(double_column(*table, "x"), (std::vector<double>{1.0, 2.0}));
    EXPECT_EQ(double_column(*table, "z"), (std::vector<double>{10.0, 20.0}));
    EXPECT_EQ(double_column(*table, "confidence"), (std::vector<double>{0.95, 0.4}));
    EXPECT_EQ(double_column(*table, "assignment_confidence"),
              (std::vector<double>{0.9, 0.75}));
    EXPECT_EQ(int32_column(*table, "cluster"), (std::vector<int32_t>{3, 0}));
    EXPECT_EQ(bool_column(*table, "is_noise"), (std::vector<bool>{false, true}));
}

TEST_F(Cov4OutputFiles, SaveSegmentedDfParquetMinimal2D) {
    baysor::MoleculeData data;
    data.x = {1.0, 2.0};
    data.y = {3.0, 4.0};
    data.gene = {0, 5};
    data.gene_names = {"A"};

    std::vector<int> assignment = {0, 2};
    const std::string out = path("molecules_min.parquet");
    baysor::save_segmented_df_parquet(data, assignment, data.gene_names, out);

    auto table = read_parquet(out);
    ASSERT_NE(table, nullptr);
    EXPECT_EQ(table->num_rows(), 2);
    // No z/confidence/cluster/ncv_color/assignment_confidence columns.
    EXPECT_EQ(field_names(*table),
              (std::vector<std::string>{"cell", "gene", "x", "y", "is_noise"}));
    // Out-of-range gene id falls back to the numeric id.
    EXPECT_EQ(string_column(*table, "gene"), (std::vector<std::string>{"0", "5"}));
    EXPECT_EQ(string_column(*table, "cell"), (std::vector<std::string>{"0", "cell_2"}));
    EXPECT_EQ(bool_column(*table, "is_noise"), (std::vector<bool>{true, false}));
}

// ============================================================================
// save_cell_stat_df (CSV + Parquet)
// ============================================================================

TEST_F(Cov4OutputFiles, SaveCellStatDfCsvWritesHeaderAndRows) {
    Eigen::MatrixXd stats(2, 3);
    stats << 1.5, 2.0, 3.0,
             4.0, 5.0, 6.5;
    std::vector<std::string> cells = {"cell_1", "cell_2"};
    std::vector<std::string> cols = {"area", "density", "elongation"};

    const std::string out = path("stats.csv");
    baysor::save_cell_stat_df(stats, cells, cols, out);
    EXPECT_EQ(read_text_file(out),
              "cell,area,density,elongation\n"
              "cell_1,1.5,2,3\n"
              "cell_2,4,5,6.5\n");

    // Empty (zero-row) table writes only the header.
    Eigen::MatrixXd empty_stats(0, 2);
    const std::string out_empty = path("stats_empty.csv");
    baysor::save_cell_stat_df(empty_stats, {}, {"area", "density"}, out_empty);
    EXPECT_EQ(read_text_file(out_empty), "cell,area,density\n");
}

TEST_F(Cov4OutputFiles, SaveCellStatDfParquetWritesColumns) {
    Eigen::MatrixXd stats(2, 2);
    stats << 10.0, 0.5,
             20.0, 1.5;
    std::vector<std::string> cells = {"cell_1", "cell_2"};
    std::vector<std::string> cols = {"area", "density"};

    const std::string out = path("cells.parquet");
    baysor::save_cell_stat_df_parquet(stats, cells, cols, out);

    auto table = read_parquet(out);
    ASSERT_NE(table, nullptr);
    EXPECT_EQ(table->num_rows(), 2);
    EXPECT_EQ(field_names(*table),
              (std::vector<std::string>{"cell", "area", "density"}));
    EXPECT_EQ(string_column(*table, "cell"), cells);
    EXPECT_EQ(double_column(*table, "area"), (std::vector<double>{10.0, 20.0}));
    EXPECT_EQ(double_column(*table, "density"), (std::vector<double>{0.5, 1.5}));
}

// ============================================================================
// save_matrix_to_tsv
// ============================================================================

TEST_F(Cov4OutputFiles, SaveMatrixToTsvWritesDenseGeneRows) {
    // 3 cells x 2 genes; cell c2 is empty (no counts at all).
    Eigen::SparseMatrix<double> matrix(3, 2);
    std::vector<Eigen::Triplet<double>> trips = {{0, 0, 2.0}, {2, 1, 5.0}};
    matrix.setFromTriplets(trips.begin(), trips.end());
    matrix.makeCompressed();

    std::vector<std::string> genes = {"G1", "G2"};
    std::vector<std::string> cells = {"c1", "c2", "c3"};

    const std::string out = path("counts.tsv");
    baysor::save_matrix_to_tsv(matrix, genes, cells, out);
    EXPECT_EQ(read_text_file(out),
              "gene\tc1\tc2\tc3\n"
              "G1\t2\t0\t0\n"
              "G2\t0\t0\t5\n");
}

TEST_F(Cov4OutputFiles, SaveMatrixToTsvRejectsMismatchedNames) {
    Eigen::SparseMatrix<double> matrix(2, 1);
    matrix.insert(0, 0) = 1.0;
    matrix.makeCompressed();

    const std::string out = path("counts_bad.tsv");
    EXPECT_THROW_MSG(baysor::save_matrix_to_tsv(matrix, {"G1", "G2"}, {"c1", "c2"}, out),
                     std::runtime_error, "save_matrix_to_tsv: gene_names length mismatch");
    EXPECT_THROW_MSG(baysor::save_matrix_to_tsv(matrix, {"G1"}, {"c1", "c2", "c3"}, out),
                     std::runtime_error, "save_matrix_to_tsv: cell_names length mismatch");
    // Unwritable target (parent directory does not exist).
    EXPECT_THROW_MSG(baysor::save_matrix_to_tsv(matrix, {"G1"}, {"c1"}, path("no_such_dir/counts.tsv")),
                     std::runtime_error, "save_matrix_to_tsv: cannot open");
}

// ============================================================================
// save_polygons_geojson
// ============================================================================

TEST_F(Cov4OutputFiles, SavePolygonsGeoJsonGeometryCollection) {
    auto coll = make_triangle_collection("cell_1");
    const std::string out = path("geom.json");
    baysor::save_polygons_geojson(coll, out, "GeometryCollection");

    auto json = nlohmann::json::parse(read_text_file(out));
    EXPECT_EQ(json["type"], "GeometryCollection");
    ASSERT_TRUE(json.contains("geometries"));
    ASSERT_EQ(json["geometries"].size(), 1u);
    EXPECT_EQ(json["geometries"][0]["type"], "Polygon");
    EXPECT_EQ(json["geometries"][0]["cell"], "cell_1");
    // Ring is closed by appending the first vertex.
    EXPECT_EQ(json["geometries"][0]["coordinates"][0].size(), 4u);
    EXPECT_FALSE(json.contains("features"));
}

TEST_F(Cov4OutputFiles, SavePolygonsGeoJsonSkipsTinyFeaturesAndRingsAreClosed) {
    baysor::PolygonCollection coll;
    Eigen::MatrixXd tiny(2, 2);  // only 2 vertices -> skipped in FeatureCollection
    tiny << 0.0, 1.0,
            0.0, 1.0;
    coll["cell_tiny"] = tiny;
    Eigen::MatrixXd open_tri(2, 3);
    open_tri << 0.0, 1.0, 0.0,
                0.0, 0.0, 1.0;
    coll["cell_open"] = open_tri;
    Eigen::MatrixXd closed(2, 4);  // last vertex equals first -> not closed again
    closed << 0.0, 1.0, 1.0, 0.0,
              0.0, 0.0, 1.0, 0.0;
    coll["cell_closed"] = closed;

    const std::string out = path("features.json");
    baysor::save_polygons_geojson(coll, out, "FeatureCollection");

    auto json = nlohmann::json::parse(read_text_file(out));
    EXPECT_EQ(json["type"], "FeatureCollection");
    ASSERT_EQ(json["features"].size(), 2u);  // tiny feature skipped
    int n_open = 0, n_closed = 0;
    for (const auto& feat : json["features"]) {
        const std::string id = feat["id"];
        const auto& ring = feat["geometry"]["coordinates"][0];
        if (id == "cell_open") {
            n_open++;
            EXPECT_EQ(ring.size(), 4u);  // 3 + appended first vertex
        } else if (id == "cell_closed") {
            n_closed++;
            EXPECT_EQ(ring.size(), 4u);  // already closed, no duplicate appended
        }
    }
    EXPECT_EQ(n_open, 1);
    EXPECT_EQ(n_closed, 1);
}

TEST_F(Cov4OutputFiles, SavePolygonsGeoJsonNoOpFormatsWriteNothing) {
    auto coll = make_triangle_collection("cell_1");
    const std::string none_out = path("none.json");
    baysor::save_polygons_geojson(coll, none_out, "none");
    EXPECT_FALSE(fs::exists(none_out));

    const std::string empty_out = path("empty.json");
    baysor::save_polygons_geojson(baysor::PolygonCollection{}, empty_out);
    EXPECT_FALSE(fs::exists(empty_out));
}

// ============================================================================
// save_polygons_geoparquet
// ============================================================================

TEST_F(Cov4OutputFiles, SavePolygonsGeoParquetHandlesNx2AndSkipsDegenerate) {
    baysor::PolygonCollection coll;
    // Nx2 layout: one (x, y) pair per row.
    Eigen::MatrixXd nx2(4, 2);
    nx2 << 0.0, 0.0,
           1.0, 0.0,
           1.0, 1.0,
           0.0, 1.0;
    coll["cell_nx2"] = nx2;
    Eigen::MatrixXd tri(2, 3);
    tri << 0.0, 1.0, 0.0,
           0.0, 0.0, 1.0;
    coll["cell_tri"] = tri;
    Eigen::MatrixXd bad(3, 3);  // neither 2xN nor Nx2 -> dropped
    bad.setZero();
    coll["cell_bad"] = bad;

    const std::string out = path("boundaries.parquet");
    baysor::save_polygons_geoparquet(coll, out, "geom");

    auto table = read_parquet(out);
    ASSERT_NE(table, nullptr);
    EXPECT_EQ(table->num_rows(), 2);
    EXPECT_EQ(field_names(*table), (std::vector<std::string>{"cell", "n_vertices", "geom"}));

    auto cells = string_column(*table, "cell");
    auto n_vertices = int32_column(*table, "n_vertices");
    std::map<std::string, int32_t> by_cell;
    for (size_t i = 0; i < cells.size(); ++i) by_cell[cells[i]] = n_vertices[i];
    ASSERT_EQ(by_cell.size(), 2u);
    EXPECT_EQ(by_cell.at("cell_nx2"), 4);  // open ring of 4 -> 4 vertices after closing
    EXPECT_EQ(by_cell.at("cell_tri"), 3);

    // WKB geometry: little-endian polygon with one ring.
    int geom_idx = table->schema()->GetFieldIndex("geom");
    auto geom = std::static_pointer_cast<arrow::BinaryArray>(
        table->column(geom_idx)->chunk(0));
    ASSERT_GE(geom->length(), 1);
    for (int64_t i = 0; i < geom->length(); ++i) {
        auto bytes = geom->GetView(i);
        ASSERT_GE(bytes.size(), 9u);
        EXPECT_EQ(static_cast<uint8_t>(bytes[0]), 1u);  // little endian
        uint32_t type = 0;
        std::memcpy(&type, bytes.data() + 1, 4);
        EXPECT_EQ(type, 3u);  // Polygon
        uint32_t n_rings = 0;
        std::memcpy(&n_rings, bytes.data() + 5, 4);
        EXPECT_EQ(n_rings, 1u);
    }

    // GeoParquet metadata names the custom geometry column.
    const std::string geo = parquet_footer_kv(out, "geo");
    EXPECT_NE(geo.find("\"primary_column\":\"geom\""), std::string::npos);
    EXPECT_NE(geo.find("\"encoding\":\"WKB\""), std::string::npos);
}

TEST_F(Cov4OutputFiles, SavePolygonsGeoParquetWritesNothingWithoutValidCells) {
    Eigen::MatrixXd bad(3, 3);
    bad.setZero();
    baysor::PolygonCollection only_bad = {{"cell_bad", bad}};

    const std::string out = path("nothing.parquet");
    baysor::save_polygons_geoparquet(only_bad, out);
    EXPECT_FALSE(fs::exists(out));

    const std::string out_empty = path("empty.parquet");
    baysor::save_polygons_geoparquet(baysor::PolygonCollection{}, out_empty);
    EXPECT_FALSE(fs::exists(out_empty));
}

// ============================================================================
// save_polygon_stack_geojson / save_polygon_stack_geoparquet
// ============================================================================

TEST_F(Cov4OutputFiles, SavePolygonStackGeoJsonWrites2dAnd3dOutputs) {
    auto coll2d = make_triangle_collection("cell_2d");
    auto coll_z = make_triangle_collection("cell_z");

    baysor::OutputPaths paths;
    paths.polygons_2d = path("stack_2d.json");
    paths.polygons_3d = path("stack_3d.json");
    baysor::PolygonStack stack = {{"2d", coll2d}, {"z0", coll_z}};
    baysor::save_polygon_stack_geojson(stack, paths);

    auto j2d = nlohmann::json::parse(read_text_file(paths.polygons_2d));
    EXPECT_EQ(j2d["type"], "FeatureCollection");
    EXPECT_EQ(j2d["features"][0]["id"], "cell_2d");

    auto j3d = nlohmann::json::parse(read_text_file(paths.polygons_3d));
    ASSERT_TRUE(j3d.contains("z0"));
    EXPECT_EQ(j3d["z0"]["type"], "FeatureCollection");
    EXPECT_EQ(j3d["z0"]["features"][0]["id"], "cell_z");
}

TEST_F(Cov4OutputFiles, SavePolygonStackGeoJsonEdgeCases) {
    auto coll2d = make_triangle_collection("cell_2d");

    // 2D-only stack: 3D file must not be created.
    baysor::OutputPaths only2d;
    only2d.polygons_2d = path("only2d.json");
    only2d.polygons_3d = path("only2d_3d.json");
    baysor::save_polygon_stack_geojson({{"2d", coll2d}}, only2d);
    EXPECT_TRUE(fs::exists(only2d.polygons_2d));
    EXPECT_FALSE(fs::exists(only2d.polygons_3d));

    // Empty stack and format "none" write nothing at all.
    baysor::OutputPaths noop;
    noop.polygons_2d = path("noop_2d.json");
    noop.polygons_3d = path("noop_3d.json");
    baysor::save_polygon_stack_geojson({}, noop);
    EXPECT_FALSE(fs::exists(noop.polygons_2d));
    EXPECT_FALSE(fs::exists(noop.polygons_3d));
    baysor::save_polygon_stack_geojson({{"2d", coll2d}}, noop, "none");
    EXPECT_FALSE(fs::exists(noop.polygons_2d));
    EXPECT_FALSE(fs::exists(noop.polygons_3d));
}

TEST_F(Cov4OutputFiles, SavePolygonStackGeoParquetSplitsLayers) {
    auto coll2d = make_triangle_collection("cell_2d");
    auto coll_z = make_triangle_collection("cell_z");
    Eigen::MatrixXd bad(3, 3);
    bad.setZero();
    coll_z["cell_bad"] = bad;

    baysor::OutputPaths paths;
    paths.polygons_2d = path("stack_2d.parquet");
    paths.polygons_3d = path("stack_3d.parquet");
    baysor::PolygonStack stack = {{"2d", coll2d}, {"z0", coll_z}};
    baysor::save_polygon_stack_geoparquet(stack, paths, "geom");

    auto t2d = read_parquet(paths.polygons_2d);
    ASSERT_NE(t2d, nullptr);
    EXPECT_EQ(t2d->num_rows(), 1);
    EXPECT_EQ(string_column(*t2d, "cell"), (std::vector<std::string>{"cell_2d"}));
    EXPECT_NE(parquet_footer_kv(paths.polygons_2d, "geo")
                  .find("\"primary_column\":\"geom\""),
              std::string::npos);

    auto t3d = read_parquet(paths.polygons_3d);
    ASSERT_NE(t3d, nullptr);
    EXPECT_EQ(t3d->num_rows(), 1);  // degenerate cell dropped
    EXPECT_EQ(field_names(*t3d),
              (std::vector<std::string>{"cell", "layer", "n_vertices", "geom"}));
    EXPECT_EQ(string_column(*t3d, "cell"), (std::vector<std::string>{"cell_z"}));
    EXPECT_EQ(string_column(*t3d, "layer"), (std::vector<std::string>{"z0"}));
}

TEST_F(Cov4OutputFiles, SavePolygonStackGeoParquetEdgeCases) {
    auto coll2d = make_triangle_collection("cell_2d");

    // 2D-only stack: no 3D parquet output.
    baysor::OutputPaths only2d;
    only2d.polygons_2d = path("p_only2d.parquet");
    only2d.polygons_3d = path("p_only2d_3d.parquet");
    baysor::save_polygon_stack_geoparquet({{"2d", coll2d}}, only2d);
    EXPECT_TRUE(fs::exists(only2d.polygons_2d));
    EXPECT_FALSE(fs::exists(only2d.polygons_3d));

    // Empty stack writes nothing.
    baysor::OutputPaths noop;
    noop.polygons_2d = path("p_noop_2d.parquet");
    noop.polygons_3d = path("p_noop_3d.parquet");
    baysor::save_polygon_stack_geoparquet({}, noop);
    EXPECT_FALSE(fs::exists(noop.polygons_2d));
    EXPECT_FALSE(fs::exists(noop.polygons_3d));
}

// ============================================================================
// save_matrix_to_loom (HDF5)
// ============================================================================

TEST_F(Cov4OutputFiles, SaveLoomWithColAttrsWritesStringAndDoubleDatasets) {
    // Row-major cells x genes; cell 2 is empty (all-zero row).
    Eigen::SparseMatrix<float, Eigen::RowMajor> matrix(2, 3);
    std::vector<Eigen::Triplet<float>> trips = {{0, 0, 2.0f}, {0, 2, 5.0f}};
    matrix.setFromTriplets(trips.begin(), trips.end());
    matrix.makeCompressed();

    std::vector<std::string> genes = {"G1", "G2", "G3"};
    std::vector<std::string> cells = {"cell_1", "cell_2"};
    baysor::LoomColAttrs col_attrs;
    col_attrs["ncv_color"] = std::vector<std::string>{"#aabbcc", "#ddeeff"};
    col_attrs["confidence"] = std::vector<double>{0.9, 0.1};

    const std::string out = path("counts.loom");
    baysor::save_matrix_to_loom(matrix, genes, cells, out, col_attrs);

    H5Handle fid(H5Fopen(out.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    ASSERT_GE(fid.id, 0);
    EXPECT_EQ(h5_read_vlen_strings(fid, "/attrs/LOOM_SPEC_VERSION"),
              std::vector<std::string>{"3.0.0"});
    EXPECT_EQ(h5_read_vlen_strings(fid, "/row_attrs/Name"), genes);
    EXPECT_EQ(h5_read_vlen_strings(fid, "/col_attrs/Name"), cells);
    EXPECT_EQ(h5_read_vlen_strings(fid, "/col_attrs/ncv_color"),
              (std::vector<std::string>{"#aabbcc", "#ddeeff"}));
    EXPECT_EQ(h5_read_doubles(fid, "/col_attrs/confidence"),
              (std::vector<double>{0.9, 0.1}));
    EXPECT_EQ(h5_read_doubles(fid, "/col_attrs/CellID"),
              (std::vector<double>{1.0, 2.0}));

    // Loom matrix is stored genes x cells; the empty cell is an all-zero column.
    const auto m = h5_read_matrix(fid);
    EXPECT_EQ(m.rows, 3u);
    EXPECT_EQ(m.cols, 2u);
    EXPECT_EQ(m.values, (std::vector<float>{2.0f, 0.0f,
                                            0.0f, 0.0f,
                                            5.0f, 0.0f}));
}

TEST_F(Cov4OutputFiles, SaveLoomColMajorOverloadWithColAttrs) {
    Eigen::SparseMatrix<float> matrix(2, 3);  // column-major cells x genes
    std::vector<Eigen::Triplet<float>> trips = {{1, 1, 7.0f}};
    matrix.setFromTriplets(trips.begin(), trips.end());
    matrix.makeCompressed();

    baysor::LoomColAttrs col_attrs;
    col_attrs["ncv_color"] = std::vector<std::string>{"x", "y"};

    const std::string out = path("counts_cm.loom");
    baysor::save_matrix_to_loom(matrix, {"G1", "G2", "G3"}, {"c1", "c2"}, out, col_attrs);

    H5Handle fid(H5Fopen(out.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose);
    ASSERT_GE(fid.id, 0);
    EXPECT_EQ(h5_read_vlen_strings(fid, "/col_attrs/ncv_color"),
              (std::vector<std::string>{"x", "y"}));
    // Gene G2 in cell c2 only.
    const auto m = h5_read_matrix(fid);
    EXPECT_EQ(m.rows, 3u);
    EXPECT_EQ(m.cols, 2u);
    EXPECT_EQ(m.values, (std::vector<float>{0.0f, 0.0f,
                                            0.0f, 7.0f,
                                            0.0f, 0.0f}));
}

TEST_F(Cov4OutputFiles, SaveLoomRejectsMismatchedNames) {
    Eigen::SparseMatrix<float> matrix(2, 3);
    matrix.insert(0, 0) = 1.0f;
    matrix.makeCompressed();

    const std::string out = path("bad.loom");
    EXPECT_THROW_MSG(baysor::save_matrix_to_loom(matrix, {"G1", "G2"}, {"c1", "c2"}, out),
                     std::runtime_error, "gene_names length mismatch");
    EXPECT_THROW_MSG(baysor::save_matrix_to_loom(matrix, {"G1", "G2", "G3"}, {"c1", "c2", "c3"}, out),
                     std::runtime_error, "cell_names length mismatch");
}

// ============================================================================
// Parquet writer error paths
// ============================================================================

TEST_F(Cov4OutputFiles, ParquetOpenFailureThrowsRuntimeError) {
    Eigen::MatrixXd stats(1, 1);
    stats << 1.0;
    EXPECT_THROW_MSG(baysor::save_cell_stat_df_parquet(stats, {"c1"}, {"area"}, path("no_such_dir/cells.parquet")),
                     std::runtime_error, "Open parquet output");
}

#ifndef _WIN32
TEST_F(Cov4OutputFiles, ParquetOpenWriteFailureThrowsRuntimeError) {
    // /dev/full accepts open() but fails the very first write (the parquet
    // magic bytes emitted by FileWriter::Open) with ENOSPC, surfacing the
    // arrow_unwrap() error path.
    if (!fs::exists("/dev/full")) {
        GTEST_SKIP() << "/dev/full is not available";
    }
    Eigen::MatrixXd stats(1, 1);
    stats << 1.0;
    EXPECT_THROW_MSG(baysor::save_cell_stat_df_parquet(stats, {"c1"}, {"area"}, "/dev/full"),
                     std::runtime_error, "Open parquet writer: IOError");
}
#endif  // !_WIN32

#ifndef _WIN32
TEST_F(Cov4OutputFiles, ParquetMidStreamWriteFailureThrowsRuntimeError) {
    // Cap the file size so the header write succeeds but a later write fails:
    // this exercises the arrow_check() error path (a failure during
    // FileWriter::Open itself would take the arrow_unwrap path instead).
    ScopedFileSizeLimit limit(4096);
    if (!limit.ok()) {
        GTEST_SKIP() << "cannot lower RLIMIT_FSIZE in this process";
    }

    // Random doubles do not compress away, so the parquet file is far larger
    // than the 4 KiB cap.
    Eigen::MatrixXd stats(5000, 1);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int r = 0; r < 5000; ++r) stats(r, 0) = dist(rng);
    std::vector<std::string> cells(5000, "cell");

    const std::string out = path("too_big.parquet");
    try {
        baysor::save_cell_stat_df_parquet(stats, cells, {"area"}, out);
        FAIL() << "expected std::runtime_error once the file-size limit is hit";
    } catch (const std::runtime_error& e) {
        const std::string msg = e.what();
        // Failure happened after the writer was opened -> arrow_check() path.
        EXPECT_EQ(msg.find("Open parquet"), std::string::npos) << msg;
        EXPECT_NE(msg.find("IOError"), std::string::npos) << msg;
        EXPECT_NE(msg.find("parquet"), std::string::npos) << msg;
    }
}
#endif  // !_WIN32
