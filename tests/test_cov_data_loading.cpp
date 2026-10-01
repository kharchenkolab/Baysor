// Tests for src/data_loading/{data,prior_segmentation}.cpp: column readers,
// CSV/Parquet loaders and the prior-segmentation loaders.

#include <gtest/gtest.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#include <tiffio.h>

#include "baysor/data_loading/data.h"
#include "baysor/utils/general.h"
#include "baysor/data_loading/prior_segmentation.h"
#include "baysor/utils/options.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "test_cov_helpers.h"

namespace {

using baysor_test::TempDir;

// ---------------------------------------------------------------------------
// Arrow array builders
// ---------------------------------------------------------------------------

template <class T> struct is_optional : std::false_type {};
template <class T> struct is_optional<std::optional<T>> : std::true_type {};

template <class Builder>
std::shared_ptr<arrow::Array> finish(Builder& b) {
    std::shared_ptr<arrow::Array> out;
    EXPECT_TRUE(b.Finish(&out).ok());
    return out;
}

// Array of `values`; std::nullopt entries become nulls.
template <class Builder, class T>
std::shared_ptr<arrow::Array> build(const std::vector<T>& values) {
    Builder b;
    for (const auto& v : values) {
        if constexpr (is_optional<T>::value) EXPECT_TRUE((v ? b.Append(*v) : b.AppendNull()).ok());
        else EXPECT_TRUE(b.Append(v).ok());
    }
    return finish(b);
}

auto arr_f64(const std::vector<double>& v) { return build<arrow::DoubleBuilder>(v); }
auto arr_f64_nulls(const std::vector<std::optional<double>>& v) { return build<arrow::DoubleBuilder>(v); }
auto arr_f32(const std::vector<float>& v) { return build<arrow::FloatBuilder>(v); }
auto arr_i64(const std::vector<int64_t>& v) { return build<arrow::Int64Builder>(v); }
auto arr_i64_nulls(const std::vector<std::optional<int64_t>>& v) { return build<arrow::Int64Builder>(v); }
auto arr_i32(const std::vector<int32_t>& v) { return build<arrow::Int32Builder>(v); }
auto arr_i16(const std::vector<int16_t>& v) { return build<arrow::Int16Builder>(v); }
auto arr_i8(const std::vector<int8_t>& v) { return build<arrow::Int8Builder>(v); }
auto arr_u64(const std::vector<uint64_t>& v) { return build<arrow::UInt64Builder>(v); }
auto arr_u32(const std::vector<uint32_t>& v) { return build<arrow::UInt32Builder>(v); }
auto arr_u16(const std::vector<uint16_t>& v) { return build<arrow::UInt16Builder>(v); }
auto arr_u8(const std::vector<uint8_t>& v) { return build<arrow::UInt8Builder>(v); }
auto arr_str(const std::vector<std::string>& v) { return build<arrow::StringBuilder>(v); }
auto arr_lstr(const std::vector<std::string>& v) { return build<arrow::LargeStringBuilder>(v); }
auto arr_bin(const std::vector<std::string>& v) { return build<arrow::BinaryBuilder>(v); }

// Dictionary-encoded array: `indices` (-1 = null) into `values`.
template <class IndexBuilder = arrow::Int8Builder>
std::shared_ptr<arrow::Array> arr_dict(const std::shared_ptr<arrow::Array>& values,
                                       const std::vector<int>& indices) {
    IndexBuilder b;
    for (int i : indices) EXPECT_TRUE((i < 0 ? b.AppendNull() : b.Append(i)).ok());
    auto result = arrow::DictionaryArray::FromArrays(finish(b), values);
    EXPECT_TRUE(result.ok()) << result.status().ToString();
    return result.ValueOrDie();
}

std::shared_ptr<arrow::Array> arr_list_i32(const std::vector<std::vector<int32_t>>& rows) {
    auto vb = std::make_shared<arrow::Int32Builder>();
    arrow::ListBuilder lb(arrow::default_memory_pool(), vb);
    for (const auto& row : rows) {
        EXPECT_TRUE(lb.Append(true).ok());
        for (int32_t v : row) EXPECT_TRUE(vb->Append(v).ok());
    }
    return finish(lb);
}

// ---------------------------------------------------------------------------
// Parquet writer (optionally storing the Arrow schema so that dictionary and
// large_string columns survive the roundtrip)
// ---------------------------------------------------------------------------

using Columns = std::vector<std::pair<std::string, std::shared_ptr<arrow::Array>>>;

std::string write_parquet(const TempDir& dir, const std::string& name, const Columns& columns,
                          bool store_schema = false) {
    const std::string path = dir.file(name);
    std::vector<std::shared_ptr<arrow::Field>> fields;
    std::vector<std::shared_ptr<arrow::Array>> arrays;
    for (const auto& [col, arr] : columns) {
        fields.push_back(arrow::field(col, arr->type()));
        arrays.push_back(arr);
    }
    auto table = arrow::Table::Make(arrow::schema(fields), arrays);
    auto sink = arrow::io::FileOutputStream::Open(path).ValueOrDie();
    auto arrow_props = store_schema
        ? parquet::ArrowWriterProperties::Builder().store_schema()->build()
        : parquet::default_arrow_writer_properties();
    auto status = parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), sink,
                                             /*row_group_size=*/3,
                                             parquet::default_writer_properties(), arrow_props);
    EXPECT_TRUE(status.ok()) << status.ToString();
    EXPECT_TRUE(sink->Close().ok());
    return path;
}


// ---------------------------------------------------------------------------
// TIFF writers
// ---------------------------------------------------------------------------

std::string write_tiff(const TempDir& dir, const std::string& name,
                       const uint8_t* row_bytes, size_t bytes_per_row,
                       uint32_t w, uint32_t h, uint16_t bps, uint16_t spp = 1,
                       bool deflate = false) {
    const std::string path = dir.file(name);
    TIFF* tif = TIFFOpen(path.c_str(), "w");
    EXPECT_NE(tif, nullptr);
    if (!tif) return path;
    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, spp);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bps);
    TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, spp == 3 ? PHOTOMETRIC_RGB : PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, h);
    if (deflate) TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_DEFLATE);
    for (uint32_t row = 0; row < h; ++row) {
        auto* row_ptr = const_cast<uint8_t*>(row_bytes + static_cast<size_t>(row) * bytes_per_row);
        EXPECT_EQ(TIFFWriteScanline(tif, row_ptr, row, 0), 1);
    }
    TIFFClose(tif);
    return path;
}

// Overwrite the compressed image data (everything between the 8-byte header
// and the IFD, which libtiff appends after the data) with 0xFF so that the
// directory still parses but every scanline decode fails.
void corrupt_tiff_pixel_data(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    ASSERT_TRUE(in.good());
    std::vector<char> bytes((std::istreambuf_iterator<char>(in)),
                            std::istreambuf_iterator<char>());
    in.close();
    ASSERT_GE(bytes.size(), 16u);
    std::uint32_t ifd_offset = 0;
    std::memcpy(&ifd_offset, bytes.data() + 4, 4); // little-endian TIFF header
    ASSERT_GT(ifd_offset, 8u);
    ASSERT_LT(ifd_offset, bytes.size());
    for (std::size_t i = 8; i < ifd_offset; ++i) bytes[i] = static_cast<char>(0xFF);
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    ASSERT_TRUE(out.good());
    out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
}

template <class T>
std::string write_tiff_mask(const TempDir& dir, const std::string& name,
                            const std::vector<T>& pixels, uint32_t w, uint32_t h) {
    return write_tiff(dir, name, reinterpret_cast<const uint8_t*>(pixels.data()), w * sizeof(T), w, h,
                      8 * sizeof(T));
}

} // namespace

// ============================================================================
// src/data_loading/data.cpp — column readers
// ============================================================================

TEST(Cov1Data_Readers, ArrowErrorOnMissingFile) {
    TempDir dir("cov1_data");
    EXPECT_THROW_MSG(baysor::read_double_column(dir.file("does_not_exist.csv"), "x"),
                     std::runtime_error, "Arrow error:");
    EXPECT_THROW_MSG(baysor::read_string_column(dir.file("does_not_exist.parquet"), "gene"),
                     std::runtime_error, "Arrow error:");
}

TEST(Cov1Data_Readers, UnsupportedFileFormat) {
    TempDir dir("cov1_data");
    const auto path = dir.write("data.txt", "x,y\n1,2\n");
    EXPECT_THROW_MSG(baysor::read_double_column(path, "x"),
                     std::runtime_error, "Unsupported file format: .txt");
}

TEST(Cov1Data_Readers, ParquetAndPqExtensions) {
    TempDir dir("cov1_data");
    auto x = arr_f64({1.5, 2.5, 3.5});
    auto y = arr_f64({4.0, 5.0, 6.0});
    auto p1 = write_parquet(dir, "a.parquet", {{"x", x}, {"y", y}});
    auto p2 = write_parquet(dir, "b.pq", {{"x", x}, {"y", y}});

    for (const auto& p : {p1, p2}) {
        auto vals = baysor::read_double_column(p, "x");
        ASSERT_EQ(vals.size(), 3u);
        EXPECT_DOUBLE_EQ(vals[0], 1.5);
        EXPECT_DOUBLE_EQ(vals[2], 3.5);
    }
}

TEST(Cov1Data_Readers, MissingColumnMessage) {
    TempDir dir("cov1_data");
    const auto path = dir.write("d.csv", "x,y\n1,2\n");
    EXPECT_THROW_MSG(baysor::read_double_column(path, "gene"), std::runtime_error,
                     "Column 'gene' not found in the data. Available columns");
}

TEST(Cov1Data_Readers, DoubleColumnNumericTypes) {
    TempDir dir("cov1_data");
    auto xf = arr_f32({1.5f, 2.5f});
    auto x32 = arr_i32({7, -3});
    auto x16 = arr_i16({5, 6});
    auto x8 = arr_u8({9, 10});
    auto x64 = arr_i64({11, 12});
    auto path = write_parquet(dir, "types.parquet",
                              {{"xf", xf}, {"x32", x32}, {"x16", x16}, {"x8", x8}, {"x64", x64}});

    const std::pair<const char*, std::vector<double>> expected[] = {
        {"xf", {1.5, 2.5}}, {"x32", {7.0, -3.0}}, {"x16", {5.0, 6.0}}, {"x8", {9.0, 10.0}}, {"x64", {11.0, 12.0}},
    };
    for (const auto& [col, values] : expected) EXPECT_EQ(baysor::read_double_column(path, col), values) << col;
}

TEST(Cov1Data_Readers, DoubleColumnCastFailureThrows) {
    TempDir dir("cov1_data");
    auto lists = arr_list_i32({{1, 2}, {3}});
    auto xs = arr_f64({1.0, 2.0});
    auto path = write_parquet(dir, "lists.parquet", {{"bad", lists}, {"x", xs}});
    EXPECT_THROW_MSG(baysor::read_double_column(path, "bad"),
                     std::runtime_error, "Cannot convert column 'bad' to double");
}

TEST(Cov1Data_Readers, StringColumnTypes) {
    TempDir dir("cov1_data");
    auto plain = arr_str({"a", "b"});
    auto large = arr_lstr({"big1", "big2"});
    auto d8 = arr_dict(arr_str({"G1", "G2"}), {0, 1});
    auto d16 = arr_dict<arrow::Int16Builder>(arr_str({"G1", "G2"}), {1, 0});
    auto d32 = arr_dict<arrow::Int32Builder>(arr_str({"G1", "G2"}), {0, 0});
    auto d64 = arr_dict<arrow::Int64Builder>(arr_str({"G1", "G2"}), {1, 1});
    auto num = arr_i64({101, 202});

    auto path = write_parquet(dir, "strings.parquet",
                              {{"plain", plain}, {"large", large}, {"d8", d8}, {"d16", d16},
                               {"d32", d32}, {"d64", d64}, {"num", num}},
                              /*store_schema=*/true);

    EXPECT_EQ(baysor::read_string_column(path, "plain"), (std::vector<std::string>{"a", "b"}));
    EXPECT_EQ(baysor::read_string_column(path, "large"), (std::vector<std::string>{"big1", "big2"}));
    EXPECT_EQ(baysor::read_string_column(path, "d8"), (std::vector<std::string>{"G1", "G2"}));
    EXPECT_EQ(baysor::read_string_column(path, "d16"), (std::vector<std::string>{"G2", "G1"}));
    EXPECT_EQ(baysor::read_string_column(path, "d32"), (std::vector<std::string>{"G1", "G1"}));
    EXPECT_EQ(baysor::read_string_column(path, "d64"), (std::vector<std::string>{"G2", "G2"}));
    EXPECT_EQ(baysor::read_string_column(path, "num"), (std::vector<std::string>{"101", "202"}));
}

// ============================================================================
// data.cpp — read_tabular_file
// ============================================================================

TEST(Cov1Data_ReadTabular, KeepsVaryingZ) {
    TempDir dir("cov1_data");
    const auto path = dir.write("z3d.csv",
        "x,y,z,gene\n"
        "1,2,10,A\n"
        "3,4,20,B\n"
        "5,6,30,A\n");
    auto raw = baysor::read_tabular_file(path, baysor::MoleculeInputOptions{});
    EXPECT_TRUE(raw.has_z);
    ASSERT_EQ(raw.z.size(), 3u);
    EXPECT_DOUBLE_EQ(raw.z[2], 30.0);
    EXPECT_EQ(raw.gene_str, (std::vector<std::string>{"A", "B", "A"}));
    ASSERT_EQ(raw.x.size(), 3u);
    EXPECT_DOUBLE_EQ(raw.x[1], 3.0);
}

TEST(Cov1Data_ReadTabular, DropsConstantZ) {
    TempDir dir("cov1_data");
    const auto path = dir.write("zconst.csv",
        "x,y,z,gene\n"
        "1,2,7,A\n"
        "3,4,7,B\n");
    auto raw = baysor::read_tabular_file(path, baysor::MoleculeInputOptions{});
    EXPECT_FALSE(raw.has_z);
    EXPECT_TRUE(raw.z.empty());
}

TEST(Cov1Data_ReadTabular, Force2DSkipsZColumn) {
    TempDir dir("cov1_data");
    const auto path = dir.write("zforce.csv",
        "x,y,z,gene\n"
        "1,2,10,A\n"
        "3,4,20,B\n");
    baysor::MoleculeInputOptions opts;
    opts.force_2d = true;
    auto raw = baysor::read_tabular_file(path, opts);
    EXPECT_FALSE(raw.has_z);
    EXPECT_TRUE(raw.z.empty());
}

TEST(Cov1Data_ReadTabular, NoZColumn) {
    TempDir dir("cov1_data");
    const auto path = dir.write("plain.csv",
        "x,y,gene\n"
        "1,2,A\n");
    auto raw = baysor::read_tabular_file(path, baysor::MoleculeInputOptions{});
    EXPECT_FALSE(raw.has_z);
    ASSERT_EQ(raw.gene_str.size(), 1u);
}

// ============================================================================
// data.cpp — load_molecules (CSV)
// ============================================================================

TEST(Cov1Data_LoadCsv, OptionalMetadataColumnsAndQvFilter) {
    TempDir dir("cov1_data");
    const auto path = dir.write("meta.csv",
        "x,y,gene,confidence,cluster,nuclei_probs,qv,transcript_id\n"
        "1,1,G1,0.9,1,0.5,30,1000\n"
        "2,2,G2,0.8,2,0.6,40,1001\n"
        "3,3,G1,0.7,1,0.4,10,1002\n");

    baysor::MoleculeInputOptions opts;
    opts.min_qv = 20.0;
    auto data = baysor::load_molecules(path, opts);

    // Third row has qv 10 < 20 and is filtered out before encoding.
    ASSERT_EQ(data.n_molecules(), 2);
    ASSERT_EQ(data.confidence.size(), 2u);
    EXPECT_DOUBLE_EQ(data.confidence[1], 0.8);
    ASSERT_EQ(data.cluster.size(), 2u);
    EXPECT_EQ(data.cluster[0], 1);
    EXPECT_EQ(data.cluster[1], 2);
    ASSERT_EQ(data.nuclei_probs.size(), 2u);
    EXPECT_DOUBLE_EQ(data.nuclei_probs[0], 0.5);
    ASSERT_EQ(data.source_transcript_id.size(), 2u);
    EXPECT_EQ(data.source_transcript_id[1], 1001u);
    EXPECT_EQ(data.gene_names.size(), 2u);
}

TEST(Cov1Data_LoadCsv, SpatialBoundsFilter) {
    TempDir dir("cov1_data");
    const auto path = dir.write("bounds.csv",
        "x,y,gene\n"
        "1,1,A\n"
        "50,1,A\n"
        "3,90,A\n"
        "2,2,B\n");
    baysor::MoleculeInputOptions opts;
    opts.x_max = 10.0;
    opts.y_max = 10.0;
    auto data = baysor::load_molecules(path, opts);
    ASSERT_EQ(data.n_molecules(), 2);
    EXPECT_DOUBLE_EQ(data.x[1], 2.0);
    EXPECT_DOUBLE_EQ(data.y[1], 2.0);
}

// ============================================================================
// data.cpp — load_molecules (Parquet)
// ============================================================================

TEST(Cov1Data_LoadParquet, AllOptionalColumns) {
    TempDir dir("cov1_data");
    auto x = arr_f64({1, 2, 3, 4});
    auto y = arr_f64({1, 1, 2, 2});
    auto z = arr_f64({10, 20, 30, 40});
    auto gene = arr_str({"A", "A", "B", "B"});
    auto qv = arr_f32({30.0f, 40.0f, 50.0f, 10.0f});
    auto conf = arr_i16({9, 8, 7, 6});
    auto clus = arr_i8({1, 1, 2, 2});
    auto nuclei = arr_u32({5, 5, 6, 6});
    auto tx = arr_i64({100, 101, 102, 103});
    auto cell = arr_str({"c1", "c1", "c2", "0"});

    auto path = write_parquet(dir, "rich.parquet",
                              {{"x", x}, {"y", y}, {"z", z}, {"gene", gene}, {"qv", qv}, {"confidence", conf},
                               {"cluster", clus}, {"nuclei_probs", nuclei}, {"transcript_id", tx}, {"cell_id", cell}});

    baysor::MoleculeInputOptions opts;
    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Column;
    prior.column_name = "cell_id";
    prior.unassigned_label = "0";
    prior.min_molecules_per_segment = 1;

    auto data = baysor::load_molecules(path, opts, prior);

    ASSERT_EQ(data.n_molecules(), 4);
    EXPECT_TRUE(data.is_3d());
    EXPECT_DOUBLE_EQ(data.z[3], 40.0);
    ASSERT_EQ(data.confidence.size(), 4u);
    EXPECT_DOUBLE_EQ(data.confidence[0], 9.0);
    ASSERT_EQ(data.cluster.size(), 4u);
    EXPECT_EQ(data.cluster[3], 2);
    ASSERT_EQ(data.nuclei_probs.size(), 4u);
    EXPECT_DOUBLE_EQ(data.nuclei_probs[2], 6.0);
    ASSERT_EQ(data.source_transcript_id.size(), 4u);
    EXPECT_EQ(data.source_transcript_id[2], 102u);
    EXPECT_EQ(data.prior_segmentation, (std::vector<int>{1, 1, 2, 0}));
}

TEST(Cov1Data_LoadParquet, CoordinateNumericTypes) {
    // Every numeric kind supported by the record-batch scan view is exercised
    // by typing the x column differently across loads.
    struct Case {
        const char* name;
        std::shared_ptr<arrow::Array> arr;
        double expect0;
    };
    TempDir dir("cov1_data");
    std::vector<Case> cases = {
        {"xf32", arr_f32({1.5f, 2.5f}), 1.5},
        {"xi64", arr_i64({7, 8}), 7.0},
        {"xi32", arr_i32({-3, 4}), -3.0},
        {"xi16", arr_i16({5, 6}), 5.0},
        {"xi8", arr_i8({-2, 3}), -2.0},
        {"xu64", arr_u64({11, 12}), 11.0},
        {"xu32", arr_u32({13, 14}), 13.0},
        {"xu16", arr_u16({15, 16}), 15.0},
        {"xu8", arr_u8({17, 18}), 17.0},
        {"xstr", arr_str({"1.25", "2.5"}), 1.25},
    };

    for (const auto& c : cases) {
        auto y = arr_f64({1, 2});
        auto gene = arr_str({"A", "A"});
        auto path = write_parquet(dir, std::string(c.name) + ".parquet",
                                  {{c.name, c.arr}, {"y", y}, {"gene", gene}});
        baysor::MoleculeInputOptions opts;
        opts.x_col = c.name;
        auto data = baysor::load_molecules(path, opts);
        ASSERT_EQ(data.n_molecules(), 2) << c.name;
        EXPECT_DOUBLE_EQ(data.x[0], c.expect0) << c.name;
    }

    // A null coordinate is NaN and fails the static row filter.
    auto xnull = arr_f64_nulls({1.0, std::nullopt});
    auto y = arr_f64({1, 2});
    auto gene = arr_str({"A", "B"});
    auto path = write_parquet(dir, "xnull.parquet", {{"x", xnull}, {"y", y}, {"gene", gene}});
    auto data = baysor::load_molecules(path, baysor::MoleculeInputOptions{});
    ASSERT_EQ(data.n_molecules(), 1);
    EXPECT_DOUBLE_EQ(data.x[0], 1.0);
}

TEST(Cov1Data_LoadParquet, TranscriptIdTypes) {
    struct Case {
        const char* name;
        std::shared_ptr<arrow::Array> arr;
        std::vector<std::uint64_t> expect;
    };
    TempDir dir("cov1_data");
    std::vector<Case> cases = {
        {"ti64", arr_i64({100, 101}), {100u, 101u}},
        {"ti32", arr_i32({200, 201}), {200u, 201u}},
        {"ti16", arr_i16({300, 301}), {300u, 301u}},
        {"ti8", arr_i8({100, 101}), {100u, 101u}},
        {"tu64", arr_u64({400, 401}), {400u, 401u}},
        {"tu32", arr_u32({500, 501}), {500u, 501u}},
        {"tu16", arr_u16({600, 601}), {600u, 601u}},
        {"tu8", arr_u8({70, 71}), {70u, 71u}},
        {"tdbl", arr_f64({800.0, 801.0}), {800u, 801u}},
        {"tflt", arr_f32({900.0f, 901.0f}), {900u, 901u}},
        {"tstr", arr_str({"1000", "1001"}), {1000u, 1001u}},
    };

    for (const auto& c : cases) {
        auto x = arr_f64({1, 2});
        auto y = arr_f64({1, 2});
        auto gene = arr_str({"A", "A"});
        auto path = write_parquet(dir, std::string(c.name) + ".parquet",
                                  {{"x", x}, {"y", y}, {"gene", gene}, {"transcript_id", c.arr}});
        auto data = baysor::load_molecules(path, baysor::MoleculeInputOptions{});
        ASSERT_EQ(data.n_molecules(), 2) << c.name;
        ASSERT_EQ(data.source_transcript_id.size(), 2u) << c.name;
        EXPECT_EQ(data.source_transcript_id, c.expect) << c.name;
    }

    // Null transcript ids become (uint64_t)-1.
    auto tnull = arr_i64_nulls({5, std::nullopt});
    auto x = arr_f64({1, 2});
    auto y = arr_f64({1, 2});
    auto gene = arr_str({"A", "A"});
    auto path = write_parquet(dir, "tnull.parquet",
                              {{"x", x}, {"y", y}, {"gene", gene}, {"transcript_id", tnull}});
    auto data = baysor::load_molecules(path, baysor::MoleculeInputOptions{});
    ASSERT_EQ(data.source_transcript_id.size(), 2u);
    EXPECT_EQ(data.source_transcript_id[0], 5u);
    EXPECT_EQ(data.source_transcript_id[1], std::uint64_t(-1));
}

TEST(Cov1Data_LoadParquet, DictionaryEncodedGene) {
    TempDir dir("cov1_data");
    auto x = arr_f64({1, 2, 3, 4, 5});
    auto y = arr_f64({1, 1, 1, 1, 1});
    // gene: [A, B, null, A, B]; cell_id: [c1, c1, null, c2, c2]
    auto gene = arr_dict(arr_str({"A", "B"}), {0, 1, -1, 0, 1});
    auto cell = arr_dict(arr_str({"c1", "c2"}), {0, 0, -1, 1, 1});

    auto path = write_parquet(dir, "dict.parquet",
                              {{"x", x}, {"y", y}, {"gene", gene}, {"cell_id", cell}},
                              /*store_schema=*/true);

    baysor::MoleculeInputOptions opts;
    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Column;
    prior.column_name = "cell_id";
    prior.unassigned_label = "";
    prior.min_molecules_per_segment = 1;

    auto data = baysor::load_molecules(path, opts, prior);

    // The null-gene row is skipped (gene_id stays 0).
    ASSERT_EQ(data.n_molecules(), 4);
    EXPECT_EQ(data.gene_names, (std::vector<std::string>{"A", "B"}));
    EXPECT_EQ(data.gene, (std::vector<int>{1, 2, 1, 2}));
    ASSERT_EQ(data.prior_segmentation.size(), 4u);
    EXPECT_EQ(data.prior_segmentation[0], 1);
    EXPECT_EQ(data.prior_segmentation[2], 2);
}

TEST(Cov1Data_LoadParquet, LargeStringGene) {
    TempDir dir("cov1_data");
    auto x = arr_f64({1, 2});
    auto y = arr_f64({1, 2});
    auto gene = arr_lstr({"GeneX", "GeneY"});
    auto path = write_parquet(dir, "lstr.parquet",
                              {{"x", x}, {"y", y}, {"gene", gene}},
                              /*store_schema=*/true);
    auto data = baysor::load_molecules(path, baysor::MoleculeInputOptions{});
    ASSERT_EQ(data.n_molecules(), 2);
    EXPECT_EQ(data.gene_names, (std::vector<std::string>{"GeneX", "GeneY"}));
}

TEST(Cov1Data_LoadParquet, DictionaryLargeStringGene) {
    TempDir dir("cov1_data");
    auto x = arr_f64({1, 2, 3});
    auto y = arr_f64({1, 1, 1});
    auto gene = arr_dict(arr_lstr({"LG1", "LG2"}), {0, 1, 0});

    auto path = write_parquet(dir, "dict_large.parquet",
                              {{"x", x}, {"y", y}, {"gene", gene}},
                              /*store_schema=*/true);
    auto data = baysor::load_molecules(path, baysor::MoleculeInputOptions{});
    ASSERT_EQ(data.n_molecules(), 3);
    EXPECT_EQ(data.gene_names, (std::vector<std::string>{"LG1", "LG2"}));
    EXPECT_EQ(data.gene, (std::vector<int>{1, 2, 1}));
}

TEST(Cov1Data_LoadParquet, BinaryDictionaryGene) {
    TempDir dir("cov1_data");
    auto x = arr_f64({1, 2});
    auto y = arr_f64({1, 2});
    auto gene = arr_dict(arr_bin({"bin1", "bin2"}), {0, 1});
    auto path = write_parquet(dir, "dict_bin.parquet",
                              {{"x", x}, {"y", y}, {"gene", gene}},
                              /*store_schema=*/true);
    auto data = baysor::load_molecules(path, baysor::MoleculeInputOptions{});
    ASSERT_EQ(data.n_molecules(), 2);
    ASSERT_EQ(data.n_genes(), 2);
    // Binary dictionary values are decoded as text (GetView), so the gene
    // names are exactly the original strings, sorted.
    EXPECT_EQ(data.gene_names, (std::vector<std::string>{"bin1", "bin2"}));
    EXPECT_EQ(data.gene, (std::vector<int>{1, 2}));
}

TEST(Cov1Data_LoadParquet, NumericDictionaryGene) {
    TempDir dir("cov1_data");
    auto x = arr_f64({1, 2});
    auto y = arr_f64({1, 2});
    // Dictionary-encoded numeric gene column. Parquet does not keep numeric
    // dictionaries, so it is read back as a plain int64 column and goes
    // through the numeric gene-name path ("7", "9").
    auto gene = arr_dict(arr_i64({7, 9}), {0, 1});
    auto path = write_parquet(dir, "dict_num.parquet",
                              {{"x", x}, {"y", y}, {"gene", gene}},
                              /*store_schema=*/true);
    auto data = baysor::load_molecules(path, baysor::MoleculeInputOptions{});
    ASSERT_EQ(data.n_molecules(), 2);
    EXPECT_EQ(data.gene_names, (std::vector<std::string>{"7", "9"}));
    EXPECT_EQ(data.gene, (std::vector<int>{1, 2}));
}

TEST(Cov1Data_LoadParquet, NumericGeneColumn) {
    TempDir dir("cov1_data");
    auto x = arr_f64({1, 2, 3});
    auto y = arr_f64({1, 1, 1});
    auto gene = arr_i64({7, 7, 9});
    auto path = write_parquet(dir, "numgene.parquet", {{"x", x}, {"y", y}, {"gene", gene}});
    auto data = baysor::load_molecules(path, baysor::MoleculeInputOptions{});
    ASSERT_EQ(data.n_molecules(), 3);
    EXPECT_EQ(data.gene_names, (std::vector<std::string>{"7", "9"}));
}

TEST(Cov1Data_LoadParquet, ExcludeGenePatternSpecialCharacters) {
    TempDir dir("cov1_data");
    auto x = arr_f64({1, 2, 3, 4});
    auto y = arr_f64({1, 1, 1, 1});
    auto gene = arr_str({"Blank-1", "GeneA", "GeneB", "MALAT1"});
    auto path = write_parquet(dir, "excl.parquet", {{"x", x}, {"y", y}, {"gene", gene}});

    baysor::MoleculeInputOptions opts;
    // One pattern containing every metacharacter the compiler escapes (none
    // of which match), plus patterns that actually match some genes.
    opts.exclude_genes = "Blank*,MALAT1,Q?rph,never[matches]^,$,x|y,(z)+,{2}\\.";
    auto data = baysor::load_molecules(path, opts);

    EXPECT_EQ(data.n_molecules(), 2);
    EXPECT_EQ(data.gene_names, (std::vector<std::string>{"GeneA", "GeneB"}));
}

TEST(Cov1Data_LoadParquet, MissingRequiredColumn) {
    TempDir dir("cov1_data");
    auto x = arr_f64({1});
    auto gene = arr_str({"A"});
    auto path = write_parquet(dir, "nocol.parquet", {{"x", x}, {"gene", gene}});
    EXPECT_THROW_MSG(baysor::load_molecules(path, baysor::MoleculeInputOptions{}),
                     std::runtime_error, "Column 'y' not found in the data");
}

// ============================================================================
// data.cpp — gene filtering helpers
// ============================================================================

TEST(Cov1Data_Genes, FilterByPatternSpecialCharacters) {
    baysor::MoleculeData data;
    baysor::encode_genes(data, {"Blank-1", "GeneA", "GeneB", "MALAT1", "MT-CO1"});
    data.x = {1, 2, 3, 4, 5};
    data.y = {1, 1, 1, 1, 1};

    // Patterns exercise every escaped metacharacter (none of which match)
    // plus three real matches.
    baysor::filter_genes_by_pattern(
        data, {"Blank*", "MALAT1", "MT*CO1", "A?B", "weird[chars]^", "$", "a|b",
               "(c)+", "{d}\\.", "GeneB"});

    EXPECT_EQ(data.n_molecules(), 1);
    EXPECT_EQ(data.gene_names, (std::vector<std::string>{"GeneA"}));
    EXPECT_EQ(data.gene, (std::vector<int>{1}));

    // A pattern list matching nothing leaves the data untouched.
    baysor::MoleculeData untouched;
    baysor::encode_genes(untouched, {"A", "B"});
    untouched.x = {1, 2};
    untouched.y = {1, 2};
    baysor::filter_genes_by_pattern(untouched, {"NoSuchGene"});
    EXPECT_EQ(untouched.n_molecules(), 2);
    EXPECT_EQ(untouched.n_genes(), 2);

    // Empty pattern list is a no-op.
    baysor::filter_genes_by_pattern(untouched, {});
    EXPECT_EQ(untouched.n_molecules(), 2);
}

// ============================================================================
// prior_segmentation.cpp — spec parsing and label encoding
// ============================================================================

TEST(Cov1Data_Prior, SpecParsingEdgeCases) {
    // Boundary detection covers .csv/.parquet/.pq; everything else is an image.
    EXPECT_EQ(baysor::parse_prior_input_spec("/x/mask.tiff").type,
              baysor::PriorInputType::Image);
    EXPECT_EQ(baysor::parse_prior_input_spec("/x/boundaries.pq").type,
              baysor::PriorInputType::Boundary);
    EXPECT_EQ(baysor::parse_prior_input_spec(":col").column_name, "col");
    EXPECT_EQ(baysor::parse_prior_input_spec("").type, baysor::PriorInputType::None);
}

TEST(Cov1Data_Prior, EncodePriorLabelsFilteringAndWarnings) {
    // Small segments are zeroed out.
    std::vector<std::string> raw = {"a", "a", "b", "0", "a", "c"};
    auto encoded = baysor::encode_prior_labels(raw, "0", /*min_molecules_per_segment=*/3);
    EXPECT_EQ(encoded, (std::vector<int>{1, 1, 0, 0, 1, 0}));

    // No unassigned values at all still encodes successfully.
    std::vector<std::string> all = {"x", "y", "x"};
    auto encoded2 = baysor::encode_prior_labels(all, "0", 1);
    EXPECT_EQ(encoded2, (std::vector<int>{1, 2, 1}));

    // Empty values count as unassigned.
    std::vector<std::string> with_empty = {"a", "", "a"};
    auto encoded3 = baysor::encode_prior_labels(with_empty, "UNASSIGNED", 0);
    EXPECT_EQ(encoded3, (std::vector<int>{1, 0, 1}));
}

TEST(Cov1Data_Prior, EstimateScaleRejectsAllUnassigned) {
    Eigen::MatrixXd pos(2, 4);
    pos << 0, 1, 2, 3,
           0, 0, 0, 0;
    std::vector<int> zeros(4, 0);
    EXPECT_THROW_MSG(baysor::estimate_scale_from_assignment(pos, zeros, 2),
                     std::runtime_error, "No assigned molecules");
}

TEST(Cov1Data_Prior, EstimateScaleFromImageAreas) {
    // Odd count works: radii = sqrt(area/pi).
    const double r = std::sqrt(100.0 / baysor::kPi);
    auto odd = baysor::estimate_scale_from_image_areas({100, 100, 100});
    EXPECT_NEAR(odd.first, r, 1e-12);
    EXPECT_NEAR(odd.second, 0.0, 1e-12);

    // Even count takes the average of the two middle radii.
    auto even = baysor::estimate_scale_from_image_areas({100, 100, 400, 400});
    const double r2 = std::sqrt(400.0 / baysor::kPi);
    EXPECT_NEAR(even.first, (r + r2) / 2.0, 1e-12);

    // Zero-area components are skipped before the n < 3 check.
    EXPECT_THROW_MSG(baysor::estimate_scale_from_image_areas({100, 100, 0, 0}),
                     std::runtime_error, "Not enough prior cells");
    EXPECT_THROW(baysor::estimate_scale_from_image_areas({}), std::runtime_error);
}

// ============================================================================
// prior_segmentation.cpp — boundary files
// ============================================================================

namespace {

baysor::MoleculeData make_molecules(const std::vector<double>& x,
                                    const std::vector<double>& y) {
    baysor::MoleculeData data;
    data.x = x;
    data.y = y;
    data.gene.assign(x.size(), 1);
    data.gene_names = {"G1"};
    return data;
}

std::string write_boundary_csv(const TempDir& dir, const std::string& name,
                               const std::string& label_col, bool far_polygon) {
    std::ostringstream ss;
    ss << "vertex_x,vertex_y," << label_col << "\n";
    auto square = [&](double x0, double y0, const std::string& lab) {
        ss << x0 << "," << y0 << "," << lab << "\n"
           << (x0 + 3) << "," << y0 << "," << lab << "\n"
           << (x0 + 3) << "," << (y0 + 3) << "," << lab << "\n"
           << x0 << "," << (y0 + 3) << "," << lab << "\n";
    };
    square(0, 0, "1");
    square(5, 0, "2");
    square(10, 0, "3");
    if (far_polygon) square(1000, 1000, "4");
    return dir.write(name, ss.str());
}

} // namespace

TEST(Cov1Data_Prior, BoundaryWithLabelIdColumn) {
    TempDir dir("cov1_data");
    const auto path = write_boundary_csv(dir, "b.csv", "label_id", /*far_polygon=*/true);

    std::vector<double> x{1, 1.5, 2, 2.5, 1.2, 6, 6.5, 7, 7.5, 6.2, 11, 11.5, 12, 12.5, 11.2};
    std::vector<double> y{1, 1.5, 2, 2.5, 1.2, 1, 1.5, 2, 2.5, 1.2, 1, 1.5, 2, 2.5, 1.2};
    auto data = make_molecules(x, y);

    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Boundary;
    prior.path = path;
    prior.min_molecules_per_segment = 2;
    prior.estimate_scale_from_prior = true;

    auto [scale, scale_std] = baysor::load_prior_segmentation(data, prior, /*min_cells=*/2);

    // The far-away polygon (outside molecule bounds) is filtered out; the
    // three overlapping squares keep their labels.
    EXPECT_EQ(data.prior_segmentation,
              (std::vector<int>{1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3}));
    EXPECT_GT(scale, 0.0);
    EXPECT_GE(scale_std, 0.0);
}

TEST(Cov1Data_Prior, BoundaryScaleEstimationFailureIsSwallowed) {
    TempDir dir("cov1_data");
    const auto path = write_boundary_csv(dir, "b2.csv", "label_id", /*far_polygon=*/false);

    std::vector<double> x{1, 1.5, 2, 6, 6.5, 7, 11, 11.5, 12};
    std::vector<double> y{1, 1.5, 2, 1, 1.5, 2, 1, 1.5, 2};
    auto data = make_molecules(x, y);

    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Boundary;
    prior.path = path;
    prior.min_molecules_per_segment = 1;

    // Only 3 molecules per segment pass min_molecules_per_cell = 5, so the
    // estimate fails and is swallowed: scale stays at the sentinel.
    auto [scale, scale_std] = baysor::load_prior_segmentation(data, prior, /*min_cells=*/5);
    EXPECT_EQ(scale, -1.0);
    EXPECT_EQ(scale_std, -1.0);
    EXPECT_EQ(data.prior_segmentation, (std::vector<int>{1, 1, 1, 2, 2, 2, 3, 3, 3}));
}

TEST(Cov1Data_Prior, ColumnPriorNotLoadedThrows) {
    auto data = make_molecules({1, 2}, {1, 2});
    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Column;
    prior.column_name = "cell_id";
    EXPECT_THROW_MSG(baysor::load_prior_segmentation(data, prior, 3),
                     std::runtime_error, "no prior labels were loaded");
}

TEST(Cov1Data_Prior, NoneReturnsSentinel) {
    auto data = make_molecules({1, 2}, {1, 2});
    baysor::PriorInputOptions prior; // type == None
    auto result = baysor::load_prior_segmentation(data, prior, 3);
    EXPECT_EQ(result.first, -1.0);
    EXPECT_EQ(result.second, -1.0);
}

// ============================================================================
// prior_segmentation.cpp — TIFF images
// ============================================================================

TEST(Cov1Data_Prior, TiffMultiChannelRejected) {
    TempDir dir("cov1_data");
    std::vector<uint8_t> rgb(8 * 8 * 3, 255);
    const auto path = write_tiff(dir, "rgb.tif", rgb.data(), 8 * 3, 8, 8, 8, /*spp=*/3);
    EXPECT_THROW_MSG(baysor::load_prior_from_image(path, {2}, {2}, 1),
                     std::runtime_error, "Only single-channel TIFF masks");
}

TEST(Cov1Data_Prior, TiffUnsupportedBitsPerSample) {
    TempDir dir("cov1_data");
    // 8 pixels wide at 4 bits/sample = 4 bytes per row.
    std::vector<uint8_t> rows(4 * 8, 0xF0);
    const auto path = write_tiff(dir, "bps4.tif", rows.data(), 4, 8, 8, /*bps=*/4);
    EXPECT_THROW_MSG(baysor::load_prior_from_image(path, {2, 3}, {2, 3}, 1),
                     std::runtime_error, "Unsupported TIFF bits/sample: 4");
}

TEST(Cov1Data_Prior, TiffCorruptScanlineThrows) {
    TempDir dir("cov1_data");
    // Deflate-compressed mask whose pixel data is corrupted after writing:
    // the directory still parses, but decoding a scanline fails.
    std::vector<uint8_t> px(8 * 8, 255);
    const auto path = write_tiff(dir, "corrupt.tif", px.data(), 8, 8, 8,
                                 /*bps=*/8, /*spp=*/1, /*deflate=*/true);
    corrupt_tiff_pixel_data(path);

    EXPECT_THROW_MSG(baysor::load_prior_from_image(path, {2, 3}, {2, 3}, 1),
                     std::runtime_error, "Error reading TIFF scanline");
}

TEST(Cov1Data_Prior, TiffNoMoleculesYieldsEmptyResult) {
    TempDir dir("cov1_data");
    std::vector<uint8_t> px(8 * 8, 255);
    const auto path = write_tiff_mask(dir, "m8.tif", px, 8, 8);

    // No molecules at all.
    auto res = baysor::load_prior_from_image(path, {}, {}, 1);
    EXPECT_TRUE(res.segment_per_molecule.empty());
    EXPECT_TRUE(res.component_pixel_areas.empty());
}

TEST(Cov1Data_Prior, TiffMoleculesOutOfBounds) {
    TempDir dir("cov1_data");
    std::vector<uint8_t> px(8 * 8, 255);
    const auto path = write_tiff_mask(dir, "m8b.tif", px, 8, 8);

    // All molecules outside the image -> nothing to load.
    auto res = baysor::load_prior_from_image(path, {100, 200}, {100, 200}, 1);
    EXPECT_EQ(res.segment_per_molecule, (std::vector<int>{0, 0}));

    // Mixed: only the in-bounds molecule is considered (and lands on mask).
    auto res2 = baysor::load_prior_from_image(path, {2, 100}, {2, 100}, 1);
    ASSERT_EQ(res2.segment_per_molecule.size(), 2u);
    EXPECT_GT(res2.segment_per_molecule[0], 0);
    EXPECT_EQ(res2.segment_per_molecule[1], 0);
}

TEST(Cov1Data_Prior, TiffFullImageWindowLogs) {
    TempDir dir("cov1_data");
    // Molecules at the extreme corners -> window equals the full image.
    std::vector<uint8_t> px(6 * 6, 255);
    const auto path = write_tiff_mask(dir, "full.tif", px, 6, 6);
    auto res = baysor::load_prior_from_image(path, {1, 6, 3}, {1, 6, 3}, 1);
    ASSERT_EQ(res.segment_per_molecule.size(), 3u);
    EXPECT_GT(res.segment_per_molecule[0], 0);
    // Binary mask -> one connected component with area 36.
    ASSERT_EQ(res.component_pixel_areas.size(), 1u);
    EXPECT_EQ(res.component_pixel_areas[0], 36u);
}

TEST(Cov1Data_Prior, Tiff16And32BitMultiLabel) {
    TempDir dir("cov1_data");
    const std::vector<double> x{2, 2, 2, 6, 6, 6};
    const std::vector<double> y{2, 3, 4, 2, 3, 4};
    auto check = [&](const std::string& path, int left, int right) {
        auto res = baysor::load_prior_from_image(path, x, y, /*min=*/1);
        EXPECT_EQ(res.segment_per_molecule, (std::vector<int>{left, left, left, right, right, right}));
        auto areas = res.component_pixel_areas;
        std::sort(areas.begin(), areas.end());
        EXPECT_EQ(areas, (std::vector<size_t>{6, 9}));
        // Molecule-count filtering zeroes every label when min is too large.
        EXPECT_EQ(baysor::load_prior_from_image(path, x, y, /*min=*/4).segment_per_molecule, std::vector<int>(6, 0));
    };
    // 8x8 masks: label `left` in columns 0-3, `right` in columns 4-7.
    auto mask = [](auto left, auto right) {
        std::vector<decltype(left)> px(8 * 8);
        for (size_t i = 0; i < px.size(); ++i) px[i] = (i % 8 < 4) ? left : right;
        return px;
    };
    check(write_tiff_mask(dir, "m16.tif", mask(uint16_t{10}, uint16_t{20}), 8, 8), 10, 20);
    check(write_tiff_mask(dir, "m32.tif", mask(uint32_t{100}, uint32_t{200}), 8, 8), 100, 200);
}

TEST(Cov1Data_Prior, ImageScaleFallbackWhenNoMoleculesInMask) {
    TempDir dir("cov1_data");
    // Valid image but molecules only sit on background pixels -> no components,
    // so the scale estimate falls back to the (all-unassigned) assignment and
    // fails, leaving the -1 sentinel.
    std::vector<uint8_t> px(8 * 8, 0);
    const auto path = write_tiff_mask(dir, "empty_mask.tif", px, 8, 8);

    auto data = make_molecules({2, 3, 4}, {2, 3, 4});
    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Image;
    prior.path = path;
    prior.min_molecules_per_segment = 1;
    prior.estimate_scale_from_prior = true;

    auto [scale, scale_std] = baysor::load_prior_segmentation(data, prior, 2);
    EXPECT_EQ(scale, -1.0);
    EXPECT_EQ(scale_std, -1.0);
    EXPECT_EQ(data.prior_segmentation, (std::vector<int>{0, 0, 0}));
}

TEST(Cov1Data_Prior, ImageScaleFromComponents) {
    TempDir dir("cov1_data");
    // Binary mask with three separate 2x2 blobs (area 4 each), two molecules
    // per blob -> the area-based estimator has 3 components and succeeds.
    std::vector<uint8_t> px(8 * 8, 0);
    auto blob = [&](uint32_t r0, uint32_t c0) {
        for (uint32_t r = r0; r <= r0 + 1; ++r)
            for (uint32_t c = c0; c <= c0 + 1; ++c) px[r * 8 + c] = 255;
    };
    blob(1, 1);  // molecules (2,2), (2,3)
    blob(5, 5);  // molecules (6,6), (6,7)
    blob(5, 1);  // molecules (2,6), (2,7)
    const auto path = write_tiff_mask(dir, "three_blobs.tif", px, 8, 8);

    auto data = make_molecules({2, 2, 6, 6, 2, 2}, {2, 3, 6, 7, 6, 7});
    baysor::PriorInputOptions prior;
    prior.type = baysor::PriorInputType::Image;
    prior.path = path;
    prior.min_molecules_per_segment = 1;
    prior.estimate_scale_from_prior = true;

    auto [scale, scale_std] = baysor::load_prior_segmentation(data, prior, 2);
    // Every component has area 4 -> radius sqrt(4/pi), zero MAD.
    EXPECT_NEAR(scale, std::sqrt(4.0 / baysor::kPi), 1e-9);
    EXPECT_NEAR(scale_std, 0.0, 1e-9);
    EXPECT_EQ(data.prior_segmentation, (std::vector<int>{1, 1, 2, 2, 3, 3}));
}
