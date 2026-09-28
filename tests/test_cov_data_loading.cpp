// Coverage tests for src/data_loading/{data,prior_segmentation}.cpp and
// include/baysor/data_loading/*.h (task COV-1).
//
// File-local helpers live in an anonymous namespace; every suite name is
// prefixed with Cov1 so it cannot clash with other coverage test files.

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
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Temporary files
// ---------------------------------------------------------------------------

class TempDir {
public:
    TempDir() {
        static std::atomic<int> counter{0};
        const auto base = std::filesystem::temp_directory_path();
        path_ = base / ("baysor_cov1_data_" + std::to_string(::getpid()) + "_" +
                        std::to_string(counter++));
        std::error_code ec;
        std::filesystem::remove_all(path_, ec);
        std::filesystem::create_directories(path_);
    }
    ~TempDir() {
        std::error_code ec;
        std::filesystem::remove_all(path_, ec);
    }
    std::string file(const std::string& name) const {
        return (path_ / name).string();
    }
    const std::filesystem::path& path() const { return path_; }

private:
    std::filesystem::path path_;
};

std::string write_csv(const TempDir& dir, const std::string& name,
                      const std::string& content) {
    const std::string p = dir.file(name);
    std::ofstream f(p);
    f << content;
    return p;
}

// ---------------------------------------------------------------------------
// Arrow array builders
// ---------------------------------------------------------------------------

template <typename Builder>
std::shared_ptr<arrow::Array> finish(Builder& b) {
    std::shared_ptr<arrow::Array> out;
    EXPECT_TRUE(b.Finish(&out).ok());
    return out;
}

std::shared_ptr<arrow::Array> arr_f64(const std::vector<double>& v) {
    arrow::DoubleBuilder b;
    b.AppendValues(v);
    return finish(b);
}

std::shared_ptr<arrow::Array> arr_f64_nulls(const std::vector<std::optional<double>>& v) {
    arrow::DoubleBuilder b;
    for (const auto& x : v) {
        if (x.has_value()) b.Append(*x);
        else b.AppendNull();
    }
    return finish(b);
}

std::shared_ptr<arrow::Array> arr_f32(const std::vector<float>& v) {
    arrow::FloatBuilder b;
    b.AppendValues(v);
    return finish(b);
}

std::shared_ptr<arrow::Array> arr_i64(const std::vector<int64_t>& v) {
    arrow::Int64Builder b;
    b.AppendValues(v);
    return finish(b);
}

std::shared_ptr<arrow::Array> arr_i64_nulls(const std::vector<std::optional<int64_t>>& v) {
    arrow::Int64Builder b;
    for (const auto& x : v) {
        if (x.has_value()) b.Append(*x);
        else b.AppendNull();
    }
    return finish(b);
}

std::shared_ptr<arrow::Array> arr_i32(const std::vector<int32_t>& v) {
    arrow::Int32Builder b;
    b.AppendValues(v);
    return finish(b);
}

std::shared_ptr<arrow::Array> arr_i16(const std::vector<int16_t>& v) {
    arrow::Int16Builder b;
    b.AppendValues(v);
    return finish(b);
}

std::shared_ptr<arrow::Array> arr_i8(const std::vector<int8_t>& v) {
    arrow::Int8Builder b;
    b.AppendValues(v);
    return finish(b);
}

std::shared_ptr<arrow::Array> arr_u64(const std::vector<uint64_t>& v) {
    arrow::UInt64Builder b;
    b.AppendValues(v);
    return finish(b);
}

std::shared_ptr<arrow::Array> arr_u32(const std::vector<uint32_t>& v) {
    arrow::UInt32Builder b;
    b.AppendValues(v);
    return finish(b);
}

std::shared_ptr<arrow::Array> arr_u16(const std::vector<uint16_t>& v) {
    arrow::UInt16Builder b;
    b.AppendValues(v);
    return finish(b);
}

std::shared_ptr<arrow::Array> arr_u8(const std::vector<uint8_t>& v) {
    arrow::UInt8Builder b;
    b.AppendValues(v);
    return finish(b);
}

std::shared_ptr<arrow::Array> arr_str(const std::vector<std::string>& v) {
    arrow::StringBuilder b;
    for (const auto& s : v) EXPECT_TRUE(b.Append(s).ok());
    return finish(b);
}

std::shared_ptr<arrow::Array> arr_lstr(const std::vector<std::string>& v) {
    arrow::LargeStringBuilder b;
    for (const auto& s : v) EXPECT_TRUE(b.Append(s).ok());
    return finish(b);
}

std::shared_ptr<arrow::Array> arr_bin(const std::vector<std::string>& v) {
    arrow::BinaryBuilder b;
    for (const auto& s : v) EXPECT_TRUE(b.Append(s).ok());
    return finish(b);
}

// Dictionary-encoded string array with explicit index type; optional nulls
// (null entries appear at the given indices).
std::shared_ptr<arrow::Array> arr_dict(const std::vector<std::string>& dict_values,
                                       const std::vector<int>& indices_with_nulls,
                                       const std::shared_ptr<arrow::DataType>& index_type) {
    std::shared_ptr<arrow::Array> dict = arr_str(dict_values);
    std::shared_ptr<arrow::Array> indices;
    if (index_type->id() == arrow::Type::INT8) {
        arrow::Int8Builder b;
        for (int i : indices_with_nulls) {
            if (i < 0) b.AppendNull();
            else b.Append(static_cast<int8_t>(i));
        }
        indices = finish(b);
    } else if (index_type->id() == arrow::Type::INT16) {
        arrow::Int16Builder b;
        for (int i : indices_with_nulls) {
            if (i < 0) b.AppendNull();
            else b.Append(static_cast<int16_t>(i));
        }
        indices = finish(b);
    } else if (index_type->id() == arrow::Type::INT32) {
        arrow::Int32Builder b;
        for (int i : indices_with_nulls) {
            if (i < 0) b.AppendNull();
            else b.Append(static_cast<int32_t>(i));
        }
        indices = finish(b);
    } else if (index_type->id() == arrow::Type::INT64) {
        arrow::Int64Builder b;
        for (int i : indices_with_nulls) {
            if (i < 0) b.AppendNull();
            else b.Append(static_cast<int64_t>(i));
        }
        indices = finish(b);
    } else {
        ADD_FAILURE() << "unsupported dictionary index type";
        return nullptr;
    }
    auto result = arrow::DictionaryArray::FromArrays(indices, dict);
    EXPECT_TRUE(result.ok()) << result.status().ToString();
    return result.ValueOrDie();
}

// Dictionary-encoded binary array (dictionary value type is NOT a string type).
std::shared_ptr<arrow::Array> arr_dict_binary(const std::vector<std::string>& dict_values,
                                              const std::vector<int>& indices) {
    std::shared_ptr<arrow::Array> dict = arr_bin(dict_values);
    arrow::Int8Builder b;
    for (int i : indices) b.Append(static_cast<int8_t>(i));
    std::shared_ptr<arrow::Array> idx = finish(b);
    auto result = arrow::DictionaryArray::FromArrays(idx, dict);
    EXPECT_TRUE(result.ok()) << result.status().ToString();
    return result.ValueOrDie();
}

std::shared_ptr<arrow::Array> arr_list_i32(const std::vector<std::vector<int32_t>>& rows) {
    auto vb = std::make_shared<arrow::Int32Builder>();
    arrow::ListBuilder lb(arrow::default_memory_pool(), vb);
    for (const auto& row : rows) {
        EXPECT_TRUE(lb.Append(true).ok());
        for (int32_t v : row) {
            EXPECT_TRUE(vb->Append(v).ok());
        }
    }
    return finish(lb);
}

// ---------------------------------------------------------------------------
// Parquet writer (optionally storing the Arrow schema so that dictionary and
// large_string columns survive the roundtrip)
// ---------------------------------------------------------------------------

std::string write_parquet(const TempDir& dir, const std::string& name,
                          const std::vector<std::shared_ptr<arrow::Field>>& fields,
                          const std::vector<std::shared_ptr<arrow::Array>>& arrays,
                          bool store_schema = false) {
    const std::string path = dir.file(name);
    auto schema = arrow::schema(fields);
    auto table = arrow::Table::Make(schema, arrays);
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

std::shared_ptr<arrow::Field> field(const std::string& name,
                                    const std::shared_ptr<arrow::Array>& arr) {
    return arrow::field(name, arr->type());
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

std::string write_tiff_u8(const TempDir& dir, const std::string& name,
                          const std::vector<uint8_t>& pixels, uint32_t w, uint32_t h) {
    return write_tiff(dir, name, pixels.data(), w, w, h, 8);
}

std::string write_tiff_u16(const TempDir& dir, const std::string& name,
                           const std::vector<uint16_t>& pixels, uint32_t w, uint32_t h) {
    return write_tiff(dir, name, reinterpret_cast<const uint8_t*>(pixels.data()),
                      static_cast<size_t>(w) * 2, w, h, 16);
}

std::string write_tiff_u32(const TempDir& dir, const std::string& name,
                           const std::vector<uint32_t>& pixels, uint32_t w, uint32_t h) {
    return write_tiff(dir, name, reinterpret_cast<const uint8_t*>(pixels.data()),
                      static_cast<size_t>(w) * 4, w, h, 32);
}

baysor::MoleculeInputOptions default_opts() {
    return baysor::MoleculeInputOptions{};
}

} // namespace

// ============================================================================
// src/data_loading/data.cpp — column readers
// ============================================================================

TEST(Cov1Data_Readers, ArrowErrorOnMissingFile) {
    TempDir dir;
    try {
        baysor::read_double_column(dir.file("does_not_exist.csv"), "x");
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("Arrow error:"), std::string::npos);
    }
    try {
        baysor::read_string_column(dir.file("does_not_exist.parquet"), "gene");
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("Arrow error:"), std::string::npos);
    }
}

TEST(Cov1Data_Readers, UnsupportedFileFormat) {
    TempDir dir;
    const auto path = write_csv(dir, "data.txt", "x,y\n1,2\n");
    try {
        baysor::read_double_column(path, "x");
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("Unsupported file format: .txt"),
                  std::string::npos);
    }
}

TEST(Cov1Data_Readers, ParquetAndPqExtensions) {
    TempDir dir;
    auto x = arr_f64({1.5, 2.5, 3.5});
    auto y = arr_f64({4.0, 5.0, 6.0});
    auto p1 = write_parquet(dir, "a.parquet", {field("x", x), field("y", y)}, {x, y});
    auto p2 = write_parquet(dir, "b.pq", {field("x", x), field("y", y)}, {x, y});

    for (const auto& p : {p1, p2}) {
        auto vals = baysor::read_double_column(p, "x");
        ASSERT_EQ(vals.size(), 3u);
        EXPECT_DOUBLE_EQ(vals[0], 1.5);
        EXPECT_DOUBLE_EQ(vals[2], 3.5);
    }
}

TEST(Cov1Data_Readers, MissingColumnMessage) {
    TempDir dir;
    const auto path = write_csv(dir, "d.csv", "x,y\n1,2\n");
    try {
        baysor::read_double_column(path, "gene");
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        const std::string msg = e.what();
        EXPECT_NE(msg.find("Column 'gene' not found in the data"), std::string::npos);
        EXPECT_NE(msg.find("Available columns"), std::string::npos);
    }
}

TEST(Cov1Data_Readers, DoubleColumnNumericTypes) {
    TempDir dir;
    auto xf = arr_f32({1.5f, 2.5f});
    auto x32 = arr_i32({7, -3});
    auto x16 = arr_i16({5, 6});
    auto x8 = arr_u8({9, 10});
    auto x64 = arr_i64({11, 12});
    auto path = write_parquet(dir, "types.parquet",
                              {field("xf", xf), field("x32", x32), field("x16", x16),
                               field("x8", x8), field("x64", x64)},
                              {xf, x32, x16, x8, x64});

    {
        auto v = baysor::read_double_column(path, "xf");
        ASSERT_EQ(v.size(), 2u);
        EXPECT_DOUBLE_EQ(v[1], 2.5);
    }
    {
        auto v = baysor::read_double_column(path, "x32");
        EXPECT_DOUBLE_EQ(v[1], -3.0);
    }
    {
        auto v = baysor::read_double_column(path, "x16");
        EXPECT_DOUBLE_EQ(v[0], 5.0);
    }
    {
        auto v = baysor::read_double_column(path, "x8");
        EXPECT_DOUBLE_EQ(v[1], 10.0);
    }
    {
        auto v = baysor::read_double_column(path, "x64");
        EXPECT_DOUBLE_EQ(v[0], 11.0);
    }
}

TEST(Cov1Data_Readers, DoubleColumnCastFailureThrows) {
    TempDir dir;
    auto lists = arr_list_i32({{1, 2}, {3}});
    auto xs = arr_f64({1.0, 2.0});
    auto path = write_parquet(dir, "lists.parquet",
                              {field("bad", lists), field("x", xs)}, {lists, xs});
    try {
        baysor::read_double_column(path, "bad");
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("Cannot convert column 'bad' to double"),
                  std::string::npos);
    }
}

TEST(Cov1Data_Readers, StringColumnTypes) {
    TempDir dir;
    auto plain = arr_str({"a", "b"});
    auto large = arr_lstr({"big1", "big2"});
    auto d8 = arr_dict({"G1", "G2"}, {0, 1}, arrow::int8());
    auto d16 = arr_dict({"G1", "G2"}, {1, 0}, arrow::int16());
    auto d32 = arr_dict({"G1", "G2"}, {0, 0}, arrow::int32());
    auto d64 = arr_dict({"G1", "G2"}, {1, 1}, arrow::int64());
    auto num = arr_i64({101, 202});

    auto path = write_parquet(dir, "strings.parquet",
                              {field("plain", plain), field("large", large),
                               field("d8", d8), field("d16", d16), field("d32", d32),
                               field("d64", d64), field("num", num)},
                              {plain, large, d8, d16, d32, d64, num},
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
    TempDir dir;
    const auto path = write_csv(dir, "z3d.csv",
        "x,y,z,gene\n"
        "1,2,10,A\n"
        "3,4,20,B\n"
        "5,6,30,A\n");
    auto raw = baysor::read_tabular_file(path, default_opts());
    EXPECT_TRUE(raw.has_z);
    ASSERT_EQ(raw.z.size(), 3u);
    EXPECT_DOUBLE_EQ(raw.z[2], 30.0);
    EXPECT_EQ(raw.gene_str, (std::vector<std::string>{"A", "B", "A"}));
    ASSERT_EQ(raw.x.size(), 3u);
    EXPECT_DOUBLE_EQ(raw.x[1], 3.0);
}

TEST(Cov1Data_ReadTabular, DropsConstantZ) {
    TempDir dir;
    const auto path = write_csv(dir, "zconst.csv",
        "x,y,z,gene\n"
        "1,2,7,A\n"
        "3,4,7,B\n");
    auto raw = baysor::read_tabular_file(path, default_opts());
    EXPECT_FALSE(raw.has_z);
    EXPECT_TRUE(raw.z.empty());
}

TEST(Cov1Data_ReadTabular, Force2DSkipsZColumn) {
    TempDir dir;
    const auto path = write_csv(dir, "zforce.csv",
        "x,y,z,gene\n"
        "1,2,10,A\n"
        "3,4,20,B\n");
    auto opts = default_opts();
    opts.force_2d = true;
    auto raw = baysor::read_tabular_file(path, opts);
    EXPECT_FALSE(raw.has_z);
    EXPECT_TRUE(raw.z.empty());
}

TEST(Cov1Data_ReadTabular, NoZColumn) {
    TempDir dir;
    const auto path = write_csv(dir, "plain.csv",
        "x,y,gene\n"
        "1,2,A\n");
    auto raw = baysor::read_tabular_file(path, default_opts());
    EXPECT_FALSE(raw.has_z);
    ASSERT_EQ(raw.gene_str.size(), 1u);
}

// ============================================================================
// data.cpp — load_molecules (CSV)
// ============================================================================

TEST(Cov1Data_LoadCsv, OptionalMetadataColumnsAndQvFilter) {
    TempDir dir;
    const auto path = write_csv(dir, "meta.csv",
        "x,y,gene,confidence,cluster,nuclei_probs,qv,transcript_id\n"
        "1,1,G1,0.9,1,0.5,30,1000\n"
        "2,2,G2,0.8,2,0.6,40,1001\n"
        "3,3,G1,0.7,1,0.4,10,1002\n");

    auto opts = default_opts();
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
    TempDir dir;
    const auto path = write_csv(dir, "bounds.csv",
        "x,y,gene\n"
        "1,1,A\n"
        "50,1,A\n"
        "3,90,A\n"
        "2,2,B\n");
    auto opts = default_opts();
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
    TempDir dir;
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
                              {field("x", x), field("y", y), field("z", z),
                               field("gene", gene), field("qv", qv),
                               field("confidence", conf), field("cluster", clus),
                               field("nuclei_probs", nuclei),
                               field("transcript_id", tx), field("cell_id", cell)},
                              {x, y, z, gene, qv, conf, clus, nuclei, tx, cell});

    auto opts = default_opts();
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
    TempDir dir;
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
                                  {arrow::field(c.name, c.arr->type()),
                                   field("y", y), field("gene", gene)},
                                  {c.arr, y, gene});
        auto opts = default_opts();
        opts.x_col = c.name;
        auto data = baysor::load_molecules(path, opts);
        ASSERT_EQ(data.n_molecules(), 2) << c.name;
        EXPECT_DOUBLE_EQ(data.x[0], c.expect0) << c.name;
    }

    // A null coordinate is NaN and fails the static row filter.
    auto xnull = arr_f64_nulls({1.0, std::nullopt});
    auto y = arr_f64({1, 2});
    auto gene = arr_str({"A", "B"});
    auto path = write_parquet(dir, "xnull.parquet",
                              {field("x", xnull), field("y", y), field("gene", gene)},
                              {xnull, y, gene});
    auto data = baysor::load_molecules(path, default_opts());
    ASSERT_EQ(data.n_molecules(), 1);
    EXPECT_DOUBLE_EQ(data.x[0], 1.0);
}

TEST(Cov1Data_LoadParquet, TranscriptIdTypes) {
    struct Case {
        const char* name;
        std::shared_ptr<arrow::Array> arr;
        std::vector<std::uint64_t> expect;
    };
    TempDir dir;
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
                                  {field("x", x), field("y", y), field("gene", gene),
                                   arrow::field("transcript_id", c.arr->type())},
                                  {x, y, gene, c.arr});
        auto data = baysor::load_molecules(path, default_opts());
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
                              {field("x", x), field("y", y), field("gene", gene),
                               field("transcript_id", tnull)},
                              {x, y, gene, tnull});
    auto data = baysor::load_molecules(path, default_opts());
    ASSERT_EQ(data.source_transcript_id.size(), 2u);
    EXPECT_EQ(data.source_transcript_id[0], 5u);
    EXPECT_EQ(data.source_transcript_id[1], std::uint64_t(-1));
}

TEST(Cov1Data_LoadParquet, DictionaryEncodedGene) {
    TempDir dir;
    auto x = arr_f64({1, 2, 3, 4, 5});
    auto y = arr_f64({1, 1, 1, 1, 1});
    // gene: [A, B, null, A, B]; cell_id: [c1, c1, null, c2, c2]
    auto gene = arr_dict({"A", "B"}, {0, 1, -1, 0, 1}, arrow::int8());
    auto cell = arr_dict({"c1", "c2"}, {0, 0, -1, 1, 1}, arrow::int8());

    auto path = write_parquet(dir, "dict.parquet",
                              {field("x", x), field("y", y),
                               field("gene", gene), field("cell_id", cell)},
                              {x, y, gene, cell},
                              /*store_schema=*/true);

    auto opts = default_opts();
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
    TempDir dir;
    auto x = arr_f64({1, 2});
    auto y = arr_f64({1, 2});
    auto gene = arr_lstr({"GeneX", "GeneY"});
    auto path = write_parquet(dir, "lstr.parquet",
                              {field("x", x), field("y", y), field("gene", gene)},
                              {x, y, gene},
                              /*store_schema=*/true);
    auto data = baysor::load_molecules(path, default_opts());
    ASSERT_EQ(data.n_molecules(), 2);
    EXPECT_EQ(data.gene_names, (std::vector<std::string>{"GeneX", "GeneY"}));
}

TEST(Cov1Data_LoadParquet, DictionaryLargeStringGene) {
    TempDir dir;
    auto x = arr_f64({1, 2, 3});
    auto y = arr_f64({1, 1, 1});
    // Build the dictionary values as large_string manually.
    auto dict_vals = arr_lstr({"LG1", "LG2"});
    arrow::Int8Builder ib;
    ib.Append(0);
    ib.Append(1);
    ib.Append(0);
    std::shared_ptr<arrow::Array> idx = finish(ib);
    auto gene_res = arrow::DictionaryArray::FromArrays(idx, dict_vals);
    ASSERT_TRUE(gene_res.ok()) << gene_res.status().ToString();
    auto gene = gene_res.ValueOrDie();

    auto path = write_parquet(dir, "dict_large.parquet",
                              {field("x", x), field("y", y), field("gene", gene)},
                              {x, y, gene},
                              /*store_schema=*/true);
    auto data = baysor::load_molecules(path, default_opts());
    ASSERT_EQ(data.n_molecules(), 3);
    EXPECT_EQ(data.gene_names, (std::vector<std::string>{"LG1", "LG2"}));
    EXPECT_EQ(data.gene, (std::vector<int>{1, 2, 1}));
}

TEST(Cov1Data_LoadParquet, BinaryDictionaryGene) {
    TempDir dir;
    auto x = arr_f64({1, 2});
    auto y = arr_f64({1, 2});
    auto gene = arr_dict_binary({"bin1", "bin2"}, {0, 1});
    auto path = write_parquet(dir, "dict_bin.parquet",
                              {field("x", x), field("y", y), field("gene", gene)},
                              {x, y, gene},
                              /*store_schema=*/true);
    auto data = baysor::load_molecules(path, default_opts());
    ASSERT_EQ(data.n_molecules(), 2);
    ASSERT_EQ(data.n_genes(), 2);
    // Binary values stringify differently from plain text, but each distinct
    // dictionary entry still becomes its own gene.
    EXPECT_NE(data.gene_names[0], data.gene_names[1]);
    EXPECT_FALSE(data.gene_names[0].empty());
    EXPECT_NE(data.gene[0], data.gene[1]);
}

TEST(Cov1Data_LoadParquet, NumericGeneColumn) {
    TempDir dir;
    auto x = arr_f64({1, 2, 3});
    auto y = arr_f64({1, 1, 1});
    auto gene = arr_i64({7, 7, 9});
    auto path = write_parquet(dir, "numgene.parquet",
                              {field("x", x), field("y", y), field("gene", gene)},
                              {x, y, gene});
    auto data = baysor::load_molecules(path, default_opts());
    ASSERT_EQ(data.n_molecules(), 3);
    EXPECT_EQ(data.gene_names, (std::vector<std::string>{"7", "9"}));
}

TEST(Cov1Data_LoadParquet, ExcludeGenePatternSpecialCharacters) {
    TempDir dir;
    auto x = arr_f64({1, 2, 3, 4});
    auto y = arr_f64({1, 1, 1, 1});
    auto gene = arr_str({"Blank-1", "GeneA", "GeneB", "MALAT1"});
    auto path = write_parquet(dir, "excl.parquet",
                              {field("x", x), field("y", y), field("gene", gene)},
                              {x, y, gene});

    auto opts = default_opts();
    // One pattern containing every metacharacter the compiler escapes (none
    // of which match), plus patterns that actually match some genes.
    opts.exclude_genes = "Blank*,MALAT1,Q?rph,never[matches]^,$,x|y,(z)+,{2}\\.";
    auto data = baysor::load_molecules(path, opts);

    EXPECT_EQ(data.n_molecules(), 2);
    EXPECT_EQ(data.gene_names, (std::vector<std::string>{"GeneA", "GeneB"}));
}

TEST(Cov1Data_LoadParquet, MissingRequiredColumn) {
    TempDir dir;
    auto x = arr_f64({1});
    auto gene = arr_str({"A"});
    auto path = write_parquet(dir, "nocol.parquet",
                              {field("x", x), field("gene", gene)}, {x, gene});
    try {
        baysor::load_molecules(path, default_opts());
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("Column 'y' not found in the data"),
                  std::string::npos);
    }
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
    try {
        baysor::estimate_scale_from_assignment(pos, zeros, 2);
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("No assigned molecules"),
                  std::string::npos);
    }
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
    try {
        baysor::estimate_scale_from_image_areas({100, 100, 0, 0});
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("Not enough prior cells"),
                  std::string::npos);
    }
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
    return write_csv(dir, name, ss.str());
}

} // namespace

TEST(Cov1Data_Prior, BoundaryWithLabelIdColumn) {
    TempDir dir;
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
    TempDir dir;
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
    try {
        baysor::load_prior_segmentation(data, prior, 3);
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("no prior labels were loaded"),
                  std::string::npos);
    }
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
    TempDir dir;
    std::vector<uint8_t> rgb(8 * 8 * 3, 255);
    const auto path = write_tiff(dir, "rgb.tif", rgb.data(), 8 * 3, 8, 8, 8, /*spp=*/3);
    try {
        baysor::load_prior_from_image(path, {2}, {2}, 1);
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("Only single-channel TIFF masks"),
                  std::string::npos);
    }
}

TEST(Cov1Data_Prior, TiffUnsupportedBitsPerSample) {
    TempDir dir;
    // 8 pixels wide at 4 bits/sample = 4 bytes per row.
    std::vector<uint8_t> rows(4 * 8, 0xF0);
    const auto path = write_tiff(dir, "bps4.tif", rows.data(), 4, 8, 8, /*bps=*/4);
    try {
        baysor::load_prior_from_image(path, {2, 3}, {2, 3}, 1);
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("Unsupported TIFF bits/sample: 4"),
                  std::string::npos);
    }
}

TEST(Cov1Data_Prior, TiffCorruptScanlineThrows) {
    TempDir dir;
    // Deflate-compressed mask whose pixel data is corrupted after writing:
    // the directory still parses, but decoding a scanline fails.
    std::vector<uint8_t> px(8 * 8, 255);
    const auto path = write_tiff(dir, "corrupt.tif", px.data(), 8, 8, 8,
                                 /*bps=*/8, /*spp=*/1, /*deflate=*/true);
    corrupt_tiff_pixel_data(path);

    try {
        baysor::load_prior_from_image(path, {2, 3}, {2, 3}, 1);
        FAIL() << "expected throw";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("Error reading TIFF scanline"),
                  std::string::npos);
    }
}

TEST(Cov1Data_Prior, TiffNoMoleculesYieldsEmptyResult) {
    TempDir dir;
    std::vector<uint8_t> px(8 * 8, 255);
    const auto path = write_tiff_u8(dir, "m8.tif", px, 8, 8);

    // No molecules at all.
    auto res = baysor::load_prior_from_image(path, {}, {}, 1);
    EXPECT_TRUE(res.segment_per_molecule.empty());
    EXPECT_TRUE(res.component_pixel_areas.empty());
}

TEST(Cov1Data_Prior, TiffMoleculesOutOfBounds) {
    TempDir dir;
    std::vector<uint8_t> px(8 * 8, 255);
    const auto path = write_tiff_u8(dir, "m8b.tif", px, 8, 8);

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
    TempDir dir;
    // Molecules at the extreme corners -> window equals the full image.
    std::vector<uint8_t> px(6 * 6, 255);
    const auto path = write_tiff_u8(dir, "full.tif", px, 6, 6);
    auto res = baysor::load_prior_from_image(path, {1, 6, 3}, {1, 6, 3}, 1);
    ASSERT_EQ(res.segment_per_molecule.size(), 3u);
    EXPECT_GT(res.segment_per_molecule[0], 0);
    // Binary mask -> one connected component with area 36.
    ASSERT_EQ(res.component_pixel_areas.size(), 1u);
    EXPECT_EQ(res.component_pixel_areas[0], 36u);
}

TEST(Cov1Data_Prior, Tiff16BitMultiLabel) {
    TempDir dir;
    std::vector<uint16_t> px(8 * 8);
    for (uint32_t r = 0; r < 8; ++r)
        for (uint32_t c = 0; c < 8; ++c)
            px[r * 8 + c] = (c < 4) ? 10 : 20;
    const auto path = write_tiff_u16(dir, "m16.tif", px, 8, 8);

    std::vector<double> x{2, 2, 2, 6, 6, 6};
    std::vector<double> y{2, 3, 4, 2, 3, 4};

    auto res = baysor::load_prior_from_image(path, x, y, /*min=*/1);
    EXPECT_EQ(res.segment_per_molecule, (std::vector<int>{10, 10, 10, 20, 20, 20}));
    ASSERT_EQ(res.component_pixel_areas.size(), 2u);
    auto areas16 = res.component_pixel_areas;
    std::sort(areas16.begin(), areas16.end());
    EXPECT_EQ(areas16, (std::vector<size_t>{6, 9}));

    // Molecule-count filtering zeroes every label when min is too large.
    auto res2 = baysor::load_prior_from_image(path, x, y, /*min=*/4);
    EXPECT_EQ(res2.segment_per_molecule, std::vector<int>(6, 0));
}

TEST(Cov1Data_Prior, Tiff32BitMultiLabel) {
    TempDir dir;
    std::vector<uint32_t> px(8 * 8);
    for (uint32_t r = 0; r < 8; ++r)
        for (uint32_t c = 0; c < 8; ++c)
            px[r * 8 + c] = (c < 4) ? 100 : 200;
    const auto path = write_tiff_u32(dir, "m32.tif", px, 8, 8);

    std::vector<double> x{2, 2, 6, 6};
    std::vector<double> y{2, 3, 2, 3};

    auto res = baysor::load_prior_from_image(path, x, y, /*min=*/1);
    EXPECT_EQ(res.segment_per_molecule, (std::vector<int>{100, 100, 200, 200}));
    ASSERT_EQ(res.component_pixel_areas.size(), 2u);
    auto areas32 = res.component_pixel_areas;
    std::sort(areas32.begin(), areas32.end());
    EXPECT_EQ(areas32, (std::vector<size_t>{4, 6}));

    auto res2 = baysor::load_prior_from_image(path, x, y, /*min=*/3);
    EXPECT_EQ(res2.segment_per_molecule, std::vector<int>(4, 0));
}

TEST(Cov1Data_Prior, ImageScaleFallbackWhenNoMoleculesInMask) {
    TempDir dir;
    // Valid image but molecules only sit on background pixels -> no components,
    // so the scale estimate falls back to the (all-unassigned) assignment and
    // fails, leaving the -1 sentinel.
    std::vector<uint8_t> px(8 * 8, 0);
    const auto path = write_tiff_u8(dir, "empty_mask.tif", px, 8, 8);

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
    TempDir dir;
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
    const auto path = write_tiff_u8(dir, "three_blobs.tif", px, 8, 8);

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
