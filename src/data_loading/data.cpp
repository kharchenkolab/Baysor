#include "baysor/data_loading/data.h"
#include "baysor/data_loading/prior_segmentation.h"
#include "baysor/utils/general.h"

#include <arrow/api.h>
#include <arrow/compute/api.h>
#include <arrow/csv/api.h>
#include <arrow/io/api.h>
#include <arrow/result.h>
#include <arrow/status.h>
#include <arrow/util/bit_util.h>
#include <parquet/arrow/reader.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <regex>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace baysor {

// ============================================================================
// Helpers: extract columns from Arrow Table
// ============================================================================

namespace {

/// Get file extension, lowercased
std::string file_extension(const std::string& path) {
    auto dot = path.rfind('.');
    if (dot == std::string::npos) return "";
    std::string ext = path.substr(dot + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext;
}

/// Macro-like helper: check Arrow status/result
#define ARROW_CHECK_OK(expr)                                                   \
    do {                                                                        \
        auto _s = (expr);                                                       \
        if (!_s.ok())                                                           \
            throw std::runtime_error(std::string("Arrow error: ") +            \
                                     _s.ToString());                            \
    } while (0)

template <typename T>
T arrow_unwrap(arrow::Result<T>&& result) {
    if (!result.ok()) {
        throw std::runtime_error(std::string("Arrow error: ") + result.status().ToString());
    }
    return std::move(*result);
}

/// Read a CSV file into an Arrow Table
std::shared_ptr<arrow::Table> read_csv_arrow(const std::string& path) {
    auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path));
    auto read_options = arrow::csv::ReadOptions::Defaults();
    auto parse_options = arrow::csv::ParseOptions::Defaults();
    auto convert_options = arrow::csv::ConvertOptions::Defaults();
    // Let Arrow auto-detect types
    auto reader = arrow_unwrap(
        arrow::csv::TableReader::Make(arrow::io::default_io_context(),
                                       input, read_options, parse_options, convert_options));
    return arrow_unwrap(reader->Read());
}

/// Read a Parquet file into an Arrow Table
std::shared_ptr<arrow::Table> read_parquet_arrow(const std::string& path) {
    auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path));
    parquet::arrow::FileReaderBuilder builder;
    ARROW_CHECK_OK(builder.Open(input));
    auto reader = arrow_unwrap(builder.Build());
    std::shared_ptr<arrow::Table> table;
    ARROW_CHECK_OK(reader->ReadTable(&table));
    return table;
}

/// Read a tabular file (detect format by extension)
std::shared_ptr<arrow::Table> read_table(const std::string& path) {
    std::string ext = file_extension(path);
    if (ext == "csv" || ext == "tsv") {
        return read_csv_arrow(path);
    } else if (ext == "parquet" || ext == "pq") {
        return read_parquet_arrow(path);
    } else {
        throw std::runtime_error("Unsupported file format: ." + ext +
                                 ". Provide a CSV or Parquet file.");
    }
}

/// Check that a column exists in the table
int find_column_index(const std::shared_ptr<arrow::Table>& table, const std::string& name) {
    auto schema = table->schema();
    int idx = schema->GetFieldIndex(name);
    if (idx < 0) {
        throw std::runtime_error("Column '" + name + "' not found in the data. "
                                 "Available columns: " +
                                 schema->ToString());
    }
    return idx;
}

/// Extract a double column from an Arrow table (handles chunked arrays, int/float types)
std::vector<double> extract_double_column(const std::shared_ptr<arrow::Table>& table,
                                           const std::string& name) {
    int idx = find_column_index(table, name);
    auto chunked = table->column(idx);
    int64_t n = chunked->length();
    std::vector<double> result(n);

    int64_t offset = 0;
    for (int c = 0; c < chunked->num_chunks(); ++c) {
        auto chunk = chunked->chunk(c);
        int64_t len = chunk->length();

        // Handle different numeric types
        if (chunk->type_id() == arrow::Type::DOUBLE) {
            auto arr = std::static_pointer_cast<arrow::DoubleArray>(chunk);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = arr->Value(i);
            }
        } else if (chunk->type_id() == arrow::Type::FLOAT) {
            auto arr = std::static_pointer_cast<arrow::FloatArray>(chunk);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = static_cast<double>(arr->Value(i));
            }
        } else if (chunk->type_id() == arrow::Type::INT64) {
            auto arr = std::static_pointer_cast<arrow::Int64Array>(chunk);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = static_cast<double>(arr->Value(i));
            }
        } else if (chunk->type_id() == arrow::Type::INT32) {
            auto arr = std::static_pointer_cast<arrow::Int32Array>(chunk);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = static_cast<double>(arr->Value(i));
            }
        } else if (chunk->type_id() == arrow::Type::INT16) {
            auto arr = std::static_pointer_cast<arrow::Int16Array>(chunk);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = static_cast<double>(arr->Value(i));
            }
        } else {
            // Try casting to double
            auto maybe_casted = arrow::compute::Cast(*chunk, arrow::TypeHolder(arrow::float64()));
            if (!maybe_casted.ok()) {
                throw std::runtime_error("Cannot convert column '" + name +
                                         "' to double: " + maybe_casted.status().ToString());
            }
            auto casted = std::static_pointer_cast<arrow::DoubleArray>(*maybe_casted);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = casted->Value(i);
            }
        }
        offset += len;
    }
    return result;
}

/// Extract a string column from an Arrow table (handles string, large_string, dict types)
std::vector<std::string> extract_string_column(const std::shared_ptr<arrow::Table>& table,
                                                const std::string& name) {
    int idx = find_column_index(table, name);
    auto chunked = table->column(idx);
    int64_t n = chunked->length();
    std::vector<std::string> result(n);

    int64_t offset = 0;
    for (int c = 0; c < chunked->num_chunks(); ++c) {
        auto chunk = chunked->chunk(c);
        int64_t len = chunk->length();

        if (chunk->type_id() == arrow::Type::STRING) {
            auto arr = std::static_pointer_cast<arrow::StringArray>(chunk);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = arr->GetString(i);
            }
        } else if (chunk->type_id() == arrow::Type::LARGE_STRING) {
            auto arr = std::static_pointer_cast<arrow::LargeStringArray>(chunk);
            for (int64_t i = 0; i < len; ++i) {
                result[offset + i] = arr->GetString(i);
            }
        } else if (chunk->type_id() == arrow::Type::DICTIONARY) {
            // Dictionary-encoded strings (common in Parquet)
            auto dict_arr = std::static_pointer_cast<arrow::DictionaryArray>(chunk);
            auto dict_values = std::static_pointer_cast<arrow::StringArray>(dict_arr->dictionary());
            auto indices = dict_arr->indices();

            // Handle different index types
            if (indices->type_id() == arrow::Type::INT32) {
                auto idx_arr = std::static_pointer_cast<arrow::Int32Array>(indices);
                for (int64_t i = 0; i < len; ++i) {
                    if (!dict_arr->IsNull(i)) {
                        result[offset + i] = dict_values->GetString(idx_arr->Value(i));
                    }
                }
            } else if (indices->type_id() == arrow::Type::INT16) {
                auto idx_arr = std::static_pointer_cast<arrow::Int16Array>(indices);
                for (int64_t i = 0; i < len; ++i) {
                    if (!dict_arr->IsNull(i)) {
                        result[offset + i] = dict_values->GetString(idx_arr->Value(i));
                    }
                }
            } else if (indices->type_id() == arrow::Type::INT8) {
                auto idx_arr = std::static_pointer_cast<arrow::Int8Array>(indices);
                for (int64_t i = 0; i < len; ++i) {
                    if (!dict_arr->IsNull(i)) {
                        result[offset + i] = dict_values->GetString(idx_arr->Value(i));
                    }
                }
            } else {
                // Fallback: cast indices
                auto idx_arr = std::static_pointer_cast<arrow::Int64Array>(
                    arrow_unwrap(arrow::compute::Cast(*indices, arrow::TypeHolder(arrow::int64()))));
                for (int64_t i = 0; i < len; ++i) {
                    if (!dict_arr->IsNull(i)) {
                        result[offset + i] = dict_values->GetString(
                            static_cast<int>(idx_arr->Value(i)));
                    }
                }
            }
        } else {
            // Numeric column used as gene — convert to string
            for (int64_t i = 0; i < len; ++i) {
                auto scalar = arrow_unwrap(chunk->GetScalar(i));
                result[offset + i] = scalar->ToString();
            }
        }
        offset += len;
    }
    return result;
}

/// Check if a column exists
bool has_column(const std::shared_ptr<arrow::Table>& table, const std::string& name) {
    return table->schema()->GetFieldIndex(name) >= 0;
}

int find_column_index(const std::shared_ptr<arrow::Schema>& schema, const std::string& name) {
    int idx = schema->GetFieldIndex(name);
    if (idx < 0) {
        throw std::runtime_error("Column '" + name + "' not found in the data. "
                                 "Available columns: " + schema->ToString());
    }
    return idx;
}

bool has_column(const std::shared_ptr<arrow::Schema>& schema, const std::string& name) {
    return schema->GetFieldIndex(name) >= 0;
}

std::shared_ptr<parquet::arrow::FileReader> open_parquet_reader(const std::string& path) {
    auto input = arrow_unwrap(arrow::io::ReadableFile::Open(path));
    parquet::arrow::FileReaderBuilder builder;
    ARROW_CHECK_OK(builder.Open(input));
    auto reader = arrow_unwrap(builder.Build());
    reader->set_use_threads(true);
    reader->set_batch_size(65536);
    return reader;
}

std::shared_ptr<arrow::Schema> get_parquet_schema(
    const std::shared_ptr<parquet::arrow::FileReader>& reader
) {
    std::shared_ptr<arrow::Schema> schema;
    ARROW_CHECK_OK(reader->GetSchema(&schema));
    return schema;
}

struct NumericArrayView {
    enum class Kind {
        Invalid,
        Double,
        Float,
        Int64,
        Int32,
        Int16,
        Int8,
        UInt64,
        UInt32,
        UInt16,
        UInt8,
        Fallback
    };

    const uint8_t* null_bitmap = nullptr;
    int64_t offset = 0;
    Kind kind = Kind::Invalid;
    const arrow::DoubleArray* double_arr = nullptr;
    const arrow::FloatArray* float_arr = nullptr;
    const arrow::Int64Array* int64_arr = nullptr;
    const arrow::Int32Array* int32_arr = nullptr;
    const arrow::Int16Array* int16_arr = nullptr;
    const arrow::Int8Array* int8_arr = nullptr;
    const arrow::UInt64Array* uint64_arr = nullptr;
    const arrow::UInt32Array* uint32_arr = nullptr;
    const arrow::UInt16Array* uint16_arr = nullptr;
    const arrow::UInt8Array* uint8_arr = nullptr;
    std::shared_ptr<arrow::Array> fallback;

    NumericArrayView() = default;
    explicit NumericArrayView(const std::shared_ptr<arrow::Array>& arr) { reset(arr); }

    void reset(const std::shared_ptr<arrow::Array>& arr) {
        null_bitmap = nullptr;
        offset = 0;
        kind = Kind::Invalid;
        double_arr = nullptr;
        float_arr = nullptr;
        int64_arr = nullptr;
        int32_arr = nullptr;
        int16_arr = nullptr;
        int8_arr = nullptr;
        uint64_arr = nullptr;
        uint32_arr = nullptr;
        uint16_arr = nullptr;
        uint8_arr = nullptr;
        fallback.reset();

        if (!arr) return;

        null_bitmap = arr->null_bitmap_data();
        offset = arr->offset();
        switch (arr->type_id()) {
            case arrow::Type::DOUBLE:
                double_arr = static_cast<const arrow::DoubleArray*>(arr.get());
                kind = Kind::Double;
                break;
            case arrow::Type::FLOAT:
                float_arr = static_cast<const arrow::FloatArray*>(arr.get());
                kind = Kind::Float;
                break;
            case arrow::Type::INT64:
                int64_arr = static_cast<const arrow::Int64Array*>(arr.get());
                kind = Kind::Int64;
                break;
            case arrow::Type::INT32:
                int32_arr = static_cast<const arrow::Int32Array*>(arr.get());
                kind = Kind::Int32;
                break;
            case arrow::Type::INT16:
                int16_arr = static_cast<const arrow::Int16Array*>(arr.get());
                kind = Kind::Int16;
                break;
            case arrow::Type::INT8:
                int8_arr = static_cast<const arrow::Int8Array*>(arr.get());
                kind = Kind::Int8;
                break;
            case arrow::Type::UINT64:
                uint64_arr = static_cast<const arrow::UInt64Array*>(arr.get());
                kind = Kind::UInt64;
                break;
            case arrow::Type::UINT32:
                uint32_arr = static_cast<const arrow::UInt32Array*>(arr.get());
                kind = Kind::UInt32;
                break;
            case arrow::Type::UINT16:
                uint16_arr = static_cast<const arrow::UInt16Array*>(arr.get());
                kind = Kind::UInt16;
                break;
            case arrow::Type::UINT8:
                uint8_arr = static_cast<const arrow::UInt8Array*>(arr.get());
                kind = Kind::UInt8;
                break;
            default:
                fallback = arr;
                kind = Kind::Fallback;
                break;
        }
    }

    bool is_valid(int64_t i) const {
        if (kind == Kind::Invalid) return false;
        return !null_bitmap || arrow::bit_util::GetBit(null_bitmap, offset + i);
    }

    double value(int64_t i) const {
        if (!is_valid(i)) return std::numeric_limits<double>::quiet_NaN();
        switch (kind) {
            case Kind::Double: return double_arr->Value(i);
            case Kind::Float: return static_cast<double>(float_arr->Value(i));
            case Kind::Int64: return static_cast<double>(int64_arr->Value(i));
            case Kind::Int32: return static_cast<double>(int32_arr->Value(i));
            case Kind::Int16: return static_cast<double>(int16_arr->Value(i));
            case Kind::Int8: return static_cast<double>(int8_arr->Value(i));
            case Kind::UInt64: return static_cast<double>(uint64_arr->Value(i));
            case Kind::UInt32: return static_cast<double>(uint32_arr->Value(i));
            case Kind::UInt16: return static_cast<double>(uint16_arr->Value(i));
            case Kind::UInt8: return static_cast<double>(uint8_arr->Value(i));
            case Kind::Fallback: {
                auto scalar = arrow_unwrap(fallback->GetScalar(i));
                return std::stod(scalar->ToString());
            }
            case Kind::Invalid: break; // GCOVR_EXCL_LINE: unreachable
        }
        return std::numeric_limits<double>::quiet_NaN(); // GCOVR_EXCL_LINE: unreachable
    }

    int64_t int64_value(int64_t i) const {
        if (!is_valid(i)) return -1;
        switch (kind) {
            case Kind::Int64: return int64_arr->Value(i);
            case Kind::Int32: return static_cast<int64_t>(int32_arr->Value(i));
            case Kind::Int16: return static_cast<int64_t>(int16_arr->Value(i));
            case Kind::Int8: return static_cast<int64_t>(int8_arr->Value(i));
            case Kind::UInt64: return static_cast<int64_t>(uint64_arr->Value(i));
            case Kind::UInt32: return static_cast<int64_t>(uint32_arr->Value(i));
            case Kind::UInt16: return static_cast<int64_t>(uint16_arr->Value(i));
            case Kind::UInt8: return static_cast<int64_t>(uint8_arr->Value(i));
            case Kind::Double: return static_cast<int64_t>(std::llround(double_arr->Value(i)));
            case Kind::Float: return static_cast<int64_t>(std::llround(float_arr->Value(i)));
            case Kind::Fallback: {
                auto scalar = arrow_unwrap(fallback->GetScalar(i));
                return static_cast<int64_t>(std::stoll(scalar->ToString()));
            }
            case Kind::Invalid: break; // GCOVR_EXCL_LINE: unreachable
        }
        return -1; // GCOVR_EXCL_LINE: unreachable
    }
};

struct StringArrayView {
    enum class Kind {
        Invalid,
        String,
        LargeString,
        DictionaryString,
        DictionaryLargeString,
        Fallback
    };

    const uint8_t* null_bitmap = nullptr;
    int64_t offset = 0;
    Kind kind = Kind::Invalid;
    const arrow::StringArray* string_arr = nullptr;
    const arrow::LargeStringArray* large_string_arr = nullptr;
    NumericArrayView dict_indices;
    // Binary base classes, so binary dictionary values decode as text too.
    const arrow::BinaryArray* dict_string_arr = nullptr;
    const arrow::LargeBinaryArray* dict_large_string_arr = nullptr;
    std::shared_ptr<arrow::Array> fallback;

    StringArrayView() = default;
    explicit StringArrayView(const std::shared_ptr<arrow::Array>& arr) { reset(arr); }

    void reset(const std::shared_ptr<arrow::Array>& arr) {
        null_bitmap = nullptr;
        offset = 0;
        kind = Kind::Invalid;
        string_arr = nullptr;
        large_string_arr = nullptr;
        dict_indices = NumericArrayView();
        dict_string_arr = nullptr;
        dict_large_string_arr = nullptr;
        fallback.reset();

        if (!arr) return;

        null_bitmap = arr->null_bitmap_data();
        offset = arr->offset();
        switch (arr->type_id()) {
            case arrow::Type::STRING:
                string_arr = static_cast<const arrow::StringArray*>(arr.get());
                kind = Kind::String;
                break;
            case arrow::Type::LARGE_STRING:
                large_string_arr = static_cast<const arrow::LargeStringArray*>(arr.get());
                kind = Kind::LargeString;
                break;
            case arrow::Type::DICTIONARY: {
                auto dict_arr = static_cast<const arrow::DictionaryArray*>(arr.get());
                dict_indices.reset(dict_arr->indices());
                auto dict = dict_arr->dictionary();
                if (dict->type_id() == arrow::Type::STRING ||
                        dict->type_id() == arrow::Type::BINARY) {
                    dict_string_arr = static_cast<const arrow::BinaryArray*>(dict.get());
                    kind = Kind::DictionaryString;
                } else if (dict->type_id() == arrow::Type::LARGE_STRING || // GCOVR_EXCL_LINE: parquet yields only string/binary dictionaries
                        dict->type_id() == arrow::Type::LARGE_BINARY) { // GCOVR_EXCL_LINE: parquet yields only string/binary dictionaries
                    dict_large_string_arr = static_cast<const arrow::LargeBinaryArray*>(dict.get()); // GCOVR_EXCL_LINE: parquet yields only string/binary dictionaries
                    kind = Kind::DictionaryLargeString; // GCOVR_EXCL_LINE: parquet yields only string/binary dictionaries
                } else {
                    fallback = arr; // GCOVR_EXCL_LINE: parquet yields only string/binary dictionaries
                    kind = Kind::Fallback; // GCOVR_EXCL_LINE: parquet yields only string/binary dictionaries
                }
                break;
            }
            default:
                fallback = arr;
                kind = Kind::Fallback;
                break;
        }
    }

    bool is_valid(int64_t i) const {
        if (kind == Kind::Invalid) return false;
        return !null_bitmap || arrow::bit_util::GetBit(null_bitmap, offset + i);
    }

    bool is_dictionary() const {
        return kind == Kind::DictionaryString || kind == Kind::DictionaryLargeString;
    }

    int64_t dictionary_length() const {
        if (kind == Kind::DictionaryString && dict_string_arr) return dict_string_arr->length();
        if (kind == Kind::DictionaryLargeString && dict_large_string_arr) return dict_large_string_arr->length(); // GCOVR_EXCL_LINE: parquet yields only string/binary dictionaries
        return 0; // GCOVR_EXCL_LINE: unreachable
    }

    int64_t dictionary_index(int64_t i) const {
        if (!is_dictionary() || !is_valid(i)) return -1;
        return dict_indices.int64_value(i);
    }

    std::string dictionary_value(int64_t dict_idx) const {
        if (dict_idx < 0) return "";
        if (kind == Kind::DictionaryString && dict_string_arr) {
            return std::string(dict_string_arr->GetView(dict_idx));
        }
        if (kind == Kind::DictionaryLargeString && dict_large_string_arr) { // GCOVR_EXCL_LINE: parquet yields only string/binary dictionaries
            return std::string(dict_large_string_arr->GetView(dict_idx)); // GCOVR_EXCL_LINE: parquet yields only string/binary dictionaries
        }
        return ""; // GCOVR_EXCL_LINE: unreachable
    }

    std::string value(int64_t i) const {
        if (!is_valid(i)) return "";
        switch (kind) {
            case Kind::String:
                return std::string(string_arr->GetView(i));
            case Kind::LargeString:
                return std::string(large_string_arr->GetView(i));
            case Kind::DictionaryString:
            case Kind::DictionaryLargeString:
                return dictionary_value(dictionary_index(i));
            case Kind::Fallback: {
                auto scalar = arrow_unwrap(fallback->GetScalar(i));
                return scalar->ToString();
            }
            case Kind::Invalid: // GCOVR_EXCL_LINE: unreachable
                break; // GCOVR_EXCL_LINE: unreachable
        }
        return ""; // GCOVR_EXCL_LINE: unreachable
    }
};

std::vector<std::regex> compile_gene_patterns(const std::vector<std::string>& patterns) {
    std::vector<std::regex> regexes;
    regexes.reserve(patterns.size());
    for (const auto& pat : patterns) {
        std::string re_str;
        for (char c : pat) {
            switch (c) {
                case '*': re_str += ".*"; break;
                case '?': re_str += "."; break;
                case '.': re_str += "\\."; break;
                case '+': re_str += "\\+"; break;
                case '(': re_str += "\\("; break;
                case ')': re_str += "\\)"; break;
                case '[': re_str += "\\["; break;
                case ']': re_str += "\\]"; break;
                case '{': re_str += "\\{"; break;
                case '}': re_str += "\\}"; break;
                case '^': re_str += "\\^"; break;
                case '$': re_str += "\\$"; break;
                case '|': re_str += "\\|"; break;
                case '\\': re_str += "\\\\"; break;
                default: re_str += c; break;
            }
        }
        regexes.emplace_back(re_str, std::regex::ECMAScript);
    }
    return regexes;
}

bool matches_any_gene_pattern(const std::string& gene, const std::vector<std::regex>& regexes) {
    for (const auto& re : regexes) {
        if (std::regex_match(gene, re)) return true;
    }
    return false;
}

struct Pass1Summary {
    std::vector<std::string> gene_names;
    std::unordered_map<std::string, int> gene_id_map;
    int64_t kept_rows = 0;
};

struct ParquetScanPlan {
    int x = -1;
    int y = -1;
    int z = -1;
    int gene = -1;
    int qv = -1;
    int transcript_id = -1;
    int prior = -1;
    int confidence = -1;
    int cluster = -1;
    int nuclei_probs = -1;
};

bool row_passes_static_filters(
    double x,
    double y,
    bool has_z,
    double z,
    bool has_qv,
    double qv,
    const MoleculeInputOptions& opts
) {
    if (!std::isfinite(x) || !std::isfinite(y)) return false;
    if (x < opts.x_min || x > opts.x_max) return false;
    if (y < opts.y_min || y > opts.y_max) return false;
    if (has_z) {
        if (!std::isfinite(z)) return false;
        if (z < opts.z_min || z > opts.z_max) return false;
    }
    if (has_qv && opts.min_qv >= 0.0) {
        if (!std::isfinite(qv) || qv < opts.min_qv) return false;
    }
    return true;
}

ParquetScanPlan make_parquet_scan_plan(
    const std::shared_ptr<arrow::Schema>& schema,
    const MoleculeInputOptions& opts,
    const std::string& prior_column_name
) {
    ParquetScanPlan plan;
    plan.x = find_column_index(schema, opts.x_col);
    plan.y = find_column_index(schema, opts.y_col);
    plan.gene = find_column_index(schema, opts.gene_col);
    if (!opts.force_2d && has_column(schema, opts.z_col)) {
        plan.z = find_column_index(schema, opts.z_col);
    }
    if (has_column(schema, opts.qv_col)) {
        plan.qv = find_column_index(schema, opts.qv_col);
    }
    if (has_column(schema, "transcript_id")) plan.transcript_id = find_column_index(schema, "transcript_id");
    if (!prior_column_name.empty()) {
        plan.prior = find_column_index(schema, prior_column_name);
    }
    if (has_column(schema, "confidence")) plan.confidence = find_column_index(schema, "confidence");
    if (has_column(schema, "cluster")) plan.cluster = find_column_index(schema, "cluster");
    if (has_column(schema, "nuclei_probs")) {
        plan.nuclei_probs = find_column_index(schema, "nuclei_probs");
    }
    return plan;
}

std::vector<int> all_row_groups(const std::shared_ptr<parquet::arrow::FileReader>& reader) {
    std::vector<int> groups(reader->num_row_groups());
    std::iota(groups.begin(), groups.end(), 0);
    return groups;
}

Pass1Summary scan_parquet_pass1(
    const std::string& path,
    const MoleculeInputOptions& opts
) {
    auto reader = open_parquet_reader(path);
    auto schema = get_parquet_schema(reader);
    auto plan = make_parquet_scan_plan(schema, opts, "");

    std::vector<int> projected = {plan.x, plan.y};
    if (plan.z >= 0) projected.push_back(plan.z);
    projected.push_back(plan.gene);
    if (plan.qv >= 0) projected.push_back(plan.qv);

    auto batch_reader = arrow_unwrap(reader->GetRecordBatchReader(all_row_groups(reader), projected));

    std::unordered_map<std::string, int64_t> gene_counts;
    std::shared_ptr<arrow::RecordBatch> batch;
    while (true) {
        ARROW_CHECK_OK(batch_reader->ReadNext(&batch));
        if (!batch) break;

        int col = 0;
        auto x_arr = batch->column(col++);
        auto y_arr = batch->column(col++);
        std::shared_ptr<arrow::Array> z_arr;
        if (plan.z >= 0) z_arr = batch->column(col++);
        auto gene_arr = batch->column(col++);
        std::shared_ptr<arrow::Array> qv_arr;
        if (plan.qv >= 0) qv_arr = batch->column(col++);

        const bool has_z = plan.z >= 0;
        const bool has_qv = plan.qv >= 0;
        NumericArrayView x_view(x_arr);
        NumericArrayView y_view(y_arr);
        NumericArrayView z_view(z_arr);
        NumericArrayView qv_view(qv_arr);
        StringArrayView gene_view(gene_arr);

        if (gene_view.is_dictionary()) {
            std::vector<int64_t> batch_gene_counts(
                static_cast<size_t>(gene_view.dictionary_length()), 0);
            for (int64_t i = 0; i < batch->num_rows(); ++i) {
                double x = x_view.value(i);
                double y = y_view.value(i);
                double z = has_z ? z_view.value(i) : 0.0;
                double qv = has_qv ? qv_view.value(i) : 0.0;
                if (!row_passes_static_filters(x, y, has_z, z, has_qv, qv, opts)) continue;

                int64_t dict_idx = gene_view.dictionary_index(i);
                if (dict_idx >= 0 && dict_idx < static_cast<int64_t>(batch_gene_counts.size())) {
                    batch_gene_counts[static_cast<size_t>(dict_idx)]++;
                }
            }

            for (size_t dict_idx = 0; dict_idx < batch_gene_counts.size(); ++dict_idx) {
                int64_t count = batch_gene_counts[dict_idx];
                if (count == 0) continue;
                gene_counts[gene_view.dictionary_value(static_cast<int64_t>(dict_idx))] += count;
            }
        } else {
            for (int64_t i = 0; i < batch->num_rows(); ++i) {
                double x = x_view.value(i);
                double y = y_view.value(i);
                double z = has_z ? z_view.value(i) : 0.0;
                double qv = has_qv ? qv_view.value(i) : 0.0;
                if (!row_passes_static_filters(x, y, has_z, z, has_qv, qv, opts)) continue;

                auto gene = gene_view.value(i);
                if (gene.empty()) continue;
                gene_counts[gene]++;
            }
        }
    }

    auto patterns = compile_gene_patterns(split_string_list(opts.exclude_genes));
    Pass1Summary summary;
    for (const auto& [gene, count] : gene_counts) {
        if (count < opts.min_molecules_per_gene) continue;
        if (matches_any_gene_pattern(gene, patterns)) continue;
        summary.gene_names.push_back(gene);
        summary.kept_rows += count;
    }
    std::sort(summary.gene_names.begin(), summary.gene_names.end());
    summary.gene_names.erase(std::unique(summary.gene_names.begin(), summary.gene_names.end()),
                             summary.gene_names.end());
    summary.gene_id_map.reserve(summary.gene_names.size());
    for (int i = 0; i < static_cast<int>(summary.gene_names.size()); ++i) {
        summary.gene_id_map[summary.gene_names[i]] = i + 1;
    }
    return summary;
}

} // anonymous namespace

// ============================================================================
// MoleculeData methods
// ============================================================================

Eigen::MatrixXd MoleculeData::position_matrix() const {
    int n = n_molecules();
    int d = n_dims();
    Eigen::MatrixXd mat(d, n);
    for (int i = 0; i < n; ++i) {
        mat(0, i) = x[i];
        mat(1, i) = y[i];
        if (d == 3) mat(2, i) = z[i];
    }
    return mat;
}

// ============================================================================
// Gene encoding
// ============================================================================

void encode_genes(MoleculeData& data, const std::vector<std::string>& gene_strings) {
    // Collect unique names and sort alphabetically
    std::set<std::string> unique_set;
    for (const auto& g : gene_strings) {
        if (!g.empty()) unique_set.insert(g);
    }

    data.gene_names.assign(unique_set.begin(), unique_set.end());

    // Build name -> 1-based ID map
    std::unordered_map<std::string, int> gene_id_map;
    gene_id_map.reserve(data.gene_names.size());
    for (int i = 0; i < static_cast<int>(data.gene_names.size()); ++i) {
        gene_id_map[data.gene_names[i]] = i + 1; // 1-based
    }

    // Encode
    data.gene.resize(gene_strings.size());
    for (size_t i = 0; i < gene_strings.size(); ++i) {
        auto it = gene_id_map.find(gene_strings[i]);
        data.gene[i] = (it != gene_id_map.end()) ? it->second : 0;
    }
}

// ============================================================================
// Gene filtering
// ============================================================================

void filter_genes_by_count(MoleculeData& data, int min_molecules_per_gene) {
    if (min_molecules_per_gene <= 1) return; // nothing to filter

    int n_genes = data.n_genes();
    int n = data.n_molecules();

    // Count molecules per gene
    std::vector<int> gene_counts(n_genes, 0);
    for (int i = 0; i < n; ++i) {
        int g = data.gene[i];
        if (g > 0 && g <= n_genes) gene_counts[g - 1]++;
    }

    // Find genes passing the threshold
    std::unordered_set<int> passing_genes;
    for (int g = 0; g < n_genes; ++g) {
        if (gene_counts[g] >= min_molecules_per_gene) {
            passing_genes.insert(g + 1); // 1-based
        }
    }

    if (static_cast<int>(passing_genes.size()) == n_genes) return; // all pass

    // Filter molecules: keep only those with passing genes
    std::vector<int> keep_indices;
    keep_indices.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (passing_genes.count(data.gene[i])) {
            keep_indices.push_back(i);
        }
    }

    // Compact data
    auto compact = [&](auto& vec) {
        if (vec.empty()) return;
        using T = typename std::decay<decltype(vec)>::type;
        T new_vec;
        new_vec.reserve(keep_indices.size());
        for (int idx : keep_indices) {
            new_vec.push_back(std::move(vec[idx]));
        }
        vec = std::move(new_vec);
    };

    compact(data.x);
    compact(data.y);
    compact(data.z);
    compact(data.confidence);
    compact(data.cluster);
    compact(data.prior_segmentation);
    compact(data.nuclei_probs);
    compact(data.source_transcript_id);

    // Re-encode genes: collect raw strings for kept molecules, then re-encode
    std::vector<std::string> kept_gene_strings;
    kept_gene_strings.reserve(keep_indices.size());
    for (int idx : keep_indices) {
        int g = data.gene[idx];
        kept_gene_strings.push_back(data.gene_names[g - 1]);
    }

    // Clear old gene data and re-encode
    data.gene.clear();
    data.gene_names.clear();
    encode_genes(data, kept_gene_strings);
}

void filter_genes_by_pattern(MoleculeData& data, const std::vector<std::string>& patterns) {
    if (patterns.empty()) return;

    int n_genes = data.n_genes();
    int n = data.n_molecules();

    // Convert wildcard patterns to regexes
    // Julia uses Regex() directly, where users write "Blank*" meaning "Blank.*"
    // We support shell-glob style: * -> .*, ? -> .
    std::vector<std::regex> regexes;
    regexes.reserve(patterns.size());
    for (const auto& pat : patterns) {
        // Escape regex special chars except * and ?
        std::string re_str;
        for (char c : pat) {
            switch (c) {
                case '*': re_str += ".*"; break;
                case '?': re_str += "."; break;
                case '.': re_str += "\\."; break;
                case '+': re_str += "\\+"; break;
                case '(': re_str += "\\("; break;
                case ')': re_str += "\\)"; break;
                case '[': re_str += "\\["; break;
                case ']': re_str += "\\]"; break;
                case '{': re_str += "\\{"; break;
                case '}': re_str += "\\}"; break;
                case '^': re_str += "\\^"; break;
                case '$': re_str += "\\$"; break;
                case '|': re_str += "\\|"; break;
                case '\\': re_str += "\\\\"; break;
                default: re_str += c; break;
            }
        }
        regexes.emplace_back(re_str, std::regex::ECMAScript);
    }

    // Find gene names matching any pattern
    std::unordered_set<int> excluded_genes; // 1-based
    std::vector<std::string> excluded_names; // GCOVR_EXCL_LINE: gcov artifact
    for (int g = 0; g < n_genes; ++g) {
        const auto& name = data.gene_names[g];
        for (const auto& re : regexes) {
            if (std::regex_match(name, re)) {
                excluded_genes.insert(g + 1);
                excluded_names.push_back(name);
                break;
            }
        }
    }

    if (excluded_genes.empty()) return;

    // Filter molecules: keep those NOT matching excluded genes
    std::vector<int> keep_indices;
    keep_indices.reserve(n);
    for (int i = 0; i < n; ++i) {
        if (!excluded_genes.count(data.gene[i])) {
            keep_indices.push_back(i);
        }
    }

    // Compact (same pattern as filter_genes_by_count)
    auto compact = [&](auto& vec) {
        if (vec.empty()) return;
        using T = typename std::decay<decltype(vec)>::type;
        T new_vec;
        new_vec.reserve(keep_indices.size());
        for (int idx : keep_indices) {
            new_vec.push_back(std::move(vec[idx]));
        }
        vec = std::move(new_vec);
    };

    compact(data.x);
    compact(data.y);
    compact(data.z);
    compact(data.confidence);
    compact(data.cluster);
    compact(data.prior_segmentation);
    compact(data.nuclei_probs);
    compact(data.source_transcript_id);

    // Re-encode
    std::vector<std::string> kept_gene_strings;
    kept_gene_strings.reserve(keep_indices.size());
    for (int idx : keep_indices) {
        int g = data.gene[idx];
        kept_gene_strings.push_back(data.gene_names[g - 1]);
    }

    data.gene.clear();
    data.gene_names.clear();
    encode_genes(data, kept_gene_strings);
}

// ============================================================================
// Read a single string column (for prior segmentation column references)
// ============================================================================

std::vector<std::string> read_string_column(const std::string& path, const std::string& col_name) {
    auto table = read_table(path);
    return extract_string_column(table, col_name);
}

std::vector<double> read_double_column(const std::string& path, const std::string& col_name) {
    auto table = read_table(path);
    return extract_double_column(table, col_name);
}

// ============================================================================
// Main loading function
// ============================================================================

RawTableData read_tabular_file(const std::string& path, const MoleculeInputOptions& opts) {
    auto table = read_table(path);
    RawTableData raw;

    // Extract x, y (required)
    raw.x = extract_double_column(table, opts.x_col);
    raw.y = extract_double_column(table, opts.y_col);

    // Extract z (optional)
    if (!opts.force_2d && has_column(table, opts.z_col)) {
        raw.z = extract_double_column(table, opts.z_col);
        raw.has_z = true;

        // Drop z if only 1 unique value
        if (!raw.z.empty()) {
            double z0 = raw.z[0];
            bool all_same = true;
            for (size_t i = 1; i < raw.z.size(); ++i) {
                if (raw.z[i] != z0) { all_same = false; break; }
            }
            if (all_same) {
                raw.z.clear();
                raw.has_z = false;
            }
        }
    }

    // Extract gene column as strings
    raw.gene_str = extract_string_column(table, opts.gene_col);

    // Stringify any non-string gene values (handles numeric gene columns)
    // Already handled in extract_string_column via scalar->ToString()

    return raw;
}

MoleculeData load_molecules(
    const std::string& path,
    const MoleculeInputOptions& opts,
    const PriorInputOptions& prior_opts
) {
    MoleculeData data;
    const bool load_prior_column = prior_opts.type == PriorInputType::Column && !prior_opts.column_name.empty();
    const std::string& prior_column_name = prior_opts.column_name;

    std::string ext = file_extension(path);
    if (ext == "parquet" || ext == "pq") {
        auto pass1 = scan_parquet_pass1(path, opts);
        data.gene_names = pass1.gene_names;

        auto reader = open_parquet_reader(path);
        auto schema = get_parquet_schema(reader);
        auto plan = make_parquet_scan_plan(schema, opts, load_prior_column ? prior_column_name : "");

        std::vector<int> projected = {plan.x, plan.y};
        if (plan.z >= 0) projected.push_back(plan.z);
        projected.push_back(plan.gene);
        if (plan.qv >= 0) projected.push_back(plan.qv);
        if (plan.transcript_id >= 0) projected.push_back(plan.transcript_id);
        if (plan.prior >= 0) projected.push_back(plan.prior);
        if (plan.confidence >= 0) projected.push_back(plan.confidence);
        if (plan.cluster >= 0) projected.push_back(plan.cluster);
        if (plan.nuclei_probs >= 0) projected.push_back(plan.nuclei_probs);

        auto batch_reader = arrow_unwrap(reader->GetRecordBatchReader(all_row_groups(reader), projected));

        data.x.resize(pass1.kept_rows);
        data.y.resize(pass1.kept_rows);
        data.gene.resize(pass1.kept_rows);
        if (plan.z >= 0) data.z.resize(pass1.kept_rows);
        if (plan.transcript_id >= 0) data.source_transcript_id.resize(pass1.kept_rows);
        if (plan.confidence >= 0) data.confidence.resize(pass1.kept_rows);
        if (plan.cluster >= 0) data.cluster.resize(pass1.kept_rows);
        if (plan.nuclei_probs >= 0) data.nuclei_probs.resize(pass1.kept_rows);

        std::vector<std::string> prior_raw;
        if (plan.prior >= 0) prior_raw.resize(pass1.kept_rows);

        int64_t out_idx = 0;
        std::shared_ptr<arrow::RecordBatch> batch;
        while (true) {
            ARROW_CHECK_OK(batch_reader->ReadNext(&batch));
            if (!batch) break;

            int col = 0;
            auto x_arr = batch->column(col++);
            auto y_arr = batch->column(col++);
            std::shared_ptr<arrow::Array> z_arr;
            if (plan.z >= 0) z_arr = batch->column(col++);
            auto gene_arr = batch->column(col++);
            std::shared_ptr<arrow::Array> qv_arr;
            if (plan.qv >= 0) qv_arr = batch->column(col++);
            std::shared_ptr<arrow::Array> transcript_id_arr;
            if (plan.transcript_id >= 0) transcript_id_arr = batch->column(col++);
            std::shared_ptr<arrow::Array> prior_arr;
            if (plan.prior >= 0) prior_arr = batch->column(col++);
            std::shared_ptr<arrow::Array> confidence_arr;
            if (plan.confidence >= 0) confidence_arr = batch->column(col++);
            std::shared_ptr<arrow::Array> cluster_arr;
            if (plan.cluster >= 0) cluster_arr = batch->column(col++);
            std::shared_ptr<arrow::Array> nuclei_probs_arr;
            if (plan.nuclei_probs >= 0) nuclei_probs_arr = batch->column(col++);

            const bool has_z = plan.z >= 0;
            const bool has_qv = plan.qv >= 0;
            NumericArrayView x_view(x_arr);
            NumericArrayView y_view(y_arr);
            NumericArrayView z_view(z_arr);
            NumericArrayView qv_view(qv_arr);
            NumericArrayView transcript_id_view(transcript_id_arr);
            NumericArrayView confidence_view(confidence_arr);
            NumericArrayView cluster_view(cluster_arr);
            NumericArrayView nuclei_probs_view(nuclei_probs_arr);
            StringArrayView gene_view(gene_arr);
            StringArrayView prior_view(prior_arr);

            std::vector<int> dict_to_gene_id;
            if (gene_view.is_dictionary()) {
                dict_to_gene_id.assign(static_cast<size_t>(gene_view.dictionary_length()), 0);
                for (size_t dict_idx = 0; dict_idx < dict_to_gene_id.size(); ++dict_idx) {
                    auto it = pass1.gene_id_map.find(
                        gene_view.dictionary_value(static_cast<int64_t>(dict_idx)));
                    if (it != pass1.gene_id_map.end()) dict_to_gene_id[dict_idx] = it->second;
                }
            }

            for (int64_t i = 0; i < batch->num_rows(); ++i) {
                double x = x_view.value(i);
                double y = y_view.value(i);
                double z = has_z ? z_view.value(i) : 0.0;
                double qv = has_qv ? qv_view.value(i) : 0.0;
                if (!row_passes_static_filters(x, y, has_z, z, has_qv, qv, opts)) continue;

                int gene_id = 0;
                if (gene_view.is_dictionary()) {
                    int64_t dict_idx = gene_view.dictionary_index(i);
                    if (dict_idx >= 0 && dict_idx < static_cast<int64_t>(dict_to_gene_id.size())) {
                        gene_id = dict_to_gene_id[static_cast<size_t>(dict_idx)];
                    }
                } else {
                    auto it = pass1.gene_id_map.find(gene_view.value(i));
                    if (it != pass1.gene_id_map.end()) gene_id = it->second;
                }
                if (gene_id == 0) continue;

                data.x[out_idx] = x;
                data.y[out_idx] = y;
                if (has_z) data.z[out_idx] = z;
                data.gene[out_idx] = gene_id;
                if (plan.transcript_id >= 0) {
                    data.source_transcript_id[out_idx] =
                        static_cast<std::uint64_t>(transcript_id_view.int64_value(i));
                }
                if (plan.prior >= 0) prior_raw[out_idx] = prior_view.value(i);
                if (plan.confidence >= 0) {
                    data.confidence[out_idx] = confidence_view.value(i);
                }
                if (plan.cluster >= 0) {
                    data.cluster[out_idx] =
                        static_cast<int>(std::llround(cluster_view.value(i)));
                }
                if (plan.nuclei_probs >= 0) {
                    data.nuclei_probs[out_idx] = nuclei_probs_view.value(i);
                }
                ++out_idx;
            }
        }

        data.x.resize(out_idx);
        data.y.resize(out_idx);
        data.gene.resize(out_idx);
        if (!data.z.empty()) data.z.resize(out_idx);
        if (!data.source_transcript_id.empty()) data.source_transcript_id.resize(out_idx);
        if (!data.confidence.empty()) data.confidence.resize(out_idx);
        if (!data.cluster.empty()) data.cluster.resize(out_idx);
        if (!data.nuclei_probs.empty()) data.nuclei_probs.resize(out_idx);

        if (!data.z.empty()) {
            bool all_same = true;
            double z0 = data.z.front();
            for (size_t i = 1; i < data.z.size(); ++i) {
                if (data.z[i] != z0) { all_same = false; break; }
            }
            if (all_same) data.z.clear();
        }

        if (load_prior_column && !prior_raw.empty()) {
            prior_raw.resize(out_idx);
            data.prior_segmentation = encode_prior_labels(
                prior_raw, prior_opts.unassigned_label, prior_opts.min_molecules_per_segment);
        }

        return data;
    }

    // Read the full table once (shared with optional column extraction below)
    auto table = read_table(path);

    // Extract required spatial columns
    data.x = extract_double_column(table, opts.x_col);
    data.y = extract_double_column(table, opts.y_col);

    // Extract z (optional)
    if (!opts.force_2d && has_column(table, opts.z_col)) {
        data.z = extract_double_column(table, opts.z_col);
        // Drop z if only 1 unique value
        if (!data.z.empty()) {
            double z0 = data.z[0];
            bool all_same = true;
            for (size_t i = 1; i < data.z.size(); ++i) {
                if (data.z[i] != z0) { all_same = false; break; }
            }
            if (all_same) {
                data.z.clear();
            }
        }
    }

    // Extract gene column
    auto gene_strings = extract_string_column(table, opts.gene_col);

    // Load optional metadata columns (if present in the input file)
    if (has_column(table, "confidence")) {
        data.confidence = extract_double_column(table, "confidence");
    }
    if (has_column(table, "cluster")) {
        // cluster column may be int or string; read as double and round
        auto cluster_dbl = extract_double_column(table, "cluster");
        data.cluster.resize(cluster_dbl.size());
        for (size_t i = 0; i < cluster_dbl.size(); ++i) {
            data.cluster[i] = static_cast<int>(std::round(cluster_dbl[i]));
        }
    }
    if (has_column(table, "nuclei_probs")) {
        data.nuclei_probs = extract_double_column(table, "nuclei_probs");
    }
    std::vector<double> qv_values;
    if (has_column(table, opts.qv_col)) {
        qv_values = extract_double_column(table, opts.qv_col);
    }
    if (has_column(table, "transcript_id")) {
        auto transcript_ids_d = extract_double_column(table, "transcript_id");
        data.source_transcript_id.resize(transcript_ids_d.size());
        for (size_t i = 0; i < transcript_ids_d.size(); ++i) {
            data.source_transcript_id[i] = static_cast<std::uint64_t>(std::llround(transcript_ids_d[i]));
        }
    }
    std::vector<std::string> prior_raw;
    if (load_prior_column) {
        prior_raw = extract_string_column(table, prior_column_name);
    }

    // Apply static row filters before gene encoding / gene-count filtering.
    std::vector<int> keep_indices;
    keep_indices.reserve(data.x.size());
    bool has_z = !data.z.empty();
    bool has_qv = !qv_values.empty();
    for (int i = 0; i < data.n_molecules(); ++i) {
        double z = has_z ? data.z[i] : 0.0;
        double qv = has_qv ? qv_values[i] : 0.0;
        if (row_passes_static_filters(data.x[i], data.y[i], has_z, z, has_qv, qv, opts)) {
            keep_indices.push_back(i);
        }
    }

    auto compact = [&](auto& vec) {
        if (vec.empty()) return;
        using T = typename std::decay<decltype(vec)>::type;
        T new_vec;
        new_vec.reserve(keep_indices.size());
        for (int idx : keep_indices) new_vec.push_back(std::move(vec[idx]));
        vec = std::move(new_vec);
    };

    compact(data.x);
    compact(data.y);
    compact(data.z);
    compact(data.confidence);
    compact(data.cluster);
    compact(data.nuclei_probs);
    compact(data.source_transcript_id);
    compact(gene_strings);
    compact(prior_raw);

    // Drop z if only 1 unique value after filtering
    if (!data.z.empty()) {
        double z0 = data.z.front();
        bool all_same = true;
        for (size_t i = 1; i < data.z.size(); ++i) {
            if (data.z[i] != z0) { all_same = false; break; }
        }
        if (all_same) data.z.clear();
    }

    // Encode genes on the filtered set
    encode_genes(data, gene_strings);

    if (load_prior_column && !prior_raw.empty()) {
        data.prior_segmentation = encode_prior_labels(
            prior_raw, prior_opts.unassigned_label, prior_opts.min_molecules_per_segment);
    }

    // Filter genes with too few molecules
    if (opts.min_molecules_per_gene > 1) {
        filter_genes_by_count(data, opts.min_molecules_per_gene);
    }

    // Filter genes matching exclude patterns
    if (!opts.exclude_genes.empty()) {
        auto patterns = split_string_list(opts.exclude_genes);
        if (!patterns.empty()) {
            filter_genes_by_pattern(data, patterns);
        }
    }

    return data;
}

} // namespace baysor
