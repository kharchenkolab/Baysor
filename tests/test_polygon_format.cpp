// Tests for the `--polygon-format` handling (kharchenkolab/Baysor#165):
// case-insensitive parsing, rejection of unknown values, and the
// `GeometryCollectionLegacy` layout (integer cell ids) used by Xenium
// Ranger 3.x. The cell sets of segmentation.csv and the polygons JSON must
// stay identical for every format.

#include <gtest/gtest.h>

#include <nlohmann/json.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "baysor/data_loading/data.h"
#include "baysor/processing/data_processing/boundary_estimation.h"
#include "baysor/reporting/output.h"

#include "test_cov_helpers.h"

namespace {

using nlohmann::json;

std::string read_text(const std::string& path) {
    std::ifstream f(path);
    EXPECT_TRUE(f.good()) << "cannot read " << path;
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

std::vector<std::string> split(const std::string& line, char sep) {
    std::vector<std::string> parts;
    std::stringstream ss(line);
    std::string part;
    while (std::getline(ss, part, sep)) parts.push_back(part);
    return parts;
}

std::set<std::string> csv_assigned_cells(const std::string& path) {
    std::ifstream f(path);
    EXPECT_TRUE(f.good()) << "cannot read " << path;
    std::string header;
    EXPECT_TRUE(std::getline(f, header));
    const auto cols = split(header, ',');
    int cell_col = -1;
    int noise_col = -1;
    for (int i = 0; i < static_cast<int>(cols.size()); ++i) {
        if (cols[i] == "cell") cell_col = i;
        if (cols[i] == "is_noise") noise_col = i;
    }
    EXPECT_GE(cell_col, 0) << header;
    EXPECT_GE(noise_col, 0) << header;

    std::set<std::string> cells;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        const auto fields = split(line, ',');
        EXPECT_GT(static_cast<int>(fields.size()), std::max(cell_col, noise_col));
        const std::string& cell = fields[cell_col];
        const std::string& noise = fields[noise_col];
        if (cell.empty() || cell == "0" || noise == "true" || noise == "1") continue;
        cells.insert("cell_" + cell.substr(cell.rfind("cell_") == 0 ? 5 : 0));
    }
    return cells;
}

/// Normalise a polygon JSON cell id to the `cell_<n>` form the CSV uses.
std::string normalize_polygon_cell(const json& cell) {
    if (cell.is_number_integer()) return "cell_" + std::to_string(cell.get<int>());
    return cell.get<std::string>();
}

std::set<std::string> json_polygon_cells(const std::string& path) {
    const json doc = json::parse(read_text(path));
    std::set<std::string> cells;
    if (doc.at("type") == "FeatureCollection") {
        for (const auto& feat : doc.at("features")) {
            cells.insert(normalize_polygon_cell(feat.at("properties").at("cell")));
        }
    } else {
        for (const auto& geom : doc.at("geometries")) {
            cells.insert(normalize_polygon_cell(geom.at("cell")));
        }
    }
    return cells;
}

}  // namespace

// ============================================================================
// parse_polygon_format / to_string
// ============================================================================

TEST(Cov165PolygonFormat, ParsingIsCaseInsensitiveAndValidated) {
    EXPECT_EQ(baysor::parse_polygon_format("FeatureCollection"),
              baysor::PolygonFormat::FeatureCollection);
    EXPECT_EQ(baysor::parse_polygon_format("featurecollection"),
              baysor::PolygonFormat::FeatureCollection);
    EXPECT_EQ(baysor::parse_polygon_format("GeometryCollection"),
              baysor::PolygonFormat::GeometryCollection);
    EXPECT_EQ(baysor::parse_polygon_format("geometrycollectionlegacy"),
              baysor::PolygonFormat::GeometryCollectionLegacy);
    EXPECT_EQ(baysor::parse_polygon_format("GeometryCollectionLegacy"),
              baysor::PolygonFormat::GeometryCollectionLegacy);
    EXPECT_EQ(baysor::parse_polygon_format("NONE"), baysor::PolygonFormat::None);

    EXPECT_EQ(baysor::to_string(baysor::PolygonFormat::FeatureCollection),
              "FeatureCollection");
    EXPECT_EQ(baysor::to_string(baysor::PolygonFormat::GeometryCollection),
              "GeometryCollection");
    EXPECT_EQ(baysor::to_string(baysor::PolygonFormat::GeometryCollectionLegacy),
              "GeometryCollectionLegacy");
    EXPECT_EQ(baysor::to_string(baysor::PolygonFormat::None), "none");

    EXPECT_THROW(baysor::parse_polygon_format("geometry"), std::invalid_argument);
    EXPECT_THROW(baysor::parse_polygon_format(""), std::invalid_argument);
}

// ============================================================================
// GeometryCollectionLegacy layout and consistency
// ============================================================================

TEST(Cov165PolygonFormat, LegacyWritesIntegerCellIdsAndKeepsCellSetConsistent) {
    baysor_test::TempDir tmp("cov165_format");
    baysor::MoleculeData data;
    data.x = {0.0, 10.0, 10.0, 0.0,
              100.0, 105.0, 110.0,
              200.0, 210.0};
    data.y = {0.0, 0.0, 10.0, 10.0,
              0.0, 0.0, 0.0,
              0.0, 0.0};
    data.gene.assign(data.x.size(), 1);
    data.gene_names = {"G1"};
    const std::vector<int> assignment = {1, 1, 1, 1, 2, 2, 2, 3, 3};
    const std::vector<std::string> cell_names = {"cell_1", "cell_2", "cell_3"};

    const std::string csv_path = tmp.file("segmentation.csv");
    baysor::save_segmented_df(data, assignment, data.gene_names, csv_path);
    const auto csv_cells = csv_assigned_cells(csv_path);

    auto [joined, stack] = baysor::boundary_polygons_auto(
        data.position_matrix(), assignment, /*estimate_per_z=*/false,
        &cell_names, /*verbose=*/false);
    (void)stack;

    for (const std::string& format :
         {"FeatureCollection", "GeometryCollection", "GeometryCollectionLegacy"}) {
        const std::string json_path = tmp.file("poly_" + format + ".json");
        baysor::save_polygons_geojson(joined, json_path, format);
        const auto poly_cells = json_polygon_cells(json_path);
        EXPECT_EQ(poly_cells, csv_cells) << "format=" << format;
    }

    // The legacy layout is a GeometryCollection with integer ids.
    const json legacy_doc = json::parse(read_text(tmp.file("poly_GeometryCollectionLegacy.json")));
    EXPECT_EQ(legacy_doc.at("type"), "GeometryCollection");
    ASSERT_EQ(legacy_doc.at("geometries").size(), 3u);
    for (const auto& geom : legacy_doc.at("geometries")) {
        EXPECT_TRUE(geom.at("cell").is_number_integer()) << geom.dump();
        EXPECT_GE(geom.at("coordinates").at(0).size(), 4u) << geom.dump();
    }
}

TEST(Cov165PolygonFormat, LegacyRejectsNonIntegerCellIds) {
    baysor::PolygonCollection coll;
    Eigen::MatrixXd tri(2, 3);
    tri << 0.0, 1.0, 0.0,
           0.0, 0.0, 1.0;
    coll["not_a_number"] = tri;

    baysor_test::TempDir tmp("cov165_legacy_bad");
    EXPECT_THROW(baysor::save_polygons_geojson(coll, tmp.file("bad.json"),
                                               "GeometryCollectionLegacy"),
                 std::runtime_error);
}

TEST(Cov165PolygonFormat, NoneFormatWritesNothingCaseInsensitively) {
    baysor::PolygonCollection coll;
    Eigen::MatrixXd tri(2, 3);
    tri << 0.0, 1.0, 0.0,
           0.0, 0.0, 1.0;
    coll["cell_1"] = tri;

    baysor_test::TempDir tmp("cov165_none");
    const std::string out = tmp.file("nothing.json");
    baysor::save_polygons_geojson(coll, out, "NONE");
    EXPECT_FALSE(std::ifstream(out).good());
}
