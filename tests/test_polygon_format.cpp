// Tests for the `--polygon-format` handling (kharchenkolab/Baysor#165):
// case-insensitive parsing, rejection of unknown values, and the
// `GeometryCollectionLegacy` layout (integer cell ids) used by Xenium
// Ranger 3.x. The cell sets of segmentation.csv and the polygons JSON must
// stay identical for every format.

#include <gtest/gtest.h>

#include <nlohmann/json.hpp>

#include <Eigen/Dense>

#include <set>
#include <string>
#include <vector>

#include "baysor/data_loading/data.h"
#include "baysor/processing/data_processing/boundary_estimation.h"
#include "baysor/reporting/output.h"

#include "test_cov_helpers.h"

namespace {

using nlohmann::json;
using baysor_test::csv_column;
using baysor_test::read_text_file;

std::set<std::string> csv_assigned_cells(const std::string& path) {
    const auto cells = csv_column(path, "cell");
    const auto noise = csv_column(path, "is_noise");
    std::set<std::string> assigned;
    for (size_t i = 0; i < cells.size(); ++i) {
        if (cells[i] != "0" && noise[i] != "true" && noise[i] != "1") assigned.insert(cells[i]);
    }
    return assigned;
}

/// Polygon cell ids in the `cell_<n>` form the CSV uses.
std::set<std::string> json_polygon_cells(const std::string& path) {
    const json doc = json::parse(read_text_file(path));
    const bool features = doc.at("type") == "FeatureCollection";
    std::set<std::string> cells;
    for (const auto& item : doc.at(features ? "features" : "geometries")) {
        const json& cell = features ? item.at("properties").at("cell") : item.at("cell");
        cells.insert(cell.is_number_integer() ? "cell_" + std::to_string(cell.get<int>())
                                              : cell.get<std::string>());
    }
    return cells;
}

baysor::PolygonCollection triangle(const std::string& cell_name) {
    Eigen::MatrixXd tri(2, 3);
    tri << 0.0, 1.0, 0.0,
           0.0, 0.0, 1.0;
    return {{cell_name, tri}};
}

}  // namespace

TEST(Cov165PolygonFormat, ParsingIsCaseInsensitiveAndValidated) {
    using baysor::PolygonFormat;
    const std::pair<const char*, PolygonFormat> cases[] = {
        {"FeatureCollection", PolygonFormat::FeatureCollection},
        {"featurecollection", PolygonFormat::FeatureCollection},
        {"GeometryCollection", PolygonFormat::GeometryCollection},
        {"geometrycollectionlegacy", PolygonFormat::GeometryCollectionLegacy},
        {"GeometryCollectionLegacy", PolygonFormat::GeometryCollectionLegacy},
        {"NONE", PolygonFormat::None},
    };
    for (const auto& [name, format] : cases) EXPECT_EQ(baysor::parse_polygon_format(name), format) << name;

    EXPECT_THROW(baysor::parse_polygon_format("geometry"), std::invalid_argument);
    EXPECT_THROW(baysor::parse_polygon_format(""), std::invalid_argument);
}

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

    const auto joined = baysor::boundary_polygons_auto(
        data.position_matrix(), assignment, /*estimate_per_z=*/false, &cell_names, /*verbose=*/false).first;

    for (const std::string format : {"FeatureCollection", "GeometryCollection", "GeometryCollectionLegacy"}) {
        const std::string json_path = tmp.file(format + ".json");
        baysor::save_polygons_geojson(joined, json_path, format);
        EXPECT_EQ(json_polygon_cells(json_path), csv_cells) << "format=" << format;
    }

    // The legacy layout is a GeometryCollection with integer ids.
    const json legacy_doc = json::parse(read_text_file(tmp.file("GeometryCollectionLegacy.json")));
    EXPECT_EQ(legacy_doc.at("type"), "GeometryCollection");
    ASSERT_EQ(legacy_doc.at("geometries").size(), 3u);
    for (const auto& geom : legacy_doc.at("geometries")) {
        EXPECT_TRUE(geom.at("cell").is_number_integer()) << geom.dump();
        EXPECT_GE(geom.at("coordinates").at(0).size(), 4u) << geom.dump();
    }
}

TEST(Cov165PolygonFormat, LegacyRejectsNonIntegerCellIds) {
    baysor_test::TempDir tmp("cov165_legacy_bad");
    EXPECT_THROW_MSG(baysor::save_polygons_geojson(triangle("not_a_number"), tmp.file("bad.json"),
                                                   "GeometryCollectionLegacy"),
                     std::runtime_error, "GeometryCollectionLegacy requires integer cell ids, got: not_a_number");
}

TEST(Cov165PolygonFormat, NoneFormatWritesNothingCaseInsensitively) {
    baysor_test::TempDir tmp("cov165_none");
    baysor::save_polygons_geojson(triangle("cell_1"), tmp.file("nothing.json"), "NONE");
    EXPECT_FALSE(std::filesystem::exists(tmp.file("nothing.json")));
}
