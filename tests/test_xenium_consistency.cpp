// Regression tests for kharchenkolab/Baysor#165.
//
// Xenium Ranger `import-segmentation` requires that segmentation.csv and the
// polygons JSON contain exactly the same cells: every cell id in the polygons
// must have at least one assigned transcript, and every transcript-assigned
// cell must have a polygon. Baysor v0.7.1 could violate both directions; this
// file pins the invariant for the GeoJSON formats and checks that a cell
// whose free-form boundary estimation fails still gets a fallback polygon.

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
#include "baysor/processing/utils/convex_hull.h"
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

/// Cell names in the segmentation CSV, excluding noise/unassigned rows.
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
        cells.insert(cell);
    }
    return cells;
}

/// Cell names in a `FeatureCollection` or string-id `GeometryCollection`.
std::set<std::string> json_polygon_cells(const std::string& path) {
    const json doc = json::parse(read_text(path));
    std::set<std::string> cells;
    if (doc.at("type") == "FeatureCollection") {
        for (const auto& feat : doc.at("features")) {
            cells.insert(feat.at("properties").at("cell").get<std::string>());
        }
    } else {
        for (const auto& geom : doc.at("geometries")) {
            cells.insert(geom.at("cell").get<std::string>());
        }
    }
    return cells;
}

/// Three cells: a square, a collinear (degenerate) cell whose Delaunay
/// triangulation has no faces, and a two-molecule cell.
baysor::MoleculeData make_consistency_data() {
    baysor::MoleculeData data;
    data.x = {0.0, 10.0, 10.0, 0.0,
              100.0, 105.0, 110.0,
              200.0, 210.0};
    data.y = {0.0, 0.0, 10.0, 10.0,
              0.0, 0.0, 0.0,
              0.0, 0.0};
    data.gene.assign(data.x.size(), 1);
    data.gene_names = {"G1"};
    return data;
}

std::vector<int> consistency_assignment() {
    return {1, 1, 1, 1, 2, 2, 2, 3, 3};
}

}  // namespace

// ============================================================================
// CSV cell set == polygons cell set, for every GeoJSON format (2D and 3D)
// ============================================================================

TEST(Cov165Consistency, CsvAndPolygonsCellSetsMatchAllFormats) {
    baysor_test::TempDir tmp("cov165_2d");
    const auto data = make_consistency_data();
    const auto assignment = consistency_assignment();
    const std::vector<std::string> cell_names = {"cell_1", "cell_2", "cell_3"};

    const std::string csv_path = tmp.file("segmentation.csv");
    baysor::save_segmented_df(data, assignment, data.gene_names, csv_path);
    const auto csv_cells = csv_assigned_cells(csv_path);
    ASSERT_EQ(csv_cells, (std::set<std::string>{"cell_1", "cell_2", "cell_3"}));

    auto [joined, stack] = baysor::boundary_polygons_auto(
        data.position_matrix(), assignment, /*estimate_per_z=*/false,
        &cell_names, /*verbose=*/false);
    (void)stack;

    const std::vector<std::pair<std::string, std::string>> formats = {
        {"FeatureCollection", "FeatureCollection"},
        {"featurecollection", "FeatureCollection"},
        {"GeometryCollection", "GeometryCollection"},
        {"geometrycollection", "GeometryCollection"},
    };
    for (const auto& [format, expected_type] : formats) {
        const std::string json_path = tmp.file("poly_" + format + ".json");
        baysor::save_polygons_geojson(joined, json_path, format);
        const auto poly_cells = json_polygon_cells(json_path);
        EXPECT_EQ(poly_cells, csv_cells) << "format=" << format;
        EXPECT_EQ(json::parse(read_text(json_path)).at("type"), expected_type)
            << "format=" << format;
    }

    // The degenerate cell must have a valid fallback polygon.
    ASSERT_EQ(joined.count("cell_2"), 1u);
    EXPECT_GE(joined.at("cell_2").cols(), 3);
    EXPECT_GT(baysor::polygon_area(joined.at("cell_2")), 0.0);
}

TEST(Cov165Consistency, CsvAndPolygonsCellSetsMatchIn3D) {
    baysor_test::TempDir tmp("cov165_3d");
    baysor::MoleculeData data = make_consistency_data();
    // Give every molecule a per-cell z so the run is 3D; the collinear cell
    // stays collinear in the pooled 2D projection.
    data.z.assign(data.x.size(), 0.0);
    for (size_t i = 0; i < data.z.size(); ++i) {
        data.z[i] = (i % 2 == 0) ? 0.0 : 5.0;
    }
    const auto assignment = consistency_assignment();
    const std::vector<std::string> cell_names = {"cell_1", "cell_2", "cell_3"};

    const std::string csv_path = tmp.file("segmentation.csv");
    baysor::save_segmented_df(data, assignment, data.gene_names, csv_path);
    const auto csv_cells = csv_assigned_cells(csv_path);

    auto [joined, stack] = baysor::boundary_polygons_auto(
        data.position_matrix(), assignment, /*estimate_per_z=*/true,
        &cell_names, /*verbose=*/false);

    baysor::OutputPaths paths;
    paths.polygons_2d = tmp.file("polygons_2d.json");
    paths.polygons_3d = tmp.file("polygons_3d.json");
    baysor::save_polygon_stack_geojson(stack, paths, "FeatureCollection");

    // The 2D handoff file must contain every assigned cell.
    EXPECT_EQ(json_polygon_cells(paths.polygons_2d), csv_cells);
    EXPECT_EQ(joined.size(), csv_cells.size());
}

// ============================================================================
// Fallback / naming helpers
// ============================================================================

TEST(Cov165Consistency, DegenerateCellGetsFallbackPolygon) {
    Eigen::MatrixXd pos(2, 5);
    pos << 100.0, 105.0, 110.0, 115.0, 120.0,
            0.0,   0.0,   0.0,   0.0,   0.0;
    std::vector<int> labels(5, 1);

    auto polys = baysor::boundary_polygons(pos, labels);
    ASSERT_EQ(polys.count("1"), 1u);
    EXPECT_GE(polys.at("1").cols(), 3);
    EXPECT_GT(baysor::polygon_area(polys.at("1")), 0.0);
}

TEST(Cov165Consistency, OutOfRangeCellNamesFallBackToCsvNames) {
    // Two separated 2x2 squares; the caller only names the first component.
    Eigen::MatrixXd pos(2, 8);
    pos << 0.0, 2.0, 2.0, 0.0, 10.0, 12.0, 12.0, 10.0,
           0.0, 0.0, 2.0, 2.0,  0.0,  0.0,  2.0,  2.0;
    std::vector<int> labels = {1, 1, 1, 1, 2, 2, 2, 2};
    std::vector<std::string> names = {"cell_1"};

    auto polys = baysor::boundary_polygons(pos, labels, &names);
    ASSERT_EQ(polys.size(), 2u);
    // Must match the CSV writer naming, not the bare-integer default.
    EXPECT_EQ(polys.count("cell_1"), 1u);
    EXPECT_EQ(polys.count("cell_2"), 1u);
    EXPECT_EQ(polys.count("2"), 0u);
}
