// Regression tests for BUG-3 (three small correctness issues):
//
//  (a) AdjList::from_edge_list overflow guards (src/processing/models/adj_list.cpp)
//  (b) Convex-hull orientation must be clockwise (behaviour pinned by
//      SquareHullIsClockwiseLikeJulia; src/processing/utils/convex_hull.cpp,
//      include/baysor/processing/utils/convex_hull.h)
//  (c) Config values that fail to parse must raise a clear error
//      (src/utils/options.cpp)
//
// Every suite name is prefixed with Bug3 so it cannot clash with the other
// test files.

#include <gtest/gtest.h>

#include "baysor/processing/models/adj_list.h"
#include "baysor/processing/utils/convex_hull.h"
#include "baysor/utils/options.h"

#include <climits>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

#include "test_cov_helpers.h"

namespace {

using TempDir = baysor_test::TempDir;

std::string write_file(const TempDir& dir, const std::string& name,
                       const std::string& content) {
    const std::string p = dir.file(name);
    std::ofstream f(p);
    f << content;
    return p;
}

// Load `content` as a TOML config; return the runtime_error message, or an
// empty string if load_config did not throw.
std::string config_error(const std::string& content) {
    TempDir dir("bug3_cfg");
    const auto path = write_file(dir, "cfg.toml", content);
    try {
        (void)baysor::load_config(path);
    } catch (const std::runtime_error& e) {
        return e.what();
    }
    return std::string();
}

// Signed shoelace area (positive = counter-clockwise, negative = clockwise).
double signed_area(const Eigen::MatrixXd& poly) {
    double s = 0.0;
    const int n = static_cast<int>(poly.cols());
    for (int i = 0; i < n; ++i) {
        const int j = (i + 1) % n;
        s += poly(0, i) * poly(1, j) - poly(1, i) * poly(0, j);
    }
    return s / 2.0;
}

} // namespace

// ============================================================================
// (a) AdjList::from_edge_list overflow guards
// ============================================================================

TEST(Bug3_AdjListGuards, NVertsAtIntMaxThrowsLengthError) {
    // n_verts + 1 overflows int at INT_MAX (undefined behaviour). The guard
    // must trigger on the arguments alone, without allocating anything big.
    try {
        baysor::AdjList::from_edge_list(nullptr, nullptr, nullptr,
                                        /*n_edges=*/0,
                                        std::numeric_limits<int>::max());
        FAIL() << "expected throw";
    } catch (const std::length_error& e) {
        EXPECT_NE(std::string(e.what()).find("n_verts"), std::string::npos)
            << e.what();
    }
}

TEST(Bug3_AdjListGuards, NEdgesPastHalfIntMaxThrowsLengthError) {
    // 2 * n_edges overflows int at n_edges = 2^30 (about 1.07e9). Without the
    // guard the unpatched code dereferences the (null) edge arrays first, so
    // this test crashes before it can pass.
    try {
        baysor::AdjList::from_edge_list(nullptr, nullptr, nullptr,
                                        /*n_edges=*/(1 << 30),
                                        /*n_verts=*/1);
        FAIL() << "expected throw";
    } catch (const std::length_error& e) {
        EXPECT_NE(std::string(e.what()).find("n_edges"), std::string::npos)
            << e.what();
    }
}

TEST(Bug3_AdjListGuards, NegativeCountsAreInvalidArguments) {
    // Negative sizes are invalid arguments; the unpatched code either returns
    // nonsense (n_verts = -1 with n_edges = 0) or dies inside
    // std::vector::resize with an unrelated message (n_edges = -1).
    try {
        baysor::AdjList::from_edge_list(nullptr, nullptr, nullptr,
                                        /*n_edges=*/0, /*n_verts=*/-1);
        FAIL() << "expected throw";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("n_verts"), std::string::npos)
            << e.what();
    }
    try {
        baysor::AdjList::from_edge_list(nullptr, nullptr, nullptr,
                                        /*n_edges=*/-1, /*n_verts=*/4);
        FAIL() << "expected throw";
    } catch (const std::invalid_argument& e) {
        EXPECT_NE(std::string(e.what()).find("n_edges"), std::string::npos)
            << e.what();
    }
}

TEST(Bug3_AdjListGuards, ValidInputStillBuilds) {
    // The guards must not disturb the normal path.
    const int src[] = {0, 1};
    const int dst[] = {1, 2};
    const double w[] = {0.5, 1.5};
    const auto adj = baysor::AdjList::from_edge_list(src, dst, w, 2, 3);
    EXPECT_EQ(adj.n_molecules(), 3);
    EXPECT_EQ(adj.nnz(), 4);
}

// ============================================================================
// (b) Convex-hull orientation
// ============================================================================

TEST(Bug3_HullOrientation, SquareHullIsClockwiseLikeJulia) {
    // origin/master:src/processing/utils/convex_hull.jl returns the hull in
    // clockwise order (its own test expects [0 0 1 1 0; 0 1 1 0 0], whose
    // signed shoelace area is -1). The C++ port must match, so the signed
    // area of a square hull must be negative (clockwise), not positive.
    Eigen::MatrixXd pts(2, 5);
    pts << 0.0, 4.0, 4.0, 0.0, 2.0,
           0.0, 0.0, 4.0, 4.0, 2.0;

    const Eigen::MatrixXd hull = baysor::convex_hull(pts);
    ASSERT_EQ(hull.rows(), 2);
    ASSERT_EQ(hull.cols(), 4);

    const double area = signed_area(hull);
    EXPECT_LT(area, 0.0) << "hull must be in clockwise order";
    EXPECT_NEAR(area, -16.0, 1e-12);
}

// ============================================================================
// (c) Config values that fail to parse raise a clear error
// ============================================================================

TEST(Bug3_ConfigErrors, UnparsableIntNamesKeyValueAndType) {
    const std::string msg = config_error(
        "[molecules]\nmin_molecules_per_cell = notanint\n");
    ASSERT_FALSE(msg.empty()) << "garbage int silently kept its default";
    EXPECT_NE(msg.find("min_molecules_per_cell"), std::string::npos) << msg;
    EXPECT_NE(msg.find("notanint"), std::string::npos) << msg;
    EXPECT_NE(msg.find("integer"), std::string::npos) << msg;
}

TEST(Bug3_ConfigErrors, PartiallyParsedIntIsRejected) {
    // std::stoi("12abc") returns 12 and stops at 'a'; Julia's TOML parser
    // rejects `12abc` outright, so the C++ port must not accept it either.
    const std::string msg = config_error(
        "[segmentation]\niters = 12abc\n");
    ASSERT_FALSE(msg.empty()) << "partially parsable int silently kept default";
    EXPECT_NE(msg.find("iters"), std::string::npos) << msg;
    EXPECT_NE(msg.find("12abc"), std::string::npos) << msg;
    EXPECT_NE(msg.find("integer"), std::string::npos) << msg;
}

TEST(Bug3_ConfigErrors, OutOfRangeIntIsRejected) {
    const std::string msg = config_error(
        "[segmentation]\nn_clusters = 99999999999\n");
    ASSERT_FALSE(msg.empty()) << "out-of-range int silently kept default";
    EXPECT_NE(msg.find("n_clusters"), std::string::npos) << msg;
    EXPECT_NE(msg.find("99999999999"), std::string::npos) << msg;
    EXPECT_NE(msg.find("integer"), std::string::npos) << msg;
}

TEST(Bug3_ConfigErrors, UnparsableDoubleNamesKeyValueAndType) {
    const std::string msg = config_error(
        "[molecules]\nmin_qv = notafloat\n");
    ASSERT_FALSE(msg.empty()) << "garbage double silently kept its default";
    EXPECT_NE(msg.find("min_qv"), std::string::npos) << msg;
    EXPECT_NE(msg.find("notafloat"), std::string::npos) << msg;
    EXPECT_NE(msg.find("number"), std::string::npos) << msg;
}

TEST(Bug3_ConfigErrors, UnparsableBoolNamesKeyValueAndType) {
    const std::string msg = config_error(
        "[segmentation]\nestimate_scale_from_centers = bogus\n");
    ASSERT_FALSE(msg.empty()) << "garbage bool silently kept its default";
    EXPECT_NE(msg.find("estimate_scale_from_centers"), std::string::npos) << msg;
    EXPECT_NE(msg.find("bogus"), std::string::npos) << msg;
    EXPECT_NE(msg.find("boolean"), std::string::npos) << msg;
}

TEST(Bug3_ConfigErrors, UnknownPriorTypeNamesKeyValueAndOptions) {
    const std::string msg = config_error(
        "[prior]\ntype = \"bogus\"\n");
    ASSERT_FALSE(msg.empty()) << "unknown prior type silently kept its default";
    EXPECT_NE(msg.find("'type'"), std::string::npos) << msg;
    EXPECT_NE(msg.find("bogus"), std::string::npos) << msg;
    EXPECT_NE(msg.find("'column'"), std::string::npos) << msg;
    EXPECT_NE(msg.find("'boundary'"), std::string::npos) << msg;
}

TEST(Bug3_ConfigErrors, ValidConfigAndMissingKeysStillUseDefaults) {
    // Keys that are absent keep their defaults (matching Julia, where missing
    // TOML keys fall back to the @option defaults), and well-formed values
    // keep parsing.
    TempDir dir("bug3_cfg");
    const auto path = write_file(dir, "ok.toml",
        "[molecules]\nmin_molecules_per_cell = 7\n");
    const auto opts = baysor::load_config(path);
    EXPECT_EQ(opts.molecules.min_molecules_per_cell, 7);
    EXPECT_DOUBLE_EQ(opts.molecules.min_qv, -1.0);  // key absent -> default
    EXPECT_EQ(opts.segmentation.iters, 500);        // key absent -> default
}
