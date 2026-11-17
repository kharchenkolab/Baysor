// Regression tests for BUG-3 (three small correctness issues):
//  (a) AdjList::from_edge_list overflow guards
//  (b) convex-hull orientation must be clockwise, as in Julia
//  (c) config values that fail to parse must raise a clear error

#include <gtest/gtest.h>

#include "baysor/processing/models/adj_list.h"
#include "baysor/processing/utils/convex_hull.h"
#include "baysor/utils/options.h"

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <type_traits>
#include <utility>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

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

// from_edge_list must throw E naming `arg` on the arguments alone, before
// allocating or touching the (null) edge arrays.
template<class E>
void expect_guard_throws(int64_t n_edges, int n_verts, const char* arg) {
    try {
        baysor::AdjList::from_edge_list(nullptr, nullptr, nullptr, n_edges, n_verts);
        FAIL() << "expected throw for n_edges=" << n_edges << ", n_verts=" << n_verts;
    } catch (const E& e) {
        EXPECT_NE(std::string(e.what()).find(arg), std::string::npos) << e.what();
    }
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

TEST(Bug3_AdjListGuards, OversizedOrNegativeCountsThrow) {
    // n_verts + 1 overflows int at INT_MAX, and negative counts are invalid.
    // (2 * n_edges is no longer limited to int: offsets are int64, see
    // Issue187_AdjList64 below.)
    expect_guard_throws<std::length_error>(0, std::numeric_limits<int>::max(), "n_verts");
    expect_guard_throws<std::invalid_argument>(0, -1, "n_verts");
    expect_guard_throws<std::invalid_argument>(-1, 4, "n_edges");
}

// Issue #187: graphs with more than 2^30 undirected edges (2D Delaunay graphs of
// more than ~358M molecules) overflowed the int32 CSR offsets.
static_assert(std::is_same_v<decltype(baysor::AdjList::indptr)::value_type, int64_t>,
              "AdjList row offsets must be 64-bit");
static_assert(std::is_same_v<decltype(std::declval<baysor::AdjList>().nnz()), int64_t>,
              "AdjList::nnz() must be 64-bit");

// Builds a graph with exactly 2^30 edges (2 * n_edges == 2^31 directed edges)
// and checks the CSR invariants. Needs ~45 GB RAM, so it only runs when
// BAYSOR_LARGE_TESTS is set.
TEST(Issue187_AdjList64, MoreThanInt32MaxDirectedEdges) {
    if (std::getenv("BAYSOR_LARGE_TESTS") == nullptr)
        GTEST_SKIP() << "set BAYSOR_LARGE_TESTS=1 to run (needs ~45 GB RAM)";
    const int64_t n_edges = int64_t{1} << 30;
    const int n_verts = 1000000;
    std::vector<int> src(n_edges), dst(n_edges);
    std::vector<double> w(n_edges, 1.0);
    for (int64_t e = 0; e < n_edges; ++e) {
        src[e] = static_cast<int>(e % n_verts);
        dst[e] = (src[e] + 1) % n_verts;
    }
    const auto adj = baysor::AdjList::from_edge_list(src.data(), dst.data(), w.data(), n_edges, n_verts);
    EXPECT_EQ(adj.nnz(), 2 * n_edges);
    EXPECT_EQ(adj.indptr.back(), 2 * n_edges);
    int64_t total = 0;
    for (int i = 0; i < n_verts; ++i) total += adj.neighbor_count(i);
    EXPECT_EQ(total, 2 * n_edges);
    ASSERT_GT(adj.neighbor_count(n_verts - 1), 0);
    EXPECT_EQ(adj.neighbor_ids(n_verts - 1)[0], 0);  // edge (n_verts-1) -- 0
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

TEST(Bug3_ConfigErrors, UnparsableValuesNameKeyValueAndType) {
    // `12abc`: std::stoi would return 12, but Julia's TOML parser rejects it.
    const struct {
        const char* toml;
        std::vector<std::string> expected;  // substrings of the error message
    } cases[] = {
        {"[molecules]\nmin_molecules_per_cell = notanint\n",
         {"min_molecules_per_cell", "notanint", "integer"}},
        {"[segmentation]\niters = 12abc\n", {"iters", "12abc", "integer"}},
        {"[segmentation]\nn_clusters = 99999999999\n", {"n_clusters", "99999999999", "integer"}},
        {"[molecules]\nmin_qv = notafloat\n", {"min_qv", "notafloat", "number"}},
        {"[segmentation]\nestimate_scale_from_centers = bogus\n",
         {"estimate_scale_from_centers", "bogus", "boolean"}},
        {"[prior]\ntype = \"bogus\"\n", {"'type'", "bogus", "'column'", "'boundary'"}},
    };
    for (const auto& c : cases) {
        const std::string msg = config_error(c.toml);
        ASSERT_FALSE(msg.empty()) << "silently kept the default: " << c.toml;
        for (const auto& part : c.expected) {
            EXPECT_NE(msg.find(part), std::string::npos) << msg;
        }
    }
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
