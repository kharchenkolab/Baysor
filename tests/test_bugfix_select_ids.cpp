// BUG-1: select_ids_uniformly divided by (n - 1) == 0 and read out of bounds
// whenever it was asked for n <= 1 centers (cell_centers_uniformly with
// n_clusters <= 1, --n-cells-init 1, tiny datasets where the inferred
// n_cells_init is 0) or when one high-confidence molecule was left after the
// n > high_conf clamp. Julia (initialization.jl, select_ids_uniformly) errors
// with "n must be > 1" in the first case and returns the single id in the
// second; the C++ port now does the same.

#include <gtest/gtest.h>

#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "baysor/processing/data_processing/initialization.h"

#include "test_cov_helpers.h"

namespace {

// Distinct coordinate sums, so the sort is deterministic.
Eigen::MatrixXd make_positions_2d(int n) {
    Eigen::MatrixXd pos(2, n);
    for (int i = 0; i < n; ++i) {
        pos(0, i) = static_cast<double>(i);
        pos(1, i) = static_cast<double>((i * 7) % 5);
    }
    return pos;
}

}  // namespace

TEST(Bug1_SelectIds, ApiAtMostOneCenterThrowsJuliaError) {
    // n_clusters = 0 is what tiny datasets infer (div(n, m) * 2); Julia passes
    // it straight into select_ids_uniformly too.
    const Eigen::MatrixXd pos = make_positions_2d(6);
    for (int n_clusters : {0, 1}) {
        try {
            auto init = baysor::cell_centers_uniformly<2>(pos, n_clusters,
                                                          /*confidences=*/nullptr,
                                                          /*scale=*/1.0);
            FAIL() << "n_clusters=" << n_clusters << ": expected an error, got "
                   << init.centers.rows() << " centers";
        } catch (const std::runtime_error& e) {
            EXPECT_STREQ(e.what(), "n must be > 1");  // same message as Julia's error()
        }
    }
}

TEST(Bug1_SelectIds, SingleHighConfidenceMoleculeYieldsSingleCenter) {
    // n >= 2 requested, but only one molecule clears the 0.25 confidence
    // threshold: the clamp drops n to 1 and Julia returns that single id.
    Eigen::MatrixXd pos = make_positions_2d(5);
    std::vector<double> confidences(5, 0.1);
    confidences[0] = 0.95;

    auto init = baysor::cell_centers_uniformly<2>(pos, /*n_clusters=*/3,
                                                  &confidences, /*scale=*/1.0);

    ASSERT_EQ(init.centers.rows(), 1);
    ASSERT_EQ(init.centers.cols(), 2);
    EXPECT_DOUBLE_EQ(init.centers(0, 0), pos(0, 0));
    EXPECT_DOUBLE_EQ(init.centers(0, 1), pos(1, 0));
    ASSERT_EQ(init.assignment.size(), 5u);
    for (int a : init.assignment) EXPECT_EQ(a, 1);
    EXPECT_EQ(init.covs.size(), 1u);
}

#if !defined(_WIN32) && defined(BAYSOR_CLI_PATH)

TEST(Bug1_Cli, AtMostOneInitialCellExitsCleanly) {
    baysor_test::TempDir tmp("bug1_cli");
    std::ostringstream csv;
    csv << "x,y,gene\n";
    const double xs[10] = {1.0, 1.5, 2.0, 2.5, 10.0, 10.5, 11.0, 11.5, 20.0, 20.5};
    const double ys[10] = {1.0, 1.2, 1.0, 1.4, 10.0, 10.2, 10.0, 10.4, 20.0, 20.2};
    for (int i = 0; i < 10; ++i) {
        csv << xs[i] << "," << ys[i] << "," << (i % 2 ? "GeneB" : "GeneA") << "\n";
    }
    const std::string path = baysor_test::cli::write_text(tmp, "tiny.csv", csv.str());

    // 10 molecules with -m 100 infer n_cells_init = div(10, 100) * 2 = 0; the
    // second run asks for one cell. Both used to crash (exit 139 in release,
    // 134 with debug assertions).
    for (const std::string args : {"-m 100 -s 2.5", "-m 10 -s 2.5 --n-cells-init 1"}) {
        auto r = baysor_test::cli::run_cli(
            tmp, "run '" + path + "' " + args + " -o '" + (tmp.path / "seg").string() + "'");

        EXPECT_EQ(r.exit_code, 1) << args << ": expected a clean error exit, not a crash"
                                  << "\n--- stdout ---\n" << r.out
                                  << "\n--- stderr ---\n" << r.err;
        EXPECT_NE(r.out.find("n must be > 1"), std::string::npos)
            << args << ": expected Julia's error message"
            << "\n--- stdout ---\n" << r.out;
    }
}

#else

TEST(Bug1_Cli, RequiresPosixAndCliPath) {
    GTEST_SKIP() << "requires POSIX and BAYSOR_CLI_PATH";
}

#endif
