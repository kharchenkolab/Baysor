// BUG-1: select_ids_uniformly crashed (NaN index -> out-of-bounds read of
// sum_ids) whenever it was asked for n == 1 centers:
//   - public API cell_centers_uniformly<N>(..., n_clusters=1, ...);
//   - --n-cells-init 1;
//   - tiny datasets where the inferred n_cells_init = div(n, m) * 2 is 0 and
//     cell_centers_uniformly passes it straight into select_ids_uniformly;
//   - exactly one high-confidence molecule left after the n > high_conf clamp.
//
// Julia (Baysor v0.7.1, src/processing/data_processing/initialization.jl,
// select_ids_uniformly) guards this with `if n <= 1 error("n must be > 1")`
// and returns the single surviving high-conf id when length(high_conf_ids)
// drops below n; the C++ port now does the same instead of dividing by
// (n - 1) == 0.
//
// The CLI tests spawn the instrumented `baysor` binary (BAYSOR_CLI_PATH,
// injected by the BAYSOR_WITH_TESTS CMake block) as a subprocess and assert a
// clean exit code 1 with Julia's error message instead of a crash (SIGSEGV
// 139 / SIGABRT 134).

#include <gtest/gtest.h>

#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "baysor/processing/data_processing/initialization.h"

#include "test_cov_helpers.h"

// ============================================================================
// API-level regressions (portable)
// ============================================================================

namespace {

// 6 well-separated 2D points; positions are irrelevant to the crash but keep
// the coordinate sums distinct so the sort is deterministic.
Eigen::MatrixXd make_positions_2d(int n) {
    Eigen::MatrixXd pos(2, n);
    for (int i = 0; i < n; ++i) {
        pos(0, i) = static_cast<double>(i);
        pos(1, i) = static_cast<double>((i * 7) % 5);
    }
    return pos;
}

}  // namespace

TEST(Bug1_SelectIds, ApiNCentersOneThrowsJuliaErrorInsteadOfReadingOutOfBounds) {
    Eigen::MatrixXd pos = make_positions_2d(6);
    // Before the fix: idx = round(0 * (hc-1) / 0) = NaN -> INT_MIN ->
    // sum_ids[INT_MIN] => segfault.
    try {
        auto init = baysor::cell_centers_uniformly<2>(pos, /*n_clusters=*/1,
                                                      /*confidences=*/nullptr,
                                                      /*scale=*/1.0);
        FAIL() << "expected std::runtime_error(\"n must be > 1\"), got "
               << init.centers.rows() << " centers";
    } catch (const std::runtime_error& e) {
        EXPECT_STREQ(e.what(), "n must be > 1");  // same message as Julia's error()
    }
}

TEST(Bug1_SelectIds, ApiNCentersZeroThrowsJuliaErrorInsteadOfReadingOutOfBounds) {
    // Tiny datasets infer n_cells_init = div(n_molecules, min_molecules_per_cell)
    // * 2 == 0. Julia passes 0 straight into select_ids_uniformly, which
    // rejects it with `error("n must be > 1")`; the C++ port now does the
    // same (there is no lower clamp), instead of crashing the way the old
    // clamp-to-1 path did.
    Eigen::MatrixXd pos = make_positions_2d(6);
    try {
        auto init = baysor::cell_centers_uniformly<2>(pos, /*n_clusters=*/0,
                                                      /*confidences=*/nullptr,
                                                      /*scale=*/1.0);
        FAIL() << "expected std::runtime_error(\"n must be > 1\"), got "
               << init.centers.rows() << " centers";
    } catch (const std::runtime_error& e) {
        EXPECT_STREQ(e.what(), "n must be > 1");
    }
}

TEST(Bug1_SelectIds, SingleHighConfidenceMoleculeYieldsSingleCenter) {
    // n >= 2 requested, but only one molecule clears the 0.25 confidence
    // threshold. The n > high_conf clamp drops n to 1, which used to divide
    // by zero; Julia returns that single id (warn + `return high_conf_ids`).
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

// ============================================================================
// CLI end-to-end regressions (POSIX subprocess runs)
// ============================================================================

#ifndef BAYSOR_CLI_PATH

TEST(Bug1_Cli, BaysorCliPathAvailable) {
    GTEST_SKIP() << "BAYSOR_CLI_PATH is not defined; CLI end-to-end tests are disabled";
}

#elif defined(_WIN32)

// The subprocess runner below shells out with sh-style quoting and decodes
// exit codes via WEXITSTATUS, so the end-to-end CLI tests are POSIX-only.
TEST(Bug1_Cli, SubprocessTestsArePosixOnly) {
    GTEST_SKIP() << "CLI subprocess tests require a POSIX shell and sys/wait.h";
}

#else  // BAYSOR_CLI_PATH && !defined(_WIN32)

namespace {

using TempDir = baysor_test::TempDir;
using baysor_test::cli::run_cli;
using baysor_test::cli::write_text;

// 10 molecules / 2 genes: fewer than the -m 100 used below.
std::string tiny_csv_content() {
    std::ostringstream ss;
    ss << "x,y,gene\n";
    const double xs[10] = {1.0, 1.5, 2.0, 2.5, 10.0, 10.5, 11.0, 11.5, 20.0, 20.5};
    const double ys[10] = {1.0, 1.2, 1.0, 1.4, 10.0, 10.2, 10.0, 10.4, 20.0, 20.2};
    for (int i = 0; i < 10; ++i) {
        ss << xs[i] << "," << ys[i] << "," << (i % 2 ? "GeneB" : "GeneA") << "\n";
    }
    return ss.str();
}

// 4 clumps x 50 molecules, 4 genes (>= the default --n-clusters 4 required by
// the default ICA molecule clustering), same shape as the COV-5 datasets.
std::string clumped_csv_content(int per_clump = 50) {
    static const char* kGenes[] = {"GeneA", "GeneB", "GeneC", "GeneD"};
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> jit(-2.0, 2.0);
    std::ostringstream ss;
    ss << "x,y,gene\n";
    int idx = 0;
    for (int c = 0; c < 4; ++c) {
        const double cx = 10.0 + 20.0 * (c % 2);
        const double cy = 10.0 + 20.0 * (c / 2);
        for (int i = 0; i < per_clump; ++i) {
            ss << (cx + jit(rng)) << "," << (cy + jit(rng)) << ","
               << kGenes[idx % 4] << "\n";
            ++idx;
        }
    }
    return ss.str();
}

}  // namespace

TEST(Bug1_Cli, FewerMoleculesThanMinMoleculesPerCellExitsCleanly) {
    TempDir tmp("bug1_tiny");
    const std::string csv = write_text(tmp, "tiny.csv", tiny_csv_content());
    const std::string out = (tmp.path / "seg").string();

    // 10 molecules < -m 100 => inferred n_cells_init = div(10, 100) * 2 = 0,
    // which used to divide by (n - 1) == 0 in select_ids_uniformly and crash
    // (exit 139 in release, 134 with debug assertions). Julia reports a clean
    // "n must be > 1" error instead.
    auto r = run_cli(tmp, "run '" + csv + "' -m 100 -s 2.5 "
                          "-o '" + out + "'");

    EXPECT_EQ(r.exit_code, 1) << "expected a clean error exit, not a crash"
                              << "\n--- stdout ---\n" << r.out
                              << "\n--- stderr ---\n" << r.err;
    EXPECT_NE(r.out.find("n must be > 1"), std::string::npos)
        << "expected Julia's error message"
        << "\n--- stdout ---\n" << r.out;
}

TEST(Bug1_Cli, NCellsInitOneExitsCleanly) {
    TempDir tmp("bug1_nci1");
    const std::string csv = write_text(tmp, "mols.csv", clumped_csv_content());
    const std::string out = (tmp.path / "seg").string();

    // --n-cells-init 1 reaches select_ids_uniformly with n == 1 through the
    // default pipeline and used to crash there (exit 139/134).
    auto r = run_cli(tmp, "run '" + csv + "' -m 10 -s 2.5 "
                          "--n-cells-init 1 -o '" + out + "'");

    EXPECT_EQ(r.exit_code, 1) << "expected a clean error exit, not a crash"
                              << "\n--- stdout ---\n" << r.out
                              << "\n--- stderr ---\n" << r.err;
    EXPECT_NE(r.out.find("n must be > 1"), std::string::npos)
        << "expected Julia's error message"
        << "\n--- stdout ---\n" << r.out;
}

#endif  // BAYSOR_CLI_PATH && !defined(_WIN32)
