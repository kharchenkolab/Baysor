// BUG-2 regression tests: ICA initialisation of molecule clustering
//   src/processing/bmm_algorithm/molecule_clustering.cpp
//   (cluster_molecules_ica / fast_ica)
//
// (a) With fewer genes than n_clusters (e.g. 3 genes and the default
//     --n-clusters 4), fast_ica used to silently clamp n_components and the
//     wrapper then indexed W.col(k) past the end of the unmixing matrix
//     (Eigen assertion abort in Debug, garbage init in Release). Julia's
//     MultivariateStats fit(ICA, X, k) rejects k > min(m, n) with an error
//     that Baysor's Julia wrapper (cluster_molecules_on_mrf, DataFrame
//     method) catches to fall back to hash/random initialisation; the C++
//     port now mirrors that: fast_ica throws std::invalid_argument and the
//     wrapper's catch falls back with exprs_init_ptr left/reset to null.
// (b) Any failure inside the ICA try block, including one after the init
//     matrix was built, must discard the partially built init so the log
//     message "falling back to hash initialization" is always true.
//
// Coverage: the catch (const std::exception&) handler (and its reset line)
// is covered by Bug2IcaFallback.FewerGenesThanClustersFallsBackToHashInit;
// catch (...) stays genuinely unreachable (no non-std exception sources).
//
// Note: the original out-of-bounds reads are caught only through Eigen
// assertions, which are compiled in for Debug builds; in Release these
// regression tests would pass even without the fix (the read would just
// return garbage instead of aborting).

#include <gtest/gtest.h>

#include "baysor/processing/bmm_algorithm/molecule_clustering.h"

#include "test_cov_helpers.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <vector>

namespace {

using baysor::AdjList;
using baysor_test::CapturingSink;
using baysor_test::LoggerGuard;

// Chain adjacency over n molecules (0-1-...-(n-1)).
AdjList bug2_chain_adj(int n) {
    std::vector<int> src, dst;
    std::vector<double> wts;
    for (int i = 0; i + 1 < n; ++i) {
        src.push_back(i);
        dst.push_back(i + 1);
        wts.push_back(1.0);
    }
    return AdjList::from_edge_list(
        src.data(), dst.data(), wts.data(), static_cast<int>(src.size()), n);
}

// 12 molecules cycling through exactly 3 genes.
constexpr int kBug2Mols = 12;
const std::vector<int> kBug2Genes = {1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3};

}  // namespace

// ============================================================================
// (a) Fewer genes than clusters: ICA refuses, wrapper falls back (BUG-2a).
// ============================================================================

TEST(Bug2IcaFallback, FewerGenesThanClustersFallsBackToHashInit) {
    auto sink = std::make_shared<CapturingSink>();
    LoggerGuard guard(sink);

    const std::vector<double> confidence(kBug2Mols, 1.0);
    auto adj = bug2_chain_adj(kBug2Mols);

    // 3 genes < 4 clusters: fast_ica would have returned a 3-column W and
    // the wrapper used to read W.col(3) out of bounds. It must instead take
    // the fallback path, exactly like Julia's wrapper on a fit() error.
    auto result = baysor::cluster_molecules_ica(
        kBug2Genes, adj, confidence, /*n_clusters=*/4,
        /*tol=*/0.0, /*mrf_weight=*/1.0, /*max_iters=*/10, /*verbose=*/true);

    const std::string logs = sink->data();
    EXPECT_NE(logs.find("falling back to hash initialization"), std::string::npos)
        << logs;
    EXPECT_EQ(logs.find("ICA initialization succeeded"), std::string::npos) << logs;
    EXPECT_NE(logs.find("k must not exceed min(m, n)"), std::string::npos) << logs;

    // Valid output on the hash fallback path: full dimensions, labels in
    // range, finite row-normalised expression profiles.
    ASSERT_EQ(result.assignment.size(), static_cast<size_t>(kBug2Mols));
    for (int a : result.assignment) {
        EXPECT_GE(a, 1);
        EXPECT_LE(a, 4);
    }
    ASSERT_EQ(result.exprs.rows(), 4);
    ASSERT_EQ(result.exprs.cols(), 3);
    for (int k = 0; k < result.exprs.rows(); ++k) {
        double row_sum = 0.0;
        for (int g = 0; g < result.exprs.cols(); ++g) {
            EXPECT_TRUE(std::isfinite(result.exprs(k, g)));
            row_sum += result.exprs(k, g);
        }
        EXPECT_NEAR(row_sum, 1.0, 1e-9);
    }
    ASSERT_EQ(result.assignment_probs.rows(), 4);
    ASSERT_EQ(result.assignment_probs.cols(), kBug2Mols);

    // The fallback must be the hash initialisation: identical to calling the
    // core EM directly with a null init pointer (deterministic init).
    auto expected = baysor::cluster_molecules_on_mrf(
        kBug2Genes, adj, confidence, /*n_clusters=*/4,
        /*tol=*/0.0, /*mrf_weight=*/1.0, /*max_iters=*/10, /*verbose=*/false,
        /*exprs_init=*/nullptr);
    EXPECT_EQ(result.assignment, expected.assignment);
    ASSERT_EQ(result.exprs.rows(), expected.exprs.rows());
    ASSERT_EQ(result.exprs.cols(), expected.exprs.cols());
    EXPECT_NEAR((result.exprs - expected.exprs).cwiseAbs().maxCoeff(), 0.0, 1e-12);
}

// Boundary: n_clusters == n_genes must still take the ICA path (the check
// rejects only k > min(m, n), not k == min(m, n)).
TEST(Bug2IcaFallback, EqualGenesAndClustersStillUsesIca) {
    auto sink = std::make_shared<CapturingSink>();
    LoggerGuard guard(sink);

    const std::vector<double> confidence(kBug2Mols, 1.0);
    auto adj = bug2_chain_adj(kBug2Mols);

    auto result = baysor::cluster_molecules_ica(
        kBug2Genes, adj, confidence, /*n_clusters=*/3,
        /*tol=*/0.0, /*mrf_weight=*/1.0, /*max_iters=*/10, /*verbose=*/true);

    const std::string logs = sink->data();
    EXPECT_NE(logs.find("ICA initialization succeeded (3 components)"),
              std::string::npos) << logs;
    EXPECT_EQ(logs.find("falling back to hash initialization"), std::string::npos)
        << logs;

    ASSERT_EQ(result.assignment.size(), static_cast<size_t>(kBug2Mols));
    ASSERT_EQ(result.exprs.rows(), 3);
    ASSERT_EQ(result.exprs.cols(), 3);
    ASSERT_EQ(result.assignment_probs.rows(), 3);
    ASSERT_EQ(result.assignment_probs.cols(), kBug2Mols);
}

// ============================================================================
// CLI end-to-end: 3-gene dataset with the default --n-clusters (4) must run
// to completion. Before the fix this aborted inside cluster_molecules_ica
// (Eigen assertion on W.col(3) of the 3x3 unmixing matrix, exit 134).
// POSIX-only, following tests/test_cov_cli.cpp.
// ============================================================================

#ifndef BAYSOR_CLI_PATH

TEST(Bug2IcaFallback, CliPathAvailable) {
    GTEST_SKIP() << "BAYSOR_CLI_PATH is not defined; CLI end-to-end test disabled";
}

#elif defined(_WIN32)

TEST(Bug2IcaFallback, CliSubprocessTestIsPosixOnly) {
    GTEST_SKIP() << "CLI subprocess tests require a POSIX shell and sys/wait.h";
}

#else  // BAYSOR_CLI_PATH && !defined(_WIN32)

namespace {

namespace fs = std::filesystem;

// 4 clumps x 30 molecules, exactly 3 genes (cycling), no prior column.
fs::path bug2_write_three_gene_csv(const baysor_test::TempDir& tmp) {
    const fs::path p = tmp.path / "mols_3gene.csv";
    std::ofstream f(p);
    f << "x,y,gene\n";
    static const char* kGenes[] = {"GeneA", "GeneB", "GeneC"};
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> jit(-2.0, 2.0);
    int idx = 0;
    for (int c = 0; c < 4; ++c) {
        const double cx = 10.0 + 20.0 * (c % 2);
        const double cy = 10.0 + 20.0 * (c / 2);
        for (int i = 0; i < 30; ++i, ++idx) {
            f << (cx + jit(rng)) << "," << (cy + jit(rng)) << ","
              << kGenes[idx % 3] << "\n";
        }
    }
    return p;
}

}  // namespace

TEST(Bug2IcaFallback, CliThreeGeneDatasetWithDefaultOptionsExitsZero) {
    baysor_test::TempDir tmp("bug2_cli_3gene");
    auto csv = bug2_write_three_gene_csv(tmp);
    const fs::path out = tmp.path / "seg";

    // Only required options (coordinates/scale); --n-clusters stays at its
    // default of 4, which exceeds the 3 genes and used to abort the run
    // inside cluster_molecules_ica before producing any output.
    auto r = baysor_test::cli::run_cli(tmp, "run '" + csv.string() +
                                       "' -m 10 -s 2.5 -o '" + out.string() + "'");
    EXPECT_EQ(r.exit_code, 0)
        << "--- stdout ---\n" << r.out << "\n--- stderr ---\n" << r.err;

    const std::string logs = r.out + r.err;
    // Default n-clusters = 4 reached the ICA wrapper ...
    EXPECT_NE(logs.find("Clustering molecules into 4 types"), std::string::npos)
        << logs;
    // ... which refused (3 genes) and fell back instead of crashing.
    EXPECT_NE(logs.find("falling back to hash initialization"), std::string::npos)
        << logs;
    EXPECT_NE(logs.find("Segmentation complete"), std::string::npos) << logs;

    const fs::path seg_csv = out / "segmentation.csv";
    ASSERT_TRUE(fs::is_regular_file(seg_csv)) << seg_csv;
    EXPECT_GT(fs::file_size(seg_csv), 0u);
}

#endif  // BAYSOR_CLI_PATH && !defined(_WIN32)
