// BUG-5: `baysor run` and `segfree` aborted on tiny datasets in the NCV colour
// embedding: the PCA truncation n_pca = min(n_pca_dims, n_components) could
// exceed the columns of the thin SVD basis (fewer anchors than n_pca_dims),
// and graph_k <= 0 left the interpolation k-NN empty. The truncation is now
// also clamped by the number of anchors and the k-NN k to >= 1. The original
// out-of-bounds reads abort only with Eigen assertions (Debug builds).

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

#include "baysor/reporting/color_utils.h"

#include "test_cov_helpers.h"

// ============================================================================
// API-level regressions (portable)
// ============================================================================

namespace {

// n_components x n_mols random gene-composition vectors; 20 matches the
// pipeline default (estimate_gene_vectors n_components=20).
Eigen::MatrixXf random_mol_vecs(int n_components, int n_mols, unsigned seed = 42) {
    Eigen::MatrixXf mol_vecs(n_components, n_mols);
    std::mt19937 rng(seed);
    std::normal_distribution<float> dist(0.0f, 1.0f);
    for (int i = 0; i < n_components; ++i)
        for (int j = 0; j < n_mols; ++j)
            mol_vecs(i, j) = dist(rng);
    return mol_vecs;
}

void expect_valid_hex_colors(const std::vector<std::string>& colors, size_t expected_n) {
    ASSERT_EQ(colors.size(), expected_n);
    for (const auto& c : colors) {
        ASSERT_EQ(c.size(), 7u) << "bad colour: " << c;
        EXPECT_EQ(c[0], '#') << "bad colour: " << c;
    }
}

// 9 molecules spread over a 20x20 window; positions only feed the neighbourhood
// graphs, but keep them well separated so kNN results are stable.
Eigen::MatrixXd tiny_positions(int n) {
    Eigen::MatrixXd pos(2, n);
    for (int i = 0; i < n; ++i) {
        pos(0, i) = 2.5 * (i % 4);
        pos(1, i) = 2.0 * (i / 4) + 0.3 * (i % 3);
    }
    return pos;
}

}  // namespace

// The reproducer shape: 20 components, 9 anchors, default n_pca_dims=10.
TEST(Bug5_NcvColor, EmbeddingWithFewerAnchorsThanPcaDimsExitsCleanly) {
    const Eigen::MatrixXf mol_vecs = random_mol_vecs(20, 9);
    const std::vector<double> confidence(9, 0.99);

    auto colors = baysor::gene_composition_color_embedding(mol_vecs, confidence);

    expect_valid_hex_colors(colors, 9u);
}

// Same input through the report entry point (include_report_umap=true), which
// additionally fits the 2D report UMAP on the same 9 anchors.
TEST(Bug5_NcvColor, ReportEmbeddingWithFewerAnchorsThanPcaDimsExitsCleanly) {
    const Eigen::MatrixXf mol_vecs = random_mol_vecs(20, 9);
    const std::vector<double> confidence(9, 0.99);

    auto res = baysor::gene_composition_report_embedding(mol_vecs, confidence);

    expect_valid_hex_colors(res.colors, 9u);
    EXPECT_GT(res.anchor_count, 1);
    ASSERT_EQ(res.sample_ids.size(), 9u);
    ASSERT_EQ(res.sample_umap_x.size(), res.sample_ids.size());
    ASSERT_EQ(res.sample_umap_y.size(), res.sample_ids.size());
    for (double v : res.sample_umap_x) EXPECT_TRUE(std::isfinite(v));
    for (double v : res.sample_umap_y) EXPECT_TRUE(std::isfinite(v));
}

// The streaming entry point used by `baysor run` / `baysor preview`: all 9
// molecules become anchors.
TEST(Bug5_NcvColor, StreamingWithFewerAnchorsThanPcaDimsExitsCleanly) {
    constexpr int n = 9;
    const Eigen::MatrixXd pos = tiny_positions(n);
    std::vector<int> genes(n);
    for (int i = 0; i < n; ++i) genes[i] = i % 5;  // 5 genes, as in tiny.csv
    const std::vector<double> confidence(n, 0.99);

    auto colors = baysor::gene_composition_color_embedding_streaming(
        pos, genes, /*n_genes=*/5, confidence,
        /*k_neighbors=*/8,
        /*basis_sample_size=*/20000,
        /*sample_size=*/20000,
        /*seed=*/42,
        /*n_pca_dims=*/10,
        /*graph_k=*/15
    );

    expect_valid_hex_colors(colors, n);
}

// Extreme case: the basis is capped at 2 anchors -> UMAP runs with
// n_neighbors=1 and the PCA basis is clamped to 2 columns.
TEST(Bug5_NcvColor, StreamingWithTwoAnchorsExitsCleanly) {
    constexpr int n = 4;
    const Eigen::MatrixXd pos = tiny_positions(n);
    std::vector<int> genes(n);
    for (int i = 0; i < n; ++i) genes[i] = i % 2;
    const std::vector<double> confidence(n, 0.99);

    auto colors = baysor::gene_composition_color_embedding_streaming(
        pos, genes, /*n_genes=*/2, confidence,
        /*k_neighbors=*/3,
        /*basis_sample_size=*/2,
        /*sample_size=*/2,
        /*seed=*/42,
        /*n_pca_dims=*/10,
        /*graph_k=*/15
    );

    expect_valid_hex_colors(colors, n);
}

// graph_k=0 must not leave the interpolation k-NN of compute_ncv_embedding empty.
TEST(Bug5_NcvColor, ZeroGraphKDoesNotReadPastEmptyKnnResults) {
    const Eigen::MatrixXf mol_vecs = random_mol_vecs(10, 30);
    const std::vector<double> confidence(30, 0.99);

    auto colors = baysor::gene_composition_color_embedding(
        mol_vecs, confidence, /*sample_size=*/20000, /*seed=*/42,
        /*n_pca_dims=*/10, /*graph_k=*/0);

    expect_valid_hex_colors(colors, 30u);
}

// ============================================================================
// CLI end-to-end regressions (POSIX subprocess runs)
// ============================================================================

#if !defined(_WIN32) && defined(BAYSOR_CLI_PATH)

namespace {

namespace fs = std::filesystem;

using TempDir = baysor_test::TempDir;
using baysor_test::cli::run_cli;
using baysor_test::cli::write_text;

// The BUG-5 reproducer table: 10 molecules / 5 genes giving 9 NCV anchors.
std::string tiny_csv_content() {
    return
        "x,y,gene\n"
        "2.687,16.949,A\n"
        "5.101,9.909,D\n"
        "9.445,7.592,B\n"
        "1.877,0.567,D\n"
        "8.655,15.246,A\n"
        "13.917,5.327,B\n"
        "11.823,2.045,C\n"
        "0.612,0.509,E\n"
        "0.184,17.625,B\n"
        "19.381,14.517,E\n";
}

}  // namespace

// The original reproducer `baysor run tiny.csv -m 2 -s 2.5 -o out`; --iters 10
// keeps every cluster method cheap.
TEST(Bug5_Cli, TinyDatasetRunExitsCleanlyForAllClusterMethods) {
    TempDir tmp("bug5_run");
    const std::string csv = write_text(tmp, "tiny.csv", tiny_csv_content());

    for (const char* method : {"mrf", "louvain", "leiden", "none"}) {
        const std::string out = (tmp.path / ("seg_" + std::string(method))).string();
        auto r = run_cli(tmp, "run '" + csv + "' -m 2 -s 2.5 --iters 10 "
                              "--cluster-method " + method + " -o '" + out + "'");
        EXPECT_EQ(r.exit_code, 0) << "cluster-method=" << method
                                  << " must exit cleanly, not crash"
                                  << "\n--- stdout ---\n" << r.out
                                  << "\n--- stderr ---\n" << r.err;
        EXPECT_TRUE(fs::exists(fs::path(out) / "segmentation.csv"))
            << "cluster-method=" << method << " produced no segmentation.csv";
    }
}

// `baysor segfree` drives the non-streaming entry point
// (gene_composition_color_embedding).
TEST(Bug5_Cli, TinyDatasetSegfreeExitsCleanly) {
    TempDir tmp("bug5_segfree");
    const std::string csv = write_text(tmp, "tiny.csv", tiny_csv_content());
    const std::string out = (tmp.path / "ncvs.loom").string();

    auto r = run_cli(tmp, "segfree '" + csv + "' -m 2 -o '" + out + "'");

    EXPECT_EQ(r.exit_code, 0) << "segfree must exit cleanly, not crash"
                              << "\n--- stdout ---\n" << r.out
                              << "\n--- stderr ---\n" << r.err;
    EXPECT_TRUE(fs::exists(out)) << "segfree produced no output loom";
}

#else

TEST(Bug5_Cli, SubprocessTestsNeedPosixAndCliPath) {
    GTEST_SKIP() << "CLI subprocess tests require POSIX and BAYSOR_CLI_PATH";
}

#endif
