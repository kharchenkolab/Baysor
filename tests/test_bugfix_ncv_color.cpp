// BUG-5: `baysor run` aborted on tiny datasets inside the NCV colour
// embedding. `fit_ncv_interpolation_model` (and its non-streaming twin
// `compute_ncv_embedding`) clamped the PCA truncation with
//   n_pca = min(n_pca_dims, n_components)
// but the thin SVD basis of an (n_components x n_anchors) matrix only has
// min(rows, cols) columns. Whenever a run produced fewer anchors/samples than
// min(n_pca_dims, n_components) -- always true for tiny datasets, e.g. 9
// anchors vs n_pca_dims=10 -- `leftCols(n_pca)` read past the end of
// matrixU(): an Eigen assertion (SIGABRT, exit 134) in Debug and a silent
// out-of-bounds read in Release.
//
// Julia (Baysor v0.7.1, src/processing/data_processing/neighborhood_composition.jl,
// gene_composition_color_embedding) has no intermediate PCA truncation at all:
// it fits UmapFit on `pca[:, sample_ids]` directly and interpolates with
// `knn_parallel(tree, x, nn_interpolate)` (umap_wrappers.jl). The C++ PCA step
// is a speed optimisation, so it now clamps to the number of columns that
// actually exist -- min(n_pca_dims, n_components, n_anchors) -- which is the
// closest match to Julia's "use everything the sample provides".
//
// The same audit found `compute_ncv_embedding` passing graph_k straight into
// knn_parallel: graph_k <= 0 made knn_parallel return no results and the
// interpolation loop read past the empty index lists. It now clamps k to >= 1,
// mirroring `fit_ncv_interpolation_model`'s `interp_k = max(1, graph_k)`.
//
// The CLI tests spawn the instrumented `baysor` binary (BAYSOR_CLI_PATH,
// injected by the BAYSOR_WITH_TESTS CMake block) as a subprocess and assert a
// clean exit code 0 instead of a crash (SIGABRT 134 in Debug, silent garbage
// colours in Release).

#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <set>
#include <sstream>
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

// The exact reproducer shape: 20 components, only 9 anchor molecules, default
// n_pca_dims=10. Before the fix n_pca = min(10, 20) = 10 > 9 thin-U columns
// => Eigen assertion (Debug) / out-of-bounds read (Release).
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

// The streaming entry point used by `baysor run` / `baysor preview`: the basis
// fit selects all 9 molecules as anchors, so fit_ncv_interpolation_model used
// to call leftCols(10) on a 9-column thin U (the reported crash site at
// color_utils.cpp:649).
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
// n_neighbors=1 and the PCA basis is clamped to 2 columns. Before the fix
// this asserted leftCols(10) on 2 columns.
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

// graph_k=0 used to produce k_interp=0; knn_parallel returns no results for
// k<=0 and the interpolation loop dereferenced the empty index lists (SIGSEGV).
// fit_ncv_interpolation_model already clamps its own k to >= 1; the direct
// interpolation path in compute_ncv_embedding now does the same.
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

#ifndef BAYSOR_CLI_PATH

TEST(Bug5_Cli, BaysorCliPathAvailable) {
    GTEST_SKIP() << "BAYSOR_CLI_PATH is not defined; CLI end-to-end tests are disabled";
}

#elif defined(_WIN32)

// The subprocess runner below shells out with sh-style quoting and decodes
// exit codes via WEXITSTATUS, so the end-to-end CLI tests are POSIX-only.
TEST(Bug5_Cli, SubprocessTestsArePosixOnly) {
    GTEST_SKIP() << "CLI subprocess tests require a POSIX shell and sys/wait.h";
}

#else  // BAYSOR_CLI_PATH && !defined(_WIN32)

#include <sys/wait.h>

namespace {

namespace fs = std::filesystem;

using TempDir = baysor_test::TempDir;

std::string read_text_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

struct CliResult {
    int exit_code = -1;  // -1 = process did not exit normally (e.g. signal)
    std::string out;     // stdout
    std::string err;     // stderr
};

CliResult run_cli(const TempDir& tmp, const std::string& args) {
    const fs::path out_p = tmp.path / "stdout.txt";
    const fs::path err_p = tmp.path / "stderr.txt";
    const std::string cmd = "'" + std::string(BAYSOR_CLI_PATH) + "' " + args +
                            " > '" + out_p.string() + "' 2> '" + err_p.string() + "'";
    const int status = std::system(cmd.c_str());

    CliResult r;
    if (status >= 0 && WIFEXITED(status)) {
        r.exit_code = WEXITSTATUS(status);
    }
    r.out = read_text_file(out_p);
    r.err = read_text_file(err_p);
    return r;
}

std::string write_text(const TempDir& tmp, const std::string& name,
                       const std::string& content) {
    const fs::path p = tmp.path / name;
    std::ofstream f(p);
    f << content;
    return p.string();
}

// The BUG-5 reproducer table: 10 molecules / 5 genes, exactly the rows that
// produced 9 NCV anchors (< n_pca_dims=10) and the Eigen assertion.
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

// The original reproducer: `baysor run tiny.csv -m 2 -s 2.5 -o out` exited 134
// (Eigen Block assertion in the NCV colour embedding). --iters 10 reaches the
// colour stage in well under a second while exercising the identical code
// path, so every cluster method stays cheap enough for the suite.
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
// (gene_composition_color_embedding) and used to hit the same assertion at
// compute_ncv_embedding's leftCols (color_utils.cpp:444).
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

#endif  // BAYSOR_CLI_PATH && !defined(_WIN32)
