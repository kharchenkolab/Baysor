// Determinism of the molecule graph (DET-RUNID): CGAL emits triangulation
// edges in an order that depends on the heap addresses of its faces, which
// shift with the length of the `-o` path and with allocator history. The
// adjacency list sorts the edges, so the graph (and the MRF summation order
// it sets) is a function of the edge set.

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "baysor/processing/data_processing/triangulation.h"

#include "test_cov_helpers.h"

namespace {

// Deterministic pseudo-random sequence in [0, 1) (no global RNG involved).
struct Lcg {
    std::uint64_t s;
    double operator()() {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return static_cast<double>((s >> 16) & 0xFFFFFF) / 16777216.0;
    }
};

Eigen::MatrixXd make_det_points(int n) {
    Eigen::MatrixXd pts(2, n);
    Lcg next{12345};
    for (int i = 0; i < n; ++i) {
        pts(0, i) = next() * 100.0;
        pts(1, i) = next() * 100.0;
    }
    return pts;
}

}  // namespace

TEST(DetGraph, AdjacencyEdgesAreCanonicallySorted) {
    const auto res = baysor::adjacency_list(make_det_points(600), /*filter=*/false);
    ASSERT_GT(res.edge_src.size(), 0u);

    std::vector<std::pair<int, int>> edges;
    for (size_t i = 0; i < res.edge_src.size(); ++i) edges.push_back(std::minmax(res.edge_src[i], res.edge_dst[i]));
    EXPECT_EQ(std::adjacent_find(edges.begin(), edges.end(), std::greater_equal<>()), edges.end())
        << "edges not strictly increasing";
}

TEST(DetGraph, AdjacencyOrderStableUnderHeapPerturbation) {
    const Eigen::MatrixXd pts = make_det_points(600);
    const auto ref = baysor::adjacency_list(pts, /*filter=*/false);
    ASSERT_GT(ref.edge_src.size(), 0u);

    // Hold differently sized chunks (spanning the CGAL block allocations,
    // ~1 KiB upwards) so that the triangulation's blocks land elsewhere.
    for (int perturb = 0; perturb < 6; ++perturb) {
        SCOPED_TRACE(perturb);
        std::vector<std::unique_ptr<char[]>> held;
        for (int k = 0; k <= perturb; ++k) {
            const size_t sz = 1024 + static_cast<size_t>(k) * 1536 + 64;
            held.push_back(std::make_unique<char[]>(sz));
            std::memset(held.back().get(), 0x5A, sz);
        }
        const auto got = baysor::adjacency_list(pts, /*filter=*/false);
        EXPECT_EQ(got.edge_src, ref.edge_src);
        EXPECT_EQ(got.edge_dst, ref.edge_dst);
        EXPECT_EQ(got.edge_dists, ref.edge_dists);
    }
}

#if !defined(_WIN32) && defined(BAYSOR_CLI_PATH)

// End to end: one thread, two output paths of different lengths.
TEST(DetCli, SameAssignmentsForDifferentOutputPathLengthsAtOneThread) {
    baysor_test::TempDir tmp("det_runid_pathlen");
    std::ostringstream csv;
    csv << "x,y,gene,cell_id\n";
    const char* genes[] = {"GeneA", "GeneB", "GeneC", "GeneD"};
    Lcg next{42};
    for (int idx = 0; idx < 400; ++idx) {
        const int c = idx / 100;
        const double x = 10.0 + 20.0 * (c % 2) + (next() * 4.0 - 2.0);
        const double y = 10.0 + 20.0 * (c / 2) + (next() * 4.0 - 2.0);
        csv << x << "," << y << "," << genes[idx % 4] << ","
            << (idx % 100 >= 5 ? ("cell" + std::to_string(c + 1)) : "0") << "\n";
    }
    const std::string csv_path = baysor_test::cli::write_text(tmp, "mols.csv", csv.str());

    std::vector<std::string> segmentations;
    for (const char* dir : {"o1", "o1_padding_16_chars_here"}) {
        const auto out = tmp.path / dir / "seg";
        const auto r = baysor_test::cli::run_cli(
            tmp, "run '" + csv_path + "' -m 10 -s 2.5 --iters 30 --cluster-method none -t 1 -o '" + out.string() + "'");
        ASSERT_EQ(r.exit_code, 0) << r.out << r.err;
        segmentations.push_back(baysor_test::cli::read_text_file(out / "segmentation.csv"));
    }
    ASSERT_GT(segmentations[0].size(), 0u);
    EXPECT_EQ(segmentations[0], segmentations[1]) << "segmentation depends on the output path length at 1 thread";
}

#else

TEST(DetCli, SameAssignmentsForDifferentOutputPathLengthsAtOneThread) {
    GTEST_SKIP() << "CLI subprocess tests require POSIX and BAYSOR_CLI_PATH";
}

#endif
