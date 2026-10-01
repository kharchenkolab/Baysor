// ICA whitening for large gene panels (molecule_clustering.cpp): above
// ica_dense_whitening_max_genes the gene co-occurrence matrix is built sparse
// and the top-k covariance eigenpairs come from a truncated SVD (irlba)
// instead of the dense eigen-decomposition. Both must give the same
// eigenvalues and eigenvectors up to sign; the truncated vectors use a fixed
// sign convention (largest-magnitude entry positive).

#include <gtest/gtest.h>

#include "baysor/processing/bmm_algorithm/molecule_clustering.h"
#include "baysor/processing/models/adj_list.h"
#include "baysor/reporting/color_utils.h"
#include "baysor/utils/thread_pool.h"

#include "test_cov_helpers.h"

#include <cmath>
#include <random>
#include <string>
#include <vector>

namespace {

using baysor_test::CapturingSink;
using baysor_test::LoggerGuard;

class PoolSizeGuard {
public:
    explicit PoolSizeGuard(int n) : old_(baysor::thread_pool_size()) {
        baysor::set_thread_pool_size(n);
    }
    ~PoolSizeGuard() { baysor::set_thread_pool_size(old_); }

private:
    int old_;
};

// Non-negative, mostly sparse "co-occurrence-like" matrix with a few dominant,
// well separated directions.
Eigen::MatrixXd structured_matrix(int n, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    Eigen::MatrixXd x = Eigen::MatrixXd::Zero(n, n);
    const int n_blocks = 6;
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c < n; ++c) {
            const bool same_block = (r * n_blocks / n) == (c * n_blocks / n);
            const double scale = same_block ? 1.0 + (r * n_blocks / n) : 0.05;
            if (unif(rng) < (same_block ? 0.3 : 0.02)) x(r, c) = scale * unif(rng);
        }
    }
    return x;
}

void expect_same_eigenpairs(const baysor::detail::IcaWhitening& dense,
                            const baysor::detail::IcaWhitening& trunc) {
    ASSERT_EQ(dense.eigenvalues.size(), trunc.eigenvalues.size());
    ASSERT_EQ(dense.eigenvectors.rows(), trunc.eigenvectors.rows());
    for (Eigen::Index i = 0; i < dense.eigenvalues.size(); ++i) {
        EXPECT_NEAR(trunc.eigenvalues(i), dense.eigenvalues(i), 1e-9 * dense.eigenvalues(0)) << i;
        const double dot = trunc.eigenvectors.col(i).dot(dense.eigenvectors.col(i));
        EXPECT_NEAR(std::abs(dot), 1.0, 1e-8) << i;
        // Sign convention of the truncated solver.
        Eigen::Index imax = 0;
        trunc.eigenvectors.col(i).cwiseAbs().maxCoeff(&imax);
        EXPECT_GT(trunc.eigenvectors(imax, i), 0.0) << i;
    }
}

// Molecules on a grid, two spatial domains with different gene sets.
struct IcaInput {
    std::vector<int> genes;
    std::vector<double> confidence;
    baysor::AdjList adj;
};

IcaInput make_domains(int side, int n_genes, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    IcaInput in;
    const int n = side * side;
    in.genes.resize(n);
    in.confidence.assign(n, 1.0);
    for (int y = 0; y < side; ++y) {
        for (int x = 0; x < side; ++x) {
            const int domain = (x < side / 2 ? 0 : 1) + (y < side / 2 ? 0 : 2);
            int g = 1 + domain * (n_genes / 4) + static_cast<int>(unif(rng) * (n_genes / 4));
            if (unif(rng) < 0.1) g = 1 + static_cast<int>(unif(rng) * n_genes);
            in.genes[y * side + x] = std::min(g, n_genes);
        }
    }
    in.genes[0] = n_genes;
    std::vector<int> src, dst;
    std::vector<double> wts;
    for (int y = 0; y < side; ++y) {
        for (int x = 0; x < side; ++x) {
            const int i = y * side + x;
            if (x + 1 < side) { src.push_back(i); dst.push_back(i + 1); wts.push_back(1.0); }
            if (y + 1 < side) { src.push_back(i); dst.push_back(i + side); wts.push_back(1.0); }
        }
    }
    in.adj = baysor::AdjList::from_edge_list(
        src.data(), dst.data(), wts.data(), static_cast<int>(src.size()), n);
    return in;
}

} // namespace

TEST(IcaTruncatedWhitening, MatchesDenseEigenDecomposition) {
    const Eigen::MatrixXd x = structured_matrix(300, 1);
    for (int k : {1, 4, 6}) {
        SCOPED_TRACE(k);
        expect_same_eigenpairs(baysor::detail::ica_whitening_dense(x, k),
                               baysor::detail::ica_whitening_truncated(x.sparseView(), k));
    }
}

TEST(IcaTruncatedWhitening, SmallMatricesUseTheExactFallback) {
    // irlba switches to an exact SVD when 2k >= min(rows, cols).
    const Eigen::MatrixXd x = structured_matrix(6, 2);
    expect_same_eigenpairs(baysor::detail::ica_whitening_dense(x, 3),
                           baysor::detail::ica_whitening_truncated(x.sparseView(), 3));
    EXPECT_THROW(baysor::detail::ica_whitening_truncated(x.sparseView(), 7), std::invalid_argument);
}

TEST(IcaTruncatedWhitening, ClusteringUsesTruncatedPathAboveThreshold) {
    const IcaInput in = make_domains(60, 40, 4);
    auto sink = std::make_shared<CapturingSink>();
    baysor::ClusteringResult ref;
    {
        LoggerGuard guard(sink);
        PoolSizeGuard pool(1);
        ref = baysor::cluster_molecules_ica(in.genes, in.adj, in.confidence, 4, 0.01, 1.0, 300,
                                            /*verbose=*/true, /*dense_whitening_max_genes=*/39);
    }
    const std::string logs = sink->data();
    EXPECT_NE(logs.find("ICA whitening: truncated eigen-decomposition (40 genes > 39, "), std::string::npos) << logs;
    EXPECT_NE(logs.find("ICA initialization succeeded (4 components)"), std::string::npos) << logs;
    EXPECT_EQ(logs.find("falling back"), std::string::npos) << logs;
    ASSERT_EQ(ref.assignment.size(), in.genes.size());

    // At or below the threshold: dense path, no truncated-whitening message.
    sink->clear();
    {
        LoggerGuard guard(sink);
        baysor::cluster_molecules_ica(in.genes, in.adj, in.confidence, 4, 0.01, 1.0, 5,
                                      /*verbose=*/true, /*dense_whitening_max_genes=*/40);
    }
    EXPECT_EQ(sink->data().find("truncated eigen-decomposition"), std::string::npos) << sink->data();

    // The truncated path is deterministic and independent of the thread count.
    for (int threads : {3, 8}) {
        PoolSizeGuard pool(threads);
        auto res = baysor::cluster_molecules_ica(in.genes, in.adj, in.confidence, 4, 0.01, 1.0, 300,
                                                 /*verbose=*/false, /*dense_whitening_max_genes=*/39);
        EXPECT_EQ(res.assignment, ref.assignment) << threads;
        EXPECT_EQ(res.diffs, ref.diffs) << threads;
    }
}

TEST(IcaTruncatedWhitening, TooFewGenesStillFallsBackToHashInit) {
    const IcaInput in = make_domains(10, 3, 5);
    auto sink = std::make_shared<CapturingSink>();
    LoggerGuard guard(sink);
    baysor::cluster_molecules_ica(in.genes, in.adj, in.confidence, 4, 0.01, 1.0, 10,
                                  /*verbose=*/true, /*dense_whitening_max_genes=*/0);
    const std::string logs = sink->data();
    EXPECT_NE(logs.find("falling back to hash initialization"), std::string::npos) << logs;
    EXPECT_NE(logs.find("k must not exceed min(m, n)"), std::string::npos) << logs;
}

TEST(IcaTruncatedWhitening, SparseCoOccurrenceMatchesDenseBuilder) {
    IcaInput in = make_domains(50, 60, 8);
    std::mt19937 rng(2);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    for (double& c : in.confidence) c = unif(rng) < 0.2 ? 0.5 : 0.97;  // threshold 0.95
    const Eigen::MatrixXd dense = baysor::pairwise_gene_spatial_cor(in.genes, in.confidence, in.adj);
    for (int threads : {1, 4}) {
        PoolSizeGuard pool(threads);
        const Eigen::SparseMatrix<double> sparse =
            baysor::detail::sparse_gene_spatial_cor(in.genes, in.confidence, in.adj);
        ASSERT_EQ(sparse.rows(), dense.rows());
        ASSERT_EQ(sparse.cols(), dense.cols());
        const Eigen::MatrixXd as_dense = Eigen::MatrixXd(sparse);
        int nnz_dense = 0;
        for (Eigen::Index c = 0; c < dense.cols(); ++c) {
            for (Eigen::Index r = 0; r < dense.rows(); ++r) {
                if (dense(r, c) != 0.0) ++nnz_dense;
                EXPECT_EQ(dense(r, c) == 0.0, as_dense(r, c) == 0.0) << r << "," << c;
                EXPECT_NEAR(as_dense(r, c), dense(r, c), 1e-14 * std::abs(dense(r, c))) << r << "," << c;
            }
        }
        EXPECT_EQ(sparse.nonZeros(), nnz_dense);
        EXPECT_GT(nnz_dense, 0);
    }
    // Whitening of both matrices agrees.
    const Eigen::SparseMatrix<double> sparse =
        baysor::detail::sparse_gene_spatial_cor(in.genes, in.confidence, in.adj);
    expect_same_eigenpairs(baysor::detail::ica_whitening_dense(dense, 4),
                           baysor::detail::ica_whitening_truncated(sparse, 4));
}
