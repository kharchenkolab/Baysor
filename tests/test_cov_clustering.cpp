// Molecule clustering: the MRF EM and its ICA initialisation
// (molecule_clustering*.cpp), the Louvain/Leiden graph backends and the
// cluster_molecules dispatcher.
#include <gtest/gtest.h>

#include "baysor/processing/bmm_algorithm/molecule_clustering.h"
#include "baysor/reporting/color_utils.h"
#include "baysor/utils/thread_pool.h"

#include "test_cov_helpers.h"

#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <vector>

namespace {

using baysor::AdjList;
using baysor::ClusterMethod;
using baysor::ClusteringOptions;
using baysor_test::CapturingSink;
using baysor_test::LoggerGuard;

class PoolSizeGuard {
public:
    explicit PoolSizeGuard(int n) : old_(baysor::thread_pool_size()) { baysor::set_thread_pool_size(n); }
    ~PoolSizeGuard() { baysor::set_thread_pool_size(old_); }

private:
    int old_;
};

// Chain 0-1-...-(n-1) with edge weights 1 + weight_step * (i % 3).
AdjList chain_adj(int n, double weight_step = 0.0) {
    std::vector<int> src, dst;
    std::vector<double> wts;
    for (int i = 0; i + 1 < n; ++i) {
        src.push_back(i);
        dst.push_back(i + 1);
        wts.push_back(1.0 + weight_step * (i % 3));
    }
    return AdjList::from_edge_list(src.data(), dst.data(), wts.data(), static_cast<int>(src.size()), n);
}

// 4-neighbour grid over side x side molecules (row-major ids).
AdjList grid_adj(int side, const std::function<double()>& weight) {
    const int n = side * side;
    std::vector<int> src, dst;
    std::vector<double> wts;
    for (int i = 0; i < n; ++i) {
        if (i % side + 1 < side) { src.push_back(i); dst.push_back(i + 1); wts.push_back(weight()); }
        if (i + side < n) { src.push_back(i); dst.push_back(i + side); wts.push_back(weight()); }
    }
    return AdjList::from_edge_list(src.data(), dst.data(), wts.data(), static_cast<int>(src.size()), n);
}

struct GridInput {
    std::vector<int> genes;  // 1-based; 0 = unknown gene
    std::vector<double> confidence;
    AdjList adj;
};

// Two spatial domains with different gene preferences, a few unknown genes,
// random confidences and edge weights.
GridInput make_mrf_grid(int side, int n_genes, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    GridInput in;
    for (int i = 0; i < side * side; ++i) {
        int g = 1 + static_cast<int>(unif(rng) * (n_genes / 2));
        if (i % side >= side / 2) g += n_genes / 2;
        if (unif(rng) < 0.2) g = 1 + static_cast<int>(unif(rng) * n_genes);
        if (unif(rng) < 0.01) g = 0;
        in.genes.push_back(std::min(g, n_genes));
        in.confidence.push_back(0.3 + 0.7 * unif(rng));
    }
    in.genes[0] = n_genes;  // make sure the largest gene id is present
    in.adj = grid_adj(side, [&] { return 0.5 + unif(rng); });
    return in;
}

// Four spatial domains with their own gene sets, unit weights and confidences.
GridInput make_ica_domains(int side, int n_genes, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    GridInput in;
    for (int i = 0; i < side * side; ++i) {
        const int domain = (i % side < side / 2 ? 0 : 1) + (i / side < side / 2 ? 0 : 2);
        int g = 1 + domain * (n_genes / 4) + static_cast<int>(unif(rng) * (n_genes / 4));
        if (unif(rng) < 0.1) g = 1 + static_cast<int>(unif(rng) * n_genes);
        in.genes.push_back(std::min(g, n_genes));
    }
    in.genes[0] = n_genes;
    in.confidence.assign(in.genes.size(), 1.0);
    in.adj = grid_adj(side, [] { return 1.0; });
    return in;
}

// Three spatial patches of 20 molecules with patch-specific gene mixtures.
struct PatchData {
    Eigen::MatrixXd pos;               // 2 x 60
    std::vector<int> genes;            // 1-based, max 4
    std::vector<double> confidence;    // all 1.0
    AdjList adj;                       // chain over all molecules
    int n = 60;
};

PatchData make_patch_data() {
    PatchData data;
    data.pos.resize(2, data.n);
    for (int i = 0; i < data.n; ++i) {
        const int patch = i / 20;
        const int local = i % 20;
        data.pos(0, i) = 5.0 * patch + 0.2 * ((local % 4) - 1.5);
        data.pos(1, i) = 0.2 * ((local / 4) - 2.0) + 0.05 * (local % 3);
        data.genes.push_back((i % 7 == 0) ? 4 : patch + 1);
        data.confidence.push_back(1.0);
    }
    data.adj = chain_adj(data.n, 0.2);
    return data;
}

int assert_labels_in_range(const std::vector<int>& labels, int lo, int hi) {
    std::set<int> unique(labels.begin(), labels.end());
    for (int v : unique) {
        EXPECT_GE(v, lo);
        EXPECT_LE(v, hi);
    }
    return static_cast<int>(unique.size());
}

void expect_bitwise_equal(const baysor::ClusteringResult& a, const baysor::ClusteringResult& b) {
    ASSERT_EQ(a.assignment, b.assignment);
    ASSERT_EQ(a.diffs, b.diffs);
    ASSERT_EQ(a.change_fracs, b.change_fracs);
    ASSERT_EQ(a.assignment_probs.rows(), b.assignment_probs.rows());
    ASSERT_EQ(a.assignment_probs.cols(), b.assignment_probs.cols());
    for (Eigen::Index i = 0; i < a.assignment_probs.size(); ++i)
        ASSERT_EQ(a.assignment_probs.data()[i], b.assignment_probs.data()[i]) << "prob " << i;
    ASSERT_EQ(a.exprs.rows(), b.exprs.rows());
    ASSERT_EQ(a.exprs.cols(), b.exprs.cols());
    for (Eigen::Index i = 0; i < a.exprs.size(); ++i)
        ASSERT_EQ(a.exprs.data()[i], b.exprs.data()[i]) << "expr " << i;
}

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

// Same eigenvalues and eigenvectors up to sign; the truncated vectors have a
// positive largest-magnitude entry.
void expect_same_eigenpairs(const baysor::detail::IcaWhitening& dense,
                            const baysor::detail::IcaWhitening& trunc) {
    ASSERT_EQ(dense.eigenvalues.size(), trunc.eigenvalues.size());
    ASSERT_EQ(dense.eigenvectors.rows(), trunc.eigenvectors.rows());
    for (Eigen::Index i = 0; i < dense.eigenvalues.size(); ++i) {
        EXPECT_NEAR(trunc.eigenvalues(i), dense.eigenvalues(i), 1e-9 * dense.eigenvalues(0)) << i;
        const double dot = trunc.eigenvectors.col(i).dot(dense.eigenvectors.col(i));
        EXPECT_NEAR(std::abs(dot), 1.0, 1e-8) << i;
        Eigen::Index imax = 0;
        trunc.eigenvectors.col(i).cwiseAbs().maxCoeff(&imax);
        EXPECT_GT(trunc.eigenvectors(imax, i), 0.0) << i;
    }
}

// 12 molecules on a chain cycling through exactly 3 genes.
const std::vector<int> kThreeGenes = {1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3};

} // namespace

// ============================================================================
// cluster_molecules_on_mrf
// ============================================================================

TEST(Cov2Clust, MrfAutoIterationsHashInitConvergesVerbosely) {
    const std::vector<int> genes = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2};
    const std::vector<double> confidence(10, 1.0);
    auto adj = chain_adj(10, 0.2);

    // max_iters = -1 selects the automatic budget; tol = 1.0 guarantees the
    // convergence check (first evaluated after 20+ iterations) triggers.
    auto result = baysor::cluster_molecules_on_mrf(
        genes, adj, confidence, /*n_clusters=*/2,
        /*tol=*/1.0, /*mrf_weight=*/1.0, /*max_iters=*/-1,
        /*verbose=*/true, /*exprs_init=*/nullptr);

    ASSERT_EQ(result.assignment.size(), 10u);
    assert_labels_in_range(result.assignment, 1, 2);
    ASSERT_EQ(result.diffs.size(), 21u);       // converged at iteration 21
    EXPECT_LT(result.diffs.back(), 1.0);
    ASSERT_EQ(result.change_fracs.size(), 21u);
    ASSERT_EQ(result.exprs.rows(), 2);
    ASSERT_EQ(result.exprs.cols(), 2);
    for (int k = 0; k < result.exprs.rows(); ++k) {
        EXPECT_NEAR(result.exprs.row(k).sum(), 1.0, 1e-9);
    }
    ASSERT_EQ(result.assignment_probs.rows(), 2);
    EXPECT_EQ(result.assignment_probs.cols(), 10);
}

TEST(Cov2Clust, MrfZeroInitGeneColumnProducesUniformProbabilities) {
    const std::vector<int> genes = {1, 1, 2, 2};
    const std::vector<double> confidence(4, 1.0);
    auto adj = chain_adj(4, 0.2);

    // Gene column 1 is all zeros: molecules carrying gene 1 get a uniform
    // initialization (col_sum == 0 branch) and later win the argmax tie by
    // landing in cluster 1, while gene 2 molecules follow the 0.9 profile.
    Eigen::MatrixXd exprs_init(2, 2);
    exprs_init <<
        0.0, 0.1,
        0.0, 0.9;

    auto result = baysor::cluster_molecules_on_mrf(
        genes, adj, confidence, /*n_clusters=*/2,
        /*tol=*/0.0, /*mrf_weight=*/0.0, /*max_iters=*/1,
        /*verbose=*/false, &exprs_init);

    EXPECT_EQ(result.assignment, (std::vector<int>{1, 1, 2, 2}));
    ASSERT_EQ(result.diffs.size(), 1u);
    // Column-normalized probabilities for the zero-profile gene are uniform.
    EXPECT_NEAR(result.assignment_probs(0, 0), 0.5, 1e-12);
    EXPECT_NEAR(result.assignment_probs(1, 0), 0.5, 1e-12);
}

// The E-step, the convergence statistics and the M-step run in parallel and
// must give bitwise identical results at any thread count.
TEST(MrfClusteringParallel, ResultsDoNotDependOnThreadCount) {
    // 80 x 80 = 6,400 molecules: 13 E-step chunks of 512 molecules.
    const GridInput in = make_mrf_grid(80, 40, 17);
    baysor::ClusteringResult ref;
    {
        PoolSizeGuard pool(1);
        ref = baysor::cluster_molecules_on_mrf(in.genes, in.adj, in.confidence, 3, 0.01, 1.0, 300, false);
    }
    ASSERT_FALSE(ref.diffs.empty());
    ASSERT_EQ(ref.diffs.size(), ref.change_fracs.size());
    for (int threads : {2, 3, 8}) {
        PoolSizeGuard pool(threads);
        auto res = baysor::cluster_molecules_on_mrf(in.genes, in.adj, in.confidence, 3, 0.01, 1.0, 300, false);
        SCOPED_TRACE(threads);
        expect_bitwise_equal(ref, res);
    }
}

TEST(MrfClusteringParallel, ConvergenceTraceIsConsistent) {
    const GridInput in = make_mrf_grid(40, 20, 3);
    PoolSizeGuard pool(4);
    auto res = baysor::cluster_molecules_on_mrf(in.genes, in.adj, in.confidence, 2, 0.01, 1.0, 500, false);
    ASSERT_FALSE(res.diffs.empty());
    for (size_t t = 0; t < res.diffs.size(); ++t) {
        EXPECT_GE(res.diffs[t], 0.0);
        EXPECT_LE(res.diffs[t], 1.0);
        EXPECT_GE(res.change_fracs[t], 0.0);
        EXPECT_LE(res.change_fracs[t], 1.0);
        // A molecule counts as changed only above 1e-7, so no change at all
        // means the maximum is at most 1e-7.
        if (res.change_fracs[t] == 0.0) EXPECT_LE(res.diffs[t], 1e-7);
    }
    // Converged: the last 21 maxima are all below tol.
    if (static_cast<int>(res.diffs.size()) < 500) {
        for (size_t t = res.diffs.size() - 21; t < res.diffs.size(); ++t)
            EXPECT_LT(res.diffs[t], 0.01);
    }
}

// The final expression profiles are the M-step without pseudocount over the
// returned probabilities; the parallel per-gene sums must equal the
// sequential molecule-order accumulation bitwise.
TEST(MrfClusteringParallel, FinalExpressionProfilesMatchSequentialMStep) {
    const GridInput in = make_mrf_grid(60, 30, 23);
    PoolSizeGuard pool(5);
    auto res = baysor::cluster_molecules_on_mrf(in.genes, in.adj, in.confidence, 4, 0.01, 1.0, 200, false);
    const auto& probs = res.assignment_probs;
    const int n_clusters = static_cast<int>(probs.rows());
    Eigen::MatrixXd exprs = Eigen::MatrixXd::Zero(n_clusters, res.exprs.cols());
    for (int i = 0; i < static_cast<int>(in.genes.size()); ++i) {
        const int g0 = in.genes[i] - 1;
        if (g0 < 0) continue;
        for (int k = 0; k < n_clusters; ++k) exprs(k, g0) += in.confidence[i] * probs(k, i);
    }
    for (int k = 0; k < n_clusters; ++k) {
        const double row_sum = exprs.row(k).sum();
        if (row_sum > 0) exprs.row(k) /= row_sum;
    }
    for (Eigen::Index i = 0; i < exprs.size(); ++i)
        ASSERT_EQ(exprs.data()[i], res.exprs.data()[i]) << "expr " << i;
}

// ============================================================================
// cluster_molecules_ica
// ============================================================================

// Like Julia's fit(ICA, X, k), the ICA rejects k > min(m, n) (here 4 clusters
// for 3 genes) on the dense and on the truncated whitening path, and the
// clustering falls back to the hash initialisation.
TEST(Bug2IcaFallback, FewerGenesThanClustersFallsBackToHashInit) {
    const std::vector<double> confidence(kThreeGenes.size(), 1.0);
    auto adj = chain_adj(static_cast<int>(kThreeGenes.size()));
    const auto expected = baysor::cluster_molecules_on_mrf(
        kThreeGenes, adj, confidence, 4, /*tol=*/0.0, 1.0, /*max_iters=*/10, /*verbose=*/false, nullptr);

    for (int dense_max_genes : {baysor::ica_dense_whitening_max_genes, 0}) {
        SCOPED_TRACE(dense_max_genes);
        auto sink = std::make_shared<CapturingSink>();
        LoggerGuard guard(sink);
        auto result = baysor::cluster_molecules_ica(
            kThreeGenes, adj, confidence, /*n_clusters=*/4,
            /*tol=*/0.0, /*mrf_weight=*/1.0, /*max_iters=*/10, /*verbose=*/true, dense_max_genes);

        const std::string logs = sink->data();
        EXPECT_NE(logs.find("falling back to hash initialization"), std::string::npos) << logs;
        EXPECT_NE(logs.find("k must not exceed min(m, n)"), std::string::npos) << logs;
        EXPECT_EQ(logs.find("ICA initialization succeeded"), std::string::npos) << logs;

        ASSERT_EQ(result.assignment.size(), kThreeGenes.size());
        assert_labels_in_range(result.assignment, 1, 4);
        ASSERT_EQ(result.exprs.rows(), 4);
        ASSERT_EQ(result.exprs.cols(), 3);
        EXPECT_TRUE(result.exprs.allFinite());
        for (int k = 0; k < 4; ++k) EXPECT_NEAR(result.exprs.row(k).sum(), 1.0, 1e-9);
        ASSERT_EQ(result.assignment_probs.rows(), 4);
        ASSERT_EQ(result.assignment_probs.cols(), static_cast<int>(kThreeGenes.size()));
        EXPECT_EQ(result.assignment, expected.assignment);
        EXPECT_NEAR((result.exprs - expected.exprs).cwiseAbs().maxCoeff(), 0.0, 1e-12);
    }
}

// n_clusters == n_genes is allowed (only k > min(m, n) is rejected).
TEST(Bug2IcaFallback, EqualGenesAndClustersStillUsesIca) {
    auto sink = std::make_shared<CapturingSink>();
    LoggerGuard guard(sink);
    const std::vector<double> confidence(kThreeGenes.size(), 1.0);

    auto result = baysor::cluster_molecules_ica(
        kThreeGenes, chain_adj(static_cast<int>(kThreeGenes.size())), confidence, /*n_clusters=*/3,
        /*tol=*/0.0, /*mrf_weight=*/1.0, /*max_iters=*/10, /*verbose=*/true);

    const std::string logs = sink->data();
    EXPECT_NE(logs.find("ICA initialization succeeded (3 components)"), std::string::npos) << logs;
    EXPECT_EQ(logs.find("falling back to hash initialization"), std::string::npos) << logs;
    ASSERT_EQ(result.assignment.size(), kThreeGenes.size());
    ASSERT_EQ(result.exprs.rows(), 3);
    ASSERT_EQ(result.exprs.cols(), 3);
    ASSERT_EQ(result.assignment_probs.rows(), 3);
    ASSERT_EQ(result.assignment_probs.cols(), static_cast<int>(kThreeGenes.size()));
}

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
    const GridInput in = make_ica_domains(60, 40, 4);
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

TEST(IcaTruncatedWhitening, SparseCoOccurrenceMatchesDenseBuilder) {
    GridInput in = make_ica_domains(50, 60, 8);
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

// ============================================================================
// cluster_molecules dispatcher
// ============================================================================

TEST(Cov2Clust, DispatcherRoutesEveryMethodAndRejectsUnknownValues) {
    auto data = make_patch_data();
    auto run = [&](ClusterMethod method, int n_clusters, int graph_k = 15) {
        ClusteringOptions opts;
        opts.method = method;
        opts.n_clusters = n_clusters;
        opts.graph_k = graph_k;
        opts.tol = 0.5;
        return baysor::cluster_molecules(data.pos, data.genes, data.adj, data.confidence, opts, /*verbose=*/false);
    };

    // MRF / ICA path.
    auto mrf = run(ClusterMethod::Mrf, 2);
    ASSERT_EQ(mrf.assignment.size(), static_cast<size_t>(data.n));
    assert_labels_in_range(mrf.assignment, 1, 2);

    // n_clusters <= 1 short-circuits to an empty result.
    auto mrf_trivial = run(ClusterMethod::Mrf, 1);
    EXPECT_TRUE(mrf_trivial.assignment.empty());
    EXPECT_EQ(mrf_trivial.exprs.size(), 0);

    for (ClusterMethod method : {ClusterMethod::Louvain, ClusterMethod::Leiden}) {
        auto graph = run(method, 3, /*graph_k=*/4);
        ASSERT_EQ(graph.assignment.size(), static_cast<size_t>(data.n));
        assert_labels_in_range(graph.assignment, 1, 3);
        ASSERT_NE(graph.ncv_projected_model, nullptr);
    }

    // Out-of-range enum values fall through the switch to an empty result.
    auto unknown = run(static_cast<ClusterMethod>(99), 4);
    EXPECT_TRUE(unknown.assignment.empty());
    EXPECT_TRUE(unknown.diffs.empty());
}

// ============================================================================
// Louvain / Leiden graph backend
// ============================================================================

TEST(Cov2Clust, LouvainBackendLogsAnchorDiagnostics) {
    auto sink = std::make_shared<CapturingSink>();
    LoggerGuard guard(sink);

    auto data = make_patch_data();
    auto empty_adj = AdjList::from_edge_list(nullptr, nullptr, nullptr, 0, data.n);

    auto result = baysor::cluster_molecules_louvain(
        data.pos, data.genes, empty_adj, data.confidence,
        /*resolution=*/1.0, /*graph_k=*/4, /*spatial_k=*/0,
        /*target_clusters=*/3, /*n_dims=*/20, /*basis_sample_size=*/100000,
        /*verbose=*/true);

    ASSERT_EQ(result.assignment.size(), static_cast<size_t>(data.n));
    assert_labels_in_range(result.assignment, 1, 3);
    ASSERT_NE(result.ncv_projected_model, nullptr);
    EXPECT_FALSE(result.diffs.empty());  // per-level move fractions

    EXPECT_NE(sink->data().find("Louvain clustering: using"), std::string::npos);
    EXPECT_NE(sink->data().find("Louvain clustering complete"), std::string::npos);
}

TEST(Cov2Clust, LeidenBackendWithZeroGraphKWarnsAboutIsolatedAnchors) {
    auto sink = std::make_shared<CapturingSink>();
    LoggerGuard guard(sink);

    auto data = make_patch_data();
    auto empty_adj = AdjList::from_edge_list(nullptr, nullptr, nullptr, 0, data.n);

    // graph_k = 0 makes build_knn_similarity_graph return an edgeless graph,
    // so every anchor is isolated and the backend logs its warning.
    auto result = baysor::cluster_molecules_leiden(
        data.pos, data.genes, empty_adj, data.confidence,
        /*resolution=*/1.0, /*graph_k=*/0, /*spatial_k=*/0,
        /*target_clusters=*/3, /*n_dims=*/20, /*basis_sample_size=*/100000,
        /*verbose=*/true);

    ASSERT_EQ(result.assignment.size(), static_cast<size_t>(data.n));
    assert_labels_in_range(result.assignment, 1, 3);
    EXPECT_NE(sink->data().find("isolated anchors"), std::string::npos);
    EXPECT_NE(sink->data().find("Leiden clustering complete"), std::string::npos);
}

TEST(Cov2Clust, GraphPartitionToTargetEdgeCases) {
    Eigen::MatrixXf vecs(2, 6);
    for (int i = 0; i < 6; ++i) {
        vecs(0, i) = static_cast<float>(i);
        vecs(1, i) = 1.0f;
    }
    const std::vector<double> confidence(6, 1.0);
    auto adj6 = chain_adj(6, 0.2);

    // target_clusters <= 1: everything is a single cluster and the summary
    // records the seed resolution.
    baysor::GraphClusteringSummary summary_one;
    auto all_one = baysor::graph_partition_to_target(
        adj6, vecs, confidence, ClusterMethod::Louvain,
        /*target_clusters=*/1, /*resolution_seed=*/1.0, /*max_passes=*/100,
        &summary_one);
    EXPECT_EQ(all_one, (std::vector<int>(6, 1)));
    EXPECT_EQ(summary_one.micro_clusters, 1);
    EXPECT_EQ(summary_one.final_clusters, 1);
    EXPECT_DOUBLE_EQ(summary_one.chosen_resolution, 1.0);

    // Non-graph methods are rejected up front.
    auto mrf_labels = baysor::graph_partition_to_target(
        adj6, vecs, confidence, ClusterMethod::Mrf,
        /*target_clusters=*/3, /*resolution_seed=*/1.0, /*max_passes=*/100,
        nullptr);
    EXPECT_TRUE(mrf_labels.empty());

    // No resolution reaches the micro target (12) on a six-node chain, but
    // several reach the final target of 2: second-choice fallback path.
    baysor::GraphClusteringSummary summary_two;
    auto two = baysor::graph_partition_to_target(
        adj6, vecs, confidence, ClusterMethod::Louvain,
        /*target_clusters=*/2, /*resolution_seed=*/1.0, /*max_passes=*/100,
        &summary_two);
    EXPECT_EQ(assert_labels_in_range(two, 1, 2), 2);
    EXPECT_EQ(summary_two.final_clusters, 2);
    EXPECT_GE(summary_two.micro_clusters, 2);

    // A single-node graph can never reach the requested count: the final
    // fallback keeps the attempt with the most communities (one).
    auto adj1 = AdjList::from_edge_list(nullptr, nullptr, nullptr, 0, 1);
    Eigen::MatrixXf vec1(2, 1);
    vec1 << 1.0f, 0.0f;
    const std::vector<double> conf1(1, 1.0);
    baysor::GraphClusteringSummary summary_solo;
    auto solo = baysor::graph_partition_to_target(
        adj1, vec1, conf1, ClusterMethod::Louvain,
        /*target_clusters=*/2, /*resolution_seed=*/1.0, /*max_passes=*/100,
        &summary_solo);
    EXPECT_EQ(solo, (std::vector<int>{1}));
    EXPECT_EQ(summary_solo.micro_clusters, 1);
    EXPECT_EQ(summary_solo.final_clusters, 1);

    // Empty graph short-circuits before any partitioning.
    auto empty_graph = AdjList::from_edge_list(nullptr, nullptr, nullptr, 0, 0);
    auto none = baysor::graph_partition_to_target(
        empty_graph, vecs, confidence, ClusterMethod::Louvain,
        /*target_clusters=*/3, /*resolution_seed=*/1.0, /*max_passes=*/100,
        nullptr);
    EXPECT_TRUE(none.empty());
}

TEST(Cov2Clust, LeidenPartitionSeparatesWeaklyConnectedCliques) {
    const int edge_src[] = {0, 1, 0, 3, 4, 3, 2};
    const int edge_dst[] = {1, 2, 2, 4, 5, 5, 3};
    const double edge_wt[] = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.05};
    auto adj = AdjList::from_edge_list(edge_src, edge_dst, edge_wt, 7, 6);

    auto membership = baysor::leiden_partition(adj, /*resolution=*/1.0, /*max_passes=*/100);

    ASSERT_EQ(membership.size(), 6u);
    EXPECT_EQ(membership[0], membership[1]);
    EXPECT_EQ(membership[1], membership[2]);
    EXPECT_EQ(membership[3], membership[4]);
    EXPECT_EQ(membership[4], membership[5]);
    EXPECT_NE(membership[2], membership[3]);
}

// ============================================================================
// CLI: a 3-gene dataset with the default --n-clusters (4) runs to completion
// on the hash-initialisation fallback.
// ============================================================================

#if !defined(_WIN32) && defined(BAYSOR_CLI_PATH)

TEST(Bug2IcaFallback, CliThreeGeneDatasetWithDefaultOptionsExitsZero) {
    baysor_test::TempDir tmp("bug2_cli_3gene");
    // 4 clumps x 30 molecules, exactly 3 genes (cycling), no prior column.
    std::ostringstream csv;
    csv << "x,y,gene\n";
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> jit(-2.0, 2.0);
    for (int idx = 0; idx < 120; ++idx) {
        const int c = idx / 30;
        csv << (10.0 + 20.0 * (c % 2) + jit(rng)) << "," << (10.0 + 20.0 * (c / 2) + jit(rng)) << ","
            << "Gene" << char('A' + idx % 3) << "\n";
    }
    const std::string csv_path = baysor_test::cli::write_text(tmp, "mols_3gene.csv", csv.str());
    const auto out = tmp.path / "seg";

    auto r = baysor_test::cli::run_cli(tmp, "run '" + csv_path + "' -m 10 -s 2.5 -o '" + out.string() + "'");
    EXPECT_EQ(r.exit_code, 0) << "--- stdout ---\n" << r.out << "\n--- stderr ---\n" << r.err;

    const std::string logs = r.out + r.err;
    EXPECT_NE(logs.find("Clustering molecules into 4 types"), std::string::npos) << logs;
    EXPECT_NE(logs.find("falling back to hash initialization"), std::string::npos) << logs;
    EXPECT_NE(logs.find("Segmentation complete"), std::string::npos) << logs;

    const auto seg_csv = out / "segmentation.csv";
    ASSERT_TRUE(std::filesystem::is_regular_file(seg_csv)) << seg_csv;
    EXPECT_GT(std::filesystem::file_size(seg_csv), 0u);
}

#else

TEST(Bug2IcaFallback, CliThreeGeneDatasetWithDefaultOptionsExitsZero) {
    GTEST_SKIP() << "CLI subprocess tests require POSIX and BAYSOR_CLI_PATH";
}

#endif
