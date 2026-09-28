// Coverage tests for molecule clustering (COV-2):
//   src/processing/bmm_algorithm/{molecule_clustering,
//                                molecule_clustering_louvain}.cpp
#include <gtest/gtest.h>

#include "baysor/processing/bmm_algorithm/molecule_clustering.h"

#include <spdlog/sinks/base_sink.h>
#include <spdlog/spdlog.h>

#include <Eigen/Dense>
#include <memory>
#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using baysor::AdjList;
using baysor::ClusterMethod;
using baysor::ClusteringOptions;

// Exception type deliberately not derived from std::exception, used to
// verify the catch (...) fallback in cluster_molecules_ica and the unwind
// cleanup paths of verbose log statements.
struct cov2_nonstd_exception {};

class cov2_capturing_sink : public spdlog::sinks::base_sink<std::mutex> {
protected:
    void sink_it_(const spdlog::details::log_msg& msg) override {
        data_.append(msg.payload.data(), msg.payload.size());
        data_.push_back('\n');
    }
    void flush_() override {}

public:
    std::string data() const { return data_; }

private:
    std::string data_;
};

// Throws a non-std exception for every payload containing `trigger`; all
// other messages are captured normally. Because spdlog only swallows
// std::exception from sinks, the non-std exception propagates out of the
// logging call — exercising the exception-cleanup line-table entries of the
// corresponding spdlog::info statements.
class cov2_throwing_sink : public spdlog::sinks::base_sink<std::mutex> {
public:
    explicit cov2_throwing_sink(std::string trigger)
        : trigger_(std::move(trigger)) {}

protected:
    void sink_it_(const spdlog::details::log_msg& msg) override {
        std::string payload(msg.payload.data(), msg.payload.size());
        if (payload.find(trigger_) != std::string::npos) {
            throw cov2_nonstd_exception{};
        }
        data_.append(payload);
        data_.push_back('\n');
    }
    void flush_() override {}

public:
    std::string data() const { return data_; }

private:
    std::string trigger_;
    std::string data_;
};

// Swaps the global default logger for the duration of a test and restores it
// afterwards (even when expectations fail).
struct cov2_logger_guard {
    std::shared_ptr<spdlog::logger> original;
    explicit cov2_logger_guard(std::shared_ptr<spdlog::sinks::sink> sink)
        : original(spdlog::default_logger()) {
        auto logger = std::make_shared<spdlog::logger>("cov2-clustering", sink);
        logger->set_level(spdlog::level::trace);
        sink->set_level(spdlog::level::trace);
        spdlog::set_default_logger(logger);
    }
    ~cov2_logger_guard() { spdlog::set_default_logger(original); }
};

AdjList cov2_chain_adj(int n) {
    std::vector<int> src, dst;
    std::vector<double> wts;
    for (int i = 0; i + 1 < n; ++i) {
        src.push_back(i);
        dst.push_back(i + 1);
        wts.push_back(1.0 + 0.2 * (i % 3));
    }
    return AdjList::from_edge_list(
        src.data(), dst.data(), wts.data(), static_cast<int>(src.size()), n);
}

AdjList cov2_two_clique_adj() {
    const int edge_src[] = {0, 1, 0, 3, 4, 3, 2};
    const int edge_dst[] = {1, 2, 2, 4, 5, 5, 3};
    const double edge_wt[] = {3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.05};
    return AdjList::from_edge_list(edge_src, edge_dst, edge_wt, 7, 6);
}

// Three spatial patches of 20 molecules with patch-specific gene mixtures.
struct cov2_patch_data {
    Eigen::MatrixXd pos;               // 2 x 60
    std::vector<int> genes;            // 1-based, max 4
    std::vector<double> confidence;    // all 1.0
    AdjList adj;                       // chain over all molecules
    int n = 60;
};

cov2_patch_data cov2_make_patch_data() {
    cov2_patch_data data;
    data.pos.resize(2, data.n);
    for (int i = 0; i < data.n; ++i) {
        const int patch = i / 20;
        const int local = i % 20;
        data.pos(0, i) = 5.0 * patch + 0.2 * ((local % 4) - 1.5);
        data.pos(1, i) = 0.2 * ((local / 4) - 2.0) + 0.05 * (local % 3);
        data.genes.push_back((i % 7 == 0) ? 4 : patch + 1);
        data.confidence.push_back(1.0);
    }
    data.adj = cov2_chain_adj(data.n);
    return data;
}

int cov2_assert_labels_in_range(const std::vector<int>& labels, int lo, int hi) {
    std::set<int> unique(labels.begin(), labels.end());
    for (int v : unique) {
        EXPECT_GE(v, lo);
        EXPECT_LE(v, hi);
    }
    return static_cast<int>(unique.size());
}

} // namespace

// ============================================================================
// cluster_molecules_on_mrf
// ============================================================================

TEST(Cov2Clust, MrfAutoIterationsHashInitConvergesVerbosely) {
    const std::vector<int> genes = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2};
    const std::vector<double> confidence(10, 1.0);
    auto adj = cov2_chain_adj(10);

    // max_iters = -1 selects the automatic budget; tol = 1.0 guarantees the
    // convergence check (first evaluated after 20+ iterations) triggers.
    auto result = baysor::cluster_molecules_on_mrf(
        genes, adj, confidence, /*n_clusters=*/2,
        /*tol=*/1.0, /*mrf_weight=*/1.0, /*max_iters=*/-1,
        /*verbose=*/true, /*exprs_init=*/nullptr);

    ASSERT_EQ(result.assignment.size(), 10u);
    cov2_assert_labels_in_range(result.assignment, 1, 2);
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
    auto adj = cov2_chain_adj(4);

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

// ============================================================================
// cluster_molecules_ica failure paths
// ============================================================================

TEST(Cov2Clust, IcaNonStdExceptionFromLoggerTriggersCatchAllFallback) {
    auto sink = std::make_shared<cov2_throwing_sink>(
        "ICA initialization succeeded");
    cov2_logger_guard guard(sink);

    const std::vector<int> genes = {1, 2, 3, 3, 4, 1, 4, 4, 3, 1, 3, 2};
    const std::vector<double> confidence(12, 1.0);
    const int edge_src[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    const int edge_dst[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    const double edge_wt[] = {1.0, 1.3, 1.6, 1.9, 1.0, 1.3, 1.6, 1.9, 1.0, 1.3, 1.6};
    auto adj = AdjList::from_edge_list(edge_src, edge_dst, edge_wt, 11, 12);

    // The sink throws a non-std exception on the success log line inside the
    // try block; the catch (...) handler must log its fallback warning and
    // the function must still return a valid clustering.
    auto result = baysor::cluster_molecules_ica(
        genes, adj, confidence, /*n_clusters=*/2,
        /*tol=*/0.5, /*mrf_weight=*/1.0, /*max_iters=*/50,
        /*verbose=*/true);

    ASSERT_EQ(result.assignment.size(), 12u);
    cov2_assert_labels_in_range(result.assignment, 1, 2);
    EXPECT_NE(sink->data().find("ICA failed, falling back to hash initialization"),
              std::string::npos);
}

// The verbose EM-loop log statements are annotated with exception-cleanup
// entries: when the sink throws a non-std exception from inside spdlog, the
// statement's cleanup runs and the exception escapes the clustering call.
TEST(Cov2Clust, MrfVerboseIterationLogPropagatesSinkFailure) {
    auto sink = std::make_shared<cov2_throwing_sink>("Clustering iter");
    cov2_logger_guard guard(sink);

    const std::vector<int> genes = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2};
    const std::vector<double> confidence(10, 1.0);
    auto adj = cov2_chain_adj(10);

    EXPECT_THROW(
        (void)baysor::cluster_molecules_on_mrf(
            genes, adj, confidence, /*n_clusters=*/2,
            /*tol=*/1.0, /*mrf_weight=*/1.0, /*max_iters=*/-1,
            /*verbose=*/true, /*exprs_init=*/nullptr),
        cov2_nonstd_exception);
}

TEST(Cov2Clust, MrfConvergenceLogPropagatesSinkFailure) {
    auto sink = std::make_shared<cov2_throwing_sink>("Clustering converged");
    cov2_logger_guard guard(sink);

    const std::vector<int> genes = {1, 1, 1, 1, 1, 2, 2, 2, 2, 2};
    const std::vector<double> confidence(10, 1.0);
    auto adj = cov2_chain_adj(10);

    // Runs the full EM loop to the convergence branch (tol = 1.0 guarantees
    // the check fires after 20+ iterations), then throws from the log line.
    EXPECT_THROW(
        (void)baysor::cluster_molecules_on_mrf(
            genes, adj, confidence, /*n_clusters=*/2,
            /*tol=*/1.0, /*mrf_weight=*/1.0, /*max_iters=*/-1,
            /*verbose=*/true, /*exprs_init=*/nullptr),
        cov2_nonstd_exception);
}

TEST(Cov2Clust, LouvainBackendAnchorLogPropagatesSinkFailure) {
    // "clustering: using" only appears on the backend's per-run anchor log
    // (the NCV fitter logs a different "NCV basis anchors" line).
    auto sink = std::make_shared<cov2_throwing_sink>("clustering: using");
    cov2_logger_guard guard(sink);

    auto data = cov2_make_patch_data();
    auto empty_adj = AdjList::from_edge_list(nullptr, nullptr, nullptr, 0, data.n);

    EXPECT_THROW(
        (void)baysor::cluster_molecules_louvain(
            data.pos, data.genes, empty_adj, data.confidence,
            /*resolution=*/1.0, /*graph_k=*/4, /*spatial_k=*/0,
            /*target_clusters=*/3, /*n_dims=*/20,
            /*basis_sample_size=*/100000, /*verbose=*/true),
        cov2_nonstd_exception);
}

// ============================================================================
// cluster_molecules dispatcher
// ============================================================================

TEST(Cov2Clust, DispatcherRoutesEveryMethodAndRejectsUnknownValues) {
    auto data = cov2_make_patch_data();

    // MRF / ICA path.
    ClusteringOptions mrf_opts;
    mrf_opts.method = ClusterMethod::Mrf;
    mrf_opts.n_clusters = 2;
    mrf_opts.tol = 0.5;
    auto mrf = baysor::cluster_molecules(
        data.pos, data.genes, data.adj, data.confidence, mrf_opts, /*verbose=*/false);
    ASSERT_EQ(mrf.assignment.size(), static_cast<size_t>(data.n));
    cov2_assert_labels_in_range(mrf.assignment, 1, 2);

    // n_clusters <= 1 short-circuits to an empty result.
    mrf_opts.n_clusters = 1;
    auto mrf_trivial = baysor::cluster_molecules(
        data.pos, data.genes, data.adj, data.confidence, mrf_opts, /*verbose=*/false);
    EXPECT_TRUE(mrf_trivial.assignment.empty());
    EXPECT_TRUE(mrf_trivial.exprs.size() == 0);

    // Louvain path.
    ClusteringOptions louvain_opts;
    louvain_opts.method = ClusterMethod::Louvain;
    louvain_opts.n_clusters = 3;
    louvain_opts.graph_k = 4;
    auto louvain = baysor::cluster_molecules(
        data.pos, data.genes, data.adj, data.confidence, louvain_opts, /*verbose=*/false);
    ASSERT_EQ(louvain.assignment.size(), static_cast<size_t>(data.n));
    cov2_assert_labels_in_range(louvain.assignment, 1, 3);
    ASSERT_NE(louvain.ncv_projected_model, nullptr);

    // Leiden path.
    ClusteringOptions leiden_opts;
    leiden_opts.method = ClusterMethod::Leiden;
    leiden_opts.n_clusters = 3;
    leiden_opts.graph_k = 4;
    auto leiden = baysor::cluster_molecules(
        data.pos, data.genes, data.adj, data.confidence, leiden_opts, /*verbose=*/false);
    ASSERT_EQ(leiden.assignment.size(), static_cast<size_t>(data.n));
    cov2_assert_labels_in_range(leiden.assignment, 1, 3);
    ASSERT_NE(leiden.ncv_projected_model, nullptr);

    // Out-of-range enum values fall through the switch to an empty result.
    ClusteringOptions unknown_opts;
    unknown_opts.method = static_cast<ClusterMethod>(99);
    auto unknown = baysor::cluster_molecules(
        data.pos, data.genes, data.adj, data.confidence, unknown_opts, /*verbose=*/false);
    EXPECT_TRUE(unknown.assignment.empty());
    EXPECT_TRUE(unknown.diffs.empty());
}

// ============================================================================
// Louvain / Leiden graph backend
// ============================================================================

TEST(Cov2Clust, LouvainBackendLogsAnchorDiagnostics) {
    auto sink = std::make_shared<cov2_capturing_sink>();
    cov2_logger_guard guard(sink);

    auto data = cov2_make_patch_data();
    auto empty_adj = AdjList::from_edge_list(nullptr, nullptr, nullptr, 0, data.n);

    auto result = baysor::cluster_molecules_louvain(
        data.pos, data.genes, empty_adj, data.confidence,
        /*resolution=*/1.0, /*graph_k=*/4, /*spatial_k=*/0,
        /*target_clusters=*/3, /*n_dims=*/20, /*basis_sample_size=*/100000,
        /*verbose=*/true);

    ASSERT_EQ(result.assignment.size(), static_cast<size_t>(data.n));
    cov2_assert_labels_in_range(result.assignment, 1, 3);
    ASSERT_NE(result.ncv_projected_model, nullptr);
    EXPECT_FALSE(result.diffs.empty());  // per-level move fractions

    EXPECT_NE(sink->data().find("Louvain clustering: using"), std::string::npos);
    EXPECT_NE(sink->data().find("Louvain clustering complete"), std::string::npos);
}

TEST(Cov2Clust, LeidenBackendWithZeroGraphKWarnsAboutIsolatedAnchors) {
    auto sink = std::make_shared<cov2_capturing_sink>();
    cov2_logger_guard guard(sink);

    auto data = cov2_make_patch_data();
    auto empty_adj = AdjList::from_edge_list(nullptr, nullptr, nullptr, 0, data.n);

    // graph_k = 0 makes build_knn_similarity_graph return an edgeless graph,
    // so every anchor is isolated and the backend logs its warning.
    auto result = baysor::cluster_molecules_leiden(
        data.pos, data.genes, empty_adj, data.confidence,
        /*resolution=*/1.0, /*graph_k=*/0, /*spatial_k=*/0,
        /*target_clusters=*/3, /*n_dims=*/20, /*basis_sample_size=*/100000,
        /*verbose=*/true);

    ASSERT_EQ(result.assignment.size(), static_cast<size_t>(data.n));
    cov2_assert_labels_in_range(result.assignment, 1, 3);
    EXPECT_NE(sink->data().find("isolated anchors"), std::string::npos);
    EXPECT_NE(sink->data().find("Leiden clustering complete"), std::string::npos);
}

// ============================================================================
// graph_partition_to_target edge cases
// ============================================================================

TEST(Cov2Clust, GraphPartitionToTargetEdgeCases) {
    Eigen::MatrixXf vecs(2, 6);
    for (int i = 0; i < 6; ++i) {
        vecs(0, i) = static_cast<float>(i);
        vecs(1, i) = 1.0f;
    }
    const std::vector<double> confidence(6, 1.0);
    auto adj6 = cov2_chain_adj(6);

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
    EXPECT_EQ(cov2_assert_labels_in_range(two, 1, 2), 2);
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
    auto adj = cov2_two_clique_adj();

    auto membership = baysor::leiden_partition(adj, /*resolution=*/1.0,
                                               /*max_passes=*/100);

    ASSERT_EQ(membership.size(), 6u);
    EXPECT_EQ(membership[0], membership[1]);
    EXPECT_EQ(membership[1], membership[2]);
    EXPECT_EQ(membership[3], membership[4]);
    EXPECT_EQ(membership[4], membership[5]);
    EXPECT_NE(membership[2], membership[3]);
}
