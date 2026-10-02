// HeapKnnResultSet must return exactly the neighbours nanoflann::KNNResultSet
// returns, ordered by (distance, index), also with many tied distances.

#include <gtest/gtest.h>

#include "baysor/processing/data_processing/heap_knn_result_set.h"

#include <third_party/nanoflann.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <limits>
#include <random>
#include <vector>

namespace {

using baysor::HeapKnnEntry;
using baysor::HeapKnnResultSet;

struct ColMajorAdaptor {
    const Eigen::MatrixXd& mat;
    explicit ColMajorAdaptor(const Eigen::MatrixXd& m) : mat(m) {}
    size_t kdtree_get_point_count() const { return static_cast<size_t>(mat.cols()); }
    double kdtree_get_pt(size_t idx, size_t dim) const {
        return mat(static_cast<Eigen::Index>(dim), static_cast<Eigen::Index>(idx));
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const { return false; }
};

using Tree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, ColMajorAdaptor>, ColMajorAdaptor, -1, int>;

using Neighbours = std::vector<std::pair<double, int>>;  // (distance, index)

// Reference: KNNResultSet, then ordered by (distance, index).
Neighbours reference_knn(const Tree& tree, const double* query, int k) {
    std::vector<int> idx(k, -1);
    std::vector<double> dist(k, 0.0);
    nanoflann::KNNResultSet<double, int> rs(k);
    rs.init(idx.data(), dist.data());
    tree.findNeighbors(rs, query, nanoflann::SearchParameters(0.0f, true));
    Neighbours out;
    for (int j = 0; j < k; ++j) out.emplace_back(dist[j], idx[j]);
    std::sort(out.begin(), out.end());
    return out;
}

Neighbours heap_knn(const Tree& tree, const double* query, int k) {
    std::vector<HeapKnnEntry> buf(k);
    HeapKnnResultSet rs(static_cast<size_t>(k), buf.data());
    tree.findNeighbors(rs, query, nanoflann::SearchParameters(0.0f, true));
    EXPECT_TRUE(rs.full());
    std::vector<int> idx(k, -1);
    std::vector<double> dist(k, 0.0);
    rs.extract_sorted(idx.data(), dist.data());
    Neighbours out;
    for (int j = 0; j < k; ++j) out.emplace_back(dist[j], idx[j]);
    return out;
}

void expect_matches_reference(const Eigen::MatrixXd& pts, const std::vector<int>& ks,
                              const std::vector<int>& queries) {
    ColMajorAdaptor adaptor(pts);
    Tree tree(static_cast<int>(pts.rows()), adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    for (int k : ks) {
        for (int q : queries) {
            ASSERT_EQ(heap_knn(tree, pts.col(q).data(), k), reference_knn(tree, pts.col(q).data(), k))
                << "k=" << k << " query=" << q;
        }
    }
}

} // namespace

TEST(HeapKnnResultSet, MatchesKnnResultSetOnTiedGridPoints) {
    // Points on a small integer grid with duplicates: almost every distance is tied.
    std::mt19937 rng(7);
    std::vector<int> queries;
    for (int q = 0; q < 2000; q += 97) queries.push_back(q);
    for (int side : {3, 8, 30}) {
        SCOPED_TRACE(side);
        std::uniform_int_distribution<int> coord(0, side - 1);
        Eigen::MatrixXd pts(2, 2000);
        for (int i = 0; i < pts.cols(); ++i) {
            pts(0, i) = coord(rng);
            pts(1, i) = coord(rng);
        }
        expect_matches_reference(pts, {1, 2, 33, 100, 517, 2000}, queries);
    }
}

TEST(HeapKnnResultSet, MatchesKnnResultSetOnContinuousPoints) {
    std::mt19937 rng(11);
    std::normal_distribution<double> nd(0.0, 1.0);
    Eigen::MatrixXd pts(3, 3000);
    for (int i = 0; i < pts.cols(); ++i)
        for (int d = 0; d < 3; ++d) pts(d, i) = nd(rng);
    // A block of exact duplicates of one point.
    for (int i = 0; i < 50; ++i) pts.col(100 + i) = pts.col(5);
    expect_matches_reference(pts, {40, 300}, {5, 100, 149, 2999});
}

// Direct addPoint sequences, including candidates that are not better than
// the current worst (nanoflann's leaf loop compares against a worst distance
// cached at leaf entry, so such candidates do reach addPoint).
TEST(HeapKnnResultSet, AddPointKeepsTheSameSetAsKnnResultSet) {
    std::mt19937 rng(3);
    std::uniform_int_distribution<int> dist_val(0, 6);
    for (int trial = 0; trial < 200; ++trial) {
        const int k = 1 + trial % 9;
        const int n = 40;
        std::vector<int> ki(k);
        std::vector<double> kd(k);
        nanoflann::KNNResultSet<double, int> ref(k);
        ref.init(ki.data(), kd.data());
        std::vector<HeapKnnEntry> buf(k);
        HeapKnnResultSet heap(static_cast<size_t>(k), buf.data());
        EXPECT_EQ(heap.worstDist(), std::numeric_limits<double>::max());
        for (int i = 0; i < n; ++i) {
            const double d = dist_val(rng);
            ref.addPoint(d, i);
            heap.addPoint(d, i);
            ASSERT_EQ(ref.worstDist(), heap.worstDist()) << "trial=" << trial << " i=" << i;
            ASSERT_EQ(ref.full(), heap.full());
            ASSERT_EQ(ref.size(), heap.size());
        }
        std::vector<std::pair<double, int>> expected;
        for (int j = 0; j < k; ++j) expected.emplace_back(kd[j], ki[j]);
        std::sort(expected.begin(), expected.end());
        std::vector<int> hi(k);
        std::vector<double> hd(k);
        heap.extract_sorted(hi.data(), hd.data());
        for (int j = 0; j < k; ++j) {
            ASSERT_EQ(expected[j].first, hd[j]) << "trial=" << trial << " j=" << j;
            ASSERT_EQ(expected[j].second, hi[j]) << "trial=" << trial << " j=" << j;
        }
    }
}
