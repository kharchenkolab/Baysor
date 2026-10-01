// HeapKnnResultSet (include/baysor/processing/data_processing/heap_knn_result_set.h)
// must return exactly the neighbours nanoflann::KNNResultSet returns, in the
// (distance, index) order the NCV k-NN callers produce with a stable sort, on
// inputs with many tied distances (duplicate points, integer grids).

#include <gtest/gtest.h>

#include "baysor/processing/data_processing/heap_knn_result_set.h"

#include <third_party/nanoflann.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <limits>
#include <numeric>
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

// Reference: KNNResultSet followed by the callers' stable sort by (distance, index).
void reference_knn(const Tree& tree, const double* query, int k,
                   std::vector<int>& idx, std::vector<double>& dist) {
    idx.assign(k, -1);
    dist.assign(k, 0.0);
    nanoflann::KNNResultSet<double, int> rs(k);
    rs.init(idx.data(), dist.data());
    tree.findNeighbors(rs, query, nanoflann::SearchParameters(0.0f, true));
    std::vector<int> order(k);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
        if (dist[a] != dist[b]) return dist[a] < dist[b];
        return idx[a] < idx[b];
    });
    std::vector<int> si(k);
    std::vector<double> sd(k);
    for (int j = 0; j < k; ++j) { si[j] = idx[order[j]]; sd[j] = dist[order[j]]; }
    idx.swap(si);
    dist.swap(sd);
}

void heap_knn(const Tree& tree, const double* query, int k,
              std::vector<int>& idx, std::vector<double>& dist) {
    std::vector<HeapKnnEntry> buf(k);
    HeapKnnResultSet rs(static_cast<size_t>(k), buf.data());
    tree.findNeighbors(rs, query, nanoflann::SearchParameters(0.0f, true));
    ASSERT_TRUE(rs.full());
    idx.assign(k, -1);
    dist.assign(k, 0.0);
    rs.extract_sorted(idx.data(), dist.data());
}

// Points on a small integer grid with duplicates: almost every distance is tied.
Eigen::MatrixXd grid_points(int n, int side, std::mt19937& rng) {
    std::uniform_int_distribution<int> coord(0, side - 1);
    Eigen::MatrixXd pts(2, n);
    for (int i = 0; i < n; ++i) {
        pts(0, i) = coord(rng);
        pts(1, i) = coord(rng);
    }
    return pts;
}

} // namespace

TEST(HeapKnnResultSet, MatchesKnnResultSetOnTiedGridPoints) {
    std::mt19937 rng(7);
    for (int side : {3, 8, 30}) {
        Eigen::MatrixXd pts = grid_points(2000, side, rng);
        ColMajorAdaptor adaptor(pts);
        Tree tree(2, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
        for (int k : {1, 2, 33, 100, 517, 2000}) {
            for (int q = 0; q < pts.cols(); q += 97) {
                std::vector<int> ri, hi;
                std::vector<double> rd, hd;
                reference_knn(tree, pts.col(q).data(), k, ri, rd);
                heap_knn(tree, pts.col(q).data(), k, hi, hd);
                ASSERT_EQ(ri, hi) << "side=" << side << " k=" << k << " query=" << q;
                ASSERT_EQ(rd, hd) << "side=" << side << " k=" << k << " query=" << q;
            }
        }
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
    ColMajorAdaptor adaptor(pts);
    Tree tree(3, adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
    for (int k : {40, 300}) {
        for (int q : {5, 100, 149, 2999}) {
            std::vector<int> ri, hi;
            std::vector<double> rd, hd;
            reference_knn(tree, pts.col(q).data(), k, ri, rd);
            heap_knn(tree, pts.col(q).data(), k, hi, hd);
            ASSERT_EQ(ri, hi) << "k=" << k << " query=" << q;
            ASSERT_EQ(rd, hd) << "k=" << k << " query=" << q;
        }
    }
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
