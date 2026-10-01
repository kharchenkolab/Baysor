#pragma once

// k-NN result set for nanoflann's KD-tree search, backed by a binary max-heap.
//
// nanoflann::KNNResultSet keeps the k best candidates in a sorted array and
// inserts by shifting, O(k) per accepted candidate. The KD-tree visits
// candidates roughly nearest-first, so while the set fills almost every
// candidate is accepted and the search costs ~k^2 per query. With the NCV
// neighbourhood size k = n_genes / 10 (hundreds to thousands of neighbours on
// large panels) that shift dominated the whole clustering phase. The heap
// makes each accepted candidate O(log k).
//
// The result is the same neighbour *set* KNNResultSet returns:
//  - a candidate is accepted iff the set is not full, or it is strictly closer
//    than the current worst (KNNResultSet::addPoint drops a candidate whose
//    distance equals the worst one when the set is full);
//  - on overflow the worst entry is evicted and, among entries with equal
//    worst distance, the one inserted last (KNNResultSet keeps equal
//    distances in insertion order and drops its last slot).
// worstDist() also matches (the k-th best distance once full, the maximum
// representable distance before), so the tree traversal and its pruning are
// unchanged. extract_sorted() returns the neighbours ordered by
// (distance, index), the order the callers produce with a stable sort of
// KNNResultSet's output.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace baysor {

struct HeapKnnEntry {
    double dist;
    std::uint32_t seq;  // insertion order, breaks ties between equal distances
    int idx;
};

class HeapKnnResultSet {
public:
    using DistanceType = double;
    using IndexType = int;
    using CountType = std::size_t;

    /// `buffer` must hold at least `capacity` entries and outlive the set.
    HeapKnnResultSet(std::size_t capacity, HeapKnnEntry* buffer)
        : heap_(buffer), capacity_(capacity) {}

    std::size_t size() const { return count_; }
    bool empty() const { return count_ == 0; }
    bool full() const { return count_ == capacity_; }

    double worstDist() const {
        return full() ? heap_[0].dist : std::numeric_limits<double>::max();
    }

    bool addPoint(double dist, int idx) {
        if (count_ < capacity_) {
            heap_[count_++] = HeapKnnEntry{dist, seq_++, idx};
            std::push_heap(heap_, heap_ + count_, WorseFirst{});
            return true;
        }
        if (!(dist < heap_[0].dist)) return true;
        std::pop_heap(heap_, heap_ + count_, WorseFirst{});
        heap_[count_ - 1] = HeapKnnEntry{dist, seq_++, idx};
        std::push_heap(heap_, heap_ + count_, WorseFirst{});
        return true;
    }

    /// Called by nanoflann for sorted searches; ordering is done in extract_sorted().
    void sort() {}

    /// Writes the size() neighbours ordered by (distance, index). Destroys the
    /// heap order, so the set must not receive further points afterwards.
    void extract_sorted(int* idx_out, double* dist_out) {
        std::sort(heap_, heap_ + count_, [](const HeapKnnEntry& a, const HeapKnnEntry& b) {
            if (a.dist != b.dist) return a.dist < b.dist;
            return a.idx < b.idx;
        });
        for (std::size_t j = 0; j < count_; ++j) {
            idx_out[j] = heap_[j].idx;
            dist_out[j] = heap_[j].dist;
        }
    }

private:
    // Max-heap order: the top is the entry evicted first, i.e. the largest
    // distance and, among equal distances, the latest insertion.
    struct WorseFirst {
        bool operator()(const HeapKnnEntry& a, const HeapKnnEntry& b) const {
            return a.dist < b.dist || (a.dist == b.dist && a.seq < b.seq);
        }
    };

    HeapKnnEntry* heap_;
    std::size_t capacity_;
    std::size_t count_ = 0;
    std::uint32_t seq_ = 0;
};

} // namespace baysor
