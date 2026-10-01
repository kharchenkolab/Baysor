#pragma once

// k-NN result set for nanoflann's KD-tree search, backed by a binary max-heap:
// O(log k) per accepted candidate instead of KNNResultSet's O(k) insertion
// shift, which dominates for the large NCV neighbourhoods (k = n_genes / 10).
//
// It keeps the neighbour set KNNResultSet keeps: a candidate is accepted iff
// the set is not full or it is strictly closer than the current worst, and on
// overflow the latest-inserted of the worst entries is evicted. worstDist()
// matches too, so the tree traversal and its pruning are unchanged.

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

    /// `buffer` must hold at least `capacity` entries and outlive the set.
    HeapKnnResultSet(std::size_t capacity, HeapKnnEntry* buffer)
        : heap_(buffer), capacity_(capacity) {}

    std::size_t size() const { return count_; }
    bool full() const { return count_ == capacity_; }

    double worstDist() const {
        return full() ? heap_[0].dist : std::numeric_limits<double>::max();
    }

    bool addPoint(double dist, int idx) {
        if (full()) {
            if (!(dist < heap_[0].dist)) return true;
            std::pop_heap(heap_, heap_ + count_, WorseFirst{});
            --count_;
        }
        heap_[count_++] = HeapKnnEntry{dist, seq_++, idx};
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
