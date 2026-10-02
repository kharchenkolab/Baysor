#pragma once

#include "baysor/utils/thread_pool.h"

#include <cstddef>
#include <cstdint>
#include <deque>
#include <stdexcept>
#include <vector>

namespace baysor {

/// History of per-molecule assignments (global component GUIDs, 0 = noise),
/// oldest entry first, stored as the newest entry in full plus backward
/// deltas: for every pair of consecutive entries (t, t+1), the molecules
/// whose value differs, with their value in entry t, in ascending molecule
/// order. In the BMM loop about 15 % of the molecules change per iteration,
/// so this takes about a third of the memory of full copies.
class AssignmentHistory {
public:
    struct Change {
        std::int32_t mol;
        std::int32_t value;  // value in the older of the two entries
    };

    AssignmentHistory() = default;
    explicit AssignmentHistory(const std::vector<std::vector<int>>& rows) {
        for (const auto& row : rows) push_back(row);
    }

    std::size_t size() const { return has_newest_ ? deltas_.size() + 1 : 0; }
    bool empty() const { return !has_newest_; }
    /// Newest entry (stored in full).
    const std::vector<int>& back() const { return newest_; }
    /// Molecules that differ between entry t and t + 1 (t in [0, size() - 1)),
    /// with their values in entry t, ascending by molecule.
    const std::vector<Change>& changes_after(std::size_t t) const { return deltas_[t]; }
    /// All entries, reconstructed (for tests and small data).
    std::vector<std::vector<int>> rows() const;

    /// Append a new newest entry. All rows must have the same length.
    void push_back(const std::vector<int>& row);

    /// push_back of the entry with value `value_of(i)` for molecule i in
    /// [0, n), in three steps so that callers can place the parallel one in
    /// their own region: begin_push (serial; drops oldest entries so that at
    /// most `max_size` remain after the push), push_rows (work-shared; every
    /// participant must call it), end_push (serial). The result does not
    /// depend on the number of threads.
    void begin_push(int n, std::size_t max_size);
    template <class F>
    void push_rows(ParallelRegion& region, int n, F&& value_of);
    void end_push();

    /// Drop the oldest entry.
    void pop_front();

private:
    std::vector<int> newest_;
    bool has_newest_ = false;
    std::deque<std::vector<Change>> deltas_;   // deltas_[t]: entry t vs t + 1
    std::vector<Change> spare_;                // recycled buffer of a dropped delta
    std::vector<std::vector<Change>> block_changes_;  // per-block scratch of push_rows
    bool first_push_ = false;
    static constexpr std::int64_t kPushBlock = 16384;
};

inline void AssignmentHistory::begin_push(int n, std::size_t max_size) {
    if (has_newest_ && n != static_cast<int>(newest_.size())) {
        throw std::invalid_argument("AssignmentHistory: all entries must have the same length");
    }
    while (!empty() && size() >= max_size) pop_front();
    first_push_ = !has_newest_;
    if (first_push_) {
        newest_.resize(n);
    } else {
        const std::size_t n_blocks = (static_cast<std::size_t>(n) + kPushBlock - 1) / kPushBlock;
        if (block_changes_.size() < n_blocks) block_changes_.resize(n_blocks);
    }
}

template <class F>
void AssignmentHistory::push_rows(ParallelRegion& region, int n, F&& value_of) {
    if (first_push_) {
        region.for_each(0, n, kPushBlock, [&](std::int64_t i) {
            newest_[static_cast<std::size_t>(i)] = value_of(static_cast<int>(i));
        });
        return;
    }
    region.for_chunks(0, n, kPushBlock, Scheduling::Dynamic,
        [&](std::int64_t b, std::int64_t e, int) {
        auto& out = block_changes_[static_cast<std::size_t>(b / kPushBlock)];
        out.clear();
        for (std::int64_t i = b; i < e; ++i) {
            const int v = value_of(static_cast<int>(i));
            int& cur = newest_[static_cast<std::size_t>(i)];
            if (v != cur) {
                out.push_back({static_cast<std::int32_t>(i), cur});
                cur = v;
            }
        }
    });
}

inline void AssignmentHistory::end_push() {
    if (first_push_) {
        has_newest_ = true;
        return;
    }
    const std::size_t n_blocks = (newest_.size() + kPushBlock - 1) / kPushBlock;
    std::size_t total = 0;
    for (std::size_t k = 0; k < n_blocks; ++k) total += block_changes_[k].size();
    std::vector<Change> delta = std::move(spare_);
    delta.clear();
    delta.reserve(total);
    for (std::size_t k = 0; k < n_blocks; ++k) {
        delta.insert(delta.end(), block_changes_[k].begin(), block_changes_[k].end());
    }
    deltas_.push_back(std::move(delta));
}

inline void AssignmentHistory::push_back(const std::vector<int>& row) {
    const int n = static_cast<int>(row.size());
    begin_push(n, static_cast<std::size_t>(-1));
    parallel_region([&](ParallelRegion& region) {
        push_rows(region, n, [&row](int i) { return row[i]; });
    });
    end_push();
}

inline std::vector<std::vector<int>> AssignmentHistory::rows() const {
    std::vector<std::vector<int>> out(size());
    if (out.empty()) return out;
    out.back() = newest_;
    for (std::size_t k = deltas_.size(); k-- > 0;) {
        out[k] = out[k + 1];
        for (const Change& c : deltas_[k]) out[k][c.mol] = c.value;
    }
    return out;
}

inline void AssignmentHistory::pop_front() {
    if (deltas_.empty()) {
        newest_.clear();
        has_newest_ = false;
        return;
    }
    spare_ = std::move(deltas_.front());
    deltas_.pop_front();
}

} // namespace baysor
