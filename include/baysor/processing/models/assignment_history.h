#pragma once

#include "baysor/utils/thread_pool.h"

#include <cstddef>
#include <cstdint>
#include <deque>
#include <initializer_list>
#include <iterator>
#include <stdexcept>
#include <vector>

namespace baysor {

/// History of per-molecule assignments (global component GUIDs, 0 = noise),
/// oldest entry first, stored as the newest entry in full plus backward
/// deltas: for every pair of consecutive entries (t, t+1), the molecules
/// whose value differs, with their value in entry t, in ascending molecule
/// order.
///
/// Lossless. In the BMM loop about 15 % of the molecules change per
/// iteration, so a depth-50 history takes about 4 + 49 * 0.15 * 8 = 63 bytes
/// per molecule instead of 200 for 50 full copies.
///
/// The container interface (size, empty, operator[], front, back,
/// iteration, construction from rows) reconstructs full rows on demand; it
/// is meant for tests and small data. The BMM code reads the deltas.
class AssignmentHistory {
public:
    struct Change {
        std::int32_t mol;
        std::int32_t value;  // value in the older of the two entries
    };

    AssignmentHistory() = default;
    AssignmentHistory(std::initializer_list<std::vector<int>> rows);
    explicit AssignmentHistory(const std::vector<std::vector<int>>& rows);

    std::size_t size() const { return has_newest_ ? deltas_.size() + 1 : 0; }
    bool empty() const { return !has_newest_; }
    int n_molecules() const { return static_cast<int>(newest_.size()); }
    void clear();

    /// Entry t (0 = oldest), reconstructed.
    std::vector<int> operator[](std::size_t t) const;
    std::vector<int> front() const { return (*this)[0]; }
    const std::vector<int>& back() const { return newest_; }
    std::vector<std::vector<int>> rows() const;

    /// Newest entry (stored in full).
    const std::vector<int>& newest() const { return newest_; }
    /// Molecules that differ between entry t and t + 1 (t in [0, size() - 1)),
    /// with their values in entry t, ascending by molecule.
    const std::vector<Change>& changes_after(std::size_t t) const { return deltas_[t]; }

    /// Append a new newest entry. All rows must have the same length.
    void push_back(const std::vector<int>& row);

    /// Append a new newest entry with value `value_of(i)` for molecule i in
    /// [0, n). The change list is built in parallel over fixed blocks and
    /// concatenated in block order (the result does not depend on threads).
    template <class F>
    void push_back_generated(int n, F&& value_of);

    /// Drop the oldest entry.
    void pop_front();

    /// Heap bytes held (capacity of the stored arrays).
    std::size_t memory_bytes() const;

    class const_iterator {
    public:
        using iterator_category = std::input_iterator_tag;
        using value_type = std::vector<int>;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using reference = std::vector<int>;

        const_iterator(const AssignmentHistory* h, std::size_t t) : h_(h), t_(t) {}
        std::vector<int> operator*() const { return (*h_)[t_]; }
        const_iterator& operator++() { ++t_; return *this; }
        bool operator==(const const_iterator& o) const { return t_ == o.t_ && h_ == o.h_; }
        bool operator!=(const const_iterator& o) const { return !(*this == o); }

    private:
        const AssignmentHistory* h_;
        std::size_t t_;
    };
    const_iterator begin() const { return {this, 0}; }
    const_iterator end() const { return {this, size()}; }

private:
    std::vector<int> newest_;
    bool has_newest_ = false;
    std::deque<std::vector<Change>> deltas_;   // deltas_[t]: entry t vs t + 1
    std::vector<Change> spare_;                // recycled buffer of a dropped delta
    std::vector<std::vector<Change>> block_changes_;  // scratch of push_back_generated

    std::vector<Change> take_spare();
    [[noreturn]] static void throw_row_length_mismatch();
};

template <class F>
void AssignmentHistory::push_back_generated(int n, F&& value_of) {
    if (has_newest_ && n != n_molecules()) {
        throw_row_length_mismatch();
    }
    if (!has_newest_) {
        newest_.resize(n);
        parallel_for(0, n, 16384, [&](int i) { newest_[i] = value_of(i); });
        has_newest_ = true;
        return;
    }

    constexpr std::int64_t kBlock = 16384;
    const std::int64_t n_blocks = (static_cast<std::int64_t>(n) + kBlock - 1) / kBlock;
    if (static_cast<std::int64_t>(block_changes_.size()) < n_blocks) {
        block_changes_.resize(static_cast<std::size_t>(n_blocks));
    }
    run_parallel_chunks(0, n, kBlock, Scheduling::Dynamic,
        [&](std::int64_t b, std::int64_t e, int) {
        auto& out = block_changes_[static_cast<std::size_t>(b / kBlock)];
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

    std::size_t total = 0;
    for (std::int64_t k = 0; k < n_blocks; ++k) total += block_changes_[static_cast<std::size_t>(k)].size();
    std::vector<Change> delta = take_spare();
    delta.clear();
    delta.reserve(total);
    for (std::int64_t k = 0; k < n_blocks; ++k) {
        const auto& bc = block_changes_[static_cast<std::size_t>(k)];
        delta.insert(delta.end(), bc.begin(), bc.end());
    }
    deltas_.push_back(std::move(delta));
}

inline AssignmentHistory::AssignmentHistory(std::initializer_list<std::vector<int>> rows) {
    for (const auto& row : rows) push_back(row);
}

inline AssignmentHistory::AssignmentHistory(const std::vector<std::vector<int>>& rows) {
    for (const auto& row : rows) push_back(row);
}

inline void AssignmentHistory::clear() {
    newest_.clear();
    has_newest_ = false;
    deltas_.clear();
}

inline std::vector<int> AssignmentHistory::operator[](std::size_t t) const {
    if (t >= size()) {
        throw std::out_of_range("AssignmentHistory: entry index out of range");
    }
    std::vector<int> row = newest_;
    for (std::size_t k = deltas_.size(); k-- > t;) {
        for (const Change& c : deltas_[k]) row[c.mol] = c.value;
    }
    return row;
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

inline void AssignmentHistory::push_back(const std::vector<int>& row) {
    push_back_generated(static_cast<int>(row.size()), [&row](int i) { return row[i]; });
}

inline void AssignmentHistory::pop_front() {
    if (!has_newest_) return;
    if (deltas_.empty()) {
        clear();
        return;
    }
    spare_ = std::move(deltas_.front());
    deltas_.pop_front();
}

inline std::vector<AssignmentHistory::Change> AssignmentHistory::take_spare() {
    return std::move(spare_);
}

inline std::size_t AssignmentHistory::memory_bytes() const {
    std::size_t bytes = newest_.capacity() * sizeof(int) + spare_.capacity() * sizeof(Change);
    for (const auto& d : deltas_) bytes += d.capacity() * sizeof(Change);
    for (const auto& d : block_changes_) bytes += d.capacity() * sizeof(Change);
    return bytes;
}

inline void AssignmentHistory::throw_row_length_mismatch() {
    throw std::invalid_argument("AssignmentHistory: all entries must have the same length");
}

} // namespace baysor
