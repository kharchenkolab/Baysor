#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

namespace baysor {

// Minimal Dict{Int,Float64}-like accumulator for Julia parity.
// It preserves Julia's integer hash, linear probing, and slot-order iteration.
// TODO(parity): This exists because the stochastic E-step is sensitive to
// candidate ordering. Revisit whether we still want to mirror Julia's Dict
// iteration order once parity is no longer the only goal.
class JuliaIntDoubleDict {
public:
    JuliaIntDoubleDict() { reset(); }

    /// Back to the state of a new dict (16 empty slots), keeping the capacity.
    void reset() { reset_storage(16); }

    void clear() {
        std::fill(slots_.begin(), slots_.end(), 0);
        count_ = 0;
        idxfloor_ = static_cast<int>(slots_.size());
    }

    void add(int key, double value) {
        for (;;) {
            const int sz = static_cast<int>(slots_.size());
            int index = hash_index(key, sz);

            for (int iter = 0; iter < sz; ++iter) {
                if (slots_[index] == 0) {
                    slots_[index] = short_hash(key);
                    keys_[index] = key;
                    vals_[index] = value;
                    ++count_;
                    if (index < idxfloor_) idxfloor_ = index;
                    return;
                }

                if (slots_[index] == short_hash(key) && keys_[index] == key) {
                    vals_[index] += value;
                    return;
                }

                index = (index + 1) & (sz - 1);
            }

            grow();
        }
    }

    template <class Fn>
    void for_each(Fn&& fn) const {
        const int sz = static_cast<int>(slots_.size());
        for (int i = idxfloor_; i < sz; ++i) {
            if ((slots_[i] & 0x80) != 0) {
                fn(keys_[i], vals_[i]);
            }
        }
    }

    int size() const { return count_; }

private:
    std::vector<std::uint8_t> slots_;
    std::vector<int> keys_;
    std::vector<double> vals_;
    int count_ = 0;
    int idxfloor_ = 0;

    static std::uint64_t hash_64_64(std::uint64_t n) {
        std::uint64_t a = n;
        a = ~a + (a << 21);
        a ^= (a >> 24);
        a = a + (a << 3) + (a << 8);
        a ^= (a >> 14);
        a = a + (a << 2) + (a << 4);
        a ^= (a >> 28);
        a = a + (a << 31);
        return a;
    }

    static std::uint64_t julia_hash_int(int key) {
        return hash_64_64(static_cast<std::uint64_t>(static_cast<std::int64_t>(key)));
    }

    static std::uint8_t short_hash(int key) {
        return static_cast<std::uint8_t>((julia_hash_int(key) >> 57) & 0x7fU) | 0x80U;
    }

    static int hash_index(int key, int sz) {
        return static_cast<int>(julia_hash_int(key) & static_cast<std::uint64_t>(sz - 1));
    }

    void reset_storage(int sz) {
        slots_.assign(sz, 0);
        keys_.assign(sz, 0);
        vals_.assign(sz, 0.0);
        count_ = 0;
        idxfloor_ = sz;
    }

    void grow() {
        const int old_sz = static_cast<int>(slots_.size());
        const int new_sz = (count_ > 64000) ? old_sz * 2 : old_sz * 4;
        auto old_slots = std::move(slots_);
        auto old_keys = std::move(keys_);
        auto old_vals = std::move(vals_);
        reset_storage(new_sz);
        for (int i = 0; i < old_sz; ++i) {
            if ((old_slots[i] & 0x80) != 0) {
                add(old_keys[i], old_vals[i]);
            }
        }
    }
};

} // namespace baysor
