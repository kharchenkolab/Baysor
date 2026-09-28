#include "baysor/utils/general.h"

#include <cassert>
#include <cmath>

namespace baysor {

namespace {
Xoshiro256pp& global_rng_storage() {
    static Xoshiro256pp rng(1);
    return rng;
}
}

// --- Counting ---

std::vector<int> count_array(const std::vector<int>& values, int max_value, bool drop_zero) {
    if (values.empty()) return {};

    if (max_value <= 0) {
        max_value = *std::max_element(values.begin(), values.end());
    }
    if (max_value <= 0) return {};

    std::vector<int> counts(max_value, 0);
    bool has_zero = false;
    for (int v : values) {
        if (v == 0) {
            has_zero = true;
            continue;
        }
        if (v > 0 && v <= max_value) {
            counts[v - 1]++;
        }
    }
    // Note: Julia's count_array warns on zero values if !drop_zero.
    // We skip that for now (caller handles).
    (void)has_zero;
    (void)drop_zero;
    return counts;
}

std::vector<double> count_array_weighted(const std::vector<int>& values,
                                          const std::vector<double>& weights,
                                          int max_value, bool drop_zero) {
    if (values.empty()) return {};

    if (max_value <= 0) {
        max_value = *std::max_element(values.begin(), values.end());
    }
    if (max_value <= 0) return {};

    std::vector<double> counts(max_value, 0.0);
    for (size_t i = 0; i < values.size(); ++i) {
        int v = values[i];
        if (v == 0) continue;
        if (v > 0 && v <= max_value) {
            counts[v - 1] += weights[i];
        }
    }
    return counts;
}

// --- Group-by / split ---

std::vector<std::vector<int>> split_ids(const std::vector<int>& factor,
                                         int max_factor, bool drop_zero) {
    if (factor.empty()) return {};

    if (max_factor <= 0) {
        max_factor = *std::max_element(factor.begin(), factor.end());
    }
    if (max_factor <= 0) return {};

    // First pass: count elements per group
    std::vector<int> counts(max_factor, 0);
    for (int f : factor) {
        if (f == 0 && drop_zero) continue;
        if (f > 0 && f <= max_factor) {
            counts[f - 1]++;
        }
    }

    // Allocate
    std::vector<std::vector<int>> result(max_factor);
    for (int k = 0; k < max_factor; ++k) {
        result[k].reserve(counts[k]);
    }

    // Second pass: fill
    for (int i = 0; i < static_cast<int>(factor.size()); ++i) {
        int f = factor[i];
        if (f == 0 && drop_zero) continue;
        if (f > 0 && f <= max_factor) {
            result[f - 1].push_back(i);
        }
    }

    return result;
}

// --- Weighted sampling ---

int fsample(const double* weights, int n, std::mt19937& rng) {
    if (n == 0) return -1;

    double total = 0.0;
    for (int i = 0; i < n; ++i) total += weights[i];

    std::uniform_real_distribution<double> dist(0.0, total);
    double t = dist(rng);

    double cw = weights[0];
    int i = 0;
    while (cw < t && i < n - 1) {
        ++i;
        cw += weights[i];
    }
    return i;
}

int fsample(const double* weights, int n, Xoshiro256pp& rng) {
    if (n == 0) return -1;

    double total = 0.0;
    for (int i = 0; i < n; ++i) total += weights[i];

    double t = rng.rand_float64() * total;

    double cw = weights[0];
    int i = 0;
    while (cw < t && i < n - 1) {
        ++i;
        cw += weights[i];
    }
    return i;
}

Xoshiro256pp& global_xoshiro_rng() {
    return global_rng_storage();
}

void reset_global_xoshiro_rng(std::uint64_t seed) {
    global_rng_storage().seed(seed);
}

// --- Weighted statistics ---

double wmean(const double* values, const double* weights, int n) {
    double s = 0.0, ws = 0.0;
    for (int i = 0; i < n; ++i) {
        s += values[i] * weights[i];
        ws += weights[i];
    }
    return ws > 0.0 ? s / ws : 0.0;
}

std::pair<double, double> wmean_std(const double* values, const double* weights, int n) {
    double w = 0.0, wx = 0.0, wxx = 0.0;
    for (int i = 0; i < n; ++i) {
        double wi = weights[i];
        double xi = values[i];
        w += wi;
        wx += wi * xi;
        wxx += wi * xi * xi;
    }
    if (w <= 0.0) return {0.0, 0.0};
    double m = wx / w;
    double s2 = (wxx / w) - m * m;
    return {m, std::sqrt(std::max(s2, 0.0))};
}

// --- Convergence ---

ConvergenceResult estimate_difference_l0(
    const double* m1, const double* m2, int rows, int cols,
    const double* col_weights, double change_threshold
) {
    double max_diff = 0.0;
    int n_changed = 0;

    for (int ci = 0; ci < cols; ++ci) {
        double c_max = 0.0;
        for (int ri = 0; ri < rows; ++ri) {
            int idx = ri + ci * rows; // column-major
            double d = std::abs(m1[idx] - m2[idx]);
            if (d > c_max) c_max = d;
        }
        if (col_weights) c_max *= col_weights[ci];
        if (c_max > change_threshold) n_changed++;
        if (c_max > max_diff) max_diff = c_max;
    }

    return {max_diff, cols > 0 ? static_cast<double>(n_changed) / cols : 0.0};
}

// --- String utilities ---

std::string strip(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

std::vector<std::string> split_string_list(const std::string& s, char sep) {
    std::vector<std::string> result;
    if (s.empty()) return result;

    size_t start = 0;
    while (start < s.size()) {
        size_t end = s.find(sep, start);
        if (end == std::string::npos) end = s.size();
        std::string token = strip(s.substr(start, end - start));
        if (!token.empty()) {
            result.push_back(std::move(token));
        }
        start = end + 1;
    }
    return result;
}

} // namespace baysor
