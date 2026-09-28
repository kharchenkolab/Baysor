#include "baysor/processing/utils/convex_hull.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace baysor {

// Cross product z-component of vectors (a→b) and (a→c).
// > 0: counter-clockwise turn, < 0: clockwise, == 0: collinear
static double cross(double ax, double ay,
                    double bx, double by,
                    double cx, double cy) {
    return ax * (by - cy) + bx * (cy - ay) + cx * (ay - by);
}

// Andrew's monotone chain convex hull (2D).
// Input:  points  — 2 × n matrix (row 0 = x, row 1 = y)
// Output: 2 × hull_size matrix of hull vertices in clockwise order
//         (matching origin/master:src/processing/utils/convex_hull.jl; its
//         test in test/runtests.jl expects [0 0 1 1 0; 0 1 1 0 0] for the
//         unit square — i.e. the signed shoelace area is negative).
// Degenerate cases: 0 points → empty; 1–2 points → return as-is.
Eigen::MatrixXd convex_hull(const Eigen::MatrixXd& points) {
    int n = static_cast<int>(points.cols());
    if (n == 0) return {};
    if (n <= 2) return points;

    // Sort indices by (x ascending, then y ascending)
    std::vector<int> idx(n);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {
        double xa = points(0, a), xb = points(0, b);
        if (xa != xb) return xa < xb;
        return points(1, a) < points(1, b);
    });

    std::vector<int> hull;
    hull.reserve(2 * n);

    // Upper hull (clockwise turns)
    for (int i = 0; i < n; ++i) {
        int p = idx[i];
        while (hull.size() >= 2) {
            int q = hull[hull.size() - 1];
            int r = hull[hull.size() - 2];
            // Remove q if not a clockwise turn (cross >= 0 means left or collinear)
            if (cross(points(0, r), points(1, r),
                      points(0, q), points(1, q),
                      points(0, p), points(1, p)) >= 0) {
                hull.pop_back();
            } else break;
        }
        hull.push_back(p);
    }

    // Lower hull (also keeping clockwise turns; scanned in reverse order)
    int upper_size = static_cast<int>(hull.size()) + 1;
    for (int i = n - 2; i >= 0; --i) {
        int p = idx[i];
        while (static_cast<int>(hull.size()) >= upper_size) {
            int q = hull[hull.size() - 1];
            int r = hull[hull.size() - 2];
            if (cross(points(0, r), points(1, r),
                      points(0, q), points(1, q),
                      points(0, p), points(1, p)) >= 0) {
                hull.pop_back();
            } else break;
        }
        hull.push_back(p);
    }

    // Last point duplicates the first (remove it)
    if (!hull.empty()) hull.pop_back();

    int h = static_cast<int>(hull.size());
    Eigen::MatrixXd result(2, h);
    for (int i = 0; i < h; ++i) {
        result(0, i) = points(0, hull[i]);
        result(1, i) = points(1, hull[i]);
    }
    return result;
}

// Shoelace formula for polygon area.
// Input: 2 × n matrix of vertices (columns), wraps around.
double polygon_area(const Eigen::MatrixXd& polygon) {
    int n = static_cast<int>(polygon.cols());
    if (n < 3) return 0.0;

    double area = 0.0;
    for (int i = 0; i < n; ++i) {
        int j = (i + 1) % n;
        area += polygon(0, i) * polygon(1, j)
              - polygon(1, i) * polygon(0, j);
    }
    return std::abs(area) / 2.0;
}

} // namespace baysor
