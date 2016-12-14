#include "util.h"

#include <cstdint>
#include <cstdlib>
#include <cmath>

/**
AlmostEqualDouble returns true if \p a and \p b are almost equal.

This functions relies on the fact that two adjacent doubles have adjacent
integer representations.

Two double are considered equal if the distance between their integer
representations is less then 7.
 */
bool Math::AlmostEqualDouble(const double a, const double b) {
    if (isnan(a) || isnan(b)) {
        return false;
    }
    int64_t _a = *((int64_t *)&a);
    int64_t _b = *((int64_t *)&b);
    if (_a < 0) {
        _a = 0x8000000000000000 - _a;
    }
    if (_b < 0) {
        _b = 0x8000000000000000 - _b;
    }
    return llabs(_a - _b) <= 6;
}
