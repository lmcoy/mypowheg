#include "powheg/pickelement.h"

#include <cassert>
#include <cfloat>

int pick_element(int n, const double *x, double random) {
    assert(random >= 0.0 && random <= 1.0);
    if (random == 1.0) {
        random -= 4.0 * DBL_EPSILON;
    }
    double sum = 0.0;
    for (int i = 0; i < n; i++) {
        sum += x[i];
    }
    if (sum == 0.0) {
        return -1;
    }
    double x_min = 0.0;
    for (int i = 0; i < n; i++) {
        if (x[i] == 0.0) {
            continue;
        }
        double p = x[i] / sum;
        double x_max = x_min + p;
        if (random >= x_min && random < x_max) {
            return i;
        }
        x_min = x_max;
        if (x_min >= 1.0) {
            break;
        }
    }
    assert(0 && "unreachable");
    return -1;
}
