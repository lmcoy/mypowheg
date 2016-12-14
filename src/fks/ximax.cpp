#include "fks/ximax.h"

namespace FKS {

double XiMaxISR(double x1_b, double x2_b, double y) {
    double x1_2 = x1_b * x1_b;
    double x2_2 = x2_b * x2_b;
    double xi_max_tmp1 =
        (2.0 * x1_2 * (1.0 + y)) /
        ((1.0 - x1_2) * (1.0 - y) +
         sqrt((1.0 + x1_2) * (1.0 + x1_2) * (1 - y) * (1.0 - y) +
              16.0 * x1_2 * y));
    double xi_max_tmp2 =
        (2.0 * x2_2 * (1.0 - y)) /
        ((1.0 - x2_2) * (1.0 + y) +
         sqrt((1.0 + x2_2) * (1.0 + x2_2) * (1 + y) * (1.0 + y) -
              16.0 * x2_2 * y));
    double xi_max = 1.0 - fmax(xi_max_tmp1, xi_max_tmp2);
    return xi_max;
}

double XiMaxCollinearISR(double x1_b, double x2_b, int y) {
    assert(y == 1 || y == -1);
    if (y == 1) {
        return 1.0 - x1_b;
    }
    return 1.0 - x2_b;
}

double XiMaxFSR(double sqrts, double len_kj_born) {
    return 2.0 * len_kj_born / sqrts;
}

} // namespace FKS
