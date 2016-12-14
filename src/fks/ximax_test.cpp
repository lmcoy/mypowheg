#include <cmath>

#include "gtest/gtest.h"

#include "fks/ximax.h"

using namespace FKS;

TEST(XiMaxISR, Limit) {
    double x1 = 0.3;
    double x2 = 0.8;

    double ximax1 = XiMaxISR(x1, x2, 1.0);
    double ximax2 = XiMaxISR(x1, x2, -1.0);

    double ximax1_lim = XiMaxCollinearISR(x1, x2, 1);
    double ximax2_lim = XiMaxCollinearISR(x1, x2, -1);

    EXPECT_DOUBLE_EQ(ximax1, ximax1_lim);
    EXPECT_DOUBLE_EQ(ximax2, ximax2_lim);
}

TEST(XiMaxISR, Perpendicular) {
    double x[][2] = { { 0.1, 0.3 }, { 0.7, 0.9 }, { 0.6, 0.1 }, { 0.5, 0.5 } };
    for (int i = 0; i < 4; i++) {
        double x1 = x[i][0];
        double x2 = x[i][1];
        double ximax = XiMaxISR(x1, x2, 0.0);
        double pre = 1.0 - fmax(x1 * x1, x2 * x2);
        EXPECT_DOUBLE_EQ(ximax, pre);
    }
}
