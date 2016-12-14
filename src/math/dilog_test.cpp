#include "gtest/gtest.h"

#include <cmath>

#include "math/dilog.h"

double diff(double a, double b) {
    return fabs(a-b)/b;
}

TEST(DILOG, Value1) {
    EXPECT_DOUBLE_EQ(Math::Dilog(1.0), M_PI*M_PI/6.0);
}

TEST(DILOG, Value0) {
    EXPECT_DOUBLE_EQ(Math::Dilog(0.0), 0.0);
}

TEST(DILOG, Values) {
    EXPECT_NEAR(diff(Math::Dilog(0.1), 0.1026177910993911311138374), 0.0,
                1e-11);
    EXPECT_NEAR(diff(Math::Dilog(0.2), 0.2110037754397047726111851), 0.0,
                1e-11);
    EXPECT_NEAR(diff(Math::Dilog(0.3), 0.3261295100754760695300357), 0.0,
                1e-11);
    EXPECT_NEAR(diff(Math::Dilog(0.4), 0.4492829744712816644647334), 0.0,
                1e-11);
    EXPECT_NEAR(
        diff(Math::Dilog(0.5), M_PI * M_PI / 12.0 - log(2.0) * log(2.0) / 2.0),
        0.0, 1e-15);
    EXPECT_NEAR(diff(Math::Dilog(0.6), 0.7275863077163333895135363), 0.0,
                1e-11);
    EXPECT_NEAR(diff(Math::Dilog(0.7), 0.8893776242860387386010063), 0.0,
                1e-11);
    EXPECT_NEAR(diff(Math::Dilog(0.8), 1.074794600008248359395452), 0.0, 1e-11);
    EXPECT_NEAR(diff(Math::Dilog(0.9), 1.299714723004958725171060), 0.0, 1e-11);
}
