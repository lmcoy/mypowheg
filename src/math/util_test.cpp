#include "gtest/gtest.h"

#include <limits>
#include <cmath>

#include "math/util.h"

using namespace Math;

TEST(AlmostEqualDouble, Self) {
    double a = 10.0;
    EXPECT_TRUE(AlmostEqualDouble(a, a));
}

TEST(AlmostEqualDouble, Zeros) {
    EXPECT_TRUE(AlmostEqualDouble(0.0, -0.0));
}

TEST(AlmostEqualDouble, False) {
    EXPECT_FALSE(AlmostEqualDouble(0.1, -0.1));
}

TEST(AlmostEqualDouble, True) {
    double a = 1.0;
    EXPECT_TRUE(AlmostEqualDouble(1.0, a));
    a += std::numeric_limits<double>::epsilon();
    EXPECT_TRUE(AlmostEqualDouble(1.0, a));
    a += std::numeric_limits<double>::epsilon();
    EXPECT_TRUE(AlmostEqualDouble(1.0, a));
    a += std::numeric_limits<double>::epsilon();
    EXPECT_TRUE(AlmostEqualDouble(1.0, a));
    a += std::numeric_limits<double>::epsilon();
    EXPECT_TRUE(AlmostEqualDouble(1.0, a));
    a += std::numeric_limits<double>::epsilon();
    EXPECT_TRUE(AlmostEqualDouble(1.0, a));
    a += std::numeric_limits<double>::epsilon();
    EXPECT_TRUE(AlmostEqualDouble(1.0, a));
    a += std::numeric_limits<double>::epsilon();
    EXPECT_FALSE(AlmostEqualDouble(1.0, a));
}

TEST(AlmostEqualDouble, NaN) {
    double NaN = nan("");
    EXPECT_FALSE(AlmostEqualDouble(1.0, NaN));
    EXPECT_FALSE(AlmostEqualDouble(NaN, 1.0));
    EXPECT_FALSE(AlmostEqualDouble(NaN, NaN));
}
