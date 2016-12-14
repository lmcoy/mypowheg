#include "gtest/gtest.h"

#include "math/fourmomentum.h"
#include "math/math.h"

#include <cmath>
#include <iostream>

using namespace Math;

// TEST if the elements set & get methods work
TEST(FourMomentumTest, SetGetConsistency) {
    FourMomentum p1;
    double value1 = 3.1415, value2 = -0.423, value3 = 17.45, value4 = -2.34;
    p1.SetE(value1);
    EXPECT_DOUBLE_EQ(value1, p1.E());
    p1.SetPX(value2);
    EXPECT_DOUBLE_EQ(value2, p1.PX());
    p1.SetPY(value3);
    EXPECT_DOUBLE_EQ(value3, p1.PY());
    p1.SetPZ(value4);
    EXPECT_DOUBLE_EQ(value4, p1.PZ());

    value1 = 0.45;
    value2 = 4.35;
    value3 = -100.3;
    p1.SetP(value1, value2, value3);
    EXPECT_DOUBLE_EQ(value1, p1.PX());
    EXPECT_DOUBLE_EQ(value2, p1.PY());
    EXPECT_DOUBLE_EQ(value3, p1.PZ());
}

TEST(FourMomentumTest, DefaultConstructorTest) {
    FourMomentum p;
    EXPECT_DOUBLE_EQ(p.E(), 0.0);
    EXPECT_DOUBLE_EQ(p.PX(), 0.0);
    EXPECT_DOUBLE_EQ(p.PY(), 0.0);
    EXPECT_DOUBLE_EQ(p.PZ(), 0.0);
}

TEST(FourMomentumTest, ElementAccessByIndex) {
    FourMomentum p;
    double value1 = 4.1315, value2 = 4.423, value3 = 1.45, value4 = -2.34;
    p.SetE(value1);
    p.SetP(value2, value3, value4);
    EXPECT_DOUBLE_EQ(value1, p.At(0));
    EXPECT_DOUBLE_EQ(value2, p.At(1));
    EXPECT_DOUBLE_EQ(value3, p.At(2));
    EXPECT_DOUBLE_EQ(value4, p.At(3));
}

TEST(Dot, Square) {
    FourMomentum p1(0.0, 0.0, 0.0, 0.0);
    EXPECT_DOUBLE_EQ(p1.Dot(p1), 0.0);
    p1.Set(3.0, 4.0, 5.0, 6.0);
    EXPECT_DOUBLE_EQ(p1.Dot(p1), -68.0);

    double xi = 0.6;
    double phi = 0.4;
    double theta = 0.3;
    p1.Set(xi, xi * sin(theta) * cos(phi), xi * sin(theta) * sin(phi),
           xi * cos(theta));
    EXPECT_NEAR(p1.Dot(p1), 0.0, 1e-14);

    double p = 3.0;
    double m = 4.0;
    theta = 0.1;
    phi = 4.1;
    p1.Set(sqrt(m * m + p * p), p * sin(theta) * cos(phi),
           p * sin(theta) * sin(phi), p * cos(theta));
    EXPECT_NEAR(p1.Dot(p1), m * m, 1e-14);
}

TEST(Operator, Plus) {
    FourMomentum p1(2.0, 3.0, 4.0, 5.0);
    FourMomentum p2(-2.0, -3.0, -4.0, -5.0);
    FourMomentum result = p1.Plus(p2);
    EXPECT_DOUBLE_EQ(result.E(), 0.0);
    EXPECT_DOUBLE_EQ(result.PX(), 0.0);
    EXPECT_DOUBLE_EQ(result.PY(), 0.0);
    EXPECT_DOUBLE_EQ(result.PZ(), 0.0);

    p1.Set(2.0, -4.0, 3.0, 5.0);
    p2.Set(-.314, 3.14, 4.23, -3.0);
    result = p1.Plus(p2);
    EXPECT_DOUBLE_EQ(result.E(), p1.E() + p2.E());
    EXPECT_DOUBLE_EQ(result.PX(), p1.PX() + p2.PX());
    EXPECT_DOUBLE_EQ(result.PY(), p1.PY() + p2.PY());
    EXPECT_DOUBLE_EQ(result.PZ(), p1.PZ() + p2.PZ());
}

TEST(Operator, Minus) {
    FourMomentum p1(2.0, 3.0, 4.0, 5.0);
    FourMomentum result = p1.Minus(p1);
    EXPECT_DOUBLE_EQ(result.E(), 0.0);
    EXPECT_DOUBLE_EQ(result.PX(), 0.0);
    EXPECT_DOUBLE_EQ(result.PY(), 0.0);
    EXPECT_DOUBLE_EQ(result.PZ(), 0.0);

    p1.Set(2.0, -4.0, 3.0, 5.0);
    FourMomentum p2(-.314, 3.14, 4.23, -3.0);
    result = p1.Minus(p2);
    EXPECT_DOUBLE_EQ(result.E(), p1.E() - p2.E());
    EXPECT_DOUBLE_EQ(result.PX(), p1.PX() - p2.PX());
    EXPECT_DOUBLE_EQ(result.PY(), p1.PY() - p2.PY());
    EXPECT_DOUBLE_EQ(result.PZ(), p1.PZ() - p2.PZ());
}

TEST(Operator, Scale) {
    FourMomentum p(1.0, 2.0, 3.0, 4.0);
    double scale = 3.1415;
    p.Scale(scale);
    EXPECT_DOUBLE_EQ(p.E(), scale);
    EXPECT_DOUBLE_EQ(p.PX(), 2.0 * scale);
    EXPECT_DOUBLE_EQ(p.PY(), 3.0 * scale);
    EXPECT_DOUBLE_EQ(p.PZ(), 4.0 * scale);
}

TEST(Equals, Self) {
    FourMomentum p(1.0, 2.0, 3.0, 4.0);
    EXPECT_TRUE(p.Equals(p));
}

TEST(Equals, Two) {
    FourMomentum p1(2.0, 4.0, 6.0, 9.0);
    FourMomentum p2(2.0, 4.0, 6.0, 9.0);
    EXPECT_TRUE(p1.Equals(p2));
    EXPECT_TRUE(p2.Equals(p1));
}

TEST(Equals, NotEqual) {
    FourMomentum p1(4.0, 2.0, 5.0, 2.0);
    FourMomentum p2(-4.0, 17.0, -1.0, 3.0);
    EXPECT_FALSE(p1.Equals(p2));
    EXPECT_FALSE(p2.Equals(p1));
}

TEST(MomentumMagnitude, MomentumMagnitude) {
    FourMomentum p(0.0, 0.0, 0.0, 0.0);
    EXPECT_DOUBLE_EQ(p.MomentumMagnitude(), 0.0);

    p.Set(3.0, 1.0, 0.0, 0.0);
    EXPECT_DOUBLE_EQ(p.MomentumMagnitude(), 1.0);

    p.Set(4.0, 0.0, 1.0, 0.0);
    EXPECT_DOUBLE_EQ(p.MomentumMagnitude(), 1.0);

    p.Set(5.0, 0.0, 0.0, 1.0);
    EXPECT_DOUBLE_EQ(p.MomentumMagnitude(), 1.0);

    p.Set(9.0, 3.0, 4.0, 0.0);
    EXPECT_DOUBLE_EQ(p.MomentumMagnitude(), 5.0);

    p.Set(1.5, 0.0, 4.0, 3.0);
    EXPECT_DOUBLE_EQ(p.MomentumMagnitude(), 5.0);

    p.Set(7.0, 6.0, 8.0, 0.0);
    EXPECT_DOUBLE_EQ(p.MomentumMagnitude(), 10.0);
}

TEST(CosTheta, Normal) {
    FourMomentum p1(1.0, 0.0, 0.0, 1.0);
    EXPECT_DOUBLE_EQ(p1.CosTheta(), 1.0);

    p1.Set(1.0, 0.0, 1.0, 0.0);
    EXPECT_DOUBLE_EQ(p1.CosTheta(), 0.0);

    p1.Set(1.0, 1.0, 1.0, 0.0);
    EXPECT_DOUBLE_EQ(p1.CosTheta(), 0.0);

    p1.Set(1.0, 0.0, 0.0, -1.0);
    EXPECT_DOUBLE_EQ(p1.CosTheta(), -1.0);

    p1.Set(1.0, 1.0, 0.0, 1.0);
    EXPECT_DOUBLE_EQ(p1.CosTheta(), 1.0 / sqrt(2.0));

    p1.Set(1.0, 0.0, 0.5, 0.5 * sqrt(3.0));
    EXPECT_DOUBLE_EQ(p1.CosTheta(), 0.5 * sqrt(3.0));
}

TEST(CosTheta, New) {
    const double theta = 0.3;
    FourMomentum p = FourMomentum::NewMasslessFromPThetaPhi(1.2, theta, 2.1);
    EXPECT_DOUBLE_EQ(p.CosTheta(), cos(theta));
}

TEST(Phi, Normal) {
    FourMomentum p(2.0, 0.0, 0.0, 1.0);
    EXPECT_DOUBLE_EQ(p.Phi(), 0.0);

    p.Set(5.0, 0.0, 0.0, -1.0);
    EXPECT_DOUBLE_EQ(p.Phi(), 0.0);

    p.Set(4.0, 1.0, 0.0, 0.0);
    EXPECT_DOUBLE_EQ(p.Phi(), 0.0);

    p.Set(0.0, 0.0, 5.0, 0.0);
    EXPECT_DOUBLE_EQ(p.Phi(), Math::Pi / 2.0);

    p.Set(1.0, 0.0, -1.0, 4.0);
    EXPECT_DOUBLE_EQ(p.Phi(), 1.5 * Math::Pi);

    p.Set(1.0, -10.0, 0.0, -4.0);
    EXPECT_DOUBLE_EQ(p.Phi(), Math::Pi);

    const int max = 10;
    for (int i = 0; i < max; i++) {
        double phi = 2.0 * Math::Pi / ((double)max) * ((double)i);
        p.Set(2.0, cos(phi), sin(phi), 0.0);
        EXPECT_DOUBLE_EQ(p.Phi(), phi);
    }
}

TEST(DeltaR, Self) {
    FourMomentum p1(1.0, 2.0, 3.0, 4.0);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p1, p1), 0.0);

    p1.Set(3.0, -4.0, 5.0, 60.0);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p1, p1), 0.0);
}

TEST(Eta, PhiZero) {
    double cost1 = 0.0; // theta = pi/2
    double sint1 = sqrt(1.0 - cost1 * cost1);
    FourMomentum p1(1.0, sint1, 0.0, cost1);
    EXPECT_DOUBLE_EQ(p1.Eta(), 0.0);

    double cost2 = 0.5; // theta = pi/3
    double sint2 = sqrt(1.0 - cost2 * cost2);
    FourMomentum p2(2.0, 0.0, sint2, cost2);
    EXPECT_DOUBLE_EQ(p2.Eta(), log(3.0) / 2.0);

    double cost3 = -0.5; // theta = 2/3*pi
    double sint3 = sqrt(1.0 - cost3 * cost3);
    FourMomentum p3(3.0, sint3, 0.0, cost3);
    EXPECT_DOUBLE_EQ(p3.Eta(), -log(3.0) / 2.0);
}

TEST(Eta, PhiNotZero) {
    double r1 = 3.0;
    double phi1 = 0.2;
    double cost1 = 0.0; // theta = pi/2
    double sint1 = sqrt(1.0 - cost1 * cost1);
    FourMomentum p1(1.0, r1 * sint1 * cos(phi1), r1 * sint1 * sin(phi1),
                    r1 * cost1);
    EXPECT_DOUBLE_EQ(p1.Eta(), 0.0);

    double r2 = 4.0;
    double phi2 = 0.4;
    double cost2 = 0.5; // theta = pi/3
    double sint2 = sqrt(1.0 - cost2 * cost2);
    FourMomentum p2(2.0, r2 * sint2 * cos(phi2), r2 * sint2 * sin(phi2),
                    r2 * cost2);
    EXPECT_DOUBLE_EQ(p2.Eta(), log(3.0) / 2.0);

    double r3 = 5.0;
    double phi3 = 3.7;
    double cost3 = -0.5; // theta = 2/3*pi
    double sint3 = sqrt(1.0 - cost3 * cost3);
    FourMomentum p3(3.0, r3 * sint3 * cos(phi3), r3 * sint3 * sin(phi3),
                    r3 * cost3);
    EXPECT_DOUBLE_EQ(p3.Eta(), -log(3.0) / 2.0);
}

TEST(Eta, LimitPhiZero) {
    double cost1 = 1.0; // theta = 0
    FourMomentum p1(1.0, 0.0, 0.0, cost1);
    EXPECT_DOUBLE_EQ(p1.Eta(), std::numeric_limits<double>::max());

    double cost2 = -1.0; // theta = pi
    FourMomentum p2(2.0, 0.0, 0.0, cost2);
    EXPECT_DOUBLE_EQ(p2.Eta(), std::numeric_limits<double>::lowest());
}

TEST(DeltaR, Consistency) {
    FourMomentum p1(2.0, 4.0, 3.0, 1.3);
    FourMomentum p2(2.0, -0.04, 3.1415, -5.3);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p1, p2), FourMomentum::DeltaR(p2, p1))
        << "angular distance cannot depend on ordering";

    p1.Set(3.0, 4.1, 5.0, 6.0);
    p2.Set(2.0, -4.5, 3.4, -6.0);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p1, p2), FourMomentum::DeltaR(p2, p1))
        << "angular distance cannot depend on ordering";
}

TEST(DeltaR, Azimuthal) {
    double phi1 = Math::Pi / 4.0;
    double phi2 = 3.0 * Math::Pi / 4.0;
    double phi3 = 5.0 * Math::Pi / 4.0;
    double phi4 = 7.0 * Math::Pi / 4.0;
    FourMomentum p1(1.0, cos(phi1), sin(phi1), 0.0);
    FourMomentum p2(2.0, cos(phi2), sin(phi2), 0.0);
    FourMomentum p3(3.0, cos(phi3), sin(phi3), 0.0);
    FourMomentum p4(4.0, cos(phi4), sin(phi4), 0.0);

    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p1, p2), Math::Pi / 2.0);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p2, p1), Math::Pi / 2.0);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p2, p3), Math::Pi / 2.0);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p3, p2), Math::Pi / 2.0);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p3, p4), Math::Pi / 2.0);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p4, p3), Math::Pi / 2.0);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p4, p1), Math::Pi / 2.0);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p1, p4), Math::Pi / 2.0);

    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p1, p3), Math::Pi);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p3, p1), Math::Pi);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p2, p4), Math::Pi);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p4, p2), Math::Pi);
}

TEST(DeltaR, Eta) {
    double cost1 = 0.0; // theta = pi/2
    double sint1 = sqrt(1.0 - cost1 * cost1);
    FourMomentum p1(1.0, sint1, 0.0, cost1);

    double cost2 = 0.5; // theta = pi/3
    double sint2 = sqrt(1.0 - cost2 * cost2);
    FourMomentum p2(2.0, sint2, 0.0, cost2);

    double cost3 = -0.5; // theta = 2/3*pi
    double sint3 = sqrt(1.0 - cost3 * cost3);
    FourMomentum p3(3.0, sint3, 0.0, cost3);

    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p1, p2), log(3.0) / 2.0);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p1, p3), log(3.0) / 2.0);
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p2, p3), log(3.0));
}

TEST(DeltaR, Normal) {
    double r1 = 3.0;
    double phi1 = 0.2;
    double cost1 = 0.0; // theta = pi/2
    double sint1 = sqrt(1.0 - cost1 * cost1);
    FourMomentum p1(1.0, r1 * sint1 * cos(phi1), r1 * sint1 * sin(phi1),
                    r1 * cost1);

    double r2 = 4.0;
    double phi2 = 0.4;
    double cost2 = 0.5; // theta = pi/3
    double sint2 = sqrt(1.0 - cost2 * cost2);
    FourMomentum p2(2.0, r2 * sint2 * cos(phi2), r2 * sint2 * sin(phi2),
                    r2 * cost2);

    double r3 = 5.0;
    double phi3 = 3.7;
    double cost3 = -0.5; // theta = 2/3*pi
    double sint3 = sqrt(1.0 - cost3 * cost3);
    FourMomentum p3(3.0, r3 * sint3 * cos(phi3), r3 * sint3 * sin(phi3),
                    r3 * cost3);

    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p1, p2),
                     sqrt(log(3.0) * log(3.0) / 4.0 + 0.2 * 0.2));
    double dphi = 2.0 * Math::Pi - 3.5;
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p1, p3),
                     sqrt(log(3.0) * log(3.0) / 4.0 + dphi * dphi));
    dphi = 2.0 * Math::Pi - 3.3;
    EXPECT_DOUBLE_EQ(FourMomentum::DeltaR(p3, p2),
                     sqrt(log(3.0) * log(3.0) + dphi * dphi));
}
