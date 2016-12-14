#include "gtest/gtest.h"

#include <cmath>

#include "phasespace/phasespace.h"
#include "phasespace/twoparticlegenerator.h"
#include "math/math.h"

TEST(BoostToLabFrame, FullEnergy) {
    double SqrtS = 8000.0;
    Phasespace::Phasespace<2> phps(SqrtS * SqrtS, { { 0.0, 0.0 } });
    Phasespace::TwoParticleGenerator gen;
    double x[] = { 1.0, 1.0, 0.0, 0.0 };
    gen(&phps, 4, x);

    Math::FourMomentum initial1_cms = phps.Momenta[0];
    Math::FourMomentum initial2_cms = phps.Momenta[1];
    Math::FourMomentum final1_cms = phps.Momenta[2];
    Math::FourMomentum final2_cms = phps.Momenta[3];

    Phasespace::Phasespace<2> phps_lab(SqrtS * SqrtS, { { 0.0, 0.0 } });
    phps_lab.SetToLabFromCMS(&phps);

    EXPECT_TRUE(initial1_cms.Equals(phps_lab.Momenta[0]))
        << "Lab frame equals parton CMS for x1 = 1, x2 = 1";
    EXPECT_TRUE(initial2_cms.Equals(phps_lab.Momenta[1]))
        << "Lab frame equals parton CMS for x1 = 1, x2 = 1";
    EXPECT_TRUE(final1_cms.Equals(phps_lab.Momenta[2]))
        << "Lab frame equals parton CMS for x1 = 1, x2 = 1";
    EXPECT_TRUE(final2_cms.Equals(phps_lab.Momenta[3]))
        << "Lab frame equals parton CMS for x1 = 1, x2 = 1";
}

TEST(BoostToLabFrame, Massless) {
    double SqrtS = 8000.0;
    double x1 = 0.3;
    double x2 = 0.7;
    double cost = Math::Pi / 4.0;
    double phi = Math::Pi / 4.0;

    Phasespace::Phasespace<2> phps(SqrtS * SqrtS, { { 0.0, 0.0 } });
    Phasespace::TwoParticleGenerator gen;
    double x[] = { x1, x2, (1. + cost) / 2., phi / (2.0 * Math::Pi) };
    gen(&phps, 4, x);

    Math::FourMomentum final1_cms = phps.Momenta[2];
    Math::FourMomentum final2_cms = phps.Momenta[3];

    Phasespace::Phasespace<2> phps_lab(SqrtS * SqrtS, { { 0.0, 0.0 } });
    phps_lab.SetToLabFromCMS(&phps);

    // check initial state momentum vectors
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[0].E(), x1 * SqrtS / 2.0);
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[0].PX(), 0.0);
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[0].PY(), 0.0);
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[0].PZ(), x1 * SqrtS / 2.0);

    EXPECT_DOUBLE_EQ(phps_lab.Momenta[1].E(), x2 * SqrtS / 2.0);
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[1].PX(), 0.0);
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[1].PY(), 0.0);
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[1].PZ(), -x2 * SqrtS / 2.0);

    // test if energy and momentum is still conserved
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[2].E() + phps_lab.Momenta[3].E(),
                     phps_lab.Momenta[0].E() + phps_lab.Momenta[1].E())
        << "energy conservation violated";
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[2].PX() + phps_lab.Momenta[3].PX(),
                     phps_lab.Momenta[0].PX() + phps_lab.Momenta[1].PX())
        << "momentum conservation violated";
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[2].PY() + phps_lab.Momenta[3].PY(),
                     phps_lab.Momenta[0].PY() + phps_lab.Momenta[1].PY())
        << "momentum conservation violated";
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[2].PZ() + phps_lab.Momenta[3].PZ(),
                     phps_lab.Momenta[0].PZ() + phps_lab.Momenta[1].PZ())
        << "momentum conservation violated";

    EXPECT_DOUBLE_EQ(phps_lab.Momenta[2].E(),
                     SqrtS / 4.0 * ((x1 + x2) + (x1 - x2) * cost));
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[2].PX(), final1_cms.PX())
        << "px component of momenta doesn't change.";
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[2].PY(), final1_cms.PY())
        << "py component of momenta doesn't change.";
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[2].PZ(),
                     SqrtS / 4.0 * ((x1 + x2) * cost + (x1 - x2)));

    EXPECT_DOUBLE_EQ(phps_lab.Momenta[3].E(),
                     SqrtS / 4.0 * ((x1 + x2) - (x1 - x2) * cost));
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[3].PX(), final2_cms.PX())
        << "px component of momenta doesn't change.";
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[3].PY(), final2_cms.PY())
        << "py component of momenta doesn't change.";
    EXPECT_DOUBLE_EQ(phps_lab.Momenta[3].PZ(),
                     SqrtS / 4.0 * (-(x1 + x2) * cost + (x1 - x2)));
}
