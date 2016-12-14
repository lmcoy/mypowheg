#include "gtest/gtest.h"

#include <cmath>

#include "phasespace/phasespace.h"
#include "phasespace/twoparticlegenerator.h"
#include "math/math.h"


TEST(Init, MasslessFinalFullEnergy) {
    Phasespace::Phasespace phps;
    Phasespace::TwoParticleGenerator gen;
    // cos(theta) has to be transformed to be in [0,1): x = 0.5*(1+cos(theta))
    double x[] = {1.0,1.0,0.5,0.0};
    double masses[2] = {0.0};
    gen(&phps, 4, x, 8000.0 * 8000.0, 2, masses);

    EXPECT_DOUBLE_EQ(phps.Momenta[0].E(), 4000.0)
        << "E of parton 1 should be SqrtS/2";
    EXPECT_DOUBLE_EQ(phps.Momenta[0].PX(), 0.0)
        << "parton 1 should fly in z direction";
    EXPECT_DOUBLE_EQ(phps.Momenta[0].PY(), 0.0)
        << "parton 1 should fly in z direction";
    EXPECT_DOUBLE_EQ(phps.Momenta[0].PZ(), 4000.0)
        << "parton 1 momentum should be SqrtS/2";

    EXPECT_DOUBLE_EQ(phps.Momenta[1].E(), 4000.0)
        << "E of parton 2 should be SqrtS/2";
    EXPECT_DOUBLE_EQ(phps.Momenta[1].PX(), 0.0)
        << "parton 2 should fly in z direction";
    EXPECT_DOUBLE_EQ(phps.Momenta[1].PY(), 0.0)
        << "parton 2 should fly in z direction";
    EXPECT_DOUBLE_EQ(phps.Momenta[1].PZ(), -4000.0)
        << "parton 2 momentum should be -SqrtS/2";

    EXPECT_DOUBLE_EQ(phps.Momenta[2].E(), 4000.0)
        << "E of final 1 should be SqrtS/2 for massless particles";
    EXPECT_DOUBLE_EQ(phps.Momenta[2].PX(), 4000.0)
        << "final 1 should fly in x direction with momentum = SqrtS/2";
    EXPECT_DOUBLE_EQ(phps.Momenta[2].PY(), 0.0)
        << "final 1 should fly in x direction for cos(theta) = 0 and phi = 0.";
    EXPECT_DOUBLE_EQ(phps.Momenta[2].PZ(), 0.0)
        << "final 1 should fly in x direction for cos(theta) = 0 and phi = 0.";

    EXPECT_DOUBLE_EQ(phps.Momenta[3].E(), 4000.0)
        << "E of final 2 should be SqrtS/2 for massless particles";
    EXPECT_DOUBLE_EQ(phps.Momenta[3].PX(), -4000.0)
        << "final 2 should fly in x direction with momentum = SqrtS/2";
    EXPECT_DOUBLE_EQ(phps.Momenta[3].PY(), 0.0)
        << "final 2 should fly in x direction for cos(theta) = 0 and phi = 0.";
    EXPECT_DOUBLE_EQ(phps.Momenta[3].PZ(), 0.0)
        << "final 2 should fly in x direction for cos(theta) = 0 and phi = 0.";
}

TEST(Init, Massless) {
    double SqrtS = 7000.0;
    double x1 = 0.3;
    double x2 = 0.7;
    double cost = Math::Pi / 4.0;
    double sint = sqrt(1 - cost * cost);
    double phi = Math::Pi / 4.0;
    Phasespace::Phasespace phps;
    Phasespace::TwoParticleGenerator gen;
    double x[] = { x1, x2, (cost + 1.0) * 0.5, phi / (2.0 * Math::Pi) };
    double masses[2] = { 0.0 };
    gen(&phps, 4, x, SqrtS * SqrtS, 2, masses);

    double sqrts = SqrtS * sqrt(x1 * x2);

    EXPECT_DOUBLE_EQ(phps.Momenta[0].E(), sqrts / 2.0)
        << "E of parton 1 should be SqrtS/2*sqrt(x1*x2)";
    EXPECT_DOUBLE_EQ(phps.Momenta[0].PX(), 0.0)
        << "parton 1 should fly in z direction";
    EXPECT_DOUBLE_EQ(phps.Momenta[0].PY(), 0.0)
        << "parton 1 should fly in z direction";
    EXPECT_DOUBLE_EQ(phps.Momenta[0].PZ(), sqrts / 2.0)
        << "parton 1 momentum should be SqrtS/2*sqrt(x1*x2)";

    EXPECT_DOUBLE_EQ(phps.Momenta[1].E(), sqrts / 2.0)
        << "E of parton 2 should be SqrtS/2*sqrt(x1*x2)";
    EXPECT_DOUBLE_EQ(phps.Momenta[1].PX(), 0.0)
        << "parton 2 should fly in z direction";
    EXPECT_DOUBLE_EQ(phps.Momenta[1].PY(), 0.0)
        << "parton 2 should fly in z direction";
    EXPECT_DOUBLE_EQ(phps.Momenta[1].PZ(), -sqrts / 2.0)
        << "parton 2 momentum should be -SqrtS/2*sqrt(x1*x2)";

    EXPECT_DOUBLE_EQ(phps.Momenta[2].E(), sqrts / 2.0)
        << "E of final 1 should be SqrtS/2*sqrt(x1*x2) for massless particles";
    EXPECT_DOUBLE_EQ(phps.Momenta[2].PX(), sqrts / 2.0 * sint * cos(phi));
    EXPECT_DOUBLE_EQ(phps.Momenta[2].PY(), sqrts / 2.0 * sint * cos(phi));
    EXPECT_DOUBLE_EQ(phps.Momenta[2].PZ(), sqrts / 2.0 * cost);

    EXPECT_DOUBLE_EQ(phps.Momenta[3].E(), sqrts / 2.0)
        << "E of final 1 should be SqrtS/2*sqrt(x1*x2) for massless particles";
    EXPECT_DOUBLE_EQ(phps.Momenta[3].PX(), -sqrts / 2.0 * sint * cos(phi));
    EXPECT_DOUBLE_EQ(phps.Momenta[3].PY(), -sqrts / 2.0 * sint * cos(phi));
    EXPECT_DOUBLE_EQ(phps.Momenta[3].PZ(), -sqrts / 2.0 * cost);
}


