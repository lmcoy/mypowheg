#include "gtest/gtest.h"

#include <iostream>
#include <cstdio>
#include <array>

#include "fks/sfunctions.h"
#include "fks/process.h"

#include "phasespace/phasespace.h"
#include "phasespace/realphasespace.h"
#include "phasespace/twoparticlegenerator.h"

using ::testing::Values;

// ---------------------------------------------------------------------------
// helper structs
// ---------------------------------------------------------------------------

struct FKSParams {
    FKSParams(double xi, double y, double phi) : Xi(xi), Y(y), Phi(phi) {}
    double Xi;
    double Y;
    double Phi;
};

// teach gtest how to print FKSParams
::std::ostream &operator<<(::std::ostream &os, const FKSParams &fks) {
    return os << "{xi = " << fks.Xi << ", y = " << fks.Y
              << ", phi = " << fks.Phi << "}";
}

struct BornPsParams {
    BornPsParams(double x1, double x2, double y, double phi)
        : X1(x1), X2(x2), Y(y), Phi(phi) {}
    double X1;
    double X2;
    double Y;
    double Phi;
};

// teach gtest how to print BornPsParams
::std::ostream &operator<<(::std::ostream &os, const BornPsParams &ps) {
    return os << "{x1 = " << ps.X1 << ", x2 = " << ps.X2 << ", y = " << ps.Y
              << ", phi = " << ps.Phi << "}";
}

struct Params {
    Params(const BornPsParams &ps, const FKSParams &fks, int mother)
        : PS(ps), FKS(fks), Mother(mother) {}
    BornPsParams PS;
    FKSParams FKS;
    int Mother;
};

// teach gtest how to print BornPsParams
::std::ostream &operator<<(::std::ostream &os, const Params &p) {
    return os << "{BornPS = " << p.PS << ", FKS = " << p.FKS << "}";
}

// ---------------------------------------------------------------------------
// test if sum_i S_i + sum_{ij} S_ij = 1
// ---------------------------------------------------------------------------

class SumSFunctionsTest : public ::testing::TestWithParam<Params> {};

// test if sum_i S_i + sum_ij S_ij = 1 for AB > CDE if there are the regions 42,
// 43.
TEST_P(SumSFunctionsTest, FinalStateOnlyFSR) {
    double SqrtS = 14000.0;

    Phasespace::Phasespace ps;
    Params params = GetParam();
    double xt[] = { params.PS.X1, params.PS.X2, params.PS.Y, params.PS.Phi };
    Phasespace::TwoParticleGenerator gen;
    double masses[2] = { 0.0 };
    gen(&ps, 4, xt, SqrtS * SqrtS, 2, masses);

    double xi = params.FKS.Xi;
    double y = params.FKS.Y;
    double phi = params.FKS.Phi;

    int mother = params.Mother;

    FKS::Real_t real;
    real.Regions.push_back(FKS::Region(4, 2));
    real.Regions.push_back(FKS::Region(4, 3));
    Phasespace::Phasespace ps_real;
    if (mother == 3 || mother == 2) {
        Phasespace::GenRealPhasespaceFSR(&ps_real, &ps, mother, xi, y, phi);
    }
    if (mother == 0 || mother == 1) {
        SUCCEED() << "no initial state radiation present.";
        return;
    }
    double S_42 = FKS::SFunction(ps_real, FKS::Region(4, 2), real, false);
    double S_43 = FKS::SFunction(ps_real, FKS::Region(4, 3), real, false);
    double sum = S_42 + S_43;
    EXPECT_DOUBLE_EQ(sum, 1.0);
}

// test if sum_i S_i + sum_ij S_ij = 1 for AB > CDE if there are the regions 4,
// 42, 43.
TEST_P(SumSFunctionsTest, DrellYanEW) {
    double SqrtS = 14000.0;

    Phasespace::Phasespace ps;
    Params params = GetParam();
    double xt[] = { params.PS.X1, params.PS.X2, params.PS.Y, params.PS.Phi };
    Phasespace::TwoParticleGenerator gen;
    double masses[2] = { 0.0 };
    gen(&ps, 4, xt, SqrtS * SqrtS, 2, masses);

    double xi = params.FKS.Xi;
    double y = params.FKS.Y;
    double phi = params.FKS.Phi;

    int mother = params.Mother;
    FKS::Real_t real;
    real.Regions.push_back(FKS::Region(4, 2));
    real.Regions.push_back(FKS::Region(4, 3));
    real.Regions.push_back(FKS::Region(4, 0));
    Phasespace::Phasespace ps_real;
    if (mother == 3 || mother == 2) {
        Phasespace::GenRealPhasespaceFSR(&ps_real, &ps, mother, xi, y, phi);
    }
    if (mother == 0 || mother == 1) {
        Phasespace::GenRealPhasespaceISR(&ps_real, &ps, xi, y, phi);
    }
    double sum = 0.0;
    for (int j = 2; j < 4; j++) {
        double Si = FKS::SFunction(ps_real, FKS::Region(4, j), real, false);
        sum += Si;
    }
    sum += FKS::SFunction(ps_real, FKS::Region(4, 0), real, false);
    EXPECT_DOUBLE_EQ(sum, 1.0);
}

INSTANTIATE_TEST_CASE_P(
    RealTestParametersFSR, SumSFunctionsTest,
    Values(
        Params(BornPsParams(0.5, 0.5, 0.4, 0.3), FKSParams(0.6, 0.3, 0.4), 3),
        Params(BornPsParams(0.2, 0.6, 0.9, 1.3), FKSParams(0.9, -0.4, 5.1), 3),
        Params(BornPsParams(0.9, 0.1, 0.6, 0.0), FKSParams(0.01, 0.99, 0.1), 3),
        Params(BornPsParams(0.5, 0.5, 0.4, 0.3), FKSParams(0.6, 0.3, 0.4), 2),
        Params(BornPsParams(0.2, 0.6, 0.9, 1.3), FKSParams(0.9, -0.4, 5.1), 2),
        Params(BornPsParams(0.9, 0.1, 0.6, 0.0), FKSParams(0.01, 0.99, 0.1),
               2)));
INSTANTIATE_TEST_CASE_P(RealTestParametersISR, SumSFunctionsTest,
                        Values(Params(BornPsParams(0.5, 0.5, 0.4, 0.3),
                                      FKSParams(0.6, 0.3, 0.4), 0),
                               Params(BornPsParams(0.2, 0.6, 0.9, 1.3),
                                      FKSParams(0.9, -0.4, 5.1), 0),
                               Params(BornPsParams(0.9, 0.1, 0.6, 0.0),
                                      FKSParams(0.01, 0.99, 0.1), 0)));
INSTANTIATE_TEST_CASE_P(SoftFSR, SumSFunctionsTest,
                        Values(Params(BornPsParams(0.5, 0.5, 0.4, 0.3),
                                      FKSParams(1e-8, 0.8, 3.1415), 3),
                               Params(BornPsParams(0.5, 0.2, 0.8, 0.0),
                                      FKSParams(1e-10, 0.6, 4.1), 3),
                               Params(BornPsParams(0.1, 0.9, 0.4, 0.3),
                                      FKSParams(1e-16, 0.2, 0.3), 3)));
INSTANTIATE_TEST_CASE_P(SoftISR, SumSFunctionsTest,
                        Values(Params(BornPsParams(0.5, 0.5, 0.4, 0.3),
                                      FKSParams(1e-8, 0.8, 3.1415), 0),
                               Params(BornPsParams(0.5, 0.2, 0.8, 0.0),
                                      FKSParams(1e-10, 0.6, 4.1), 0),
                               Params(BornPsParams(0.1, 0.9, 0.4, 0.3),
                                      FKSParams(1e-16, 0.2, 0.3), 0)));
INSTANTIATE_TEST_CASE_P(CollinearFSR, SumSFunctionsTest,
                        Values(Params(BornPsParams(0.1, 0.6, 0.8, 0.2),
                                      FKSParams(0.8, 0.99999, 3.1415), 3),
                               Params(BornPsParams(0.3, 0.6, 0.2, 0.1),
                                      FKSParams(0.3, 0.9999999, 4.1), 3),
                               Params(BornPsParams(0.7, 0.5, 0.4, 0.7),
                                      FKSParams(0.4, 0.999999999999, 0.3), 3)));
INSTANTIATE_TEST_CASE_P(CollinearISR, SumSFunctionsTest,
                        Values(Params(BornPsParams(0.1, 0.6, 0.8, 0.2),
                                      FKSParams(0.8, 0.99999, 3.1415), 0),
                               Params(BornPsParams(0.3, 0.6, 0.2, 0.1),
                                      FKSParams(0.3, 0.9999999, 4.1), 0),
                               Params(BornPsParams(0.7, 0.5, 0.4, 0.7),
                                      FKSParams(0.4, 0.999999999999, 0.3), 0)));
INSTANTIATE_TEST_CASE_P(SoftCollinearFSR, SumSFunctionsTest,
                        Values(Params(BornPsParams(0.8, 0.1, 0.5, 0.2),
                                      FKSParams(1e-8, 0.99999, 3.1415), 3),
                               Params(BornPsParams(0.5, 0.4, 0.2, 0.5),
                                      FKSParams(1e-10, 0.9999999, 4.1), 3),
                               Params(BornPsParams(0.1, 0.1, 0.9, 0.9),
                                      FKSParams(1e-15, 0.999999999999, 0.3),
                                      3)));
INSTANTIATE_TEST_CASE_P(SoftCollinearISR, SumSFunctionsTest,
                        Values(Params(BornPsParams(0.8, 0.1, 0.5, 0.2),
                                      FKSParams(1e-8, 0.99999, 3.1415), 0),
                               Params(BornPsParams(0.5, 0.4, 0.2, 0.5),
                                      FKSParams(1e-10, 0.9999999, 4.1), 0),
                               Params(BornPsParams(0.1, 0.1, 0.9, 0.9),
                                      FKSParams(1e-15, 0.999999999999, 0.3),
                                      0)));

// ---------------------------------------------------------------------------
// test if the properties of S functions are fulfilled
// ---------------------------------------------------------------------------

class SPropertiesTest : public ::testing::TestWithParam<Params> {};

// test that one S function is 1 while the others are 0 in regions where the FKS
// parton is collinear. The regions are S_4, S_42, S_43 for a 2 -> 3 process.
TEST_P(SPropertiesTest, DY_EW_CollinearInitial1) {
    double SqrtS = 14000.0;

    Phasespace::Phasespace ps;
    Params params = GetParam();
    double xt[] = { params.PS.X1, params.PS.X2, params.PS.Y, params.PS.Phi };
    Phasespace::TwoParticleGenerator gen;
    double masses[2] = {0.0};
    gen(&ps, 4, xt, SqrtS * SqrtS, 2, masses);

    double xi = params.FKS.Xi;
    double y = 1.0;
    double phi = params.FKS.Phi;
    int mother = params.Mother;

    FKS::Real_t real;
    real.Regions.push_back(FKS::Region(4, 2));
    real.Regions.push_back(FKS::Region(4, 3));
    real.Regions.push_back(FKS::Region(4, 0));
    Phasespace::Phasespace ps_real;
    if (mother == 0 || mother == 1) {
        Phasespace::GenRealPhasespaceISR(&ps_real, &ps, xi, y, phi);
    }
    if (mother == 2 || mother == 3) {
        Phasespace::GenRealPhasespaceFSR(&ps_real, &ps, mother, xi, y, phi);
    }

    double S_4 = FKS::SFunction(ps_real, FKS::Region(4, 0), real, false);
    double S_42 = FKS::SFunction(ps_real, FKS::Region(4, 2), real, false);
    double S_43 = FKS::SFunction(ps_real, FKS::Region(4, 3), real, false);
    if (mother == 0 || mother == 1) {
        EXPECT_DOUBLE_EQ(S_4, 1.0);
        EXPECT_DOUBLE_EQ(S_42, 0.0);
        EXPECT_DOUBLE_EQ(S_43, 0.0);
    }
    if (mother == 2) {
        EXPECT_NEAR(S_4, 0.0, 5e-16);
        EXPECT_NEAR(S_42, 1.0, 5e-16);
        EXPECT_NEAR(S_43, 0.0, 5e-16);
    }
    if (mother == 3) {
        EXPECT_NEAR(S_4, 0.0, 5e-16);
        EXPECT_NEAR(S_42, 0.0, 5e-16);
        EXPECT_NEAR(S_43, 1.0, 5e-16);
    }
}

// test ISR
INSTANTIATE_TEST_CASE_P(RealTestParametersISR, SPropertiesTest,
                        Values(Params(BornPsParams(0.5, 0.5, 0.4, 0.3),
                                      FKSParams(0.6, 0.3, 0.4), 0),
                               Params(BornPsParams(0.2, 0.6, 0.9, 1.3),
                                      FKSParams(0.9, -0.4, 5.1), 0),
                               Params(BornPsParams(0.9, 0.1, 0.6, 0.0),
                                      FKSParams(0.01, 0.99, 0.1), 0)));

// test FSR with j = 2
INSTANTIATE_TEST_CASE_P(RealTestParametersFSR2, SPropertiesTest,
                        Values(Params(BornPsParams(0.5, 0.5, 0.4, 0.3),
                                      FKSParams(0.6, 0.3, 0.4), 2),
                               Params(BornPsParams(0.2, 0.6, 0.9, 1.3),
                                      FKSParams(0.9, -0.4, 5.1), 2),
                               Params(BornPsParams(0.9, 0.1, 0.6, 0.0),
                                      FKSParams(0.01, 0.99, 0.1), 2)));

// test FSR with j = 3
INSTANTIATE_TEST_CASE_P(RealTestParametersFSR3, SPropertiesTest,
                        Values(Params(BornPsParams(0.5, 0.5, 0.4, 0.3),
                                      FKSParams(0.6, 0.3, 0.4), 3),
                               Params(BornPsParams(0.2, 0.6, 0.9, 1.3),
                                      FKSParams(0.9, -0.4, 5.1), 3),
                               Params(BornPsParams(0.9, 0.1, 0.6, 0.0),
                                      FKSParams(0.01, 0.99, 0.1), 3)));
