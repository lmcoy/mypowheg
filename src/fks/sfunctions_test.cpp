#include "gtest/gtest.h"

#include <array>
#include <cstdio>
#include <iostream>

#include "fks/process.h"
#include "fks/sfunctions.h"
#include "fks/ximax.h"

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

static Math::FourMomentum boost_from_rest(const Math::FourMomentum &q,
                                          const Math::FourMomentum &ref) {
    double b[] = {ref.At(1) / ref.At(0), ref.At(2) / ref.At(0),
                  ref.At(3) / ref.At(0)};
    double b2 = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
    double gamma = 1 / sqrt(1.0 - b2);
    double bp = b[0] * q.At(1) + b[1] * q.At(2) + b[2] * q.At(3);
    double gamma2 = (gamma - 1) / b2;

    double t = gamma * (q.At(0) + bp);
    double x = q.At(1) + gamma2 * bp * b[0] + gamma * b[0] * q.At(0);
    double y = q.At(2) + gamma2 * bp * b[1] + gamma * b[1] * q.At(0);
    double z = q.At(3) + gamma2 * bp * b[2] + gamma * b[2] * q.At(0);
    return Math::FourMomentum(t, x, y, z);
}

static Phasespace::Phasespace GenPhasespace3(double S,
                                             const std::array<double, 7> &v) {
    assert(v[0] >= 0.0 && v[0] <= 1.0);
    assert(v[1] >= 0.0 && v[1] <= 1.0);
    assert(v[2] >= 0.0 && v[2] <= 1.0);
    assert(v[3] >= 0.0 && v[3] <= 1.0);
    assert(v[4] >= 0.0 && v[4] <= 1.0);
    assert(v[5] >= 0.0 && v[5] <= 1.0);
    assert(v[6] >= 0.0 && v[6] <= 1.0);
    double x1 = v[0];
    double x2 = v[1];
    double s = x1 * x2 * S;
    double p = sqrt(s) / 2.0;

    Math::FourMomentum k1(p, 0.0, 0.0, p);
    Math::FourMomentum k2(p, 0.0, 0.0, -p);

    double Q2 = v[2] * s;

    double E1 = (s - Q2) / sqrt(s) / 2.0;
    double costheta1 = (2 * v[3] - 1.0);
    double sintheta1 = sqrt(1 - costheta1 * costheta1);
    double phi1 = 2. * M_PI * v[4];
    double costheta2 = (2 * v[5] - 1.0);
    double sintheta2 = sqrt(1 - costheta2 * costheta2);
    double phi2 = 2. * M_PI * v[6];

    Math::FourMomentum p1(E1, E1 * sintheta1 * cos(phi1),
                          E1 * sintheta1 * sin(phi1), E1 * costheta1);
    Math::FourMomentum q(2. * p - p1.At(0), -p1.At(1), -p1.At(2), -p1.At(3));

    // p2, p3 in rest frame of q
    double t = q.Dot(q);
    double E2 = 0.5 * sqrt(t);
    Math::FourMomentum p2(E2, E2 * sintheta2 * cos(phi2),
                          E2 * sintheta2 * sin(phi2), E2 * costheta2);
    Math::FourMomentum p3(p2.At(0), -p2.At(1), -p2.At(2), -p2.At(3));

    // boost p2, p3 to frame of p1 and q
    p2 = boost_from_rest(p2, q);
    p3 = boost_from_rest(p3, q);

    double p2p2 = p2.Dot(p2);
    if (p2p2 < 0.0) {
        double x = sqrt(-p2p2 + p2.MomentumMagnitudeSqr());
        p2.SetE(x);
    }
    double p3p3 = p3.Dot(p3);
    if (p3p3 < 0.0) {
        double x = sqrt(-p3p3 + p3.MomentumMagnitudeSqr());
        p3.SetE(x);
    }

    double jac = (1 - v[2]) * s / (8 * M_PI * M_PI * M_PI * 16.0);

    Phasespace::Phasespace ps;
    ps.N = 3;
    ps.S = S;
    ps.X1 = x1;
    ps.X2 = x2;
    ps.Jacobian = jac;
    ps.Momenta[0] = k1;
    ps.Momenta[1] = k2;
    ps.Momenta[2] = p2;
    ps.Momenta[3] = p3;
    ps.Momenta[4] = p1;

    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;

    return ps;
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
    double xt[] = {params.PS.X1, params.PS.X2, params.PS.Y, params.PS.Phi};
    Phasespace::TwoParticleGenerator gen;
    double masses[2] = {0.0};
    gen(&ps, 4, xt, SqrtS * SqrtS, 2, masses);

    double xi = params.FKS.Xi;
    double y = params.FKS.Y;
    double phi = params.FKS.Phi;

    int mother = params.Mother;

    FKS::Real_t real;
    real.Flavours = {1, 1, 1, 1, 1, 1, 1, 1}; // no gluons!
    real.AllRegions.push_back(FKS::Region(4, 2));
    real.AllRegions.push_back(FKS::Region(4, 3));
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
    double xt[] = {params.PS.X1, params.PS.X2, params.PS.Y, params.PS.Phi};
    Phasespace::TwoParticleGenerator gen;
    double masses[2] = {0.0};
    gen(&ps, 4, xt, SqrtS * SqrtS, 2, masses);

    double xi = params.FKS.Xi;
    double y = params.FKS.Y;
    double phi = params.FKS.Phi;

    int mother = params.Mother;
    FKS::Real_t real;
    real.Flavours = {1, 1, 1, 1, 1, 1, 1, 1}; // no gluons!
    real.AllRegions.push_back(FKS::Region(4, 2));
    real.AllRegions.push_back(FKS::Region(4, 3));
    real.AllRegions.push_back(FKS::Region(4, 0));
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

TEST_P(SumSFunctionsTest, WjQCD_1) {
    double SqrtS = 14000.0;

    Params params = GetParam();
    std::array<double, 7> xt = {
        params.PS.X1, params.PS.X2, params.PS.Y, params.PS.Phi, 0.1, 0.4, 0.3};
    auto ps = GenPhasespace3(SqrtS * SqrtS, xt);

    double xi = params.FKS.Xi;
    double y = params.FKS.Y;
    double phi = params.FKS.Phi;

    int mother = params.Mother;
    if (mother >= 2) {
        mother += 2; // radiation from 4 or 5
    }
    int ifks = 5;
    if (mother == 5) {
        ifks = 4;
    }
    double ximax = 1.0;
    if (mother >= 2) {
        ximax = FKS::XiMaxFSR(sqrt(ps.X1 * ps.X2 * ps.S), ps.Momenta[4].E());
    } else {
        ximax = FKS::XiMaxISR(ps.X1, ps.X2, y);
    }
    xi *= ximax;

    FKS::Real_t real;
    real.Flavours = {2, -1, -13, 14, 21, 21};
    real.AllRegions.push_back(FKS::Region(5, 0));
    real.AllRegions.push_back(FKS::Region(4, 0));
    real.AllRegions.push_back(FKS::Region(5, 4));
    real.AllRegions.push_back(FKS::Region(4, 5));
    Phasespace::Phasespace ps_real;
    Phasespace::GenRealPhasespace(&ps_real, &ps, ifks, mother, xi, y, phi);

    double S40 = FKS::SFunction(ps_real, FKS::Region(4, 0), real, false);
    double S50 = FKS::SFunction(ps_real, FKS::Region(5, 0), real, false);
    double S54 = FKS::SFunction(ps_real, FKS::Region(5, 4), real, false);
    double S45 = FKS::SFunction(ps_real, FKS::Region(4, 5), real, false);
    double sum = S40 + S50 + S45 + S54;

    EXPECT_NEAR(sum, 1.0, 5e-8) << "S50 = " << S50 << " S40 = " << S40
                                << " S54 = " << S54 << " S45 = " << S45;
}

TEST_P(SumSFunctionsTest, WjQCD_2) {
    double SqrtS = 14000.0;

    Params params = GetParam();
    std::array<double, 7> xt = {
        params.PS.X1, params.PS.X2, params.PS.Y, params.PS.Phi, 0.1, 0.2, 0.3};
    auto ps = GenPhasespace3(SqrtS * SqrtS, xt);

    double xi = params.FKS.Xi;
    double y = params.FKS.Y;
    double phi = params.FKS.Phi;

    int mother = params.Mother;
    if (mother >= 2) {
        mother += 2; // radiation from 4 or 5
    }
    int ifks = 5;
    if (mother == 5) {
        ifks = 4;
    }
    double ximax = 1.0;
    if (mother >= 2) {
        ximax = FKS::XiMaxFSR(sqrt(ps.X1 * ps.X2 * ps.S), ps.Momenta[4].E());
    } else {
        ximax = FKS::XiMaxISR(ps.X1, ps.X2, y);
    }
    xi *= ximax;

    FKS::Real_t real;
    real.Flavours = {2, 21, -13, 14, 1, 21};
    real.AllRegions.push_back(FKS::Region(5, 0));
    real.AllRegions.push_back(FKS::Region(5, 4));
    Phasespace::Phasespace ps_real;
    Phasespace::GenRealPhasespace(&ps_real, &ps, ifks, mother, xi, y, phi);

    double sum = 0.0;
    sum += FKS::SFunction(ps_real, FKS::Region(5, 0), real, false);
    sum += FKS::SFunction(ps_real, FKS::Region(5, 4), real, false);
    EXPECT_DOUBLE_EQ(sum, 1.0);
}

INSTANTIATE_TEST_CASE_P(
    RealTestParametersFSR, SumSFunctionsTest,
    Values(
        Params(BornPsParams(0.5, 0.5, 0.4, 0.3), FKSParams(0.6, 0.3, 0.4), 3),
        Params(BornPsParams(0.2, 0.6, 0.9, 0.7), FKSParams(0.9, -0.4, 5.1), 3),
        Params(BornPsParams(0.9, 0.1, 0.6, 0.0), FKSParams(0.01, 0.99, 0.1), 3),
        Params(BornPsParams(0.5, 0.5, 0.4, 0.3), FKSParams(0.6, 0.3, 0.4), 2),
        Params(BornPsParams(0.2, 0.6, 0.9, 0.7), FKSParams(0.9, -0.4, 5.1), 2),
        Params(BornPsParams(0.9, 0.1, 0.6, 0.0), FKSParams(0.01, 0.99, 0.1),
               2)));
INSTANTIATE_TEST_CASE_P(RealTestParametersISR, SumSFunctionsTest,
                        Values(Params(BornPsParams(0.5, 0.5, 0.4, 0.3),
                                      FKSParams(0.6, 0.3, 0.4), 0),
                               Params(BornPsParams(0.2, 0.6, 0.9, 0.7),
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
    double xt[] = {params.PS.X1, params.PS.X2, params.PS.Y, params.PS.Phi};
    Phasespace::TwoParticleGenerator gen;
    double masses[2] = {0.0};
    gen(&ps, 4, xt, SqrtS * SqrtS, 2, masses);

    double xi = params.FKS.Xi;
    double y = 1.0;
    double phi = params.FKS.Phi;
    int mother = params.Mother;

    FKS::Real_t real;
    real.Flavours = {1, 1, 1, 1, 1, 1, 1, 1}; // no gluons!
    real.AllRegions.push_back(FKS::Region(4, 2));
    real.AllRegions.push_back(FKS::Region(4, 3));
    real.AllRegions.push_back(FKS::Region(4, 0));
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

TEST_P(SPropertiesTest, WjQCD_1) {
    double SqrtS = 14000.0;

    Params params = GetParam();
    std::array<double, 7> xt = {
        params.PS.X1, params.PS.X2, params.PS.Y, params.PS.Phi, 0.5, 0.2, 0.0};
    auto ps = GenPhasespace3(SqrtS * SqrtS, xt);

    double xi = params.FKS.Xi;
    double y = 0.999999999;
    double phi = params.FKS.Phi;

    int mother = params.Mother;
    if (mother >= 2) {
        mother += 2; // radiation from 4 or 5
    }
    int ifks = 5;
    if (mother == 5) {
        ifks = 4;
    }
    double ximax = 1.0;
    if (mother >= 2) {
        ximax = FKS::XiMaxFSR(sqrt(ps.X1 * ps.X2 * ps.S), ps.Momenta[4].E());
    } else {
        ximax = FKS::XiMaxISR(ps.X1, ps.X2, y);
    }
    xi *= ximax;

    FKS::Real_t real;
    real.Flavours = {2, -1, -13, 14, 21, 21};
    real.AllRegions.push_back(FKS::Region(5, 0));
    real.AllRegions.push_back(FKS::Region(4, 0));
    real.AllRegions.push_back(FKS::Region(5, 4));
    real.AllRegions.push_back(FKS::Region(4, 5));
    Phasespace::Phasespace ps_real;
    Phasespace::GenRealPhasespace(&ps_real, &ps, ifks, mother, xi, y, phi);

    double S40 = FKS::SFunction(ps_real, FKS::Region(4, 0), real, false);
    double S50 = FKS::SFunction(ps_real, FKS::Region(5, 0), real, false);
    double S54 = FKS::SFunction(ps_real, FKS::Region(5, 4), real, false);
    double S45 = FKS::SFunction(ps_real, FKS::Region(4, 5), real, false);
    if (mother == 0 && ifks == 4) {
        EXPECT_NEAR(S40, 1.0, 1e-6) << "S40 = " << S40 << " S50 = " << S50
                                    << " S54 = " << S54 << " S45 = " << S45;
    }
    if (mother == 0 && ifks == 5) {
        EXPECT_NEAR(S50, 1.0, 1e-6) << "S40 = " << S40 << " S50 = " << S50
                                    << " S54 = " << S54 << " S45 = " << S45;
    }
    if (mother == 4 && ifks == 5) {
        EXPECT_NEAR(S54 + S45, 1.0, 1e-6) << "S40 = " << S40 << " S50 = " << S50
                                          << " S54 = " << S54
                                          << " S45 = " << S45;
    }
    if (mother == 5 && ifks == 4) {
        EXPECT_NEAR(S54 + S45, 1.0, 1e-6) << "S40 = " << S40 << " S50 = " << S50
                                          << " S54 = " << S54
                                          << " S45 = " << S45;
    }
}

// test ISR
INSTANTIATE_TEST_CASE_P(RealTestParametersISR, SPropertiesTest,
                        Values(Params(BornPsParams(0.5, 0.5, 0.4, 0.3),
                                      FKSParams(0.6, 0.3, 0.4), 0),
                               Params(BornPsParams(0.2, 0.6, 0.9, 0.7),
                                      FKSParams(0.9, -0.4, 5.1), 0),
                               Params(BornPsParams(0.9, 0.1, 0.6, 0.5),
                                      FKSParams(0.01, 0.99, 0.1), 0)));

// test FSR with j = 2
INSTANTIATE_TEST_CASE_P(RealTestParametersFSR2, SPropertiesTest,
                        Values(Params(BornPsParams(0.5, 0.5, 0.4, 0.3),
                                      FKSParams(0.6, 0.3, 0.4), 2),
                               Params(BornPsParams(0.2, 0.6, 0.9, 0.7),
                                      FKSParams(0.9, -0.4, 5.1), 2),
                               Params(BornPsParams(0.9, 0.1, 0.6, 0.5),
                                      FKSParams(0.11, 0.99, 0.1), 2)));

// test FSR with j = 3
INSTANTIATE_TEST_CASE_P(RealTestParametersFSR3, SPropertiesTest,
                        Values(Params(BornPsParams(0.5, 0.5, 0.4, 0.3),
                                      FKSParams(0.6, 0.3, 0.4), 3),
                               Params(BornPsParams(0.2, 0.6, 0.9, 0.7),
                                      FKSParams(0.9, -0.4, 5.1), 3),
                               Params(BornPsParams(0.9, 0.1, 0.6, 0.5),
                                      FKSParams(0.1, 0.99, 0.1), 3)));

TEST(SumS, WjQCD_2) {

    FKS::Real_t real;
    real.Flavours = {2, -1, -13, 14, 21, 21};
    real.AllRegions.push_back(FKS::Region(5, 0));
    real.AllRegions.push_back(FKS::Region(4, 0));
    Phasespace::Phasespace ps_real;

    ps_real.X1 = 0.052623947831342881;
    ps_real.X2 = 0.014467093764705709;
    ps_real.S = 169000000;

    ps_real.Momenta[0].Set(179.3476611884779, 0, 0, 179.3476611884779);
    ps_real.Momenta[1].Set(179.3476611884779, 0, 0, -179.3476611884779);
    ps_real.Momenta[2].Set(75.224784573900919, 47.552754750211932,
                           57.572350825193595, 9.1064894617424912);
    ps_real.Momenta[3].Set(106.75125076102408, -75.687227309694265,
                           -22.442539390673975, -71.858232561841902);
    ps_real.Momenta[4].Set(51.478050427941611, -38.35280213933104,
                           13.682460442345597, -31.493531401966834);
    ps_real.Momenta[5].Set(125.24123661408927, 66.487274698813366,
                           -48.812271876865211, 94.245274502066252);
    ps_real.Jacobian = 441.05559073555139;

    double S50 = FKS::SFunction(ps_real, FKS::Region(5, 0), real, false);
    double S40 = FKS::SFunction(ps_real, FKS::Region(4, 0), real, false);

    double sum = S50 + S40;
    EXPECT_DOUBLE_EQ(sum, 1.0);
    if (HasFailure()) {

        std::cout << ps_real.Momenta[0].ToString() << "\n";
        std::cout << ps_real.Momenta[1].ToString() << "\n";
        std::cout << ps_real.Momenta[2].ToString() << "\n";
        std::cout << ps_real.Momenta[3].ToString() << "\n";
        std::cout << ps_real.Momenta[4].ToString() << "\n";
        std::cout << ps_real.Momenta[5].ToString() << "\n\n";
        printf("S50 = %.12g\n", S50);
        printf("S40 = %.12g\n", S40);
    }
}
