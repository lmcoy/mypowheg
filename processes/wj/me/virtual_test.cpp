#include "gtest/gtest.h"

#include <cfenv>

#include "virtual.h"

#include "fks/scales.h"
#include "fks/virtual.h"
#include "parameters_sm.h"
#include "phasespace/phasespace.h"
#include "util/matrix.h"

TEST(EpsPoles, QCD) {
#ifdef _GNU_SOURCE
// feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
    Parameters_sm par;
    par.Set(1.1663900000000000E-005, 132.50700000000001, 91.188000000000002,
            2.4414039999999999);
    par.SetAlphaS(0.12);

    int perm[] = {0, 1, 3, 2, 4};

    Phasespace::Phasespace ps;

    ps.Momenta[0].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, 7.500000000000000000000000e+02);
    ps.Momenta[1].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, -7.500000000000000000000000e+02);
    ps.Momenta[2].Set(
        6.878681818281603455034201e+02, 2.541798304645197390527755e+02,
        5.694804931172980104747694e+02, -2.902537119753787351328356e+02);
    ps.Momenta[3].Set(
        5.460999311052265738908318e+02, -2.749480393978777570396232e+01,
        -5.215564519790506210483727e+02, 1.595244116380622187989502e+02);
    ps.Momenta[4].Set(
        2.660318870666130806057481e+02, -2.266850265247319953232363e+02,
        -4.792404113824733968840519e+01, 1.307293003373165447555948e+02);

    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;

    double mu = 80.3;

    double born0 = 0.0;
    double virt0 = 0.0;
    Wj_InitDimReg(par, 0.0);
    Wj_MatrixElement(0, &born0, &virt0, ps.Momenta.data(), perm, par.alpha,
                     par.alphaS(), mu, Type::QCD, Scheme::Alpha0);

    double inveps1 = 3.0;
    double born1 = 0.0;
    double virt1 = 0.0;
    Wj_InitDimReg(par, inveps1);
    Wj_MatrixElement(0, &born1, &virt1, ps.Momenta.data(), perm, par.alpha,
                     par.alphaS(), mu, Type::QCD, Scheme::Alpha0);

    double inveps2 = 5.0;
    double born2 = 0.0;
    double virt2 = 0.0;
    Wj_InitDimReg(par, inveps2);
    Wj_MatrixElement(0, &born2, &virt2, ps.Momenta.data(), perm, par.alpha,
                     par.alphaS(), mu, Type::QCD, Scheme::Alpha0);

    double inveps12 = inveps1 * inveps1;
    double inveps22 = inveps2 * inveps2;

    // extract poles from
    // V[0] = virt0
    // virt1 = V[inveps1] = a * inveps1^2 + b * inveps1 + virt0
    // virt2 = V[inveps2] = a * inveps2^2 + b * inveps2 + virt0
    double b =
        (virt0 * (inveps22 - inveps12) + inveps12 * virt2 - inveps22 * virt1) /
        (inveps1 * inveps2 * (inveps1 - inveps2));

    double a = (inveps2 * virt1 - inveps1 * virt2 + virt0 * inveps1 -
                virt0 * inveps2) /
               (inveps12 * inveps2 - inveps1 * inveps22);

    // test pole extraction
    EXPECT_NEAR(inveps12 * a + inveps1 * b + virt0, virt1, 1e-15);
    EXPECT_NEAR(inveps22 * a + inveps2 * b + virt0, virt2, 1e-15);

    double CF = 4.0 / 3.0;
    double CA = 3;
    double NC = 3.0;
    EXPECT_NEAR(a / born0, -2.0 * CF - CA, 1e-13);

    int pdgs[] = {2, -1, 14, -13, 21};
    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    ColorCorr.Set(0, 1, -1 / (2 * NC) * born0);
    ColorCorr.Set(1, 0, -1 / (2 * NC) * born0);
    ColorCorr.Set(4, 1, NC / 2.0 * born0);
    ColorCorr.Set(1, 4, NC / 2.0 * born0);
    ColorCorr.Set(4, 0, NC / 2.0 * born0);
    ColorCorr.Set(0, 4, NC / 2.0 * born0);
    double epspole = FKS::QCD::EpsPole(5, pdgs, ps.Momenta.data(), born0,
                                       ColorCorr, mu * mu);
    EXPECT_NEAR(-b / epspole, 1.0, 1e-15);
}

TEST(CompareMassAndDimReg, EW) {
#ifdef _GNU_SOURCE
// feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
    Parameters_sm par;
    par.Set(1.1663900000000000E-005, 132.50700000000001, 91.188000000000002,
            2.4414039999999999);
    par.SetAlphaS(0.12);

    Phasespace::Phasespace ps;

    ps.S = 196000000;
    ps.X1 = 0.8540346744275702;
    ps.X2 = 0.477559445749165;
    ps.Jacobian = 15691.1274156559;
    ps.N = 3;
    ps.Momenta[0].Set(4470.432189703782, 0, 0, 4470.432189703782);
    ps.Momenta[1].Set(4470.432189703782, 0, 0, -4470.432189703782);
    ps.Momenta[2].Set(3123.569155464953, 1748.007090708135, -583.1801005865702,
                      2522.113488740338);
    ps.Momenta[3].Set(2334.690107086365, 1153.833353480185, -1146.097665913952,
                      -1675.084066162113);
    ps.Momenta[4].Set(3482.605116856246, -2901.840444188321, 1729.277766500523,
                      -847.0294225782245);
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;

    double mu = 1.0;
    double muF = mu;

    double born0 = 0.0;
    double virt0 = 0.0;
    int perm[] = {0, 1, 3, 2, 4};
    Wj_InitDimReg(par, 0.0);
    Wj_MatrixElement(0, &born0, &virt0, ps.Momenta.data(), perm, par.alpha,
                     par.alphaS(), mu, Type::EW, Scheme::Gmu);

    int pdgs[] = {2, -1, -11, 12, 0};
    double m[] = {par.MUQuark, par.MDQuark, par.MElectron, 0.0, 0.0};
    Util::Matrix2 unused(1, born0);
    FKS::Scales scales;
    scales.muF = muF;
    scales.mu = mu;
    scales.Q2 = scales.mu * scales.mu;
    double s = ps.X1 * ps.X2 * ps.S;

    double dimreg = FKS::QED::Virtual(5, pdgs, ps.Momenta.data(), born0, unused,
                                      virt0, par.alpha, sqrt(s), scales);

    double lambda = 1.0;
    Wj_InitMassReg(par, lambda);
    double born1 = 0.0;
    double virt1 = 0.0;
    Wj_MatrixElement(0, &born1, &virt1, ps.Momenta.data(), perm, par.alpha,
                     par.alphaS(), mu, Type::EW, Scheme::Gmu);

    double v = FKS::QED::VirtualMReg(ps, pdgs, m, lambda, muF * muF);
    double coup = par.alpha * M_1_PI * 0.5;
    double massreg = coup * (born1 * v + virt1);

    EXPECT_NEAR(dimreg / massreg, 1.0, 1e-14);
}

TEST(CheckDeltaInPaper, Check) {
    Parameters_sm par;
    par.Set(1.1663900000000000E-005, 132.50700000000001, 91.188000000000002,
            2.4414039999999999);
    par.SetAlphaS(0.12);

    Phasespace::Phasespace ps;

    ps.S = 196000000;
    ps.X1 = 0.8540346744275702;
    ps.X2 = 0.477559445749165;
    ps.Jacobian = 15691.1274156559;
    ps.N = 3;
    ps.Momenta[0].Set(4470.432189703782, 0, 0, 4470.432189703782);
    ps.Momenta[1].Set(4470.432189703782, 0, 0, -4470.432189703782);
    ps.Momenta[2].Set(3123.569155464953, 1748.007090708135, -583.1801005865702,
                      2522.113488740338);
    ps.Momenta[3].Set(2334.690107086365, 1153.833353480185, -1146.097665913952,
                      -1675.084066162113);
    ps.Momenta[4].Set(3482.605116856246, -2901.840444188321, 1729.277766500523,
                      -847.0294225782245);
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;

    double mu = 2.0;
    double muF = mu;

    double born0 = 0.0;
    double virt0 = 0.0;
    int perm[] = {0, 1, 3, 2, 4};
    Wj_InitDimReg(par, 0.0);
    Wj_MatrixElement(0, &born0, &virt0, ps.Momenta.data(), perm, par.alpha,
                     par.alphaS(), mu, Type::EW, Scheme::Gmu);

    int pdgs[] = {2, -1, -11, 12, 0};
    double m[] = {par.MUQuark, par.MDQuark, par.MElectron, 0.0, 0.0};
    Util::Matrix2 unused(1, born0);
    FKS::Scales scales;
    scales.muF = muF;
    scales.mu = mu;
    scales.Q2 = scales.mu * scales.mu;
    double s = ps.X1 * ps.X2 * ps.S;

    double lambda = 0.1;
    Wj_InitMassReg(par, lambda);
    double born1 = 0.0;
    double virt1 = 0.0;
    Wj_MatrixElement(0, &born1, &virt1, ps.Momenta.data(), perm, par.alpha,
                     par.alphaS(), mu, Type::EW, Scheme::Gmu);

    double v = FKS::QED::VirtualMReg(ps, pdgs, m, lambda, muF * muF);
    double coup = par.alpha * M_1_PI * 0.5;

    double q[] = {2.0/3.0, 1.0/3.0, 1.0, 0.0, 0.0};
    double Delta = 0.0;
    for( int i = 0; i < 5; i++) {
        double r = q[i]*q[i];
        if (m[i] == 0.0) continue;
        r *= (-2.0*log(lambda/mu)-log(m[i]/mu)+0.5 *pow( 2*log(m[i]/mu), 2) + 2 + M_PI*M_PI/6.0);
        Delta += r*born1;
    }

    double Bij[5][5] = {{0, -2.0 / 9.0 * born1, 2.0 / 3.0 * born1, 0.0, 0.0},
                        {-2.0 * born1, 0.0, 1. / 3.0 * born1, 0.0, 0.0},
                        {2.0 / 3.0 * born1, 1.0 / 3.0 * born1, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0, 0.0, 0.0}};

    for( int i = 0; i < 5; i++) {
        for(int j = i+1; j < 5; j++) {
            if (m[i] == 0.0) continue;
            if (m[j] == 0.0) continue;
            double pipj = ps.Momenta[i].Dot(ps.Momenta[j]);
            double r = 0.5 * log( m[i]*m[i]*m[j]*m[j]/4.0/pipj/pipj)*2.0*log(lambda/mu);
            Delta -= 2.0*Bij[i][j] * r;
        }
    }
    EXPECT_NEAR((virt1-Delta)/virt0, 1.0, 1e-13);
}
