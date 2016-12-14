#include "gtest/gtest.h"

#include <cstdio>
#include <iostream>

#include "virtual_qed.h"

#include "parameters_sm.h"

#include "phasespace/phasespace.h"
#include "phasespace/twoparticlegenerator.h"
#include "fks/virtual.h"
#include "fks/scales.h"

TEST(Born, UXU_MUPMUM) {
    double SqrtS = 14000.0;
    Phasespace::Phasespace ps;
    Phasespace::TwoParticleGenerator ps_gen;
    double x[4] = { 0.6, 0.2, 0.6, 0.2 };
    double masses[2] = {0.0};
    ps_gen(&ps, 4, x, SqrtS*SqrtS, 2 , masses);

    Parameters_sm param;
    param.Set(1.1663900000000000E-005, 132.50700000000001,
                        91.188000000000002, 2.4414039999999999);

    param.SetAlphaS(0.118);

    int perm[4] = { 0, 1, 2, 3 };
    double born = -1.0;
    double virt = -1.0;
    double ct = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Alpha0, 0.0, &ct);
    Vfin_qxq_lxl(&born, &virt, ps.Momenta.data(), perm, CorrectionType::None, 2,
                 param.alpha, param.alphaS(), 91.188, ct);

    EXPECT_DOUBLE_EQ(born, 0.0025453361707912514);
}

TEST(Born, DXD_MUPMUM) {
    double SqrtS = 14000.0;
    Phasespace::Phasespace ps;
    Phasespace::TwoParticleGenerator ps_gen;
    double x[4] = { 0.6, 0.2, 0.6, 0.2 };
    double masses[2] = {0.0};
    ps_gen(&ps, 4, x, SqrtS*SqrtS, 2 , masses);

    Parameters_sm param;
    param.Set(1.1663900000000000E-005, 132.50700000000001,
                        91.188000000000002, 2.4414039999999999);

    param.SetAlphaS(0.118);

    int perm[4] = { 0, 1, 2, 3 };
    double born = -1.0;
    double virt = -1.0;
    double ct = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Alpha0, 0.0, &ct);
    Vfin_qxq_lxl(&born, &virt, ps.Momenta.data(), perm, CorrectionType::None, 1,
                 param.alpha, param.alphaS(), 91.188, ct);

    EXPECT_DOUBLE_EQ(born, 0.0013592103989550372);
}

TEST(Divergence, UXU_MUPMUM) {
    double SqrtS = 14000.0;
    Phasespace::Phasespace ps;
    Phasespace::TwoParticleGenerator ps_gen;
    double x[4] = { 0.1, 0.3, 0.9, 0.5 };
    double masses[2] = {0.0};
    ps_gen(&ps, 4, x, SqrtS*SqrtS, 2 , masses);

    Parameters_sm param;
    param.Set(1.1663900000000000E-005, 132.50700000000001, 91.188000000000002,
              2.4414039999999999);
    param.SetAlphaS(0.118);

    int perm[4] = { 0, 1, 2, 3 };
    double born0 = -1.0;
    double virt0 = -1.0;
    double ct = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Alpha0, 0.0, &ct);
    Vfin_qxq_lxl(&born0, &virt0, ps.Momenta.data(), perm, CorrectionType::EW, 2,
                 param.alpha, param.alphaS(), 91.188 * 91.188, ct);

    double inveps1 = 3.0;
    double born1 = -1.0;
    double virt1 = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Alpha0, inveps1, &ct);
    Vfin_qxq_lxl(&born1, &virt1, ps.Momenta.data(), perm, CorrectionType::EW, 2,
                 param.alpha, param.alphaS(), 91.188 * 91.188, ct);

    double inveps2 = 5.0;
    double born2 = -1.0;
    double virt2 = -1.0;
    InitVfin_qxq_lxl(param,  AlphaScheme::Alpha0,inveps2, &ct);
    Vfin_qxq_lxl(&born2, &virt2, ps.Momenta.data(), perm, CorrectionType::EW, 2,
                 param.alpha, param.alphaS(), 91.188 * 91.188, ct);

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

    int pdgs[] = { -2, 2, -13, 13 };
    double V1 =
        FKS::QED::EpsPole(4, pdgs, ps.Momenta.data(), born1, 91.188 * 91.188);
    double V0 = FKS::QED::Eps2Pole(4, pdgs, born1);

    EXPECT_NEAR(b + V1, 0.0, 1e-15) << "1/eps doesn't cancel";
    EXPECT_NEAR(a + V0, 0.0, 1e-15) << "1/eps^2 doesn't cancel";
}

TEST(Divergence, DXD_MUPMUM) {
    double SqrtS = 14000.0;
    Phasespace::Phasespace ps;
    Phasespace::TwoParticleGenerator ps_gen;
    double x[4] = { 0.1, 0.3, 0.9, 0.5 };
    double masses[2] = {0.0};
    ps_gen(&ps, 4, x, SqrtS*SqrtS, 2 , masses);

    Parameters_sm param;
    param.Set(1.1663900000000000E-005, 132.50700000000001, 91.188000000000002,
              2.4414039999999999);
    param.SetAlphaS(0.118);

    int perm[4] = { 0, 1, 2, 3 };
    double born0 = -1.0;
    double virt0 = -1.0;
    double ct = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Alpha0, 0.0, &ct);
    Vfin_qxq_lxl(&born0, &virt0, ps.Momenta.data(), perm, CorrectionType::EW, 1,
                 param.alpha, param.alphaS(), 91.188 * 91.188, ct);

    double inveps1 = 3.0;
    double born1 = -1.0;
    double virt1 = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Alpha0, inveps1, &ct);
    Vfin_qxq_lxl(&born1, &virt1, ps.Momenta.data(), perm, CorrectionType::EW, 1,
                 param.alpha, param.alphaS(), 91.188 * 91.188, ct);

    double inveps2 = 5.0;
    double born2 = -1.0;
    double virt2 = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Alpha0, inveps2, &ct);
    Vfin_qxq_lxl(&born2, &virt2, ps.Momenta.data(), perm, CorrectionType::EW, 1,
                 param.alpha, param.alphaS(), 91.188 * 91.188, ct);

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

    int pdgs[] = { -1, 1, -13, 13 };
    double V1 =
        FKS::QED::EpsPole(4, pdgs, ps.Momenta.data(), born1, 91.188 * 91.188);
    double V0 = FKS::QED::Eps2Pole(4, pdgs, born1);

    EXPECT_NEAR(b + V1, 0.0, 1e-15) << "1/eps doesn't cancel";
    EXPECT_NEAR(a + V0, 0.0, 1e-16) << "1/eps^2 doesn't cancel";
}

TEST(GmuScheme, UXU_MUPMUM) {
    // the virtual correction should be invariant under fermion mass changes in
    // the Gmu scheme.
    double SqrtS = 14000.0;
    Phasespace::Phasespace ps;
    Phasespace::TwoParticleGenerator ps_gen;
    double x[4] = { 0.6, 0.2, 0.6, 0.2 };
    double masses[2] = {0.0};
    ps_gen(&ps, 4, x, SqrtS*SqrtS, 2 , masses);

    Parameters_sm param;
    param.Set(1.1663900000000000E-005, 132.50700000000001, 91.188000000000002,
              2.4414039999999999);

    param.SetAlphaS(0.118);

    int perm[4] = { 0, 1, 2, 3 };
    double born1 = -1.0;
    double virt1 = -1.0;
    double ct = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Gmu, 0.0, &ct);
    Vfin_qxq_lxl(&born1, &virt1, ps.Momenta.data(), perm, CorrectionType::EW, 2,
                 param.alpha, param.alphaS(), 91.188 * 91.188, ct);

    Parameters_sm param2;
    param2.Set(1.1663900000000000E-005, 132.50700000000001, 91.188000000000002,
               2.4414039999999999);

    param2.MElectron = 2.0;
    param2.MMuon = 4.0;
    param2.MTau = 6.0;
    param2.MUQuark = 1.0;
    param2.MDQuark = 1.0;

    double born2 = -1.0;
    double virt2 = -1.0;
    InitVfin_qxq_lxl(param2, AlphaScheme::Gmu, 0.0, &ct);
    Vfin_qxq_lxl(&born2, &virt2, ps.Momenta.data(), perm, CorrectionType::EW, 2,
                 param.alpha, param.alphaS(), 91.188 * 91.188, ct);

    EXPECT_NEAR(virt1, virt2, 1e-15);
}

TEST(GmuScheme, DXD_MUPMUM) {
    // the virtual correction should be invariant under fermion mass changes in
    // the Gmu scheme.
    double SqrtS = 14000.0;
    Phasespace::Phasespace ps;
    Phasespace::TwoParticleGenerator ps_gen;
    double x[4] = { 0.6, 0.2, 0.6, 0.2 };
    double masses[2] = {0.0};
    ps_gen(&ps, 4, x, SqrtS*SqrtS, 2 , masses);

    Parameters_sm param;
    param.Set(1.1663900000000000E-005, 132.50700000000001, 91.188000000000002,
              2.4414039999999999);

    param.SetAlphaS(0.118);

    int perm[4] = { 0, 1, 2, 3 };
    double born1 = -1.0;
    double virt1 = -1.0;
    double ct = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Gmu, 0.0, &ct);
    Vfin_qxq_lxl(&born1, &virt1, ps.Momenta.data(), perm, CorrectionType::EW, 1,
                 param.alpha, param.alphaS(), 91.188 * 91.188, ct);

    Parameters_sm param2;
    param2.Set(1.1663900000000000E-005, 132.50700000000001, 91.188000000000002,
               2.4414039999999999);

    param2.MElectron = 2.0;
    param2.MMuon = 4.0;
    param2.MTau = 6.0;
    param2.MUQuark = 1.0;
    param2.MDQuark = 1.0;
    param2.SetAlphaS(0.118);

    double born2 = -1.0;
    double virt2 = -1.0;
    InitVfin_qxq_lxl(param2,  AlphaScheme::Gmu, 0.0, &ct);
    Vfin_qxq_lxl(&born2, &virt2, ps.Momenta.data(), perm, CorrectionType::EW, 1,
                 param.alpha, param.alphaS(), 91.188 * 91.188, ct);

    EXPECT_NEAR(virt1, virt2, 1e-15);
}

TEST(MuDependence, UXU_MUPMUM) {
    double SqrtS = 14000.0;
    Phasespace::Phasespace ps;
    Phasespace::TwoParticleGenerator ps_gen;
    double x[4] = { 0.6, 0.2, 0.6, 0.2 };
    double masses[2] = {0.0};
    ps_gen(&ps, 4, x, SqrtS*SqrtS, 2 , masses);

    double sqrts = sqrt(ps.X1 * ps.X2 * ps.S);

    Parameters_sm param;
    param.Set(1.1663900000000000E-005, 132.50700000000001, 91.188000000000002,
              2.4414039999999999);

    param.SetAlphaS(0.118);

    int perm[4] = { 0, 1, 2, 3 };
    double born1 = -1.0;
    double virt1 = -1.0;
    double ct = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Gmu, 0.0, &ct);
    Vfin_qxq_lxl(&born1, &virt1, ps.Momenta.data(), perm, CorrectionType::EW, 2,
                 param.alpha, param.alphaS(), 91.188 * 91.188, ct);

    int pdgs[] = { -2, 2, -13, 13 };
    FKS::Scales scales;
    scales.muF = 91.188;
    scales.mu = 91.188;
    scales.Q2 = scales.mu * scales.mu;
    Util::Matrix2 unused(1);
    double v1 = FKS::QED::Virtual(4, pdgs, ps.Momenta.data(), born1, unused,
                                  virt1, param.alpha, sqrts, scales);

    double born2 = -1.0;
    double virt2 = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Gmu, 0.0, &ct);
    Vfin_qxq_lxl(&born2, &virt2, ps.Momenta.data(), perm, CorrectionType::EW, 2,
                 param.alpha, param.alphaS(), 191.188 * 191.188, ct);

    scales.mu = 191.188;
    scales.Q2 = scales.mu * scales.mu;
    double v2 = FKS::QED::Virtual(4, pdgs, ps.Momenta.data(), born2, unused,
                                  virt2, param.alpha, sqrts, scales);

    EXPECT_NEAR(v1, v2, 1e-15);

    double born3 = -1.0;
    double virt3 = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Gmu, 0.0, &ct);
    Vfin_qxq_lxl(&born3, &virt3, ps.Momenta.data(), perm, CorrectionType::EW, 2,
                 param.alpha, param.alphaS(), 1191.188 * 1191.188, ct);

    scales.mu = 1191.188;
    scales.Q2 = scales.mu * scales.mu;
    double v3 = FKS::QED::Virtual(4, pdgs, ps.Momenta.data(), born3, unused,
                                  virt3, param.alpha, sqrts, scales);
    EXPECT_DOUBLE_EQ(v1, v3);
}

TEST(MuDependence, DXD_MUPMUM) {
    double SqrtS = 14000.0;
    Phasespace::Phasespace ps;
    Phasespace::TwoParticleGenerator ps_gen;
    double x[4] = { 0.6, 0.2, 0.6, 0.2 };
    double masses[2] = {0.0};
    ps_gen(&ps, 4, x, SqrtS*SqrtS, 2 , masses);

    double sqrts = sqrt(ps.X1 * ps.X2 * ps.S);

    Parameters_sm param;
    param.Set(1.1663900000000000E-005, 132.50700000000001, 91.188000000000002,
              2.4414039999999999);

    param.SetAlphaS(0.118);

    int perm[4] = { 0, 1, 2, 3 };
    double born1 = -1.0;
    double virt1 = -1.0;
    double ct = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Gmu, 0.0, &ct);
    Vfin_qxq_lxl(&born1, &virt1, ps.Momenta.data(), perm, CorrectionType::EW, 1,
                 param.alpha, param.alphaS(), 91.188 * 91.188, ct);

    int pdgs[] = { -1, 1, -13, 13 };
    FKS::Scales scales;
    scales.muF = 91.188;
    scales.mu = 91.188;
    scales.Q2 = scales.mu * scales.mu;
    Util::Matrix2 unused(1);
    double v1 = FKS::QED::Virtual(4, pdgs, ps.Momenta.data(), born1, unused,
                                  virt1, param.alpha, sqrts, scales);

    double born2 = -1.0;
    double virt2 = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Gmu, 0.0, &ct);
    Vfin_qxq_lxl(&born2, &virt2, ps.Momenta.data(), perm, CorrectionType::EW, 1,
                 param.alpha, param.alphaS(), 191.188 * 191.188, ct);

    scales.mu = 191.188;
    scales.Q2 = scales.mu * scales.mu;
    double v2 = FKS::QED::Virtual(4, pdgs, ps.Momenta.data(), born2, unused,
                                  virt2, param.alpha, sqrts, scales);

    EXPECT_DOUBLE_EQ(v1, v2);

    double born3 = -1.0;
    double virt3 = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Gmu, 0.0, &ct);
    Vfin_qxq_lxl(&born3, &virt3, ps.Momenta.data(), perm, CorrectionType::EW, 1,
                 param.alpha, param.alphaS(), 1191.188 * 1191.188, ct);

    scales.mu = 1191.188;
    scales.Q2 = scales.mu * scales.mu;
    double v3 = FKS::QED::Virtual(4, pdgs, ps.Momenta.data(), born3, unused,
                                  virt3, param.alpha, sqrts, scales);
    EXPECT_DOUBLE_EQ(v1, v3);
}

TEST(Vfin_qxq_lxl, QCD) {
    double SqrtS = 14000.0;
    Phasespace::Phasespace ps;
    Phasespace::TwoParticleGenerator ps_gen;
    double x[4] = { 0.3, 0.2, 0.6, 0.2 };
    double masses[2] = {0.0};
    ps_gen(&ps, 4, x, SqrtS*SqrtS, 2 , masses);

    double sqrts = sqrt(ps.X1 * ps.X2 * ps.S);
    double s = sqrts * sqrts;

    Parameters_sm param;
    param.Set(1.1663900000000000E-005, 132.50700000000001, 91.188000000000002,
              2.4414039999999999);

    param.SetAlphaS(0.118);

    int perm[4] = { 0, 1, 2, 3 };
    double born1 = -1.0;
    double virt1 = -1.0;
    double ct = -1.0;
    InitVfin_qxq_lxl(param, AlphaScheme::Alpha0, 0.0, &ct);
    ASSERT_DOUBLE_EQ(ct, 0.0);
    Vfin_qxq_lxl(&born1, &virt1, ps.Momenta.data(), perm, CorrectionType::QCD,
                 2, param.alpha, param.alphaS(), 91.188 * 91.188, ct);

    // Q^2 is set to mu^2!
    double tmp = log(91.188 * 91.188 / s);
    double vfin = 4.0 / 3.0 * (M_PI * M_PI - 8.0 - tmp * (tmp + 3.0)) * born1;
    EXPECT_NEAR(virt1, vfin, 1e-15);
}
