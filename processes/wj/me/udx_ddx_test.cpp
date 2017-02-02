#include "gtest/gtest.h"

#include <complex>
#include <cstdio>

#include "parameters_sm.h"
#include "phasespace/phasespace.h"
#include "udx_ddx.h"

#include "math/fourmomentum.h"

TEST(UDX_DDX, PS1) {
    Parameters_sm param;
    double aew = 1.325070e+02;
    double mz = 9.1188e+01;
    param.MW = 8.04190024457561634108e+01;
    param.WidthW = 2.0476;

    double Gf = M_PI / sqrt(2.0) / aew / param.MW / param.MW *
                (1.0 / (1.0 - param.MW * param.MW / mz / mz));

    param.Set(Gf, aew, mz, 2.441404);

    param.SetAlphaS(0.11799999999999999);

    Phasespace::Phasespace ps;

    ps.Momenta[0].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, 7.500000000000000000000000e+02);
    ps.Momenta[1].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, -7.500000000000000000000000e+02);
    ps.Momenta[2].Set(
        1.328269995817544781857578e+02, -3.315103543153497867024271e+01,
        6.012052978752799248240990e+01, -1.137081464354049415987902e+02);
    ps.Momenta[3].Set(
        4.924941288406478179240366e+02, -1.557744178251844289206929e+02,
        -4.529006330843101295613451e+02, 1.147423820807488255013595e+02);
    ps.Momenta[4].Set(
        2.285371642011459698551334e+02, -1.588214394998882710297039e+02,
        -1.465644574904635533130204e+02, 7.432257784018922563973319e+01);
    ps.Momenta[5].Set(
        6.461417073764514498179778e+02, 3.477468927566076786206395e+02,
        5.393445607872457685516565e+02, -7.535681348553316638572142e+01);

    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Masses[3] = 0.0;

    double Vud = 9.74619954999999982270e-01;
    double Vus = 2.25300000000000000266e-01;
    auto Vub = std::complex<double>(1.219744e-03, -3.15100661526525631026e-03);
    double Vcd = 2.25300000000000000266e-01;
    double Vcs = 9.74619954999999982270e-01;
    double Vcb = 4.10141527200000022280e-02;

    double Vud2 = Vud * Vud;
    double Vus2 = Vus * Vus;
    double Vub2 = std::norm(Vub);
    double Vcd2 = Vcd * Vcd;
    double Vcs2 = Vcs * Vcs;
    double Vcb2 = Vcb * Vcb;

    UDX_DDX me;

    int perm1[] = { 0, 1, 2, 3, 5, 4 };
    double m2 = 0.0;

    // u b~ > mu+ vm b~ b
    m2 = me.Calculate(ps, perm1, param) * Vub2;
    ASSERT_NEAR(m2 / 5.881135636801859870528818e-17, 1.0, 1e-6);

    // u s~ > mu+ vm s~ s
    m2 = me.Calculate(ps, perm1, param) * Vus2;
    ASSERT_NEAR(m2 / 2.614845653639532519026610e-13, 1.0, 1e-6);

    // u d~ > mu+ vm d~ d
    m2 = me.Calculate(ps, perm1, param) * Vud2;
    ASSERT_NEAR(m2 / 4.893214722594411567720175e-12, 1.0, 1e-6);

    // c b~ > mu+ vm b~ b
    m2 = me.Calculate(ps, perm1, param) * Vcb2;
    ASSERT_NEAR(m2 / 8.665450861390081254413362e-15, 1.0, 1e-6);

    // c s~ > mu+ vm s~ s
    m2 = me.Calculate(ps, perm1, param) * Vcs2;
    ASSERT_NEAR(m2 / 4.893214722594411567720175e-12, 1.0, 1e-6);

    // c d~ > mu+ vm d~ d
    m2 = me.Calculate(ps, perm1, param) * Vcd2;
    ASSERT_NEAR(m2 / 2.614845653639532519026610e-13, 1.0, 1e-6);

    int perm2[] = { 1, 0, 2, 3, 5, 4 };
    // b~ u > mu+ vm b~ b
    m2 = me.Calculate(ps, perm2, param) * Vub2;
    ASSERT_NEAR(m2 / 8.589718327491816092376938e-17, 1.0, 1e-6);

    // b~ c > mu+ vm b~ b
    m2 = me.Calculate(ps, perm2, param) * Vcb2;
    ASSERT_NEAR(m2 / 1.265636208324876748893943e-14, 1.0, 1e-6);

    // s~ u > mu+ vm s~ s
    m2 = me.Calculate(ps, perm2, param) * Vus2;
    ASSERT_NEAR(m2 / 3.819124234115424835365740e-13, 1.0, 1e-6);

    // s~ c > mu+ vm s~ s
    m2 = me.Calculate(ps, perm2, param) * Vcs2;
    ASSERT_NEAR(m2 / 7.146806123634737291008950e-12, 1.0, 1e-6);

    // d~ u > mu+ vm d~ d
    m2 = me.Calculate(ps, perm2, param) * Vud2;
    ASSERT_NEAR(m2 / 7.146806123634737291008950e-12, 1.0, 1e-6);

    // d~ c > mu+ vm d~ d
    m2 = me.Calculate(ps, perm2, param) * Vcd2;
    ASSERT_NEAR(m2 / 3.819124234115424835365740e-13, 1.0, 1e-6);
}
