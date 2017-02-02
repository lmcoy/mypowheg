#include "gtest/gtest.h"

#include <complex>
#include <cstdio>

#include "dxdx_uxdx.h"
#include "parameters_sm.h"
#include "phasespace/phasespace.h"

#include "math/fourmomentum.h"

TEST(DXDX_UXDX, PS1) {
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

    DXDX_UXDX me;

    int perm1[] = { 0, 1, 2, 3, 4, 5 };
    double m2 = 0.0;

    // b~ b~ > mu+ vm c~ b~
    m2 = me.Calculate(ps, perm1, param) * Vcb2;
    ASSERT_NEAR(m2 / 3.224466795828124779053593e-15, 1.0, 1e-6);

    // b~ b~ > mu+ vm u~ b~
    m2 = me.Calculate(ps, perm1, param) * Vub2;
    ASSERT_NEAR(m2 / 2.188406222130145373509527e-17, 1.0, 1e-6);

    // s~ s~ > mu+ vm c~ s~
    m2 = me.Calculate(ps, perm1, param) * Vcs2;
    ASSERT_NEAR(m2 / 1.820794861138010305516519e-12, 1.0, 1e-6);

    // s~ s~ > mu+ vm u~ s~
    m2 = me.Calculate(ps, perm1, param) * Vus2;
    ASSERT_NEAR(m2 / 9.729999190167486299667012e-14, 1.0, 1e-6);

    // d~ d~ > mu+ vm c~ d~
    m2 = me.Calculate(ps, perm1, param) * Vcd2;
    ASSERT_NEAR(m2 / 9.729999190167486299667012e-14, 1.0, 1e-6);

    // d~ d~ > mu+ vm u~ d~
    m2 = me.Calculate(ps, perm1, param) * Vud2;
    ASSERT_NEAR(m2 / 1.820794861138010305516519e-12, 1.0, 1e-6);
}
