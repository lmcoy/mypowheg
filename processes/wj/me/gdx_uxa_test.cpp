#include "gtest/gtest.h"

#include <cstdio>

#include "gdx_uxa.h"
#include "parameters_sm.h"
#include "phasespace/phasespace.h"

#include "math/fourmomentum.h"

TEST(MatrixElementTest, PS1) {
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

    GDX_UXA me;
    // calculate matrix element more than once because Calculate determines in
    // the first calculations if a helicity can be neglected or not.
    for (int i = 0; i < 10; i++) {
        int perm[] = {0, 1, 2, 3, 4, 5};
        double m2 = me.Calculate(ps, perm, param);

        double pre = 7.592660551456415850371919e-14;

        ASSERT_NEAR(m2 / pre, 1.0, 1e-6) << "m2 = " << m2 << ", pre = " << pre;

        perm[0] = 1;
        perm[1] = 0;
        pre = 3.713548849191971800999538e-14;
        m2 = me.Calculate(ps, perm, param);

        ASSERT_NEAR(m2 / pre, 1.0, 1e-6);
    }
}