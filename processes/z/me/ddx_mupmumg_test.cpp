#include "gtest/gtest.h"

#include <cstdio>

#include "phasespace/phasespace.h"
#include "ddx_mupmumg.h"
#include "parameters_sm.h"

#include "math/fourmomentum.h"


TEST(MatrixElementTest, PS1) {
    Parameters_sm param;
    double aew =  1.325070e+02;
    double mz = 9.1188e+01;
    param.MW = 8.04190024457561634108e+01;
    param.WidthW = 2.0476;

    double Gf = M_PI / sqrt(2.0) / aew / param.MW / param.MW *
                (1.0 / (1.0 - param.MW * param.MW / mz / mz));

    param.Set(Gf, aew, mz, 2.441404);

    param.SetAlphaS(0.11799999999999999);

    Phasespace::Phasespace ps;

    ps.Momenta[0]
        .Set(7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
             0.000000000000000000000000e+00, 7.500000000000000000000000e+02);
    ps.Momenta[1]
        .Set(7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
             0.000000000000000000000000e+00, -7.500000000000000000000000e+02);
    ps.Momenta[2]
        .Set(6.878681818281603455034201e+02, 2.541798304645197390527755e+02,
             5.694804931172980104747694e+02, -2.902537119753787351328356e+02);
    ps.Momenta[3]
        .Set(5.460999311052265738908318e+02, -2.749480393978777570396232e+01,
             -5.215564519790506210483727e+02, 1.595244116380622187989502e+02);
    ps.Momenta[4]
        .Set(2.660318870666130806057481e+02, -2.266850265247319953232363e+02,
             -4.792404113824733968840519e+01, 1.307293003373165447555948e+02);

    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;

    ME_ddx_mupmumg me;
    int perm[] = { 0, 1, 2, 3, 4 };
    double m2 = me.Calculate(ps, perm, param);
    double pre = 2.7595429126450211e-07;

    EXPECT_NEAR(m2, pre, 1e-12);

    perm[0] = 1;
    perm[1] = 0;
    pre = 7.7980332583945798e-08;
    m2 = me.Calculate(ps, perm, param);

    EXPECT_NEAR(m2, pre, 1e-13);
}

