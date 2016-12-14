#include "gtest/gtest.h"

#include <cstdio>

#include "phasespace/phasespace.h"
#include "uux_mupmuma.h"
#include "parameters_sm.h"

#include "fks/limits.h"

TEST(MatrixElementTest, PS1) {
    Parameters_sm param;
    double aew =  1.325070e+02;
    double mz = 9.1188e+01;
    param.MW = 8.0419e+01;
    param.WidthW = 2.0476;

    double Gf = M_PI / sqrt(2.0) / aew / param.MW / param.MW *
                (1.0 / (1.0 - param.MW * param.MW / mz / mz));

    param.Set(Gf, aew, mz, 2.441404);

    Parameters_alphaS param_aS;
    param_aS.Set(0.11799999999999999);

    Phasespace::Phasespace ps;

    ps.Momenta[0].Set(7.500000000000e+02, 0.000000000000e+00,
                      0.000000000000e+00, 7.500000000000e+02);
    ps.Momenta[1].Set(7.500000000000e+02, 0.000000000000e+00,
                      0.000000000000e+00, -7.500000000000e+02);
    ps.Momenta[2].Set(6.878681818282e+02, 2.541798304645e+02,
                      5.694804931173e+02, -2.902537119754e+02);
    ps.Momenta[3].Set(5.460999311052e+02, -2.749480393979e+01,
                      -5.215564519791e+02, 1.595244116381e+02);
    ps.Momenta[4].Set(2.660318870666e+02, -2.266850265247e+02,
                      -4.792404113825e+01, 1.307293003373e+02);

    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;

    ME_uux_mupmuma me;
    int perm[] = { 0, 1, 2, 3, 4 };
    double m2 = me.Calculate(ps, perm, param, param_aS);
    double pre = 3.1483409516656004e-08;

    EXPECT_NEAR(m2, pre, 1e-16);

    perm[0] = 1;
    perm[1] = 0;
    pre = 7.2700219903302722e-09;
    m2 = me.Calculate(ps, perm, param, param_aS);

    EXPECT_NEAR(m2, pre, 1e-16);
}
