#include "gtest/gtest.h"

#include <cfenv>

#include "virtual.h"

#include "phasespace/phasespace.h"
#include "phasespace/twoparticlegenerator.h"
#include "parameters_sm.h"

TEST(test1, test1) {
#ifdef _GNU_SOURCE
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
  Parameters_sm par;
  par.Set(1.1663900000000000E-005, 132.50700000000001,
                        91.188000000000002, 2.4414039999999999);
  par.SetAlphaS(0.12);
  
  

  int perm[] = {0,1,2,3};
  
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
    
    double born0 = 0.0;
    double virt0 = 0.0;
    W_InitDimReg(par, 0.0);
    W_MatrixElement(&born0, &virt0,
                     ps.Momenta.data(), perm,
                    par.alpha, par.alphaS(), 80.3, Type::QCD, Scheme::Gmu);
    
    double inveps1 = 3.0;
    double born1 = 0.0;
    double virt1 = 0.0;
    W_InitDimReg(par, inveps1);
    W_MatrixElement(&born1, &virt1,
                     ps.Momenta.data(), perm,
                    par.alpha, par.alphaS(), 80.3, Type::QCD, Scheme::Gmu);
    
    double inveps2 = 5.0;
    double born2 = 0.0;
    double virt2 = 0.0;
    W_InitDimReg(par, inveps2);
    W_MatrixElement(&born2, &virt2,
                     ps.Momenta.data(), perm,
                    par.alpha, par.alphaS(), 80.3, Type::QCD, Scheme::Gmu);
    
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
    EXPECT_NEAR(a/born0, -2.0*CF, 1e-13);
    
    double s = 2.0 * ps.Momenta[0].Dot(ps.Momenta[1]);
    EXPECT_NEAR(b/born0, CF*(-3.0 + 2.0*log(s/80.3/80.3)),1e-13);
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
    W_InitDimReg(param, 0.0);
    W_MatrixElement(&born1, &virt1,
                    ps.Momenta.data(), perm,
                    param.alpha, param.alphaS(), 91.188, Type::QCD, Scheme::Alpha0);
    
    // Q^2 is set to mu^2!
    double tmp = log(91.188 * 91.188 / s);
    double vfin = 4.0 / 3.0 * (M_PI * M_PI - 8.0 - tmp * (tmp + 3.0)) * born1;
    EXPECT_NEAR(virt1, vfin, 1e-15);
}
