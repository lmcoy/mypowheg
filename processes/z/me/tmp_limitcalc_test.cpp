#include "phasespace/phasespace.h"
#include "uux_mupmuma.h"
#include "ddx_mupmuma.h"
#include "ddx_mupmum.h"
#include "uux_mupmum.h"
#include "parameters_sm.h"

#include "phasespace/realphasespace.h"
#include "phasespace/twoparticlegenerator.h"
#include "math/fourmomentum.h"

#include "fks/limits.h"

#include <cstdio>

int main() {
    Parameters_sm param;
    param.Set(1.1663900000000000E-005, 132.50700000000001, 91.188000000000002,
              2.4414039999999999);

    Parameters_alphaS param_aS;
    param_aS.Set(0.11799999999999999);

    Phasespace::Phasespace<2> ps(14000.0*14000.0, { { 0.0, 0.0 } });
    Phasespace::Phasespace<3> ps_real(14000.0*14000.0, { { 0.0, 0.0, 0.0 } });

    double params[] = { 0.3, 0.4, 0.18, 0.36 };
    Phasespace::TwoParticleGenerator ps_gen;
    ps_gen(&ps, 4, params);
    std::cout << "genps\n";

    Math::FourMomentum mll = ps.Momenta[2].Plus(ps.Momenta[3]);
    double mll2 = mll.Dot(mll);
    std::cout << "m_ll = " << sqrt(mll2) << "\n";

    double y = -0.2;
    double xi = 1e-9;
    double phi = 3.455751919;
    Phasespace::GenRealPhasespaceISR<2>(&ps_real, &ps, xi, y, phi);

    Math::FourMomentum mll_r = ps_real.Momenta[2].Plus(ps_real.Momenta[3]);
    double mll2_r = mll_r.Dot(mll_r);
    std::cout << "m_ll_r = " << sqrt(mll2_r) << "\n";

    std::cout << "born\n";
    std::cout << ps.Momenta[0].ToString() << "\n";
    std::cout << ps.Momenta[1].ToString() << "\n";
    std::cout << ps.Momenta[2].ToString() << "\n";
    std::cout << ps.Momenta[3].ToString() << "\n";

    std::cout << "real\n";
    std::cout << ps_real.Momenta[0].ToString() << "\n";
    std::cout << ps_real.Momenta[1].ToString() << "\n";
    std::cout << ps_real.Momenta[2].ToString() << "\n";
    std::cout << ps_real.Momenta[3].ToString() << "\n";
    std::cout << ps_real.Momenta[4].ToString() << "\n";

    Math::FourMomentum knp1(ps_real.Momenta[4]);
    knp1.Scale(1.0 / xi);

    std::cout << "scaled = " << knp1.ToString() << "\n";
    double y_calc = 1. - ps_real.Momenta[2].Dot(ps_real.Momenta[4]) /
                             ps_real.Momenta[2].E() / ps_real.Momenta[4].E();
    std::cout << "y_calc = " << y_calc << "\n";

    ME_uux_mupmuma me;
    int perm[] = { 1, 0, 2, 3, 4 };
    double m2 =
        (1 - y * y) * xi * xi * me.Calculate(ps_real, perm, param, param_aS);
    ME_uux_mupmum meb;
    double m2b = meb.Calculate(ps, perm, param, param_aS);
    printf("real me = %.16g\n", m2);
    printf("born me = %.16g\n", m2b);

    double s = ps.X1 * ps.X2 * ps.S;
    printf("s = %.16g\n", s);
    double s_real = ps_real.X1 * ps_real.X2 * ps_real.S;
    printf("s_r = %.16g\n", s_real);

    int pdg[] = { -2, 2, -13, 13 };
    //double limit = FKS::QED::SoftCollinearLimitISR(2, 2, s, 1.0 / 132.507, m2b);
     double limit = FKS::QED::SoftLimit(4, ps.Momenta.data(), pdg, s, 0,
     1.0/132.507, m2b, y, phi);
    printf("limit me = %.16g\n", limit);
    printf("-----------------------------\n");
    printf("limit/real = %.16g\n", limit / m2);
    printf("real/limit = %.16g\n", m2 / limit);

    ////
    double sin_phi = sin(phi);
    double cos_phi = cos(phi);
    double sqrts = sqrt(s);
    double k0 = sqrts * 0.5;
    double kx_p = k0 * sqrt(1 - y * y) * cos_phi;
    double ky_p = k0 * sqrt(1 - y * y) * sin_phi;
    double kz_p = k0 * y;
    Math::FourMomentum kip(k0, kx_p, ky_p, kz_p);
    std::cout << "nr = " << kip.ToString() << "\n";

    int j = 2;
    double len_kj_under = ps.Momenta[j].MomentumMagnitude();
    double cos_theta_j = ps.Momenta[j].PZ() / len_kj_under;
    assert(fabs(cos_theta_j) <= 1.0);
    double sin_theta_j = sqrt(1.0 - cos_theta_j * cos_theta_j);
    double phi_j = atan2(ps.Momenta[j].PY(), ps.Momenta[j].PX());
    double sin_phi_j = sin(phi_j);
    double cos_phi_j = cos(phi_j);

    // rotate k_{n+1} and k_j to the direction of the underlying born k_j
    Math::FourMomentum ki =
        rotate(kip, cos_theta_j, sin_theta_j, cos_phi_j, sin_phi_j);
    std::cout << ki.ToString() << "\n";
}
