// Implementation of template functions.
// This file has to be included in the header!

#include <cmath>

#include "math/fourmomentum.h"
#include "math/math.h"
#include "phasespace/phasespace.h"
#include "realphasespace.h"

#include <iomanip>
#include <iostream>

namespace {
double gen_x1(double x1bar, double xi, double y) {
    return x1bar *
           sqrt((2.0 - xi * (1.0 - y)) / ((2.0 - xi * (1 + y)) * (1.0 - xi)));
}

double gen_x2(double x2bar, double xi, double y) {
    return x2bar *
           sqrt((2.0 - xi * (1.0 + y)) / ((2.0 - xi * (1 - y)) * (1.0 - xi)));
}

Math::FourMomentum rotate(const Math::FourMomentum &p, double costh,
                          double sinth, double cosphi, double sinphi) {
    double px = p.PX();
    double py = p.PY();
    double pz = p.PZ();
    double px_p = cosphi * costh * px - sinphi * py + sinth * cosphi * pz;
    double py_p = sinphi * costh * px + cosphi * py + sinphi * sinth * pz;
    double pz_p = -sinth * px + costh * pz;
    return Math::FourMomentum(p.E(), px_p, py_p, pz_p);
}

} // namespace

namespace Phasespace {

/**
@brief Create FKS phase space for FSR

 */
void GenRealPhasespaceFSR(Phasespace *RESTRICT ps_real,
                          Phasespace const *const RESTRICT ps_born, int j,
                          double xi, double y, double phi) {
    assert(j >= 2 && "particle j must be a final state particle!");
    assert(j < ps_born->N + 2 && "particle index j too large");
    assert(ps_born->Masses[j - 2] == 0.0 &&
           "the particle which radiates a particle is assumed to be massless!");
    assert(xi >= 0.0 && xi <= 1.0);
    assert(y >= -1.0 && y <= 1.0);
    assert(phi >= 0.0 && phi < 2.0 * Math::Pi);
    double s = ps_born->X1 * ps_born->X2 * ps_born->S;
    double sqrts = sqrt(s);
    /*if (xi > 2.0 * ps_born->Momenta[j].MomentumMagnitude() / sqrts) {
        printf("xi =    %25.15g\n", xi);
        printf("ximax = %25.15g\n",
               2.0 * ps_born->Momenta[j].MomentumMagnitude() / sqrts);
    }*/
    // Math::FourMomentum rec =
    //     Math::FourMomentum(sqrts, 0, 0, 0).Minus(ps_born->Momenta[j]);
    // double Mrec2 = rec.Dot(rec);
    double Mrec2 = s - 2.0 * ps_born->Momenta[j].E() * sqrts;
    double k0 = sqrts * xi * 0.5;
    double len_kj =
        (s - Mrec2 - 2 * sqrts * k0) / (2 * (sqrts - k0 * (1.0 - y)));

    double len_mother = sqrt(len_kj * len_kj + k0 * k0 + 2.0 * k0 * len_kj * y);

    // angle between mother and radiated parton
    double cos_psi = (k0 + len_kj * y) / len_mother;

    if (fabs(cos_psi) > 1.0 && fabs(cos_psi) - 1.0 < 1e-6) {
        // correct numerical errors
        cos_psi = copysign(1.0, cos_psi);
    } else {
        assert(fabs(cos_psi) <= 1.0);
    }
    double sin_psi = sqrt(1.0 - cos_psi * cos_psi);
    double sin_phi = sin(phi);
    double cos_phi = cos(phi);
    // construct k_{n+1} in a frame where k points in z direction
    double kx_p = k0 * sin_psi * cos_phi;
    double ky_p = k0 * sin_psi * sin_phi;
    double kz_p = k0 * cos_psi;
    Math::FourMomentum kip(k0, kx_p, ky_p, kz_p);
    // construct vec{k}_j = vec{k} - vec{k}_{n+1} and k_j = len_kj
    Math::FourMomentum kjp(len_kj, -kx_p, -ky_p, len_mother - kz_p);

    // k points by construction in the same direction as k_j in the underlying
    // born. Get direction of k_j
    double len_kj_under = ps_born->Momenta[j].MomentumMagnitude();
    double cos_theta_j = ps_born->Momenta[j].PZ() / len_kj_under;
    assert(fabs(cos_theta_j) <= 1.0);
    double sin_theta_j = sqrt(1.0 - cos_theta_j * cos_theta_j);
    double phi_j = atan2(ps_born->Momenta[j].PY(), ps_born->Momenta[j].PX());
    double sin_phi_j = sin(phi_j);
    double cos_phi_j = cos(phi_j);

    // rotate k_{n+1} and k_j to the direction of the underlying born k_j
    Math::FourMomentum ki =
        rotate(kip, cos_theta_j, sin_theta_j, cos_phi_j, sin_phi_j);
    Math::FourMomentum kj =
        rotate(kjp, cos_theta_j, sin_theta_j, cos_phi_j, sin_phi_j);

    // calculate k_rec = q - k
    Math::FourMomentum k = ki.Plus(kj);
    Math::FourMomentum krec(sqrts - k.E(), -k.PX(), -k.PY(), -k.PZ());

    // boost velocity beta
    double t = krec.E() + krec.MomentumMagnitude();
    double beta = (s - t * t) / (s + t * t);

    double inv_abs_krec = 1.0 / krec.MomentumMagnitude();
    double t2 = beta * inv_abs_krec;
    double betax = krec.PX() * t2;
    double betay = krec.PY() * t2;
    double betaz = krec.PZ() * t2;
    double betax_s = krec.PX() * inv_abs_krec;
    double betay_s = krec.PY() * inv_abs_krec;
    double betaz_s = krec.PZ() * inv_abs_krec;
    double gamma = 1.0 / sqrt(1.0 - beta * beta);

    double Lambda[4][4];
    Lambda[0][0] = gamma;
    Lambda[0][1] = -betax * gamma;
    Lambda[0][2] = -betay * gamma;
    Lambda[0][3] = -betaz * gamma;

    Lambda[1][0] = -gamma * betax;
    Lambda[1][1] = (1.0 + (gamma - 1.0) * betax_s * betax_s);
    Lambda[1][2] = (gamma - 1.0) * betax_s * betay_s;
    Lambda[1][3] = (gamma - 1.0) * betax_s * betaz_s;

    Lambda[2][0] = -gamma * betay;
    Lambda[2][1] = (gamma - 1.0) * betay_s * betax_s;
    Lambda[2][2] = (1.0 + (gamma - 1.0) * betay_s * betay_s);
    Lambda[2][3] = (gamma - 1.0) * betay_s * betaz_s;

    Lambda[3][0] = -gamma * betaz;
    Lambda[3][1] = (gamma - 1.0) * betaz_s * betax_s;
    Lambda[3][2] = (gamma - 1.0) * betaz_s * betay_s;
    Lambda[3][3] = (1.0 + (gamma - 1.0) * betaz_s * betaz_s);

    int N = ps_born->N;
    for (int i = 2; i < N + 2; i++) {
        if (i == j) {
            continue;
        }
        double pold[] = {
            ps_born->Momenta[i].E(), ps_born->Momenta[i].PX(),
            ps_born->Momenta[i].PY(), ps_born->Momenta[i].PZ(),
        };
        double pnew[4] = {0.0};
        for (int k = 0; k < 4; k++) {
            for (int l = 0; l < 4; l++) {
                pnew[k] += Lambda[k][l] * pold[l];
            }
        }
        ps_real->Momenta[i].Set(pnew[0], pnew[1], pnew[2], pnew[3]);
    }
    ps_real->Momenta[N + 2].Set(ki.E(), ki.PX(), ki.PY(), ki.PZ());
    ps_real->Momenta[j].Set(kj.E(), kj.PX(), kj.PY(), kj.PZ());
    ps_real->Momenta[0] = ps_born->Momenta[0];
    ps_real->Momenta[1] = ps_born->Momenta[1];

    // set energy of n+1 phase space
    ps_real->S = ps_born->S;
    ps_real->X1 = ps_born->X1;
    ps_real->X2 = ps_born->X2;

    // set final state masses
    for (int i = 0; i < N; i++) {
        ps_real->Masses[i] = ps_born->Masses[i];
    }
    ps_real->Masses[N] = 0.0;

    double real_over_born = ps_real->Momenta[j].MomentumMagnitude() /
                            ps_born->Momenta[j].MomentumMagnitude();
    ps_real->Jacobian = ps_born->Jacobian * 2.0 * s * xi /
                        (64.0 * Math::Pi_Cube) / (2.0 - xi * (1 - y)) *
                        real_over_born;
    ps_real->N = ps_born->N + 1;
}

/**
@brief Create FKS phase space for ISR

GenRealPhasespaceISR creates a FKS phase space from a born phase space ps_born
and the parametrization of the FKS parton.
The born phase space ps_born has to be in its partonic CM frame.

The FKS phase space is in its CM frame is written to ps_real.

@param ps_born born phase space (in partonic CMS)
@param xi energy fraction of the FKS parton (xi in (0,1))
@param y cos theta of FKS parton in n+1 CMS (y in [-1,1])
@param phi azimuthal angle of FKS parton in n+1 CMS (phi in [0,2Pi))
@param[out] ps_real The constructed phase space is written to ps_real.
*/
void GenRealPhasespaceISR(Phasespace *RESTRICT ps_real,
                          Phasespace const *const RESTRICT ps_born, double xi,
                          double y, double phi, bool boost_to_cms) {
    // notation:
    //    n+1 momenta: k
    //
    int N = ps_born->N;
    const double x1bar = ps_born->X1;
    const double x2bar = ps_born->X2;
    const double S = ps_born->S;
    const double sqrts_bar = sqrt(x1bar * x2bar * S);
    const double x1 = gen_x1(x1bar, xi, y);
    const double x2 = gen_x2(x2bar, xi, y);
    const double sqrts = sqrt(x1 * x2 * S);
    const double t0 = 0.5 * sqrts_bar * x1 / x1bar;
    const double t1 = 0.5 * sqrts_bar * x2 / x2bar;

    // create n+1 initial state momenta from born momenta
    std::array<Math::FourMomentum, MAXMOM> k_final;
    k_final[0].Set(t0, 0.0, 0.0, t0);
    k_final[1].Set(t1, 0.0, 0.0, -t1);

    const double k_cm_0 = 0.5 * sqrts * xi;
    const double sin_t = sqrt(1.0 - y * y);
    Math::FourMomentum k_cm(k_cm_0, k_cm_0 * sin_t * cos(phi),
                            k_cm_0 * sin_t * sin(phi), k_cm_0 * y);
    // boost k_cm to frame of k1 and k2
    const double t2 = 1.0 / (2.0 * sqrt(x1 * x2 * x1bar * x2bar));
    const double ch_eta = t2 * (x1bar * x2 + x1 * x2bar);
    const double sh_eta_m = -t2 * (x1bar * x2 - x1 * x2bar);

    Math::FourMomentum k(k_cm.E() * ch_eta + sh_eta_m * k_cm.PZ(), k_cm.PX(),
                         k_cm.PY(), k_cm.E() * sh_eta_m + k_cm.PZ() * ch_eta);

#ifdef SEPARATEBOOSTS
    // Here we construct a separte boost for every component of k_tot.
    // boost B1 to remove x component of k_tot = k1 + k2 - k
    double sinh_n1 =
        -xi * cos(phi) *
        sqrt((y * y - 1) / (xi * xi * (y * y - 1.0) * sin(phi) * sin(phi) +
                            4.0 * (xi - 1.0)));
    double cosh_n1 = sqrt(1.0 + sinh_n1 * sinh_n1);

    // boost to remoe y component of B1*k_tot
    double sinh_n2 = 0.5 * xi * sin(phi) * sqrt((-y * y + 1.0) / (1.0 - xi));
    double cosh_n2 = sqrt(1.0 + sinh_n2 * sinh_n2);
#else
    // uses the boost formular for an arbitrary direction.

    // construct a boost for the final state momenta
    // beta_i = ktot_i/ktot_0
    double beta_abs = -xi * sqrt((1.0 - y * y) /
                                 (xi * xi * (1.0 - y * y) + 4.0 * (1.0 - xi)));

    double betax = cos(phi) * beta_abs;
    double betay = sin(phi) * beta_abs;
    // betaz = 0
    double beta_sqr = betax * betax + betay * betay;
    double gamma = 1.0 / sqrt(1.0 - beta_sqr);
    double gamma2 = beta_sqr > 0 ? (gamma - 1.0) / beta_sqr : 0.0;
#endif

    // boost final state momenta rest frame of k_tot & set masses
    std::array<double, MAXMOM - 2> masses;
    for (int i = 2; i < N + 2; i++) {
        double E_old = ps_born->Momenta[i].E();
        double px_old = ps_born->Momenta[i].PX();
        double py_old = ps_born->Momenta[i].PY();
#ifdef SEPARATEBOOSTS
        double E = E_old * cosh_n1 * cosh_n2 + px_old * sinh_n1 -
                   py_old * cosh_n1 * sinh_n2;
        double px = px_old * cosh_n1 + E_old * cosh_n2 * sinh_n1 -
                    py_old * sinh_n1 * sinh_n2;
        double py = py_old * cosh_n2 - E_old * sinh_n2;
        double pz = ps_born->Momenta[i].PZ();
#else
        double bp = betax * px_old + betay * py_old;
        double px = px_old + gamma2 * bp * betax + gamma * betax * E_old;
        double py = py_old + gamma2 * bp * betay + gamma * betay * E_old;
        double pz = ps_born->Momenta[i].PZ();
        double E = gamma * (E_old + bp);
#endif
        k_final[i].Set(E, px, py, pz);

        masses[i - 2] = ps_born->Masses[i - 2];
    }
    // set mass & momenta of fks parton
    masses[N] = 0.0;
    k_final[N + 2] = k;

    // boost n+1 to cms
    if (boost_to_cms) {
        double ch = t2 * (x1bar * x2 + x1 * x2bar);
        double sh = t2 * (x1bar * x2 - x1 * x2bar);
        for (int i = 0; i < N + 3; i++) {
            Math::FourMomentum k = k_final[i];
            double E = k.E() * ch + sh * k.PZ();
            double px = k.PX();
            double py = k.PY();
            double pz = k.E() * sh + k.PZ() * ch;
            k_final[i].Set(E, px, py, pz);
        }
    }

    // write n+1 phasespace to ps_real
    ps_real->S = ps_born->S;
    ps_real->X1 = x1;
    ps_real->X2 = x2;
    ps_real->Masses = masses;
    ps_real->Momenta = k_final;
    ps_real->Jacobian = ps_born->Jacobian * x1 * x2 * ps_real->S * xi /
                        (64.0 * Math::Pi_Cube * (1.0 - xi));
    ps_real->N = ps_born->N + 1;
}

void GenRealPhasespace(Phasespace *RESTRICT ps_real,
                       Phasespace const *const RESTRICT ps_born, int i, int j,
                       double xi, double y, double phi, bool boost_to_cms) {
    assert(i >= 2);
    int mother = (i < j) ? i : j;
    if (mother >= 2) {
        GenRealPhasespaceFSR(ps_real, ps_born, mother, xi, y, phi);
    } else {
        GenRealPhasespaceISR(ps_real, ps_born, xi, y, phi, boost_to_cms);
    }

    // move last particle to position i
    int N = ps_real->N;
    auto p_rad = ps_real->Momenta[N + 1];
    for (int n = N + 1; n > i; n--) {
        ps_real->Momenta[n] = ps_real->Momenta[n - 1];
    }
    ps_real->Momenta[i] = p_rad;
    if (i < j && mother + 1 != j) {
        // move mother particle to j (note: mother is now at position mother+1)
        auto pm = ps_real->Momenta[mother + 1];
        if (mother + 1 > j) {
            for (int n = mother + 1; n > j; n--) {
                ps_real->Momenta[n] = ps_real->Momenta[n - 1];
            }
        }
        if (mother + 1 < j) {
            for (int n = mother + 1; n < j; n++) {
                ps_real->Momenta[n] = ps_real->Momenta[n + 1];
            }
        }
        ps_real->Momenta[j] = pm;
    }
}

} // end namespace Phasespace
