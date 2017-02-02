#include "fks/limits.h"

#include <cassert>
#include <cstdlib>
#include <iostream>

#include "fks/splitting.h"
#include "math/fourmomentum.h"
#include "math/math.h"
#include "phasespace/realphasespace.h"
#include "physics/pdgcode.h"

using namespace Physics;

namespace {

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

Math::FourMomentum rotateFourMomentum(const Math::FourMomentum &p,
                                      const Math::FourMomentum &mother) {
    double len_mother = mother.MomentumMagnitude();
    double cos_theta_m = mother.PZ() / len_mother;
    assert(fabs(cos_theta_m) <= 1.0);
    double sin_theta_m = sqrt(1.0 - cos_theta_m * cos_theta_m);
    double phi_m = atan2(mother.PY(), mother.PX());
    double sin_phi_m = sin(phi_m);
    double cos_phi_m = cos(phi_m);
    return rotate(p, cos_theta_m, sin_theta_m, cos_phi_m, sin_phi_m);
}

double sigma_f(int i, int pdg) {
    // incoming anti fermion
    if (i < 2 && pdg < 0) {
        return -1.0;
    }
    // outgoing fermion
    if (i >= 2 && pdg > 0) {
        return -1.0;
    }
    return 1.0;
}

bool is_valid_splitting_QED(int realpdg, int bornpdg) {
    if (PDG::IsChargedFermion(realpdg)) {
        if (PDG::IsChargedFermion(bornpdg) && realpdg == bornpdg) {
            return true;
        }
        if (PDG::IsPhoton(bornpdg)) {
            return true;
        }
    }
    if (PDG::IsPhoton(realpdg) && PDG::IsChargedFermion(bornpdg)) {
        return true;
    }
    return false;
}

bool is_valid_splitting_QCD(int realpdg, int bornpdg) {
    if (PDG::IsQuark(realpdg)) {
        if (PDG::IsQuark(bornpdg) && realpdg == bornpdg) {
            return true;
        }
        if (PDG::IsGluon(bornpdg)) {
            return true;
        }
    }
    if (PDG::IsGluon(realpdg) && PDG::IsQuark(bornpdg)) {
        return true;
    }
    if (PDG::IsGluon(realpdg) && PDG::IsGluon(bornpdg)) {
        return true;
    }
    return false;
}

/**
 * @brief get_splitting_ISR returns the pdgs for the splitting functions
 *
 * A real matrix element can be approximated in the coll. limit by P*B where P
 * is the splitting function and B is the born matrix element. If the matrix
 * element has a collinear singularity, it is canceled by the FKS prefactor
 * (1-y^2). If there is no singularity, the prefactor evaluates to 0. This
 * function returns if there is a singularity for y == +1 or y == -1 and it
 * stores the pdgs for the splitting function in bornpdg and realpdg.
 *
 * Example 1:
 * bornpdgs = u u~ -> X
 * realpdgs = u u~ -> X + g
 * There is a singularity for y == +1 (collinear to u) and y == -1 (collinear to
 * u~). The splitting function to connect the born and the real process is
 * P_{qq}.
 *
 * Example 2:
 * bornpdgs = u u~ -> X
 * realpdgs = g u~ -> X u~
 * There is a singularity for y == +1 (collinear to g) but no singularity
 * for y == -1. The splitting function is P_{q<-g}.
 *
 * @param bornpdgs pdgs of the born process
 * @param realpdgs pdgs of the real process
 * @param [out] bornpdg born pdg of the splitting
 * @param [out] realpdg real pdg of the splitting
 * @param y direction of the radiated parton y == +1 or y == -1
 *
 * @return returns if there is a singularity for y.
 */
template <bool (*is_valid_splitting)(int, int)>
bool get_splitting_ISR(const int *bornpdgs, const int *realpdgs, int *bornpdg,
                       int *realpdg, int y) {
    assert(y == 1 || y == -1);
    if (realpdgs[0] == bornpdgs[0] && realpdgs[1] == bornpdgs[1]) {
        if (y == 1) {
            if (!is_valid_splitting(realpdgs[0], bornpdgs[0])) {
                // no radiation from this leg possible => no singularity
                return false;
            }
            *realpdg = realpdgs[0];
            *bornpdg = bornpdgs[0];
            return true;
        }
        if (y == -1) {
            if (!is_valid_splitting(realpdgs[1], bornpdgs[1])) {
                // no radiation from this leg possible => no singularity
                return false;
            }
            *realpdg = realpdgs[1];
            *bornpdg = bornpdgs[1];
            return true;
        }
    }

    if (realpdgs[0] == bornpdgs[0] && realpdgs[1] != bornpdgs[1] && y == -1) {
        assert(is_valid_splitting(realpdgs[1], bornpdgs[1]));
        *realpdg = realpdgs[1];
        *bornpdg = bornpdgs[1];
        return true;
    }

    if (realpdgs[0] != bornpdgs[0] && realpdgs[1] == bornpdgs[1] && y == 1) {
        assert(is_valid_splitting(realpdgs[0], bornpdgs[0]));
        *realpdg = realpdgs[0];
        *bornpdg = bornpdgs[0];
        return true;
    }
    // it is not possible that both initial state partons change their flavour.
    assert(realpdgs[0] == bornpdgs[0] || realpdgs[1] == bornpdgs[1]);
    return false;
}

bool get_splitting_ISR_QCD(const int *bornpdgs, const int *realpdgs,
                           int *bornpdg, int *realpdg, int y) {
    return get_splitting_ISR<is_valid_splitting_QCD>(bornpdgs, realpdgs,
                                                     bornpdg, realpdg, y);
}

bool get_splitting_ISR_QED(const int *bornpdgs, const int *realpdgs,
                           int *bornpdg, int *realpdg, int y) {
    return get_splitting_ISR<is_valid_splitting_QED>(bornpdgs, realpdgs,
                                                     bornpdg, realpdg, y);
}

double limit_with_spincorr_qcd_isr(double sb,
                                   const UserProcess::SpinCorrelated &spincorr,
                                   double xi, double phi, int realpdg) {
    // double y = 0.9999999999;
    // Phasespace::Phasespace ps_real;
    // Phasespace::GenRealPhasespace(&ps_real, &ps, 0, xi, y, phi);

    Math::FourMomentum kT(0.0, cos(phi), sin(phi), 0.0);
    kT.Scale(1.0 / kT.MomentumMagnitude());

    double z = 1.0 - xi;

    double splitting1 = 0.0;
    double charge_factor = 1.0;
    if (realpdg == 0 || realpdg == 21) {
        splitting1 = 2.0 * spincorr.Born() * (z / (1 - z) + z * (1 - z));
        charge_factor = 3.0; // CA
    } else {
        splitting1 = z * spincorr.Born();
        charge_factor = 4.0 / 3.0; // CF
    }

    double splitting2 = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            splitting2 += 4 * (1 - z) / z * kT.At(mu) * kT.At(nu) *
                          std::real(spincorr[mu][nu]);
        }
    }

    // k2 = sb/2 *xi*(1-y);
    // FKS xi^2(1-y^2) cancel with k2 -> 2 * xi remaining
    return (splitting1 + splitting2) * charge_factor / sb * 32 * M_PI * xi;
}

double
limit_with_spincorr_qcd_isr_soft(double sb,
                                 const UserProcess::SpinCorrelated &spincorr,
                                 double phi, int realpdg) {
    double splitting1 = 0.0;
    double charge_factor = 1.0;
    if (realpdg == 0 || realpdg == 21) {
        splitting1 = 2.0 * spincorr.Born();
        charge_factor = 3.0; // CA
    }

    return splitting1 * charge_factor / sb * 32 * M_PI;
}

double limit_with_spincorr_qcd_fsr(const Phasespace::Phasespace &ps,
                                   const UserProcess::SpinCorrelated &spincorr,
                                   double xi, double phi, int realpdg,
                                   int splitting_j) {
    Phasespace::Phasespace ps_real_c;
    double y = 0.9999999;
    Phasespace::GenRealPhasespace(&ps_real_c, &ps, ps.N + 2, splitting_j, xi, y,
                                  phi);

    const auto &pj = ps_real_c.Momenta[splitting_j];
    const auto &pi = ps_real_c.Momenta[ps_real_c.N + 1]; // radiated particle is
                                                         // always last particle
    double z = ps_real_c.Momenta[splitting_j].E() / ps.Momenta[splitting_j].E();

    double splitting1 = 0.0;
    double splitting2_sign = 0.0;
    double charge_factor = 1.0;
    if (realpdg == 0 || realpdg == 21) {
        splitting1 = 2 * ((1 - z) / z + z / (1 - z));
        splitting2_sign = +1.0;
        charge_factor = 3.0 / 2.0; // CA/2 (identical particles)
    } else {
        splitting2_sign = -1.0;
        splitting1 = 1.0;
        charge_factor = 0.5; // TF
    }
    double splitting2 = 0.0;
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            splitting2 += pj.At(mu) * pj.At(nu) * std::real(spincorr[mu][nu]);
        }
    }

    splitting2 *= 2 / pi.Dot(pj) * splitting2_sign;
    double result = splitting1 * spincorr.Born() + splitting2;
    // k2 is actually proportional to (1-y) but it cancels with the FKS
    // prefactor.
    double k2 = 2 * ps.Momenta[4].E() * ps.Momenta[4].E() * z * (1 - z);

    result *= 8.0 * M_PI / k2;
    // returns xi^2*(1-y)*ME -> (1-y) canceled with k2
    result *= xi * xi;
    return charge_factor * result;
}

double
limit_with_spincorr_qcd_fsr_soft(const Phasespace::Phasespace &ps,
                                 const UserProcess::SpinCorrelated &spincorr,
                                 double phi, int realpdg, int splitting_j) {

    double s = ps.X1 * ps.X2 * ps.S;
    double splitting1 = 0.0;
    double charge_factor = 1.0;
    if (realpdg == 0 || realpdg == 21) {
        splitting1 = 4.0 / s;
        charge_factor = 3.0 / 2.0; // CA/2 (identical particles)
    }

    double result = splitting1 * spincorr.Born();

    return 8.0 * M_PI * charge_factor * result;
}

} // end namespace

namespace FKS {

namespace QED {

double CollinearLimitFSR(int realpdg, const int bornpdg,
                         const Phasespace::Phasespace &ps, int splitting_j,
                         double xi, double phi, double alpha, double bornme,
                         const UserProcess::SpinCorrelated &spincorr) {
    assert(splitting_j >= 2);
    assert(splitting_j < ps.N + 2);
    if (bornpdg != 22) {
        // normal splitting functions

        // z is defined as p_j.E()/p.E() where p is the born momentum and p_j
        // the n+1 momentum. xi is the FKS xi and  defined as follows: p_j.E() =
        // xi*sqrt(s)/2. For more than two particles there is an upper bound for
        // xi.
        // Therefore, z != 1 - xi.
        double z =
            1 -
            xi * sqrt(ps.X1 * ps.X2 * ps.S) / 2.0 / ps.Momenta[splitting_j].E();
        // the  xi here is not the FKS xi from the function parameters! It's the
        // splitting function's xi = 1-z.
        double splittingXxi =
            FKS::QED::splittingTimesXi(realpdg, bornpdg, 1 - z);
        // k2 is actually proportional to (1-y) but it cancels with the FKS
        // prefactor.
        double k2 = 2 * ps.Momenta[splitting_j].E() *
                    ps.Momenta[splitting_j].E() * z * (1 - z);
        return xi * xi * 8 * M_PI * alpha / (1 - z) / k2 * splittingXxi *
               bornme;
    }
    assert(0 && "not implemented");
    return 0.0;
}

double CollinearLimitISR(int Nborn, const int *realpdgs, const int *bornpdgs,
                         double xi, int y, double phi, double sb, double alpha,
                         double bornme,
                         const UserProcess::SpinCorrelated &spincorr) {
    assert(y == -1 || y == 1);
    int bornpdg = 0xffffff;
    int realpdg = 0xffffff;

    bool has_singularity =
        get_splitting_ISR_QED(bornpdgs, realpdgs, &bornpdg, &realpdg, y);
    if (!has_singularity) {
        // If the real matrix element has no singularity in y == +1 or y == -1,
        // the prefactor (1-y)(1+y) is not canceled. Therefore, the prefactor
        // evaluates to 0.0 and the overall result is 0.0.
        return 0.0;
    };

    if (bornpdg == 22) {
        // spin correlated
        assert(0 && "not implemented");
    }

    // normal splitting functions
    double sp = FKS::QED::splittingTimesXi(realpdg, bornpdg, xi);
    return 32.0 * Math::Pi * alpha / sb * sp * bornme;
}

/*
* SoftLimit computes the soft limit for photon radiation.
*
* @note this function is only needed if the emitted particle is a photon. For
* charged particles the limit is 0.
*
* @param N number of momenta in the phase space
* @param born_momenta phase space of the born configuration
* @param pdg pdg numbers for all momenta
* @param s partonic center of mass energy for the born phase space
* @param jmother jmother is the index of the final state particle which is used
*                to parametrize the n+1 phase space in terms of the FKS
*                variables. (Should correspondend to j for the generation of the
*                n+1 final state FKS phase space.)
* @param alpha electro magnetic alpha
* @param born_me born matrix element for born_momenta
* @param y FKS y
* @param phi FKS phi
*/
double SoftLimit(int N, const Math::FourMomentum *born_momenta, const int *pdg,
                 int realpdg, double s, int jmother, double alpha,
                 const Util::Matrix2 &Born, double y, double phi) {
    assert(jmother >= 0);
    assert(jmother < N);
    assert(Born.GetLen() == 1);
    if (!PDG::IsPhoton(realpdg)) {
        return 0.0;
    }
    double born_me = Born.Get(0, 0);
    double limit = 0.0;
    double sqrts = sqrt(s);
    double p0 = sqrts * 0.5; // energy of radiated photon without xi

    double sinf = sin(phi);
    double cosf = cos(phi);
    double sint = sqrt(1.0 - y * y);
    // photon momentum in frame where mother particle points in z direction.
    // note that the mother particle has the same direction as the k_j of the
    // final state since a soft photon doesn't affect the mother particle.
    Math::FourMomentum p(p0, p0 * sint * cosf, p0 * sint * sinf, p0 * y);

    double y_pre = 0.0;
    if (jmother >= 2) {
        // rotate the photon momentum to the frame of the mother particle.
        p = rotateFourMomentum(p, born_momenta[jmother]);
        y_pre = 1.0 - y;
    } else { // initial state radiation
        y_pre = 1.0 - y * y;
    }

    for (int i = 0; i < N; i++) {
        double chargei = Physics::PDG::Charge(abs(pdg[i]));
        double sigma_i = sigma_f(i, pdg[i]);
        if (fabs(chargei) < 1e-6) {
            continue;
        }
        assert(fabs(born_momenta[i].Dot(born_momenta[i])) < 1e-4 &&
               "massive particles are not implemented");

        for (int j = i + 1; j < N; j++) {
            double chargej = Physics::PDG::Charge(abs(pdg[j]));
            double sigma_j = sigma_f(j, pdg[j]);
            if (fabs(chargej) < 1e-6) {
                continue;
            }
            Math::FourMomentum ki = born_momenta[i];
            Math::FourMomentum kj = born_momenta[j];
            double eikonal = ki.Dot(kj) / (ki.Dot(p) * kj.Dot(p));
            limit -= eikonal * chargej * chargei * sigma_i * sigma_j;
        }
    }
    // the factor 2 corrects for the j > i for loop instead of j = 0..N-1
    return 2.0 * y_pre * 4.0 * Math::Pi * alpha * born_me * limit;
}

double SoftCollinearLimitFSR(int realpdg, const int bornpdg,
                             const Phasespace::Phasespace &ps, int splitting_j,
                             double phi, double alpha, double bornme,
                             const UserProcess::SpinCorrelated &spincorr) {
    assert(splitting_j >= 2);
    if (bornpdg != 22) {
        // normal splitting functions
        double splittingXxi = FKS::QED::splittingTimesXiSoft(realpdg, bornpdg);
        double s = ps.X1 * ps.X2 * ps.S;
        return 16 * M_PI * alpha * splittingXxi * bornme / s;
    }
    assert(0 && "not implemented");
    return 0.0;
}

double SoftCollinearLimitISR(int Nborn, const int *realpdgs,
                             const int *bornpdgs, int y, double phi, double sb,
                             double alpha, double bornme,
                             const UserProcess::SpinCorrelated &spincorr) {
    assert(y == -1 || y == 1);
    int bornpdg = 0xffffff;
    int realpdg = 0xffffff;

    bool has_singularity =
        get_splitting_ISR_QED(bornpdgs, realpdgs, &bornpdg, &realpdg, y);
    if (!has_singularity) {
        // If the real matrix element has no singularity in y == +1 or y == -1,
        // the prefactor (1-y)(1+y) is not canceled. Therefore, the prefactor
        // evaluates to 0.0 and the overall result is 0.0.
        return 0.0;
    };

    if (bornpdg == 22) {
        // spin correlated
        assert(0 && "not implemented");
    }

    // normal splitting functions
    double sp = FKS::QED::splittingTimesXiSoft(realpdg, bornpdg);
    return 32.0 * Math::Pi * alpha / sb * sp * bornme;
}

} // end namespace QED

namespace QCD {

double CollinearLimitISR(int Nborn, const int *realpdgs, const int *bornpdgs,
                         double xi, int y, double phi, double sb, double alpha,
                         double bornme,
                         const UserProcess::SpinCorrelated &spincorr) {
    assert(y == -1 || y == 1);
    int bornpdg = 0xffffff;
    int realpdg = 0xffffff;

    bool has_singularity =
        get_splitting_ISR_QCD(bornpdgs, realpdgs, &bornpdg, &realpdg, y);

    if (!has_singularity) {
        // If the real matrix element has no singularity in y == +1 or y == -1,
        // the prefactor (1-y)(1+y) is not canceled. Therefore, the prefactor
        // evaluates to 0.0 and the overall result is 0.0.
        return 0.0;
    };

    if (bornpdg == 0 || bornpdg == 21) {
        int radpdg = 21;
        if (!(realpdg == 21 || realpdg == 0)) {
            radpdg = realpdg;
        }
        int n = 1;
        for (int i = 2; i < Nborn; i++) {
            if (bornpdgs[i] == radpdg) {
                n += 1;
            }
            if (radpdg == 21 && bornpdgs[i] == 0) {
                n += 1;
            }
        }
        double sym = 1.0 / ((double)n);
        // spin correlated
        return alpha * sym *
               limit_with_spincorr_qcd_isr(sb, spincorr, xi, phi, realpdg);
    }

    int radpdg = 21; // gluon
    if (realpdg == 0 || realpdg == 21) {
        radpdg = -bornpdg;
    }
    int n = 1;
    for (int i = 2; i < Nborn; i++) {
        if (bornpdgs[i] == radpdg) {
            n += 1;
        }
        if (radpdg == 21 && bornpdgs[i] == 0) {
            n += 1;
        }
    }
    double sym = 1.0 / ((double)n);

    // normal splitting functions
    double sp = FKS::QCD::splittingTimesXi(realpdg, bornpdg, xi);
    return 32.0 * Math::Pi * alpha / sb * sp * bornme * sym;
}

double SoftCollinearLimitISR(int Nborn, const int *realpdgs,
                             const int *bornpdgs, int y, double phi, double sb,
                             double alpha, double bornme,
                             const UserProcess::SpinCorrelated &spincorr) {
    assert(y == -1 || y == 1);
    int bornpdg = 0xffffff;
    int realpdg = 0xffffff;

    bool has_singularity =
        get_splitting_ISR_QCD(bornpdgs, realpdgs, &bornpdg, &realpdg, y);
    if (!has_singularity) {
        // If the real matrix element has no singularity in y == +1 or y == -1,
        // the prefactor (1-y)(1+y) is not canceled. Therefore, the prefactor
        // evaluates to 0.0 and the overall result is 0.0.
        return 0.0;
    };

    if (bornpdg == 0 || bornpdg == 21) {
        // spin correlated
        return alpha *
               limit_with_spincorr_qcd_isr_soft(sb, spincorr, phi, realpdg);
    }

    int radpdg = 21; // gluon
    if (realpdg == 0 || realpdg == 21) {
        radpdg = -bornpdg;
    }
    int n = 1;
    for (int i = 2; i < Nborn; i++) {
        if (bornpdgs[i] == radpdg) {
            n += 1;
        }
        if (radpdg == 21 && bornpdgs[i] == 0) {
            n += 1;
        }
    }
    double sym = 1.0 / ((double)n);

    // normal splitting functions
    double sp = FKS::QCD::splittingTimesXiSoft(realpdg, bornpdg);
    return 32.0 * Math::Pi * alpha / sb * sp * bornme * sym;
}

double CollinearLimitFSR(int realpdg, const int bornpdg,
                         const Phasespace::Phasespace &ps, int splitting_j,
                         double xi, double phi, double alpha, double bornme,
                         const UserProcess::SpinCorrelated &spincorr) {
    assert(splitting_j >= 2);
    assert(splitting_j < ps.N + 2);
    if (bornpdg != 0 && bornpdg != 21) {
        // normal splitting functions

        // z is defined as p_j.E()/p.E() where p is the born momentum and p_j
        // the n+1 momentum. xi is the FKS xi and  defined as follows: p_j.E() =
        // xi*sqrt(s)/2. For more than two particles there is an upper bound for
        // xi.
        // Therefore, z != 1 - xi.
        double z =
            1 -
            xi * sqrt(ps.X1 * ps.X2 * ps.S) / 2.0 / ps.Momenta[splitting_j].E();
        // the  xi here is not the FKS xi from the function parameters! It's the
        // splitting function's xi = 1-z.
        double splittingXxi =
            FKS::QCD::splittingTimesXi(realpdg, bornpdg, 1 - z);
        // k2 is actually proportional to (1-y) but it cancels with the FKS
        // prefactor.
        double k2 = 2 * ps.Momenta[splitting_j].E() *
                    ps.Momenta[splitting_j].E() * z * (1 - z);
        return xi * xi * 8 * M_PI * alpha / (1 - z) / k2 * splittingXxi *
               bornme;
    }
    return alpha * limit_with_spincorr_qcd_fsr(ps, spincorr, xi, phi, realpdg,
                                               splitting_j);
}

double SoftCollinearLimitFSR(int realpdg, const int bornpdg,
                             const Phasespace::Phasespace &ps, int splitting_j,
                             double phi, double alpha, double bornme,
                             const UserProcess::SpinCorrelated &spincorr) {
    assert(splitting_j >= 2);
    assert(splitting_j < ps.N + 2);
    if (bornpdg != 0 && bornpdg != 21) {
        // normal splitting functions
        double splittingXxi = FKS::QCD::splittingTimesXiSoft(realpdg, bornpdg);
        double s = ps.X1 * ps.X2 * ps.S;
        return 16 * M_PI * alpha * splittingXxi * bornme / s;
    }
    return alpha * limit_with_spincorr_qcd_fsr_soft(ps, spincorr, phi, realpdg,
                                                    splitting_j);
}

double SoftLimit(int N, const Math::FourMomentum *born_momenta, const int *pdg,
                 int realpdg, double s, int jmother, double alpha_s,
                 const Util::Matrix2 &ColorCorrelatedBorn, double y,
                 double phi) {
    assert(jmother >= 0);
    assert(jmother < N);
    assert(ColorCorrelatedBorn.GetLen() == N);
    if (!PDG::IsGluon(realpdg)) {
        // only gluon has soft limit different from 0.0. The divergene of an
        // emitted quark cancels.
        return 0.0;
    }
    double limit = 0.0;

    double sinf = sin(phi);
    double cosf = cos(phi);
    double sint = sqrt(1.0 - y * y);
    // gluon momentum in frame where mother particle points in z direction.
    // note that the mother particle has the same direction as the k_j of the
    // final state since a soft gluon doesn't affect the mother particle.
    // The energy is set to 1.0. We add it at the end.
    Math::FourMomentum p(1.0, sint * cosf, sint * sinf, y);

    double y_pre = 0.0;
    if (jmother >= 2) {
        // rotate the photon momentum to the frame of the mother particle.
        p = rotateFourMomentum(p, born_momenta[jmother]);
        y_pre = 1.0 - y;
    } else { // initial state radiation
        y_pre = 1.0 - y * y;
    }

    for (int i = 0; i < N; i++) {
        assert(fabs(born_momenta[i].Dot(born_momenta[i])) < 1e-4 &&
               "massive particles are not implemented");

        int pdgi = std::abs(pdg[i]);
        if (pdgi > 6 && pdgi != 21) {
            continue;
        }
        for (int j = i + 1; j < N; j++) {
            int pdgj = std::abs(pdg[j]);
            if (pdgj > 6 && pdgj != 21) {
                continue;
            }
            Math::FourMomentum ki = born_momenta[i];
            Math::FourMomentum kj = born_momenta[j];
            double eikonal = ki.Dot(kj) / (ki.Dot(p) * kj.Dot(p));
            limit += eikonal * ColorCorrelatedBorn.Get(i, j);
            limit += eikonal * ColorCorrelatedBorn.Get(j, i);
        }
    }

    int ngluon = 1;
    for (int i = 2; i < N; i++) {
        if (PDG::IsGluon(pdg[i])) {
            ngluon += 1;
        }
    }
    double sym = 1.0 / ((double)ngluon);
    // 4/s is the factor gluon energy times xi^2 in eikonal.
    return y_pre * 4.0 * Math::Pi * alpha_s * limit / s * 4.0 * sym;
}

} // end namespace QCD

} // end namespace FKS
