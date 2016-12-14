#include "fks/limits.h"

#include <cstdlib>
#include <cassert>
#include <iostream>

#include "math/math.h"
#include "math/fourmomentum.h"

#include "fks/splitting.h"
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
bool get_splitting_ISR(const int *bornpdgs, const int *realpdgs,
                              int *bornpdg, int *realpdg, int y) {
    assert(y == 1 || y == -1);
    if (realpdgs[0] == bornpdgs[0] && realpdgs[1] == bornpdgs[1]) {
        if (y == 1) {
            if(!is_valid_splitting(realpdgs[0], bornpdgs[0])) {
                // no radiation from this leg possible => no singularity
                return false;
            }
            *realpdg = realpdgs[0];
            *bornpdg = bornpdgs[0];
            return true;
        }
        if (y == -1) {
            if(!is_valid_splitting(realpdgs[1], bornpdgs[1])) {
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
    return get_splitting_ISR<is_valid_splitting_QCD>(bornpdgs, realpdgs, bornpdg, realpdg, y);
}

bool get_splitting_ISR_QED(const int *bornpdgs, const int *realpdgs,
                              int *bornpdg, int *realpdg, int y) {
    return get_splitting_ISR<is_valid_splitting_QED>(bornpdgs, realpdgs, bornpdg, realpdg, y);
}

} // end namespace

namespace FKS {

namespace QED {

double CollinearLimitFSR(int realpdg, int bornpdg, double xi, double s,
                         double alpha, double born_me) {
    double sp = FKS::QED::splittingTimesXi(realpdg, bornpdg, xi);
    return 16.0 * Math::Pi * alpha / s / (1.0 - xi) * sp * born_me;
}

double CollinearLimitISR(const int *realpdgs, const int *bornpdgs, double xi,
                         int y, double s, double alpha, double born_me) {
    int realpdg = 0xffffff;
    int bornpdg = 0xffffff;

    bool has_singularity = get_splitting_ISR_QCD(bornpdgs, realpdgs, &bornpdg, &realpdg, y);

    if (!has_singularity) {
        // If the real matrix element has no singularity in y == +1 or y == -1,
        // the prefactor (1-y)(1+y) is not canceled. Therefore, the prefactor
        // evaluates to 0.0 and the overall result is 0.0.
        return 0.0; 
    };

    double sp = FKS::QED::splittingTimesXi(realpdg, bornpdg, xi);
    return 32.0 * Math::Pi * alpha / s * sp * born_me;
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
    double sint = sqrt(1.0 - y*y);
    // photon momentum in frame where mother particle points in z direction.
    // note that the mother particle has the same direction as the k_j of the
    // final state since a soft photon doesn't affect the mother particle.
    Math::FourMomentum p(p0, p0*sint*cosf, p0*sint*sinf, p0*y);

    double y_pre = 0.0;
    if (jmother >= 2) {
        // rotate the photon momentum to the frame of the mother particle.
        p = rotateFourMomentum(p, born_momenta[jmother]);
        y_pre = 1.0 - y;
    } else { // initial state radiation
        y_pre = 1.0 - y * y;
    }

    for (int i = 0; i < N; i++) {
        double chargei = Physics::PDG::Charge( abs(pdg[i]));
        double sigma_i = sigma_f(i, pdg[i]);
        if(fabs(chargei) < 1e-6) {
            continue;
        }
        assert(fabs(born_momenta[i].Dot(born_momenta[i])) < 1e-4 &&
               "massive particles are not implemented");

        for (int j = i + 1; j < N ; j++) {
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

double SoftCollinearLimitFSR(int realpdg, int bornpdg, double s,
                         double alpha, double born_me) {
    double sp = FKS::QED::splittingTimesXiSoft(realpdg, bornpdg);
    return 16.0 * Math::Pi * alpha / s * sp * born_me;
}

double SoftCollinearLimitISR(const int *realpdgs, const int *bornpdgs, int y,
                             double s, double alpha, double born_me) {
    int realpdg = 0xffffff;
    int bornpdg = 0xffffff;

    bool has_singularity = get_splitting_ISR_QCD(bornpdgs, realpdgs, &bornpdg, &realpdg, y);

    if (!has_singularity) {
        // If the real matrix element has no singularity in y == +1 or y == -1,
        // the prefactor (1-y)(1+y) is not canceled. Therefore, the prefactor
        // evaluates to 0.0 and the overall result is 0.0.
        return 0.0; 
    };

    double sp = FKS::QED::splittingTimesXiSoft(realpdg, bornpdg);
    return 32.0 * Math::Pi * alpha / s * sp * born_me;
}

} // end namespace QED

namespace QCD {

double CollinearLimitISR(const int *realpdgs, const int *bornpdgs,
                         double xi, int y, double s, double alpha_s, double born_me) {
    int realpdg = 0xffffff;
    int bornpdg = 0xffffff;

    bool has_singularity = get_splitting_ISR_QCD(bornpdgs, realpdgs, &bornpdg, &realpdg, y);

    if (!has_singularity) {
        // If the real matrix element has no singularity in y == +1 or y == -1,
        // the prefactor (1-y)(1+y) is not canceled. Therefore, the prefactor
        // evaluates to 0.0 and the overall result is 0.0.
        return 0.0; 
    };

    double sp = FKS::QCD::splittingTimesXi(realpdg, bornpdg, xi);
    return 32.0 * Math::Pi * alpha_s / s * sp * born_me;
}

double SoftCollinearLimitISR(const int *realpdgs, const int *bornpdgs, int y,
                             double s, double alpha_s, double born_me) {
    int realpdg = 0xffffff;
    int bornpdg = 0xffffff;

    bool has_singularity = get_splitting_ISR_QCD(bornpdgs, realpdgs, &bornpdg, &realpdg, y);

    if (!has_singularity) {
        // If the real matrix element has no singularity in y == +1 or y == -1,
        // the prefactor (1-y)(1+y) is not canceled. Therefore, the prefactor
        // evaluates to 0.0 and the overall result is 0.0.
        return 0.0; 
    };

    double sp = FKS::QCD::splittingTimesXiSoft(realpdg, bornpdg);
    return 32.0 * Math::Pi * alpha_s / s * sp * born_me;
}

double CollinearLimitFSR(int realpdg, int bornpdg, double xi, double s,
                         double alpha, double born_me) {
    assert(0 && "not implemented");
    return 0.0;
}

double SoftCollinearLimitFSR(int realpdg, int bornpdg, double s, double alpha,
                             double born_me) {
    assert(0 && "not implemented");
    return 0.0;
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
        assert(fabs(born_momenta[i].Dot(born_momenta[i])) < 1e-4 &&
               "massive particles are not implemented");

        int pdgi = std::abs(pdg[i]);
        if (pdgi == 0 || (pdgi > 6 && pdgi != 21)) {
            continue;
        }
        for (int j = i + 1; j < N ; j++) {
            int pdgj = std::abs(pdg[j]);
            if (pdgj == 0 || (pdgj > 6 && pdgj != 21)) {
                 continue;
            }
            Math::FourMomentum ki = born_momenta[i];
            Math::FourMomentum kj = born_momenta[j];
            double eikonal = ki.Dot(kj) / (ki.Dot(p) * kj.Dot(p));
            limit += eikonal * ColorCorrelatedBorn.Get(i,j);
            limit += eikonal * ColorCorrelatedBorn.Get(j,i);
        }
    }
    return y_pre * 4.0 * Math::Pi * alpha_s * limit;
}

} // end namespace QCD

} // end namespace FKS
