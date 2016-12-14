#include "phasespace/recombination.h"

#include "physics/pdgcode.h"

namespace Phasespace {
Phasespace Recombine(const int *pdgs, const Phasespace &ps_real,
                                 double dR) {
    Phasespace ps_recombined;
    ps_recombined = ps_real;
    int i_photon = ps_real.N + 1;
    if (!Physics::PDG::IsPhoton(pdgs[i_photon])) {
        return ps_recombined;
    }
    if (ps_real.Momenta[i_photon].MomentumMagnitudeSqr() > 1e-16) {
        for (int i = 2; i < i_photon; i++) {
            if( !Physics::PDG::IsChargedLepton(pdgs[i]) ) {
                continue;
            }
            if (ps_real.Momenta[i].MomentumMagnitudeSqr() > 1e-16) {
                double dr = Math::FourMomentum::DeltaR(
                    ps_real.Momenta[i], ps_real.Momenta[i_photon]);
                if (dr < dR) {
                    ps_recombined.Momenta[i] =
                        ps_real.Momenta[i].Plus(ps_real.Momenta[i_photon]);
                    break;
                }
            }
        }
    }

    return ps_recombined;
}

}
