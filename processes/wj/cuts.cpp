#include "cuts.h"
#include "process/cuts.h"

#include <cmath>

#include "math/fourmomentum.h"

static bool isphoton(int pdg) {
    return pdg == 22;
}

static bool isquark(int pdg) {
    return pdg > -6 && pdg < 6 && pdg != 0;
}

bool WjCuts::ApplyCuts(int n, const int *pdgs,
                       const Math::FourMomentum *momenta) const {
    double pT = momenta[4].PT();

    if (momenta[5].E() < 1e-22) {
        if (pT < PTCUT) {
            return false;
        }
        return true;
    } else {
        double dR = Math::FourMomentum::DeltaR(momenta[4], momenta[5]);
        if (isquark(pdgs[4]) && isphoton(pdgs[5]) && dR < photonDR) {
            // recombine
            auto p = momenta[4].Plus(momenta[5]);
            if (p.PT() > PTCUT) {
                return true;
            }
            return false;
        }
        if (dR < jetDR) {
            // recombine jet
            auto p = momenta[4].Plus(momenta[5]);
            if (isphoton(pdgs[5])) {
                double Ea = momenta[5].E();
                double Ejet = p.E();
                if (Ea / Ejet < z_thr && p.PT() > PTCUT) {
                    return true;
                }
                return false;
            }
            if (p.PT() > PTCUT) {
                return true;
            }
            return false;
        }
        if (momenta[4].PT() > PTCUT) {
            // at least one QCD jet
            return true;
        }
        if (!isphoton(pdgs[5]) && momenta[5].PT() > PTCUT) {
            return true;
        }
        return false;
    }
    return true;
}
