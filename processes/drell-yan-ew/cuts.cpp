#include "cuts.h"
#include "process/cuts.h"

#include <cmath>

#include "math/fourmomentum.h"
bool ApplyCuts(int n, const int *pdgs, const Math::FourMomentum *momenta,
               const Cuts &cuts) {
    // cut on inv. mass of leptons
    Math::FourMomentum mll = momenta[2].Plus(momenta[3]);
    if (mll.Dot(mll) < cuts.mllmin * cuts.mllmin) {
        return false;
    }
    if (mll.Dot(mll) > cuts.mllmax * cuts.mllmax) {
        return false;
    }

    // cut on lepton pseudo rapidity
    if ( momenta[2].Eta() > cuts.EtaMax) {
        return false;
    }
    if (momenta[2].Eta() < cuts.EtaMin) {
        return false;
    }
    if ( momenta[3].Eta() > cuts.EtaMax) {
        return false;
    }
    if (momenta[3].Eta() < cuts.EtaMin) {
        return false;
    }

    // cut on lepton rapidity
    if ( momenta[2].Rapidity() > cuts.Ymax) {
        return false;
    }
    if (momenta[2].Rapidity() < cuts.Ymin) {
        return false;
    }
    if ( momenta[3].Rapidity() > cuts.Ymax) {
        return false;
    }
    if (momenta[3].Rapidity() < cuts.Ymin) {
        return false;
    }

    // cut on lepton pT
    if (momenta[2].PT() < cuts.pTmin) {
        return false;
    }
    if (momenta[3].PT() < cuts.pTmin) {
        return false;
    }
    if (momenta[2].PT() > cuts.pTmax) {
        return false;
    }
    if (momenta[3].PT() > cuts.pTmax) {
        return false;
    }
    return true;
}
