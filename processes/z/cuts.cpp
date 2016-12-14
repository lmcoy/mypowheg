#include "cuts.h"
#include "process/cuts.h"

#include <cmath>

#include "math/fourmomentum.h"

bool DrellYanCuts::ApplyCuts(int n, const int *pdgs,
                             const Math::FourMomentum *momenta) const {
    // cut on inv. mass of leptons
    Math::FourMomentum mll = momenta[2].Plus(momenta[3]);
    if (mll.Dot(mll) < mllmin * mllmin) {
        return false;
    }
    if (mll.Dot(mll) > mllmax * mllmax) {
        return false;
    }

    // cut on lepton pseudo rapidity
    if (momenta[2].Eta() > EtaMax) {
        return false;
    }
    if (momenta[2].Eta() < EtaMin) {
        return false;
    }
    if (momenta[3].Eta() > EtaMax) {
        return false;
    }
    if (momenta[3].Eta() < EtaMin) {
        return false;
    }

    // cut on lepton rapidity
    if (momenta[2].Rapidity() > Ymax) {
        return false;
    }
    if (momenta[2].Rapidity() < Ymin) {
        return false;
    }
    if (momenta[3].Rapidity() > Ymax) {
        return false;
    }
    if (momenta[3].Rapidity() < Ymin) {
        return false;
    }

    // cut on lepton pT
    if (momenta[2].PT() < pTmin) {
        return false;
    }
    if (momenta[3].PT() < pTmin) {
        return false;
    }
    if (momenta[2].PT() > pTmax) {
        return false;
    }
    if (momenta[3].PT() > pTmax) {
        return false;
    }
    return true;
}
