#include "cuts.h"
#include "process/cuts.h"

#include <cmath>

#include "math/fourmomentum.h"

bool WCuts::ApplyCuts(int n, const int *pdgs,
                      const Math::FourMomentum *momenta) const {
    double mllmin = 1.0;
    // cut on inv. mass of leptons
    Math::FourMomentum mll = momenta[2].Plus(momenta[3]);
    if (mll.Dot(mll) < mllmin * mllmin) {
        return false;
    }
    return true;
}
