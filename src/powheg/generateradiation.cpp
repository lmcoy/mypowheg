#include "powheg/generateradiation.h"

#include <cassert>

#include "powheg/radiationgenerator.h"
#include "powheg/qcdisr.h"
#include "powheg/qedisr.h"
#include "powheg/qcdfsr.h"
#include "powheg/qedfsr.h"

#include "fks/process.h"
#include "random/rnd.h"
#include "phasespace/phasespace.h"
#include "process/data.h"

namespace Powheg {

namespace QCD {

RadiationType GenerateRadiation(int pdf, double B,
                                const FKS::RadiationRegion &radreg,
                                const Phasespace::Phasespace &ps, double pTmin2,
                                UserProcess::Data *userdata, Random::RNG *rng,
                                Powheg::Radiation *rad) {

    RadiationAlphaS alphas(userdata->AlphaS);
    if (radreg.Region.J == 0) {
        QCDISR qcdisr(alphas);
        RadiationGenerator<QCDISR> generator(pTmin2);
        return generator.Generate(qcdisr, pdf, B, radreg, ps, userdata, rng, rad);
    }
    QCDFSR qcdfsr(alphas);
    RadiationGenerator<QCDFSR> generator(pTmin2);

    return generator.Generate(qcdfsr, pdf, B, radreg, ps, userdata, rng, rad);
}

} /* QCD */

namespace QED {

RadiationType GenerateRadiation(int pdf, double B,
                                const FKS::RadiationRegion &radreg,
                                const Phasespace::Phasespace &ps, double pTmin2,
                                UserProcess::Data *userdata, Random::RNG *rng,
                                Powheg::Radiation *rad) {

    if (radreg.Region.J == 0) {
        QEDISR qedisr;
        RadiationGenerator<QEDISR> generator(pTmin2);
        return generator.Generate(qedisr, pdf, B, radreg, ps, userdata, rng, rad);
    }
    QEDFSR qedfsr;
    RadiationGenerator<QEDFSR> generator(pTmin2);

    return generator.Generate(qedfsr, pdf, B, radreg, ps, userdata, rng, rad);
}

} /* QED */

} /* Powheg */
