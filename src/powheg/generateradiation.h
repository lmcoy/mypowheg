#ifndef POWHEG_RADIATIONGENERATOR_H
#define POWHEG_RADIATIONGENERATOR_H

#include "libconfig.h"

namespace Phasespace { class Phasespace; }
namespace FKS { class RadiationRegion; class Real_t; }
namespace Random { class RNG; }
namespace UserProcess { class Data; }

namespace Powheg {

struct LIB_LOCAL Radiation {
    double xi = 0.0;
    double y = 0.0;
    double phi = 0.0;
    const FKS::Real_t * Real;
    double kT2 = 0.0;
    int j = 0;
};

enum class RadiationType {
    BORN,   ///< born like event
    REAL,   ///< real like event
    ERROR,  ///< error
    ENORM,  ///< normalization for upper bounding function is to small
    ERADVAR ///< radiation variables are out of bounds
};

namespace QCD {

LIB_LOCAL RadiationType
GenerateRadiation(int pdf, double B, const FKS::RadiationRegion &radreg,
                  const Phasespace::Phasespace &ps, double,
                  UserProcess::Data *userdata, Random::RNG *rng,
                  Powheg::Radiation *rad);

} /* QCD */

namespace QED {

LIB_LOCAL RadiationType
GenerateRadiation(int pdf, double B, const FKS::RadiationRegion &radreg,
                  const Phasespace::Phasespace &ps, double,
                  UserProcess::Data *userdata, Random::RNG *rng,
                  Powheg::Radiation *rad);

} /* QED */

} // end namespace Powheg


#endif /* end of include guard: POWHEG_RADIATIONGENERATOR_H */ 
