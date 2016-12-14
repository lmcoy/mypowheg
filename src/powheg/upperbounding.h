#ifndef POWHEG_UPPERBOUNDING_H_
#define POWHEG_UPPERBOUNDING_H_

#include "libconfig.h"
#include "fks/process.h"

namespace FKS { class RadiationRegion; }
namespace UserProcess { class Data; }
namespace Phasespace { class Phasespace; }

namespace Powheg {

LIB_PUBLIC void GetNormForBoundingFcts(FKS::Type_t type, bool,
                                       FKS::RadiationRegion &radregion,
                                       const Phasespace::Phasespace &ps,
                                       double x1, double x2, double x3,
                                       UserProcess::Data *userdata);

} // end namespace Powheg

#endif
