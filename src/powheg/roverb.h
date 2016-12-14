#ifndef POWHEG_ROVERB_H_
#define POWHEG_ROVERB_H_

#include <vector>
#include <fks/regions.h>
#include <fks/process.h>
#include <array>

#include "util/staticmatrix.h"

namespace Phasespace { class Phasespace; }
namespace UserProcess { class Data; }
namespace FKS { class PDF; class RadiationRegion; }

namespace Powheg {

std::array<double, 8> RoverB(double B, const FKS::RadiationRegion &radreg,
                             const Phasespace::Phasespace &ps,
                             const Phasespace::Phasespace &ps_real,
                             double rad_alpha, UserProcess::Data *userdata);

void LumiRatio(const FKS::RadiationRegion &radreg,
               const Phasespace::Phasespace &ps_born,
               const Phasespace::Phasespace &ps_real,
               UserProcess::Data *userdata, double scale, int usepdf,
               Util::StaticMatrix32 *output);

} // end namespace powheg

#endif
