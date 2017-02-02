#ifndef POWHEG_ROVERB_H_
#define POWHEG_ROVERB_H_

#include <array>
#include <vector>

#include <fks/process.h>
#include <fks/regions.h>

#include "util/staticmatrix.h"

namespace Phasespace {
class Phasespace;
}
namespace UserProcess {
struct Data;
}
namespace FKS {
class PDF;
struct RadiationRegion;
}

namespace Powheg {

static constexpr size_t MaxR = 12;

std::array<double, MaxR> RoverB(double B, const FKS::RadiationRegion &radreg,
                                const Phasespace::Phasespace &ps,
                                const Phasespace::Phasespace &ps_real,
                                double rad_alpha, UserProcess::Data *userdata);

void LumiRatio(const FKS::RadiationRegion &radreg,
               const Phasespace::Phasespace &ps_born,
               const Phasespace::Phasespace &ps_real,
               UserProcess::Data *userdata, double scale, int usepdf,
               Util::StaticMatrix128 *output);

} // end namespace powheg

#endif
