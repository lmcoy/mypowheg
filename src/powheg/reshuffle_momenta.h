#ifndef RESHUFFLE_MOMENTA_H_0HKNZ816
#define RESHUFFLE_MOMENTA_H_0HKNZ816

#include "libconfig.h"

namespace LHE {
class Event;
}
namespace UserProcess {
class Data;
}

namespace Powheg {

class Resonance;

void LIB_LOCAL reshuffle_momenta(LHE::Event *event,
                                 const Powheg::Resonance &resonance,
                                 const UserProcess::Data *userdata);
} /* Powheg */

#endif /* end of include guard: RESHUFFLE_MOMENTA_H_0HKNZ816 */
