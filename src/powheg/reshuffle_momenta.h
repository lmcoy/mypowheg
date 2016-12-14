#ifndef RESHUFFLE_MOMENTA_H_0HKNZ816
#define RESHUFFLE_MOMENTA_H_0HKNZ816

#include "libconfig.h"

namespace LHE {
struct Event;
}
namespace UserProcess {
struct Data;
}

namespace Powheg {

struct Resonance;

void LIB_LOCAL reshuffle_momenta(LHE::Event *event,
                                 const Powheg::Resonance &resonance,
                                 const UserProcess::Data *userdata);
} /* Powheg */

#endif /* end of include guard: RESHUFFLE_MOMENTA_H_0HKNZ816 */
