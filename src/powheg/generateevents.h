#ifndef POWHEG_GENERATEEVENTS_H_ 
#define POWHEG_GENERATEEVENTS_H_

#include "libconfig.h"

namespace Phasespace {
class Phasespace;
}
namespace UserProcess {
struct Data;
}

namespace Powheg {

LIB_PUBLIC int GenerateEvents(const Phasespace::Phasespace &, double, double,
                              double, double, double *, UserProcess::Data *);

} // end namespace Powheg

#endif
