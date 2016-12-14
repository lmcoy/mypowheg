#ifndef PHASESPACE_RECOMBINATION_H_
#define PHASESPACE_RECOMBINATION_H_

#include "phasespace/phasespace.h"

namespace Phasespace {
Phasespace Recombine(const int * pdgs, const Phasespace &ps_real, double dR);
}

#endif
