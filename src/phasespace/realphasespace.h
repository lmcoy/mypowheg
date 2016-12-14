#ifndef PHASESPACE_REALPHASESPACE_H_
#define PHASESPACE_REALPHASESPACE_H_

#include "phasespace/phasespace.h"

#define RESTRICT __restrict__

namespace Phasespace {

void GenRealPhasespaceISR(Phasespace * RESTRICT ps_real,
                               Phasespace const *const RESTRICT ps_born, double xi,
                               double y, double phi, bool boost_to_cms = true);

void GenRealPhasespaceFSR(Phasespace * RESTRICT ps_real,
                          Phasespace const *const RESTRICT ps_born, int j, double xi,
                          double y, double phi);

void GenRealPhasespace(Phasespace * RESTRICT ps_real,
                          Phasespace const *const RESTRICT ps_born, int j, double xi,
                          double y, double phi);
} // end namespace Phasespace

#endif
