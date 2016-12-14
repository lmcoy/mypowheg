#ifndef PHASESPACE_TWOPARTICLEGENERATOR_H_
#define PHASESPACE_TWOPARTICLEGENERATOR_H_

#include "phasespace/phasespace.h"

namespace Phasespace {

/**
Generator for a two particle phase space.
*/
class TwoParticleGenerator {
  public:
    void operator()(Phasespace *ps, int n, double const *const x, double S,
                    int m, double *masses) const;
};

} // end namspace phasespace
#endif
