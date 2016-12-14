#ifndef LHE_EVENT_H
#define LHE_EVENT_H

#include "math/fourmomentum.h"

namespace LHE {

struct Particle {
    int PDG = 0;
    int Status = 0;
    int Mother1 = 0;
    int Mother2 = 0;
    int Color1 = 0;
    int Color2 = 0;
    Math::FourMomentum Momentum;
    double LifeTime = 0.0;
    double Spin = 9.0;
};

struct Event {
    int ID;
    double Alpha;
    double AlphaS;
    double Scale;
    double Weight;
    int N;
    Particle Particles[10];
};

} /* LHE  */ 



#endif /* end of include guard: LHE_EVENT_H */ 
