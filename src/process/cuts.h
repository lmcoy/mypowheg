#ifndef PROCESS_CUTS_H_
#define PROCESS_CUTS_H_

class Cuts;
namespace Math {
class FourMomentum;
}

bool ApplyCuts(int n, const int *pdgs, const Math::FourMomentum *momenta,
               const Cuts &cuts);
#endif

