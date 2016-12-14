#ifndef CUTS_H_SDFI3FSS
#define CUTS_H_SDFI3FSS

#include "process/cuts.h"

class WCuts : public UserProcess::ICuts {
  public:
    virtual bool ApplyCuts(int n, const int *pdgs,
                           const Math::FourMomentum *momenta) const;
};

#endif /* end of include guard: CUTS_H_ */

