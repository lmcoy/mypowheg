#ifndef CUTS_H_JSD7NSO3
#define CUTS_H_JSD7NSO3

#include "process/cuts.h"

class WjCuts : public UserProcess::ICuts {
  public:
    virtual bool ApplyCuts(int n, const int *pdgs,
                           const Math::FourMomentum *momenta) const;
    double PTCUT = 5.0;
    double jetDR = 0.4;
    double z_thr = 0.8;
    double photonDR = 0.1;
};

#endif /* end of include guard: CUTS_H_ */
