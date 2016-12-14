#ifndef CUTS_H_
#define CUTS_H_

#include "process/cuts.h"

class DrellYanCuts : public UserProcess::ICuts {
  public:
    virtual bool ApplyCuts(int n, const int *pdgs,
                           const Math::FourMomentum *momenta) const;

    double EtaMin = -1e10;
    double EtaMax = 1e10;
    double pTmin = 1.0;
    double pTmax = 1e10;
    double mllmin = 30;
    double mllmax = 1e10;
    double Ymin = -1e10;
    double Ymax = 1e10;
};

#endif /* end of include guard: CUTS_H_ */

