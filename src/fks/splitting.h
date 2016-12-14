#ifndef FKS_SPLITTING_H_
#define FKS_SPLITTING_H_

#include "libconfig.h"

namespace FKS {

struct Splitting {
    Splitting(int realpdg, int bornpdg) : BornPDG(bornpdg), RealPDG(realpdg) {}
    int BornPDG;
    int RealPDG;
};

namespace QCD {
LIB_LOCAL double splittingEps(int real_parton, int born_parton, double xi);
LIB_LOCAL double splittingTimesXi(int real_parton, int born_parton, double xi);
LIB_LOCAL double splittingTimesXiSoft(int real_parton, int born_parton);
} // QCD

namespace QED {
LIB_LOCAL double splittingTimesXi(int realpdg, int bornpdg, double xi);
LIB_LOCAL double splittingTimesXiSoft(int realpdg, int bornpdg);
LIB_LOCAL double splittingEps(int realpdg, int bornpdg, double xi);
} // QED 

} // end namespace FKS
#endif
