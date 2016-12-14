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
/**
 * @brief order epsilon QCD splitting fct.
 *
 * splittingEps returns the order(epsilon) part of the QCD splitting functions.
 *
 * @param real_parton real parton
 * @param born_parton this parton splits into the real parton
 * @param xi          energy fraction (xi=1-z)
 */
LIB_LOCAL double splittingEps(int real_parton, int born_parton, double xi);
    
/**
 * @brief QCD splitting function times xi
 * 
 * splittingTimesXi returns the QCD splitting function times xi.
 *
 * @param real_parton real parton
 * @param born_parton this parton splits into the real parton
 * @param xi          energy fraction (xi=1-z)
 */
LIB_LOCAL double splittingTimesXi(int real_parton, int born_parton, double xi);

/**
 * @brief QCD splitting fct. times xi for xi=0
 *
 * splittingTimesXiSoft returns the QCD splitting function times xi
 * in the limit xi = 0.
 *
 * @param real_parton real parton
 * @param born_parton this parton splits into the real parton
 */
LIB_LOCAL double splittingTimesXiSoft(int real_parton, int born_parton);
} // QCD

namespace QED {
LIB_LOCAL double splittingTimesXi(int realpdg, int bornpdg, double xi);
LIB_LOCAL double splittingTimesXiSoft(int realpdg, int bornpdg);
LIB_LOCAL double splittingEps(int realpdg, int bornpdg, double xi);
} // QED 

} // end namespace FKS
#endif
