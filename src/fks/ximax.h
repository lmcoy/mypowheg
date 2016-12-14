#ifndef FKS_XIMAX_H_
#define FKS_XIMAX_H_

#include <cmath>
#include <cassert>

#include "libconfig.h"

namespace FKS {

/**
 * @brief Maximal energy fractions
 *
 * @see FKS::XiMaxISR() 
 * @see FKS::XiMaxFSR()
 * @see FKS::XiMaxCollinearISR()
 */
struct LIB_LOCAL Xi {
    Xi() : Max(0.0), Max_Coll1(0.0), Max_Coll2(0.0) {}
    double Max;
    double Max_Coll1;
    double Max_Coll2;
};

/**
 * @brief Max. energy fraction for ISR
 *
 * XiMax returns the maximal value for xi.
 *
 * @param x1_b momentum fraction of the underlying born
 * @param x2_b momentum fraction of the underlying born
 * @param y FKS y of the emitted particle
 */
LIB_LOCAL double XiMaxISR(double x1_b, double x2_b, double y);

/**
 * @brief Max. energy fraction for FSR
 *
 * XiMaxFSR returns the maximal value for xi for FSR radiation.
 *
 * @param sqrts \f$\sqrt{s}\f$ of the partonic process
 * @param len_kj_born length of momentum of the radiating born particle.
 */
LIB_LOCAL double XiMaxFSR(double sqrts, double len_kj_born);

/**
 * @brief Max. energy fraction for ISR in collinear limit
 *
 * XiMaxCollinearSR returns the maximal value for xi in the collinear limit.
 *
 * @param x1_b momentum fraction of the underlying born
 * @param x2_b momentum fraction of the underlying born
 * @param y FKS y of the emitted particle (-1 or 1)
 */
LIB_LOCAL double XiMaxCollinearISR(double x1_b, double x2_b, int y);

} // namespace FKS

#endif
