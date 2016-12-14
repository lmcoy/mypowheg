#ifndef G_FSR_H_
#define G_FSR_H_

#include "libconfig.h"

namespace Phasespace {
class Phasespace;
}

namespace UserProcess {
struct Data;
}

namespace Util {
class Matrix2;
}

namespace FKS {

class Param;
class Param_as;
struct Xi;
struct Real_t;
struct Region;
struct MatrixElement;
struct Phasespaces;

namespace QED {

/**
 * @brief limits of FKS subtracted real EW matrix element
 *
 * Limits returns the soft, collinear and soft/collinear limits of SxG().
 *
 * @see QED::SxG
 */
LIB_LOCAL FKS::MatrixElement
Limits(const Real_t &real, const FKS::Phasespaces &PS, double bornme,
       const Util::Matrix2 &collcorr, const Region &region, double x,
       const FKS::Xi &Xi, double y, double phi, const UserProcess::Data *);

/**
 * @brief FKS subtracted real EW matrix element
 *
 * SxG returns S * R where S is the FKS S function and R is the matrix element
 * multiplied by xi^2*(1-y) for FSR and and xi^2*(1-y^2) for ISR.
 * The result has to be integrated over the n+1 phase-space and the FKS plus
 * distributions in the FKS subtraction scheme.
 *
 * This function uses the real matrix element for EW corrections.
 */
LIB_LOCAL double SxG(const Real_t &real, const FKS::Phasespaces &PS,
                     double bornme, const Util::Matrix2 &collcorr,
                     const Region &region, double xi, double y, double phi,
                     const UserProcess::Data *);

} // namespace QED

namespace QCD {

/**
 * @brief limits of FKS subtracted real QCD matrix element
 *
 * Limits returns the soft, collinear and soft/collinear limits of SxG().
 * 
 * @see QCD::SxG
 */
LIB_LOCAL FKS::MatrixElement
Limits(const Real_t &real, const FKS::Phasespaces &PS, double bornme,
       const Util::Matrix2 &collcorr, const Region &region, double x,
       const FKS::Xi &Xi, double y, double phi, const UserProcess::Data *);

/**
 * @brief FKS subtracted real QCD matrix element
 *
 * SxG returns S * R where S is the FKS S function and R is the matrix element
 * multiplied by xi^2*(1-y) for FSR and and xi^2*(1-y^2) for ISR.
 * The result has to be integrated over the n+1 phase-space and the FKS plus
 * distributions in the FKS subtraction scheme.
 *
 * This function uses the real matrix element for QCD corrections.
 */
LIB_LOCAL double SxG(const Real_t &real, const FKS::Phasespaces &PS,
                     double bornme, const Util::Matrix2 &collcorr,
                     const Region &region, double xi, double y, double phi,
                     const UserProcess::Data *);

} // namespace QCD

} // end namespace FKS

#endif /* end of include guard: G_FSR_H_ */

