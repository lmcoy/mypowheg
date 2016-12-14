#ifndef FKS_XSEC_H_
#define FKS_XSEC_H_

#include <array>

#include "libconfig.h"

namespace UserProcess {
struct Data;
}

namespace Phasespace {
class Phasespace;
}

namespace FKS {

constexpr size_t MAXPDF = 6;
typedef std::array<double, MAXPDF> Result;

/**
 * @brief XSecFullByPDF returns the cross section integrand.
 *
 * @return
 * This function returns the integrand of
 * \sigma = \int dps dx_1 dx_2 dx_3 f(ps,x_1, x_2, x_3)
 * The elements of the returned value represent the different pdfs which are
 * connected to the flavour configs. The jacobian of the born phase-space is
 * included in this function.
 *
 * @param ps  born phase space
 * @param x1  value between [0,1]. Represents the energy of n+1 particle.
 * @param x2  value between [0,1]. Represents cos theta of n+1 particle.
 * @param x3  value between [0,1]. Represents phi of n+1 particle.
 * @param wgt weight of the event
 * @param params process definition
 */
LIB_PUBLIC Result XSecFullByPDF(const Phasespace::Phasespace &ps, double x1, double x2,
                     double x3, double wgt, UserProcess::Data *params);
/**
 * @brief Born part of XSecFullByPDF()
 */
LIB_PUBLIC Result XSecBornByPDF(const Phasespace::Phasespace &ps, double x1, double x2,
                     double x3, double wgt, UserProcess::Data *params);

/**
 * @brief Real part of XSecFullByPDF()
 */
LIB_PUBLIC Result XSecRealByPDF(const Phasespace::Phasespace &ps, double x1, double x2,
                     double x3, double wgt, UserProcess::Data *params);

/**
 * @brief Virtual part of XSecFullByPDF()
 */
LIB_PUBLIC Result XSecVirtualByPDF(const Phasespace::Phasespace &ps, double x1, double x2,
                     double x3, double wgt, UserProcess::Data *params);

/**
 * @brief PDF renorm. remnant of XSecFullByPDF()
 */
LIB_PUBLIC Result XSecRemnantByPDF(const Phasespace::Phasespace &ps, double x1, double x2,
                     double x3, double wgt, UserProcess::Data *params);
/**
 * @brief XSecFull returns the cross section integrand.
 *
 * XSecFull returns the sum of XSecFullByPDF().
 */
LIB_PUBLIC int XSecFull(const Phasespace::Phasespace &, double, double, double, double,
             double *, UserProcess::Data *);

/**
 * @brief Real part of XSecFull().
 */
LIB_PUBLIC int XSecReal(const Phasespace::Phasespace &, double, double, double, double,
             double *, UserProcess::Data *);

/**
 * @brief Virtual part of XSecFull().
 */
LIB_PUBLIC int XSecVirtual(const Phasespace::Phasespace &, double, double, double, double,
                double *, UserProcess::Data *);

/**
 * @brief Born part of XSecFull().
 */
LIB_PUBLIC int XSecBorn(const Phasespace::Phasespace &, double, double, double, double,
             double *, UserProcess::Data *);

/**
 * @brief Collinear remnant from pdf renormalization of XSecFull()
 */
LIB_PUBLIC int XSecRemnant(const Phasespace::Phasespace &, double, double,
                           double, double, double *, UserProcess::Data *);

/**
 * @brief Virtual part + collinear remnant of XSecFull().
 */
LIB_PUBLIC int XSecVirtualPlusRemnant(const Phasespace::Phasespace &, double, double,
                           double, double, double *, UserProcess::Data *);

} // end namespace FKS

#endif

