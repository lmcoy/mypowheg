#ifndef FKS_REMNANTS_H_
#define FKS_REMNANTS_H_

#include "fks/splitting.h"

namespace FKS {

// forward declarations
struct Scales;
struct LumRemnants;

enum class PDFRenorm {
    MSbar,
    DIS
};

namespace QED {

/**
 * @brief PDF renormalization
 *
 * Remnant is the integrand of the pdf renormalization term which is integrated
 * in xi in (0,1-x).
 */
double Remnant(const Scales &scales, const Splitting &sp, double xi,
               double xi_max, double s_b, double alpha,
               const LumRemnants &lumir, double born, PDFRenorm pdfren);
} // QED

namespace QCD {

/**
 * @brief PDF renormalization
 *
 * Remnant is the integrand of the pdf renormalization term which is integrated
 * in xi in (0,1-x).
 */
double Remnant(const Scales &scales, const Splitting &sp, double xi,
               double xi_max, double s_b, double alpha_s,
               const LumRemnants &lumir, double born, PDFRenorm pdfren);
} // QCD

} // FKS

#endif /* end of include guard: FKS_REMNANTS_H_ */

