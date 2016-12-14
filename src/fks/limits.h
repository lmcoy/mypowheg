#ifndef FKS_LIMITS_H_
#define FKS_LIMITS_H_

#include <array>

#include "libconfig.h"

#include "phasespace/phasespace.h"
#include "util/matrix.h"

namespace FKS {

struct LIB_LOCAL MatrixElement {
    MatrixElement()
        : Soft(0.0), Collinear1(0.0), Collinear2(0.0), SoftCollinear1(0.0),
          SoftCollinear2(0.0), Real(0.0) {
    }
    double Soft;
    double Collinear1;
    double Collinear2;
    double SoftCollinear1;
    double SoftCollinear2;
    double Real;
    void SetLimits(const MatrixElement& lim) {
        Soft = lim.Soft;
        Collinear1 = lim.Collinear1;
        Collinear2 = lim.Collinear2;
        SoftCollinear1 = lim.SoftCollinear1;
        SoftCollinear2 = lim.SoftCollinear2;
    }
};

namespace QED {

LIB_LOCAL double CollinearLimitFSR(int realpdg, int bornpdg, double xi,
                                   double s, double alpha, double born_me);

LIB_LOCAL double CollinearLimitISR(const int *realpdgs, const int *bornpdgs,
                                   double xi, int y, double s, double alpha,
                                   double born_me);

LIB_LOCAL double SoftLimit(int N, const Math::FourMomentum *born_momenta,
                           const int *pdg, int realpdg, double s, int jmother,
                           double alpha, const Util::Matrix2 &Born, double y,
                           double phi);

LIB_LOCAL double SoftCollinearLimitFSR(int realpdg, int bornpdg, double s,
                                       double alpha, double born_me);

LIB_LOCAL double SoftCollinearLimitISR(const int *realpdgs, const int *bornpdgs,
                                       int y, double s, double alpha,
                                       double born_me);
} // end namespace QED

namespace QCD {

LIB_LOCAL double CollinearLimitFSR(int realpdg, int bornpdg, double xi,
                                   double s, double alpha, double born_me);

LIB_LOCAL double CollinearLimitISR(const int *realpdgs, const int *bornpdgs,
                                   double xi, int y, double s, double alpha_s,
                                   double born_me);

LIB_LOCAL double SoftLimit(int N, const Math::FourMomentum *born_momenta,
                           const int *pdg, int realpdg, double s, int jmother,
                           double alpha_s,
                           const Util::Matrix2 &ColorCorrelatedBorn, double y,
                           double phi);

LIB_LOCAL double SoftCollinearLimitFSR(int realpdg, int bornpdg, double s,
                                       double alpha, double born_me);

LIB_LOCAL double SoftCollinearLimitISR(const int *realpdgs, const int *bornpdgs,
                                       int y, double s, double alpha_s,
                                       double born_me);
} // end namespace QCD

} // end namespace FKS

#endif
