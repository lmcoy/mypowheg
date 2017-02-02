#ifndef FKS_LIMITS_H_
#define FKS_LIMITS_H_

#include <array>

#include "libconfig.h"

#include "phasespace/phasespace.h"
#include "process/matrixelement.h"
#include "util/matrix.h"

namespace FKS {

struct LIB_LOCAL MatrixElement {
    MatrixElement()
        : Soft(0.0), Collinear1(0.0), Collinear2(0.0), SoftCollinear1(0.0),
          SoftCollinear2(0.0), Real(0.0) {}
    double Soft;
    double Collinear1;
    double Collinear2;
    double SoftCollinear1;
    double SoftCollinear2;
    double Real;
    void SetLimits(const MatrixElement &lim) {
        Soft = lim.Soft;
        Collinear1 = lim.Collinear1;
        Collinear2 = lim.Collinear2;
        SoftCollinear1 = lim.SoftCollinear1;
        SoftCollinear2 = lim.SoftCollinear2;
    }
};

namespace QED {

LIB_LOCAL double CollinearLimitFSR(int realpdg, const int bornpdg,
                                   const Phasespace::Phasespace &ps,
                                   int splitting_j, double xi, double phi,
                                   double alpha, double bornme,
                                   const UserProcess::SpinCorrelated &spincorr);

LIB_LOCAL double CollinearLimitISR(int Nborn, const int *realpdgs,
                                   const int *bornpdgs, double xi, int y,
                                   double phi, double sb, double alpha,
                                   double bornme,
                                   const UserProcess::SpinCorrelated &spincorr);

LIB_LOCAL double SoftLimit(int N, const Math::FourMomentum *born_momenta,
                           const int *pdg, int realpdg, double s, int jmother,
                           double alpha, const Util::Matrix2 &Born, double y,
                           double phi);

LIB_LOCAL double
SoftCollinearLimitFSR(int realpdg, const int bornpdg,
                      const Phasespace::Phasespace &ps, int splitting_j,
                      double phi, double alpha, double bornme,
                      const UserProcess::SpinCorrelated &spincorr);

LIB_LOCAL double
SoftCollinearLimitISR(int Nborn, const int *realpdgs, const int *bornpdgs,
                      int y, double phi, double sb, double alpha, double bornme,
                      const UserProcess::SpinCorrelated &spincorr);
} // end namespace QED

namespace QCD {

LIB_LOCAL double CollinearLimitFSR(int realpdg, const int bornpdg,
                                   const Phasespace::Phasespace &ps,
                                   int splitting_j, double xi, double phi,
                                   double alpha, double bornme,
                                   const UserProcess::SpinCorrelated &spincorr);

LIB_LOCAL double CollinearLimitISR(int Nborn, const int *realpdgs,
                                   const int *bornpdgs, double xi, int y,
                                   double phi, double sb, double alpha,
                                   double bornme,
                                   const UserProcess::SpinCorrelated &spincorr);

LIB_LOCAL double SoftLimit(int N, const Math::FourMomentum *born_momenta,
                           const int *pdg, int realpdg, double s, int jmother,
                           double alpha_s,
                           const Util::Matrix2 &ColorCorrelatedBorn, double y,
                           double phi);

LIB_LOCAL double
SoftCollinearLimitFSR(int realpdg, const int bornpdg,
                      const Phasespace::Phasespace &ps, int splitting_j,
                      double phi, double alpha, double bornme,
                      const UserProcess::SpinCorrelated &spincorr);

LIB_LOCAL double
SoftCollinearLimitISR(int Nborn, const int *realpdgs, const int *bornpdgs,
                      int y, double phi, double sb, double alpha, double bornme,
                      const UserProcess::SpinCorrelated &spincorr);
} // end namespace QCD

} // end namespace FKS

#endif
