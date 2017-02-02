#ifndef FKS_VIRTUAL_H_
#define FKS_VIRTUAL_H_

// forward declaration
namespace Math {
class FourMomentum;
}

namespace Util {
class Matrix2;
}

namespace Phasespace {
class Phasespace;
}

namespace FKS {

// forward declaration
struct Scales;

namespace QCD {

double Virtual(int n, int bornpdgs[], const Math::FourMomentum *momenta,
               double born, const Util::Matrix2 &ColorCorrelated, double Vfin,
               double alpha_s, double sqrts, const FKS::Scales &scales);

double Eps2Pole(int n, int bornpdgs[], double bornme);

double EpsPole(int n, int bornpdgs[], const Math::FourMomentum *momenta,
               double born, const Util::Matrix2 &ColorCorrelatedBorn,
               double Q2);

} // namespace QCD

namespace QED {

double Virtual(int n, int bornpdgs[], const Math::FourMomentum *momenta,
               double born, const Util::Matrix2 &unused, double Vfin,
               double alpha, double sqrts, const FKS::Scales &scales);

double Eps2Pole(int n, int bornpdgs[], double bornme);

double EpsPole(int n, int bornpdgs[], const Math::FourMomentum *momenta,
               double born, double Q2);

double VirtualMReg(const Phasespace::Phasespace &ps, int *pdgs, double *m,
                   double lambda, double muF2);

} // namespace QED

} // namespace FKS

#endif
