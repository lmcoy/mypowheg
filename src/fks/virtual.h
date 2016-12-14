#ifndef FKS_VIRTUAL_H_
#define FKS_VIRTUAL_H_

// forward declaration
namespace Math {
class FourMomentum;
}

namespace Util {
class Matrix2;
}

namespace FKS {

// forward declaration
class Scales;

namespace QCD {

double Virtual(int n, int bornpdgs[], const Math::FourMomentum *momenta,
               double born, const Util::Matrix2 &ColorCorrelated, double Vfin,
               double alpha_s, double sqrts, const FKS::Scales &scales);

} // namespace QCD 

namespace QED {

double Virtual(int n, int bornpdgs[], const Math::FourMomentum *momenta,
               double born, const Util::Matrix2 &unused, double Vfin,
               double alpha, double sqrts, const FKS::Scales &scales);

double Eps2Pole(int n, int bornpdgs[], double bornme);

double EpsPole(int n, int bornpdgs[], const Math::FourMomentum *momenta,
               double born, double Q2);

} // namespace QED

} // namespace FKS

#endif

