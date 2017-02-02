#ifndef MATRIXELEMENT_H_SHDNAK6A
#define MATRIXELEMENT_H_SHDNAK6A

#include "process/matrixelement.h"

namespace FKS {
class Param;
class Param_as;
}

namespace Phasespace {
class Phasespace;
}

enum class AlphaScheme { Alpha0, Gmu };

class MatrixElementWj : public UserProcess::IMatrixElement {
  public:
    void Init(const FKS::Param *param, AlphaScheme scheme, double inveps);

    virtual Result Born(int process, const Phasespace::Phasespace &ps, double Q,
                        const FKS::Param *param, bool QCD, bool EW);

    virtual double Real(int process, const Phasespace::Phasespace &ps,
                        const FKS::Param *param,
                        Diagrams diagrams = Diagrams::ALL);
};

#endif