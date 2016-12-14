#ifndef MATRIXELEMENT_H_DAHDA872
#define MATRIXELEMENT_H_DAHDA872

#include "process/matrixelement.h"

namespace FKS {
class Param;
class Param_as;
}

namespace Phasespace {
class Phasespace;
}

enum class AlphaScheme {
    Alpha0,
    Gmu
};

class MatrixElementW : public UserProcess::IMatrixElement {
  public:
    struct SubProcesses {
        // QCD
        static constexpr int QXG_MUPNU_QX = 0;
        static constexpr int QXQ_MUPNU_G = 1;
        static constexpr int GQX_MUPNU_QX = 2;
        static constexpr int GQ_MUPNU_Q = 3;
        static constexpr int QQX_MUPNU_G = 4;
        static constexpr int QG_MUPNU_Q = 5;
        // EW
        static constexpr int QXQ_MUPNU_A = 6;
        static constexpr int QQX_MUPNU_A = 7;
        // Born
        static constexpr int QXQ_MUPNU = 8;
        static constexpr int QQX_MUPNU = 9;
    };

    void Init(const FKS::Param *param, AlphaScheme scheme, double inveps);

    virtual Result Born(int process, const Phasespace::Phasespace &ps, double Q,
                        const FKS::Param *param, bool QCD, bool EW);

    virtual double Real(int process, const Phasespace::Phasespace &ps,
                        const FKS::Param *param,
                        Diagrams diagrams = Diagrams::ALL);
};

#endif

