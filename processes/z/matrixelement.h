#ifndef MATRIXELEMENT_H_K3XX23LA4
#define MATRIXELEMENT_H_K3XX23LA4

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

class DrellYanME : public UserProcess::IMatrixElement {
  public:

    struct SubProcesses {
        // QCD
        static constexpr int UXU_MUXMU_G = 0;
        static constexpr int UXG_MUXMU_UX = 1;
        static constexpr int GU_MUXMU_U = 2;
        static constexpr int UUX_MUXMU_G = 3;
        static constexpr int UG_MUXMU_U = 4;
        static constexpr int GUX_MUXMU_UX = 5;
        static constexpr int DXD_MUXMU_G = 6;
        static constexpr int DXG_MUXMU_DX = 7;
        static constexpr int GD_MUXMU_D = 8;
        static constexpr int DDX_MUXMU_G = 9;
        static constexpr int DG_MUXMU_D = 10;
        static constexpr int GDX_MUXMU_DX = 11;
        // EW
        static constexpr int UXU_MUXMU_A = 12;
        static constexpr int UUX_MUXMU_A = 13;
        static constexpr int DXD_MUXMU_A = 14;
        static constexpr int DDX_MUXMU_A = 15;
    };

    void Init(const FKS::Param *param, AlphaScheme scheme, double inveps);

    virtual Result Born(int process, const Phasespace::Phasespace &ps, double Q,
                        const FKS::Param *param, bool QCD, bool EW);

    virtual double Real(int process, const Phasespace::Phasespace &ps,
                        const FKS::Param *param,
                        Diagrams diagrams = Diagrams::ALL);

  private:
    double counterterm = 0.0;
    bool initialized = false;
};

#endif

