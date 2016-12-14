#ifndef PROCESS_MATRIXELEMENT_H_
#define PROCESS_MATRIXELEMENT_H_

#include "fks/regions.h"
#include "phasespace/phasespace.h"
#include "util/matrix.h"

namespace FKS {
class Param;
class Param_as;
}

enum class CorrectionType {
    None,
    QCD,
    EW
};

enum class AlphaScheme {
    Alpha0,
    Gmu
};

enum class Diagrams {
    ALL,
    ONLYFSR,
    ONLYISR
};

struct BornMEOut {
    double M2;
    double VfinQCD;
    double VfinEW;
    Util::Matrix2 ColorCorr;
};

void InitBornME(const FKS::Param *param, const FKS::Param_as *param_as,
                AlphaScheme, double mu2, double inveps, double *counterterm);

int BornME(int process, const Phasespace::Phasespace &ps, double Q,
           const FKS::Param *param, const FKS::Param_as *param_aS,
           BornMEOut *me, bool QCD, bool EW, double counterterm);

double RealME(int process, const Phasespace::Phasespace &ps,
              const FKS::Param *param, const FKS::Param_as *param_aS,
              Diagrams diagrams = Diagrams::ALL);

#endif

