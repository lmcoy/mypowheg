#ifndef PROCESS_MATRIXELEMENT_H_
#define PROCESS_MATRIXELEMENT_H_

#include <memory>

#include "fks/regions.h"
#include "phasespace/phasespace.h"
#include "util/matrix.h"

namespace FKS {
class Param;
}

namespace UserProcess {

class IMatrixElement {
  public:
    struct Result {
        double M2;
        double VfinQCD;
        double VfinEW;
        Util::Matrix2 ColorCorr;
    };
    enum class Diagrams {
        ALL,
        ONLYFSR,
        ONLYISR
    };
    virtual Result Born(int process, const Phasespace::Phasespace &ps, double Q,
                        const FKS::Param *param, bool QCD, bool EW) = 0;

    virtual double Real(int process, const Phasespace::Phasespace &ps,
                        const FKS::Param *param,
                        Diagrams diagrams = Diagrams::ALL) = 0;
};

using IMatrixElementPtr = std::shared_ptr<IMatrixElement>;

} // end namespace UserProcess

#endif

