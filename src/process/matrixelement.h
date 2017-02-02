#ifndef PROCESS_MATRIXELEMENT_H_
#define PROCESS_MATRIXELEMENT_H_

#include <complex>
#include <memory>

#include "fks/regions.h"
#include "phasespace/phasespace.h"
#include "util/matrix.h"

namespace FKS {
class Param;
}

namespace UserProcess {

class SpinCorrelated {
  public:
    SpinCorrelated() {
        for (int i = 0; i < 16; i++) {
            data[i] = std::complex<double>(0.0, 0.0);
        }
    }
    std::complex<double> *operator[](int row) { return &data[4 * row]; }
    const std::complex<double> *operator[](int row) const {
        return &data[4 * row];
    }

    double Born() const {
        return std::real(-(*this)[0][0] + (*this)[1][1] + (*this)[2][2] +
                         (*this)[3][3]);
    }

    std::complex<double> At(int i, int j) const { return data[4 * i + j]; }

    void Set(int i, int j, std::complex<double> v) { data[4 * i + j] = v; }

  private:
    std::complex<double> data[16];
};

class IMatrixElement {
  public:
    struct Result {
        double M2 = 0.0;
        double VfinQCD = 0.0;
        double VfinEW = 0.0;
        Util::Matrix2 ColorCorr;
        SpinCorrelated SpinCorr;
    };
    enum class Diagrams { ALL, ONLYFSR, ONLYISR };
    virtual Result Born(int process, const Phasespace::Phasespace &ps, double Q,
                        const FKS::Param *param, bool QCD, bool EW) = 0;

    virtual double Real(int process, const Phasespace::Phasespace &ps,
                        const FKS::Param *param,
                        Diagrams diagrams = Diagrams::ALL) = 0;
};

using IMatrixElementPtr = std::shared_ptr<IMatrixElement>;

} // end namespace UserProcess

#endif
