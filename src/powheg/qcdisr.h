#ifndef QCDISR_H_RFP9MTVQ
#define QCDISR_H_RFP9MTVQ

#include <cmath>

#include "powheg/solverqcdisr.h"
#include "powheg/ptgenerator.h"

#include "libconfig.h"
#include "fks/param.h"
#include "process/data.h"
#include "powheg/alphas.h"
#include "powheg/isr.h"

namespace {

class LIB_LOCAL QCDISR : public ISR {
  public:
    explicit QCDISR(RadiationAlphaS as) : alphas(as) {}

    double GenPT2(double pT2min, double pT2max, double kT2tilde, double ximax,
                  double sb, double N, Random::RNG *rng) const {
        SolverQCDISR solver(pT2min, N, kT2tilde, sb, alphas);
        PTGenerator<SolverQCDISR> generator(pT2min);
        return generator.GenPT2(solver, pT2max, kT2tilde, ximax, sb, rng);
    }

    double GetUpperboundingCoupling(double kT2,
                                    UserProcess::Data *userdata) const {
        return alphas.SimpleAlphaS(kT2);
    }

    double GetCoupling(double kT2, UserProcess::Data *userdata) const {
        return alphas.AlphaSNLL(kT2);
    }

    static double PDFScale(double kT2, UserProcess::Data *userdata) {
        return sqrt(kT2);
    }

  private:
    RadiationAlphaS alphas;
};

} // end namespace

#endif /* end of include guard: QCDISR_H_RFP9MTVQ */
