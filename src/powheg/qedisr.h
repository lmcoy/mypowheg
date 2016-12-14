#ifndef POWHEG_QEDISR_H_
#define POWHEG_QEDISR_H_

#include <cmath>

#include "libconfig.h"
#include "powheg/ptgenerator.h"
#include "powheg/solverqedisr.h"
#include "powheg/alphas.h"
#include "powheg/isr.h"

namespace {

class LIB_LOCAL QEDISR : public ISR {
  public:
    static constexpr double alpha = 1.0 / 137.03599;
    double GenPT2(double pT2min, double pT2max, double kT2tilde, double ximax,
                  double sb, double N, Random::RNG *rng) const {
        SolverQEDISR solver(pT2min, N, kT2tilde, sb, alpha);
        PTGenerator<SolverQEDISR> generator(pT2min);
        return generator.GenPT2(solver, pT2max, kT2tilde, ximax, sb, rng);
    }

    double GetUpperboundingCoupling(double kT2,
                                    UserProcess::Data *userdata) const {
        return alpha;
    }

    double GetCoupling(double kT2, UserProcess::Data *userdata) const {
        return alpha;
    }

    static double PDFScale(double kT2, UserProcess::Data *userdata) {
        return sqrt(kT2);
    }
};

} // end namespace

#endif /* end of include guard: POWHEG_QEDISR_H_ */
