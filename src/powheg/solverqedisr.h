#ifndef SOLVERQEDISR_H_Z0K8KBI9
#define SOLVERQEDISR_H_Z0K8KBI9

#include <cmath>
#include <cassert>
#include <boost/math/tools/roots.hpp>

#include "libconfig.h"

namespace {

class LIB_LOCAL SolverQEDISR {
  public:
    SolverQEDISR(double pT2min, double N, double kT2tilde, double sb,
                 double alpha) {
        pre1_ = kT2tilde + sb;
        pre2_ = 2.0 / (M_PI * N * alpha);
        pT2min_ = pT2min;
    }

    double solve(double p2, double u) const {
        double tmp1 = log(p2 / pre1_);
        double tmp2 = pre2_ * log(u);
        double Exp = exp(-sqrt(tmp1 * tmp1 - tmp2));
        double pT2 = pre1_ * Exp;
        if (pT2 < pT2min_) {
            return 0.0;
        }
        return pT2;
    }

  private:
    double pre1_;
    double pre2_;
    double pT2min_;
};

} // end namespace

#endif /* end of include guard: SOLVERQEDISR_H_Z0K8KBI9 */ 
