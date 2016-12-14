#ifndef SOLVERQCDISR_H_POFIYGQX
#define SOLVERQCDISR_H_POFIYGQX

#include <cmath>
#include <cassert>
#include <boost/math/tools/roots.hpp>

#include "libconfig.h"
#include "powheg/alphas.h"

namespace {

class LIB_LOCAL SolverQCDISR {
  public:
    SolverQCDISR(double pT2min, double N, double kT2tilde, double sb,
                 RadiationAlphaS alphas) {
        invLambda2_ = alphas.SimpleInvLambda2();
        assert(invLambda2_ * pT2min > 1.0);
        double b0 = RadiationAlphaS::b0;
        pre1_ = M_PI * N / b0;
        pre2_ = log((kT2tilde + sb) * invLambda2_);
        pT2min_ = pT2min;
    }

    double solve(double p2, double u) const {
        double pre1 = pre1_;
        double pre2 = pre2_;
        double invLambda2 = invLambda2_;
        auto eqn = [pre1, pre2, p2, u, invLambda2](double x) {
            return pre1 *
                       (pre2 * log(log(x * invLambda2) / log(p2 * invLambda2)) +
                        log(p2 / x)) -
                   log(u);
        };

        double fa = eqn(pT2min_);
        double fb = eqn(p2);
        if (fa * fb > 0.0) {
            // born
            return 0.0;
        }

        boost::uintmax_t max_iter = 100;

        // relative difference: 2^(1-15) = 6e-5
        boost::math::tools::eps_tolerance<double> tol(15);
        auto bracket = boost::math::tools::toms748_solve(eqn, pT2min_, p2, fa,
                                                         fb, tol, max_iter);
        if (max_iter == 100) {
            fprintf(stderr, "not converged: result in [%.16g, %.16g]\n",
                    bracket.first, bracket.second);
        }

        return bracket.first;
    }

  private:
    double pre1_;
    double pre2_;
    double invLambda2_;
    double pT2min_;
};

} // end namespace

#endif /* end of include guard: SOLVERQCDISR_H_POFIYGQX */
