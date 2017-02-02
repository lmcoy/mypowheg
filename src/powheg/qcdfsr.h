#ifndef QCDFSR_H_8QZSWIML
#define QCDFSR_H_8QZSWIML

#include <cmath>
#include <cassert>

#include <boost/math/tools/roots.hpp>

#include "libconfig.h"
#include "process/data.h"
#include "fks/param.h"
#include "powheg/alphas.h"
#include "phasespace/phasespace.h"
#include "powheg/fsr.h"

namespace {

class LIB_LOCAL QCDFSR : public FSR {
  public:
    explicit QCDFSR(RadiationAlphaS as) : alphas(as) {}

    double GenPT2(double pT2min, double pT2max, double kT2tilde, double ximax,
                  double sb, double N, Random::RNG *rng) const {
        double pre1 = M_PI * N / RadiationAlphaS::b0;
        double invLambda2 = alphas.SimpleInvLambda2();
        double pre2 = log(ximax * ximax * sb * invLambda2);
        double u = rng->Random();

        auto eqn = [pre1, pre2, invLambda2, u, pT2max](double x) {
            return pre1 * (pre2 * log(log(x * invLambda2) /
                                      log(pT2max * invLambda2)) +
                           log(pT2max / x)) -
                   log(u);
        };
        if (pT2min >= pT2max) {
            return 0.0;
        }

        double fa = eqn(pT2min);
        double fb = eqn(pT2max);
        if (fa * fb > 0.0) {
            return 0.0;
        }
        boost::uintmax_t max_iter = 100;

        // relative difference: 2^(1-15) = 6e-5
        boost::math::tools::eps_tolerance<double> tol(15);
        auto bracket = boost::math::tools::toms748_solve(eqn, pT2min, pT2max,
                                                         fa, fb, tol, max_iter);
        if (max_iter == 100) {
            fprintf(stderr, "not converged: result in [%.16g, %.16g]\n",
                    bracket.first, bracket.second);
        }

        return bracket.first;
    }

    double GetUpperboundingCoupling(double kT2,
                                    UserProcess::Data *userdata) const {
        return alphas.SimpleAlphaS(kT2);
    }

    double GetCoupling(double kT2, UserProcess::Data *userdata) const {
        return alphas.AlphaSNLL(kT2);
    }

    static double PDFScale(double kT2, UserProcess::Data *userdata) {
        return 50.0; // pdfs cancel for FSR, use arbitrary pdf scale
    }

  private:
    RadiationAlphaS alphas;
};

} // end namespace

#endif /* end of include guard: QCDFSR_H_8QZSWIML */
