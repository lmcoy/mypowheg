#ifndef PTGENERATOR_H_2TSAGG6N
#define PTGENERATOR_H_2TSAGG6N

#include <cmath>
#include <cstdio>
#include "random/rnd.h"

namespace {

template <typename Solver> class PTGenerator {
  public:
    explicit PTGenerator(double pT2min) : pT2min_(pT2min) {}
    double GenPT2(const Solver &solver, double pT2max, double kT2tilde,
                  double ximax, double sb, Random::RNG *rng) const {
        if (pT2min_ >= pT2max) {
            fprintf(stderr, "%s:%d: pT2min = %g  >= %g = pT2max\n", __FILE__,
                    __LINE__, pT2min_, pT2max);
            return 0.0;
        }
        double p2 = pT2max;
        while (p2 >= pT2min_) {
            double u = rng->Random();
            double pT2 = solver.solve(p2, u);
            if (pT2 < pT2min_) {
                return 0.0;
            }
            double v = rng->Random();
            double IoverIbar = I(pT2, ximax, sb) / Ibar(pT2, kT2tilde, sb);
            if (v < IoverIbar) {
                return pT2;
            }
            p2 = pT2;
        }
        return 0.0;
    }

  private:
    double pT2min_;

    static double I(double kT2, double ximax, double sb) {
        double xi1 = -2.0 / sb * kT2 + sqrt(4.0 * kT2 / sb * (1.0 + kT2 / sb));
        double xi2 = -2.0 / sb * kT2 - sqrt(4.0 * kT2 / sb * (1.0 + kT2 / sb));

        double arg1 = sqrt(ximax - xi1) + sqrt(ximax - xi2);
        double arg2 = xi1 - xi2;
        return 2.0 * log(arg1) - log(arg2);
    }

    static double Ibar(double kT2, double kT2tilde, double sb) {
        return 0.5 * log((kT2tilde + sb) / kT2);
    }
};

} // end namespace

#endif /* end of include guard: PTGENERATOR_H_2TSAGG6N */ 
