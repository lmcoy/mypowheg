#ifndef POWHEG_QEDFSR_H_
#define POWHEG_QEDFSR_H_

#include <cmath>

#include "libconfig.h"
#include "powheg/ptgenerator.h"
#include "powheg/alphas.h"
#include "powheg/fsr.h"

namespace {

class LIB_LOCAL QEDFSR : public FSR {
  public:
    static constexpr double alpha = 1.0 / 137.03599;
    double GenPT2(double pT2min, double pT2max, double kT2tilde, double ximax,
                  double sb, double N, Random::RNG *rng) const {
        double u = rng->Random();
        double tmp = log((ximax * ximax * sb) / pT2max);
        double arg =
            -sqrt(tmp * tmp - 2.0 / (N * M_PI * alpha) * log(u));
        double Exp = exp(arg);
        double pT2 = ximax * ximax * sb * Exp;
        assert(pT2 <= pT2max);
        if (pT2 < pT2min) {
            return 0.0;
        }
        return pT2;
    }

    double GetUpperboundingCoupling(double kT2,
                                    UserProcess::Data *userdata) const {
        return alpha;
    }

    double GetCoupling(double kT2, UserProcess::Data *userdata) const {
        return alpha;
    }

    static double PDFScale(double kT2, UserProcess::Data *userdata) {
        return 50.0; // pdfs cancel for FSR, use arbitrary pdf scale
    }
};

} // end namespace

#endif /* end of include guard: POWHEG_QEDFSR_H_ */
