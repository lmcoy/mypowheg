#include <cassert>

#include "physics/alphas.h"

using namespace Physics;

AlphaSRunning::AlphaSRunning(double alphas_at_mz,
                             AlphaSRunning::AlphaSOrder order, int nf) {
    switch (order) {
    case AlphaSRunning::AlphaSOrder::LO:
    case AlphaSRunning::AlphaSOrder::NLO:
        LambdaQCD_NL(alphas_at_mz, nf);
        break;
    case AlphaSRunning::AlphaSOrder::NNLO:
        assert(0 && "not implemented");
        break;
    }
}

double AlphaSRunning::AlphaS(double scale2) const {
    if (scale2 > M_B2) {
        return alphas_2loop(scale2, invLambda2, 5);
    }
    if (scale2 >= M_C2) {
        double t1 = 1.0 / alphas_2loop(scale2, invLambda2, 4);
        double t2 = 1.0 / alphas_2loop(M_B2, invLambda2, 5);
        double t3 = 1.0 / alphas_2loop(M_B2, invLambda2, 4);
        double t = t1 + t2 - t3;
        return 1.0 / t;
    }
    double t1 = 1.0 / alphas_2loop(scale2, invLambda2, 3);
    double t2 = 1.0 / AlphaS(M_C2 + 1e-6);
    double t3 = 1.0 / alphas_2loop(M_C2, invLambda2, 3);
    double t = t1 + t2 - t3;
    return 1.0 / t;
}

void AlphaSRunning::LambdaQCD_NL(double as_at_mz, int nf) {
    constexpr double mz = 91.1876;
    assert(nf >= 3 && nf <= 5);
    double b = (33.0 - 2. * nf) / (12.0 * M_PI);
    double bp = (153. - 19. * nf) / (2 * M_PI * (33. - 2. * nf));
    double t = 1.0 / (b * as_at_mz);
    double ot = t * 2.0; // just make sure that ot != t
    while (fabs(ot - t) / ot > 1e-8) {
        ot = t;
        double logt = log(t);
        double as0 = 1.0 / (b * t) - bp * logt / (b * b * t * t);
        double as1 =
            -1.0 / (b * t * t) - bp / (b * b) * (1. - 2. * logt) / (t * t * t);
        t += (as_at_mz - as0) / as1;
    }
    invLambda2 = exp(t) / mz / mz;
}

double AlphaSRunning::alphas_2loop(double scale2, double invLambda2, int nf) {
    double t1 = log(scale2 * invLambda2);
    double denom = (33.0 - 2.0 * (double)nf);
    double alpha_1loop = 12.0 * M_PI / (denom * t1);
    return alpha_1loop *
           (1.0 - alpha_1loop / (M_PI * 2.0) * (153. - 19. * (double)nf) /
                      denom * log(t1));
}
