#ifndef FSR_H_6QHYAOSZ
#define FSR_H_6QHYAOSZ

#include "libconfig.h"
#include "phasespace/phasespace.h"
#include "random/rnd.h"

namespace {

class LIB_LOCAL FSR {
  public:
    static double XiMax(const Phasespace::Phasespace &ps_born, int j) {
        assert(j >= 2);
        double s = ps_born.X1 * ps_born.X2 * ps_born.S;
        double Mrec2 = s - 2.0 * ps_born.Momenta[j].E() * sqrt(s);
        return 1 - Mrec2 / s;
    }

    static double PT2Max(const Phasespace::Phasespace &ps_born, double ximax, double s) {
        return ximax * ximax * s;
    }

    static double UpperBounding(double N, double xi, double y, double alpha) {
        return N / (xi * (1.0 - y)) * alpha;
    }

    static bool GenRadiationVariables(double pT2, double s, double xi_max,
                                      Random::RNG *rng, double *xi, double *y,
                                      double *phi) {
        double u = rng->Random();
        double xi_ = exp(u * log(xi_max) + (u - 1.0) * 0.5 * log(s / pT2));
        if(!(xi_ >= 0.0 && xi_ <= 1.0)) {
            return false;
        }
        double y_ = 1.0 - 2.0 / (xi_ * xi_) * pT2 / s;
        if( !(y_ >= -1.0 && y_ <= 1.0) ) {
            return false;
        }
        *phi = rng->Random() * 2.0 * M_PI;
        *xi = xi_;
        *y = y_;
        return true;
    }
};

} // end namespace

#endif /* end of include guard: FSR_H_6QHYAOSZ */
