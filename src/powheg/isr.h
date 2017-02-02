#ifndef ISR_H_DNHKLBQW
#define ISR_H_DNHKLBQW

#include "libconfig.h"
#include "phasespace/phasespace.h"
#include "random/rnd.h"

namespace {

class LIB_LOCAL ISR {
  public:
    static double XiMax(const Phasespace::Phasespace &ps_born, int j) {
        assert(j < 2);
        return 1.0 - ps_born.X1 * ps_born.X2;
    }

    static double PT2Max(const Phasespace::Phasespace &ps_born, double ximax, double sb) {
        double x1 = ps_born.X1;
        double x2 = ps_born.X2;
        double f = (1.0 - x1 * x1) * (1.0 - x2 * x2) / (x1 + x2) / (x1 + x2);
        return sb * f;
    }

    static double UpperBounding(double N, double xi, double y, double alpha) {
        return N / (xi * (1.0 - y * y)) * alpha;
    }

    static bool GenRadiationVariables(double pT2, double sb, double xi_max,
                                      Random::RNG *rng, double *xi, double *y,
                                      double *phi) {
        double xi_ = gen_xi(pT2, sb, xi_max, rng->Random());
        if (!(xi_ >= 0.0 && xi_ <= 1.0)) {
            return false;
        }
        double y_ = sqrt(1.0 - 4.0 * (1.0 - xi_) / (xi_ * xi_) * pT2 / sb);
        if (!(y_ >= 0.0 && y_ <= 1.0)) {
            return false;
        }
        *phi = rng->Random() * 2.0 * M_PI;
        if (rng->Random() < 0.5) {
            y_ = -y_;
        }
        *xi = xi_;
        *y = y_;
        return true;
    }

  private:
    static double gen_xi(double pT2, double sb, double xi_max, double u) {
        double t1 = pT2 / sb;
        auto f = [t1](double xi) {
            return sqrt(4.0 * t1 * (xi - 1.0) + xi * xi) + 2.0 * t1 + xi;
        };
        double f_max = f(xi_max);
        double f_min = 2.0 * sqrt(t1 * (1.0 + t1));
        double F = exp(u * log(f_max / f_min)) * f_min;
        return (F * F - 4.0 * t1 * F + 4.0 * t1 * t1 + 4.0 * t1) / (2.0 * F);
    }
};

} // end namespace

#endif /* end of include guard: ISR_H_DNHKLBQW */
