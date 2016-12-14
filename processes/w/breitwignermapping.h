#ifndef BREITWIGNERMAPPING_H_DJ983NF8
#define BREITWIGNERMAPPING_H_DJ983NF8

#include "powheggen/integraltransformation.h"
#include "process/data.h"

#include "cuts.h"

class BreitWignerMapping : public IntegralTransformation {
  public:
    virtual double Transform(int n, const double *x_in, double *x_out,
                             UserProcess::Data *userdata);

  private:
    double clip(double x) {
        if (x >= 1.0 && x < 1.0 + 1e-9) {
            return 1.0 - 1e-12;
        }
        if (x <= 0.0 && x > -1e-9) {
            return 1e-12;
        }
        return x;
    }
};

double BreitWignerMapping::Transform(int n, const double *x_in, double *x_out,
                                     UserProcess::Data *userdata) {
    double S = userdata->SqrtS * userdata->SqrtS;
    double Mll = 1.0;
    constexpr double M = 80.385;
    constexpr double G = 2.085;
    double Delta =
        atan((S - M * M) / (M * G)) - atan((Mll * Mll - M * M) / (M * G));
    double s = x_in[0];
    double z = x_in[1];
    double ts = Delta * s + atan((Mll * Mll - M * M) / (M * G));
    double tau = clip((G * tan(ts) + M) * M / S);
    LIB_ASSERT(tau >= 0.0 && tau <= 1.0, "tau = %.16g\n", tau);
    double x1 = (1.0 - tau) * z + tau;
    x_out[0] = clip(x1);
    x_out[1] = clip(tau / x1);

    for (int i = 2; i < n; i++) {
        x_out[i] = clip(x_in[i]);
    }

    double tmp = tau * S - M * M;
    double inv_breit_wigner = (M * M * G * G + tmp * tmp) / (M * G * S);
    return Delta * (1.0 - tau) / x1 * inv_breit_wigner;
}

#endif
