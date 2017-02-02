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

    constexpr double M = 80.385;
    constexpr double G = 2.085;
    const double atan_M_over_G = atan(M / G);
   
    for (int i = 0; i < n; i++) {
        x_out[i] = x_in[i];
    }
    // double x1 = x_in[0];
    // double x2 = x_in[1];
    // double t = x_in[2];
    // double s = x1 * x2 * S;
    // double Deltap = atan_M_over_G - atan((M * M - s) / G / M);
    // double taup = Deltap * t - atan_M_over_G;
    // double z = (G * M * tan(taup) + M * M) / s;
    // double costaup = cos(taup);
    // double inv_breit_wigner = M * G / S / (costaup * costaup);
    // double jac = Deltap / (x1 * x2) * inv_breit_wigner;
    // x_out[2] = clip(z);
    
    constexpr double a = 1.0;
    double tau = a * x_in[0];
    double x1p = x_in[1];
    double t = x_in[2];
    double s = tau * S;
    double Deltap = atan_M_over_G - atan((M * M - s) / G / M);
    double taup = Deltap * t - atan_M_over_G;
    double z = (G * M * tan(taup) + M * M) / s;
    double costaup = cos(taup);
    double inv_breit_wigner = M * G / S / (costaup * costaup);
    double jac = Deltap / tau * inv_breit_wigner;
    double x1 = (1.0 - tau) * x1p + tau;
    jac *= (1.0 - tau) / x1 * a;
    x_out[0] = clip(x1);
    x_out[1] = clip(tau / x1);
    x_out[2] = clip(z);
    return jac;
}

#endif
