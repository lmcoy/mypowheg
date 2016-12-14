#ifndef INTEGRALTRANSFORMATION_H_HZ7DGH89
#define INTEGRALTRANSFORMATION_H_HZ7DGH89

#include <memory>

#include "process/data.h"

class IntegralTransformation {
  public:
    virtual double Transform(int n, const double *x_in, double *x_out,
                             UserProcess::Data *userdata) = 0;
};
typedef std::shared_ptr<IntegralTransformation> IntegralTransformationPtr;

class NoTransformation {
  public:
    virtual double Transform(int n, const double *x_in, double *x_out,
                             UserProcess::Data *userdata) {
        for (int i = 0; i < n; i++) {
            x_out[i] = x_in[i];
        }
        return 1.0;
    }
};

#endif
