#ifndef SCALES_H_SEFKO3S9
#define SCALES_H_SEFKO3S9

#include "phasespace/phasespace.h"
#include "math/fourmomentum.h"

class FixedScales : public UserProcess::IScales {
  public:
    FixedScales(double mu, double muF) : mu_(mu), muF_(muF) {}
    virtual double Factorization(const Phasespace::Phasespace &ps) const {
        return mu_;
    }
    virtual double Renorm(const Phasespace::Phasespace &ps) const { return muF_; }

  private:
    double mu_;
    double muF_;
};

#endif
