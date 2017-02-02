#ifndef SCALES_H_JD8SM27S
#define SCALES_H_JD8SM27S

#include "math/fourmomentum.h"
#include "phasespace/phasespace.h"

class FixedScales : public UserProcess::IScales {
  public:
    FixedScales(double mu, double muF) : mu_(mu), muF_(muF) {
    }
    virtual double Factorization(const Phasespace::Phasespace &ps) const {
        return mu_;
    }
    virtual double Renorm(const Phasespace::Phasespace &ps) const {
        return muF_;
    }

  private:
    double mu_;
    double muF_;
};

class PTWScales : public UserProcess::IScales {
  public:
    virtual double Factorization(const Phasespace::Phasespace &ps) const {
        return ps.Momenta[4].PT();
    }
    virtual double Renorm(const Phasespace::Phasespace &ps) const {
        return ps.Momenta[4].PT();
    }

};

#endif
