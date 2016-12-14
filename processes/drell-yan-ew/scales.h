#ifndef SCALES_H_
#define SCALES_H_

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

class InvMassScales : public UserProcess::IScales {
  public:
    virtual double Factorization(const Phasespace::Phasespace &ps) const {
        return invMassLeptons(ps);
    }

    virtual double Renorm(const Phasespace::Phasespace &ps) const {
        return invMassLeptons(ps);
    }

  private:
    double invMassLeptons(const Phasespace::Phasespace &ps) const {
        auto pll = ps.Momenta[2].Plus(ps.Momenta[3]);
        double mll2 = pll.Dot(pll);
        return sqrt(mll2);
    }
};

#endif
