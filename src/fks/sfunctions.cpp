#include "fks/sfunctions.h"

#include "fks/regions.h"
#include "fks/process.h"
#include "math/fourmomentum.h"
#include "phasespace/phasespace.h"
#include "physics/pdgcode.h"

//#define SF1

namespace {

static const double a = 1.0;
static const double b = 1.0;

double h(double z, double c) {
    double t = pow(1.0-z,c);
    return t/( pow(z,c) + t);
}

double Energy(const Math::FourMomentum &initial1, const Math::FourMomentum initial2,
          const Math::FourMomentum &ki, double sqrts) {
    Math::FourMomentum q = initial1.Plus(initial2);
    return 1.0 / sqrts * q.Dot(ki);
}


double di(const Math::FourMomentum &initial1, const Math::FourMomentum initial2,
          const Math::FourMomentum &ki, double sqrts) {
    double Ei = Energy(initial1, initial2, ki, sqrts);

    double costhi = 1.0 - 2.0*ki.Dot(initial1)/Ei/sqrts;
    double pre = (1.0 - costhi) * (1.0 + costhi);

#ifdef SF1
    double ret = sqrts * 0.5 * Ei * pre;
#else
    double ret = Ei * Ei * pre;
#endif
    
    if (ret < 0.0) {
        // ret < 0.0 means numerical errors in costhi calculation
        return 0.0;
    }
    return ret;
}

double dij(const Math::FourMomentum &initial1,
           const Math::FourMomentum initial2, const Math::FourMomentum &ki,
           const Math::FourMomentum &kj, double sqrts) {
    double Ei = Energy(initial1, initial2, ki, sqrts);
    double Ej = Energy(initial1, initial2, kj, sqrts);
    double tmp = ki.Dot(kj)/(Ei*Ej);

    // assert that cos = 1 - tmp < 1 + 1e-10
    assert(tmp > -1e-10);

#ifdef SF1
    double ret = Ei * Ej * tmp;
#else
    Math::FourMomentum ki_dir = ki;
    ki_dir.Scale(1.0/Ei);
    double ct = ki_dir.Dot(kj);
    double ret = 2.0 * Ei * Ei * Ej / (Ei + Ej) / (Ei + Ej) * ct;
#endif

    if (ret < 0.0) {
        return 0.0;
    }
    return ret;
}

double Charge(const FKS::Region &region, const FKS::Real_t &real) {
    double Q = Physics::PDG::Charge(real.Flavours[region.J]);
    if (Q != 0.0) {
        return Q;
    }
    return Physics::PDG::Charge(real.Flavours[region.I]);
}

} // namespace

namespace FKS {

double SFunction(const Phasespace::Phasespace &ps, const Region &region,
                 const Real_t &real, int resonance) {
    bool use_charge = real.Type == Type_t::EW && resonance > 1;
    bool UseResonances = resonance > 0;
    double sqrts = sqrt(ps.X1 * ps.X2 * ps.S);
    double d_pre = 0.0;
    int j = region.J;
    int i = region.I;
    if (j < 2) {
        d_pre = di(ps.Momenta[0], ps.Momenta[1], ps.Momenta[i], sqrts);
    } else {
        d_pre = dij(ps.Momenta[0], ps.Momenta[1], ps.Momenta[i], ps.Momenta[j],
                    sqrts);
    }

    double P_pre = 1.0;
    if (UseResonances) {
        P_pre = real.Resonances[region.ResonanceID].BreitWigner(ps, use_charge);
    }

    assert(P_pre > 0.0);
    double denom = 0.0;

    size_t lenf = real.Resonances.size();
    for (size_t f = 0; f < lenf; f++) {
        if (f > 0 && !UseResonances) {
            break;
        }
        double P = 1.0;
        if (f > 0 && UseResonances) {
            P = real.Resonances[f].BreitWigner(ps, use_charge);
        }
        for (const auto &reg : real.Regions) {
            if (!UseResonances || reg.ResonanceID == f) {
                double d = 0.0;
                int k = reg.I;
                int l = reg.J;
                if (reg.J < 2) {
                    d = di(ps.Momenta[0], ps.Momenta[1], ps.Momenta[k], sqrts);
                } else {
                    d = dij(ps.Momenta[0], ps.Momenta[1], ps.Momenta[k],
                            ps.Momenta[l], sqrts);
                }
                if(d_pre == 0.0 && d == 0.0) {
                    if(i == k && l == j) {
                        return 1.0;
                    } else {
                        assert(0 && "FKS parton is coll to initial state and "
                                    "coll to a final state "
                                    "=> 2 regions overlap. I don't now what to "
                                    "do with those "
                                    "points."
                                    "=> introduce cuts to remove this region "
                                    "(it's not physical "
                                    "anyway)");
                    }
                }
                assert(d >= 0.0);
                if (d == 0.0) {
                    return 0.0;
                }

                denom += d_pre / d * P / P_pre;
            }
        }
    }

    assert(denom > 0.0);
    return 1.0 / denom;
}

} // namespace FKS

