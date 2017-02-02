#include "fks/sfunctions.h"

#include <iostream>

#include "fks/process.h"
#include "fks/regions.h"
#include "math/fourmomentum.h"
#include "phasespace/phasespace.h"
#include "physics/pdgcode.h"

//#define SF1

namespace {

double Energy(const Math::FourMomentum &initial1,
              const Math::FourMomentum initial2, const Math::FourMomentum &ki,
              double sqrts) {
    Math::FourMomentum q = initial1.Plus(initial2);
    return 1.0 / sqrts * q.Dot(ki);
}

double dip(const Math::FourMomentum &initial1,
           const Math::FourMomentum initial2, const Math::FourMomentum &ki,
           double sqrts, int j) {
    double Ei = Energy(initial1, initial2, ki, sqrts);
    double costhi = 1.0 - 2.0 * ki.Dot(initial1) / Ei / sqrts;

    double sign = 1.0;
    if (j == -2) {
        sign = -1.0;
    }
#ifdef SF1
    double d = 0.0;
    assert(0 && "not implemented");
#else
    double d = 2.0 * Ei * Ei * (1.0 - sign * costhi);
#endif
    if (d < 0.0) {
        return 0.0;
    }
    return d;
}

double d0(const Math::FourMomentum &initial1, const Math::FourMomentum initial2,
          const Math::FourMomentum &ki, double sqrts) {
    double Ei = Energy(initial1, initial2, ki, sqrts);

    double costhi = 1.0 - 2.0 * ki.Dot(initial1) / Ei / sqrts;
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

double di(const Math::FourMomentum &initial1, const Math::FourMomentum initial2,
          const Math::FourMomentum &ki, double sqrts, int j) {
    assert(j < 2);
    assert(j >= -2);
    if (j < 0) {
        return dip(initial1, initial2, ki, sqrts, j);
    }
    return d0(initial1, initial2, ki, sqrts);
}

double dij(const Math::FourMomentum &initial1,
           const Math::FourMomentum initial2, const Math::FourMomentum &ki,
           const Math::FourMomentum &kj, double sqrts) {
    double Ei = Energy(initial1, initial2, ki, sqrts);
    double Ej = Energy(initial1, initial2, kj, sqrts);
//    double tmp = ki.Dot(kj) / (Ei * Ej);

// assert that cos = 1 - tmp < 1 + 1e-10
// assert(tmp > -1e-10);

#ifdef SF1
    double ret = Ei * Ej * tmp;
#else
    double ct = ki.Dot(kj);
    double ret = 2.0 * Ei * Ej / (Ei + Ej) / (Ei + Ej) * ct;
#endif

    if (ret < 0.0) {
        return 0.0;
    }
    return ret;
}

double collinearityISR(const Math::FourMomentum &initial1,
                       const Math::FourMomentum initial2,
                       const Math::FourMomentum &ki, double sqrts) {
    double Ei = Energy(initial1, initial2, ki, sqrts);

    return fabs(1.0 - 2.0 * ki.Dot(initial1) / Ei / sqrts);
}

double collinearityFSR(const Math::FourMomentum &initial1,
                       const Math::FourMomentum initial2,
                       const Math::FourMomentum &ki,
                       const Math::FourMomentum &kj, double sqrts) {
    double Ei = Energy(initial1, initial2, ki, sqrts);
    double Ej = Energy(initial1, initial2, kj, sqrts);
    return 1.0 - ki.Dot(kj) / (Ei * Ej);
}

double h(int i, int j, int pdg_i, int pdg_j, const Phasespace::Phasespace &ps,
         double sqrts) {
    if (j < 2) {
        return 1.0;
    }
    if (pdg_i != 0 && pdg_i != 21) {
        return 1.0;
    }
    if (pdg_j != 0 && pdg_j != 21) {
        return 1.0;
    }
    double Ei = Energy(ps.Momenta[0], ps.Momenta[1], ps.Momenta[i], sqrts);
    double Ej = Energy(ps.Momenta[0], ps.Momenta[1], ps.Momenta[j], sqrts);
    return Ej / (Ej + Ei); // actually: Ej^c/(Ej^c+Ei^c) with c = 1
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
    double coll_pre = 0.0;
    if (j < 2) {
        d_pre = di(ps.Momenta[0], ps.Momenta[1], ps.Momenta[i], sqrts, j);
        coll_pre =
            collinearityISR(ps.Momenta[0], ps.Momenta[1], ps.Momenta[i], sqrts);
    } else {
        d_pre = dij(ps.Momenta[0], ps.Momenta[1], ps.Momenta[i], ps.Momenta[j],
                    sqrts);
        coll_pre = collinearityFSR(ps.Momenta[0], ps.Momenta[1], ps.Momenta[i],
                                   ps.Momenta[j], sqrts);
    }
    int pdgj = 0;
    if (j >= 0) {
        pdgj = real.Flavours[j];
    }
    double h_pre = h(i, j, real.Flavours[i], pdgj, ps, sqrts);
    if (h_pre == 0.0) {
        return 0.0;
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
        for (const auto &reg : real.AllRegions) {
            if (!UseResonances || reg.ResonanceID == f) {
                double d = 0.0;
                int k = reg.I;
                int l = reg.J;
                if (i == k && j == l) {
                    denom += P / P_pre;
                    continue;
                }
                double coll = 0.0;
                if (reg.J < 2) {
                    d = di(ps.Momenta[0], ps.Momenta[1], ps.Momenta[k], sqrts,
                           reg.J);
                    coll = collinearityISR(ps.Momenta[0], ps.Momenta[1],
                                           ps.Momenta[k], sqrts);
                } else {
                    d = dij(ps.Momenta[0], ps.Momenta[1], ps.Momenta[k],
                            ps.Momenta[l], sqrts);
                    coll = collinearityFSR(ps.Momenta[0], ps.Momenta[1],
                                           ps.Momenta[k], ps.Momenta[l], sqrts);
                }
                int pdgl = (l >= 0) ? real.Flavours[l] : 0;
                double h_kl = h(k, l, real.Flavours[k], pdgl, ps, sqrts);
                if (i == k && j == l) {
                    denom += P / P_pre;
                    continue;
                }
                if (i == l && j == k) {
                    denom += P / P_pre * h_kl / h_pre;
                    continue;
                }
                assert(d >= 0.0);
                if (1.0 - coll_pre < 1e-8 && 1.0 - coll < 1e-8) {
                    if (i != k || j != l) {
                        std::cout << "region overlap d = " << d << ", ";
                        std::cout << "d_pre = " << d_pre << ", ";
                        std::cout << "i = " << i << ", j = " << j << ", ";
                        std::cout << "k = " << k << ", l = " << l << "\n";
                        std::cout << ps.Momenta[0].ToString() << "\n";
                        std::cout << ps.Momenta[1].ToString() << "\n";
                        std::cout << ps.Momenta[2].ToString() << "\n";
                        std::cout << ps.Momenta[3].ToString() << "\n";
                        std::cout << ps.Momenta[4].ToString() << "\n";
                        std::cout << ps.Momenta[5].ToString() << "\n";
                    }
                }
                if (d == 0.0) {
                    return 0.0;
                }

                denom += d_pre / d * P / P_pre * h_kl / h_pre;
            }
        }
    }

    assert(denom > 0.0);
    return 1.0 / denom;
}

} // namespace FKS
