#include "fks/virtual.h"
#include "fks/scales.h"
#include "math/fourmomentum.h"
#include "math/dilog.h"
#include "math/math.h"
#include "physics/pdgcode.h"
#include "util/matrix.h"

namespace {

enum class Type {
    QED,
    QCD
};

template <Type type> bool particle_contributes(int pdg);

template <> bool particle_contributes<Type::QED>(int pdg) {
    int apdg = abs(pdg);
    assert(apdg != 22 && "photon in born inital state not implemented");
    static const bool lookup[] = { false, // no particle
                                   true,  // d
                                   true,  // u
                                   true,  // s
                                   true,  // c
                                   true,  // b
                                   false, // t
                                   false, // b' not implemented
                                   false, // t' not implemented
                                   false, // no particle
                                   false, // no particle
                                   true,  // electron
                                   false, // nu_e
                                   true,  // mu
                                   false, // nu_mu
                                   false, // tau
                                   false  // nu_tau
    };
    if (apdg < 17) {
        return lookup[apdg];
    }
    return false;
}

template <> bool particle_contributes<Type::QCD>(int pdg) {
    int apdg = abs(pdg);
    assert(apdg != 0 && apdg != 21 &&
           "gluon in born inital state not implemented");
    if (apdg > 0 && apdg < 6) {
        return true;
    }
    return false;
}

const double gamma[] = { 0.0,       // no particle
                         1.0 / 6.0, // d
                         2.0 / 3.0, // u
                         1.0 / 6.0, // s
                         2.0 / 3.0, // c
                         1.0 / 6.0, // b
                         2.0 / 3.0, // t
                         0.0,       // b' not implemented
                         0.0,       // t' not implemented
                         0.0,       // no particle
                         0.0,       // no particle
                         1.5,       // electron
                         0.0,       // nu_e
                         1.5,       // mu
                         0.0,       // nu_mu
                         1.5,       // tau
                         0.0        // nu_tau
};

static double CF = 4.0 / 3.0;
static double nf = 5.0;
static double TF = 0.5;
static double CA = 3.0;

template <Type type> double gamma_i(int pdg);

template <> double gamma_i<Type::QCD>(int pdg) {
    int apdg = abs(pdg);
    if (apdg > 0 && apdg < 6) {
        return (3.0 / 2.0) * CF;
    }
    if (apdg == 0 || apdg == 21) {
        return (11.0 / 6.0) * CA - (2.0 / 3.0) * nf * TF;
    }
    return 0;
}

template <> double gamma_i<Type::QED>(int pdg) {
    int apdg = abs(pdg);
    assert(apdg != 24 && "radiation from W (not impl.)");
    assert(apdg != 22 && "radiation from photon (not impl.)");
    if (apdg < 16) {
        return gamma[apdg];
    }
    return 0.0;
}

template <Type type> double gammap_i(int);

template <> double gammap_i<Type::QCD>(int pdg) {
    static const double gamma_prime_gluon =
        ((67.0 / 9.0) - (2.0 / 3.0) * Math::Pi * Math::Pi) * CA -
        (23.0 / 9.0) * TF * nf;
    static const double gamma_prime_quark =
        (13.0 / 2.0 - 2.0 * Math::Pi * Math::Pi / 3.0) * CF;
    if (pdg == 0 || pdg == 21) {
        return gamma_prime_gluon;
    }
    if (pdg > -6 && pdg < 6) {
        return gamma_prime_quark;
    }

    return 0.0;
}

template <> double gammap_i<Type::QED>(int pdg) {
    static double pre = 6.5 - 2.0 / 3.0 * M_PI * M_PI;
    static double gammap[] = { 0.0,             // no particle
                               pre / 9.0,       // d
                               pre * 4.0 / 9.0, // u
                               pre / 9.0,       // s
                               pre * 4.0 / 9.0, // c
                               pre / 9.0,       // b
                               pre * 2.0 / 3.0, // t
                               0.0,             // b' not implemented
                               0.0,             // t' not implemented
                               0.0,             // no particle
                               0.0,             // no particle
                               pre,             // electron
                               0.0,             // nu_e
                               pre,             // mu
                               0.0,             // nu_mu
                               pre,             // tau
                               0.0              // nu_tau
    };
    assert(pdg != 24 && "radiation from W (not impl.)");
    assert(pdg != 22 && "radiation from photon (not impl.)");
    int apdg = abs(pdg);
    if (apdg < 16) {
        return gammap[apdg];
    }
    return 0.0;
}

template <Type type> double gamma0_i(int);

template <> double gamma0_i<Type::QCD>(int pdg) {
    if (pdg == 0 || pdg == 21) {
        return 2.0 * CA;
    }
    if (pdg > -6 && pdg < 6) {
        return 2.0 * CF;
    }
    return 0.0;
}

template <> double gamma0_i<Type::QED>(int pdg) {
    static double gamma0[] = { 0.0,       // no particle
                               2.0 / 9.0, // d
                               8.0 / 9.0, // u
                               2.0 / 9.0, // s
                               8.0 / 9.0, // c
                               2.0 / 9.0, // b
                               8.0 / 9.0, // t
                               0.0,       // b' not implemented
                               0.0,       // t' not implemented
                               0.0,       // no particle
                               0.0,       // no particle
                               2.0,       // electron
                               0.0,       // nu_e
                               2.0,       // mu
                               0.0,       // nu_mu
                               2.0,       // tau
                               0.0        // nu_tau
    };
    assert(pdg != 24 && "radiation from W (not impl.)");
    assert(pdg != 22 && "radiation from photon (not impl.)");
    int apdg = abs(pdg);
    if (apdg < 16) {
        return gamma0[apdg];
    }
    return 0.0;
}

template <Type type>
double Qfin(int n, const int bornpdgs[], const Math::FourMomentum *momenta,
            double sqrts, const FKS::Scales &scales) {
    double res = 0.0;
    double tmp1 = 2.0 / sqrts;
    double tmp2 = log(sqrts * sqrts / scales.Q2);
    Math::FourMomentum in(momenta[0].Plus(momenta[1]));
    // final state
    for (int i = 2; i < n; i++) {
        if (!particle_contributes<type>(bornpdgs[i])) {
            continue;
        }
        int pdg = bornpdgs[i];
        double E = in.Dot(momenta[i]) / sqrts;
        double L = log(E * tmp1);
        double g = gamma_i<type>(pdg);
        double gp = gammap_i<type>(pdg);
        double g0 = gamma0_i<type>(pdg);
        res += gp - tmp2 * (g - g0 * L) + g0 * L * L - 2.0 * g * L;
    }

    // initial state
    double tmp3 = log(scales.muF * scales.muF / scales.Q2);
    int pdg0 = bornpdgs[0];
    int pdg1 = bornpdgs[1];
    res -= tmp3 * (gamma_i<type>(pdg0) + gamma_i<type>(pdg1));

    return res;
}

double sigma_f(int i, int pdg) {
    // incoming anti fermion
    if (i < 2 && pdg < 0) {
        return -1.0;
    }
    // outgoing fermion
    if (i >= 2 && pdg > 0) {
        return -1.0;
    }
    return 1.0;
}

/**
 * k_i, k_j must be in the partonic center of mass frame!
 */
double FiniteI(const Math::FourMomentum &k_i, const Math::FourMomentum &k_j,
               double log_s_over_Q2) {
    double t0 = log_s_over_Q2;
    double t1 = k_i.Dot(k_j) / (2.0 * k_i.E() * k_j.E());
    if (t1 > 1.0 && 1.0 - t1 < 1e-15) {
        t1 = 1.0;
    }
    assert(t1 > 0.0 && "negative argument of log");
    double log_t1 = log(t1);
    double t2 = 0.0;
    if (t1 > 6e-17 && t1 < 1.0 - 1.0e-15) {
        t2 = log(1.0 - t1) * log_t1;
    }
    return 0.5 * t0 * t0 + t0 * log_t1 - Math::Dilog(t1) +
           0.5 * log_t1 * log_t1 - t2;
}

template <Type type>
double VirtualI(int n, int bornpdgs[], const Math::FourMomentum *momenta,
                const Util::Matrix2 &, double sqrts, double Q2);
template <>
double VirtualI<Type::QED>(int n, int bornpdgs[],
                           const Math::FourMomentum *momenta,
                           const Util::Matrix2 &born, double sqrts, double Q2) {
    assert(born.GetLen() == 1);
    // check that momenta are in partonic cms
    assert(momenta[0].Plus(momenta[1]).Equals(
        Math::FourMomentum(sqrts, 0.0, 0.0, 0.0)));
    double ret = 0.0;
    double Log = log(sqrts * sqrts / Q2);
#ifndef NDEBUG
    double charge = 0.0;
#endif
    for (int i = 0; i < n; i++) {
        if (bornpdgs[i] == 22) {
            assert(0 && "photon in born not implemented");
        }
        if (!particle_contributes<Type::QED>(bornpdgs[i])) {
            continue;
        }
        double sigma_i = sigma_f(i, bornpdgs[i]);
        double charge_i = Physics::PDG::Charge(abs(bornpdgs[i]));
#ifndef NDEBUG
        charge += charge_i * sigma_i;
#endif
        for (int j = i + 1; j < n; j++) {
            if (bornpdgs[j] == 22) {
                assert(0 && "photon in born not implemented");
            }
            if (!particle_contributes<Type::QED>(bornpdgs[j])) {
                continue;
            }
            double sigma_j = sigma_f(j, bornpdgs[j]);
            double charge_j = Physics::PDG::Charge(abs(bornpdgs[j]));
            double I_ij = FiniteI(momenta[i], momenta[j], Log);
            ret -= I_ij * sigma_i * sigma_j * charge_i * charge_j;
        }
    }
    assert(fabs(charge) < 1e-15 && "charge conservation violated");

    // factor 2 is there because the sum is over i,j, i != j and the for
    // loops are only over j > i.
    return 2.0 * ret * born.Get(0, 0);
}

template <>
double VirtualI<Type::QCD>(int n, int bornpdgs[],
                           const Math::FourMomentum *momenta,
                           const Util::Matrix2 &ColorCorrelatedBorn,
                           double sqrts, double Q2) {
    assert(ColorCorrelatedBorn.GetLen() == n);
    // check that momenta are in partonic cms
    assert(momenta[0].Plus(momenta[1]).Equals(
        Math::FourMomentum(sqrts, 0.0, 0.0, 0.0)));
    double ret = 0.0;
    double Log = log(sqrts * sqrts / Q2);
    for (int i = 0; i < n; i++) {
        if (bornpdgs[i] == 0 || bornpdgs[i] == 21) {
            assert(0 && "gluon in born not implemented");
        }
        if (!(bornpdgs[i] > -6 && bornpdgs[i] < 6)) {
            continue;
        }
        for (int j = i + 1; j < n; j++) {
            if (bornpdgs[i] == 0 || bornpdgs[i] == 21) {
                assert(0 && "gluon in born not implemented");
            }
            if (!(bornpdgs[j] > -6 && bornpdgs[j] < 6)) {
                continue;
            }
            double I_ij = FiniteI(momenta[i], momenta[j], Log);
            ret += I_ij * ColorCorrelatedBorn.Get(i, j);
            ret += I_ij * ColorCorrelatedBorn.Get(j, i);
        }
    }
    return ret;
}

} // namespace

namespace FKS {

namespace QCD {

double Virtual(int n, int bornpdgs[], const Math::FourMomentum *momenta,
               double born, const Util::Matrix2 &ColorCorrelated, double Vfin,
               double alpha_s, double sqrts, const FKS::Scales &scales) {
    double Q = Qfin<Type::QCD>(n, bornpdgs, momenta, sqrts, scales);
    double IB = VirtualI<Type::QCD>(n, bornpdgs, momenta, ColorCorrelated,
                                    sqrts, scales.Q2);

    return alpha_s / (2.0 * Math::Pi) * (Q * born + IB + Vfin);
}

/**
 * @brief 1/eps^2 pole
 *
 * Eps2Pole returns the 1/eps^2 pole of the real correction. The correction is
 * normalized to (4pi*mu^2/Q^2)^eps * alpha_s/2pi, i.e.
 * R = (4pi*mu^2/Q^2)^eps * alpha_s/2pi ( A/eps^2 + B/eps + C + O(eps) )
 * where A is the number which is returned.
 */
double Eps2Pole(int n, int bornpdgs[], double bornme) {
    double res = 0.0;
    for (int i = 0; i < n; ++i)
    {
        double C = 0.0;
        if (Physics::PDG::IsQuark(bornpdgs[i])) {
            C = CF;
        }
        if (Physics::PDG::IsGluon(bornpdgs[i])) {
            C = CA;
        }
        res += C * bornme;
    }

    return res;
}

/**
 * @brief 1/eps pole
 *
 * EpsPole returns the 1/eps pole of the real correction. The correction is
 * normalized to (4pi*mu^2/Q^2)^eps * alpha_s/2pi, i.e.
 * R = (4pi*mu^2/Q^2)^eps * alpha_s/2pi ( A/eps^2 + B/eps + C + O(eps) )
 * where B is the number which is returned.
 */
double EpsPole(int n, int bornpdgs[], const Math::FourMomentum *momenta,
               double born, const Util::Matrix2 &ColorCorrelatedBorn,
               double Q2) {
    assert(Q2 > 0.0);
    
    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        res -= gamma_i<Type::QCD>(bornpdgs[i]) * born;
        for (int j = i+1; j < n; ++j) {
            double Cij = ColorCorrelatedBorn.Get(i, j);
            double Cji = ColorCorrelatedBorn.Get(j, i);
            double pre = log(2.0 * momenta[i].Dot(momenta[j]) / Q2);
            res += Cij * pre;
            res += Cji * pre;
        }
    }
    return -res;
}

} // namespace QCD

namespace QED {

double Virtual(int n, int bornpdgs[], const Math::FourMomentum *momenta,
               double born, const Util::Matrix2 &unused, double Vfin,
               double alpha, double sqrts, const FKS::Scales &scales) {
    double Q = Qfin<Type::QED>(n, bornpdgs, momenta, sqrts, scales);
    double IB = VirtualI<Type::QED>(n, bornpdgs, momenta,
                                    Util::Matrix2(1, born), sqrts, scales.Q2);

    return alpha / (2.0 * Math::Pi) * (Q * born + IB + Vfin);
}

double Eps2Pole(int n, int bornpdgs[], double bornme) {
    double res = 0.0;
    for (int i = 0; i < n; ++i) {
        double q = Physics::PDG::Charge(bornpdgs[i]);
        res += q * q;
    }
    return res * bornme;
}

double EpsPole(int n, int bornpdgs[], const Math::FourMomentum *momenta,
               double born, double Q2) {
    double sum = 0.0;
#ifndef NDEBUG
    double charge = 0.0;
#endif
    for (int i = 0; i < n; i++) {
        if (bornpdgs[i] == 22) {
            assert(0 && "photon in born not implemented");
        }
        if (!particle_contributes<Type::QED>(bornpdgs[i])) {
            continue;
        }
        double sigma_i = sigma_f(i, bornpdgs[i]);
        double charge_i = Physics::PDG::Charge(abs(bornpdgs[i]));

        sum += (3.0 / 2.0) * charge_i * charge_i * born;

#ifndef NDEBUG
        charge += charge_i * sigma_i;
#endif
        for (int j = i + 1; j < n; j++) {
            if (bornpdgs[j] == 22) {
                assert(0 && "photon in born not implemented");
            }
            if (!particle_contributes<Type::QED>(bornpdgs[j])) {
                continue;
            }
            double sigma_j = sigma_f(j, bornpdgs[j]);
            double charge_j = Physics::PDG::Charge(abs(bornpdgs[j]));
            double Log = log(2.0 * momenta[i].Dot(momenta[j]) / Q2);
            // factor 2 is there because the sum is over i,j, i != j and the for
            // loops are only over j > i.
            sum += 2.0 * Log * sigma_i * sigma_j * charge_i * charge_j * born;
        }
    }
    assert(fabs(charge) < 1e-15 && "charge conservation violated");

    return sum;
}

} // namespace QED

} // namespace FKS
