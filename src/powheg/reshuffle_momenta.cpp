
#include <iostream>

#include <boost/math/tools/roots.hpp>

#include "powheg/reshuffle_momenta.h"

#include "libconfig.h"
#include "lhe/event.h"
#include "powheg/resonance.h"
#include "math/fourmomentum.h"
#include "process/data.h"
#include "fks/param.h"

namespace Powheg {

namespace {

class Boost {
  public:
    Boost() {
        data_[index(0, 0)] = 1.0;
        data_[index(1, 1)] = 1.0;
        data_[index(2, 2)] = 1.0;
        data_[index(3, 3)] = 1.0;
    }
    double Get(int i, int j) const { return data_[index(i, j)]; }
    void Set(int i, int j, double v) { data_[index(i, j)] = v; }
    void SetRow(int i, double a, double b, double c, double d) {
        data_[index(i, 0)] = a;
        data_[index(i, 1)] = b;
        data_[index(i, 2)] = c;
        data_[index(i, 3)] = d;
    }
    Math::FourMomentum Apply(const Math::FourMomentum &p) {
        double E = 0.0;
        double X = 0.0;
        double Y = 0.0;
        double Z = 0.0;
        for (int i = 0; i < 4; i++) {
            double x = p.At(i);
            E += Get(0, i) * x;
            X += Get(1, i) * x;
            Y += Get(2, i) * x;
            Z += Get(3, i) * x;
        }
        return Math::FourMomentum(E, X, Y, Z);
    }

  private:
    int index(int i, int j) const { return 4 * i + j; }
    double data_[16] = { 0.0 };
};

Boost gen_boost_to_rest(const Math::FourMomentum & p) {
    double t = p.E();
    double x = p.PX();
    double y = p.PY();
    double z = p.PZ();
    double l = sqrt(x * x + y * y + z * z);
    if (l < 1e-12) {
        return Boost();
    }
    double nx = x / l;
    double ny = y / l;
    double nz = z / l;
    double nx2 = nx * nx;
    double ny2 = ny * ny;
    double nz2 = nz * nz;
    double b = l / t;
    double g = 1.0 / sqrt(1.0 - b * b);
    double gm = g - 1.0;
    double gb = g * b;
    Boost boost;
    boost.SetRow(0, g, -gb * nx, -gb * ny, -gb * nz);
    boost.SetRow(1, -gb * nx, 1 + gm * nx2, gm * nx * ny, gm * nx * nz);
    boost.SetRow(2, -gb * ny, gm * ny * nx, 1 + gm * ny2, gm * ny * nz);
    boost.SetRow(3, -gb * nz, gm * nz * nx, gm * nz * ny, 1 + gm * nz2);
    return boost;
}

Boost invert_boost(const Boost & boost) {
    Boost ret = boost;
    ret.Set(0,1, -boost.Get(0,1));
    ret.Set(0,2, -boost.Get(0,2));
    ret.Set(0,3, -boost.Get(0,3));

    ret.Set(1,0, -boost.Get(1,0));
    ret.Set(2,0, -boost.Get(2,0));
    ret.Set(3,0, -boost.Get(3,0));
    return ret;
}

double scale_factor_2(double m1, double m2, const Math::FourMomentum &p ) {
    double m12 = m1 * m1;
    double m22 = m2 * m2;
    double p2 = p.E() * p.E();
    double p4 = p2 * p2;
    double c2 =
        1.0 - 0.5 * (m12 + m22) / p2 + (m12 - m22) * (m12 - m22) / 16.0 / p4;
    return sqrt(c2);
}

double scale_factor_3(double m1, double m2, double m3,
                      const Math::FourMomentum &p1, Math::FourMomentum &p2,
                      Math::FourMomentum &p3) {
    auto eqn = [m1, m2, m3, p1, p2, p3](double c2) {
        double p12 = p1.MomentumMagnitudeSqr();
        double p22 = p2.MomentumMagnitudeSqr();
        double p32 = p3.MomentumMagnitudeSqr();
        double t1 = sqrt(m1*m1 + c2*p12);
        double t2 = sqrt(m2*m2 + c2*p22);
        double t3 = sqrt(m3 * m3 + c2 * p32);
        double E = p1.E() + p2.E() + p3.E();
        return t1 + t2 + t3 - E;
    };

    double fa = eqn(1e-12);
    double fb = eqn(1.0);
    if ( fa * fb > 0.0 ) {
        // no solution
        std::cerr << "no solution for momenta reshuffle\n";
        return 1.0;
    }
    boost::uintmax_t max_iter = 100;
    // relative difference: 2^(1-15) = 6e-5
    boost::math::tools::eps_tolerance<double> tol(15);
    auto bracket = boost::math::tools::toms748_solve(eqn, 1e-12, 1.0, fa, fb,
                                                     tol, max_iter);
    return sqrt(bracket.first);
}

} // namespace

void LIB_LOCAL reshuffle_momenta(LHE::Event *event,
                                 const Powheg::Resonance &resonance,
                                 const UserProcess::Data *userdata) {
    int mother_id = 2;
    // we have to shift by one because the resonance is already inserted in
    // event.
    int daughter1_id = resonance.ID[0] + 1;
    int daughter2_id = resonance.ID[1] + 1;
    int daughter3_id = resonance.ID[2] + 1;

    int num_daughter = 2;
    if (daughter3_id > 2) {
        num_daughter = 3;
    }

    Boost torestframe = gen_boost_to_rest(event->Particles[mother_id].Momentum);
    auto boost_back = invert_boost(torestframe);
 
    auto Daughter1_r =
        torestframe.Apply(event->Particles[daughter1_id].Momentum);

    double m1 = userdata->Params->Mass(event->Particles[daughter1_id].PDG);
    double m2 = userdata->Params->Mass(event->Particles[daughter2_id].PDG);
    if (num_daughter == 2) {
        double c = scale_factor_2(m1, m2, Daughter1_r);
        Daughter1_r.Scale(c);
        double p2 = Daughter1_r.MomentumMagnitudeSqr();
        Daughter1_r.SetE(sqrt(m1 * m1 + p2));
        Math::FourMomentum Daughter2_r = Daughter1_r;
        Daughter2_r.Scale(-1.0);
        Daughter2_r.SetE(sqrt(m2 * m2 + p2));

        event->Particles[daughter1_id].Momentum = boost_back.Apply(Daughter1_r);
        event->Particles[daughter2_id].Momentum = boost_back.Apply(Daughter2_r);

        return;
    }

    double m3 = userdata->Params->Mass(event->Particles[daughter3_id].PDG);

    auto Daughter2_r =
        torestframe.Apply(event->Particles[daughter2_id].Momentum);
    auto Daughter3_r =
        torestframe.Apply(event->Particles[daughter3_id].Momentum);
    double c =
        scale_factor_3(m1, m2, m3, Daughter1_r, Daughter2_r, Daughter3_r);
    Daughter1_r.Scale(c);
    double p12 = Daughter1_r.MomentumMagnitudeSqr();
    Daughter1_r.SetE(sqrt(m1 * m1 + p12));

    Daughter2_r.Scale(c);
    double p22 = Daughter2_r.MomentumMagnitudeSqr();
    Daughter2_r.SetE(sqrt(m2 * m2 + p22));

    Daughter3_r.Scale(c);
    double p32 = Daughter3_r.MomentumMagnitudeSqr();
    Daughter3_r.SetE(sqrt(m3 * m3 + p32));

    event->Particles[daughter1_id].Momentum = boost_back.Apply(Daughter1_r);
    event->Particles[daughter2_id].Momentum = boost_back.Apply(Daughter2_r);
    event->Particles[daughter3_id].Momentum = boost_back.Apply(Daughter3_r);
}

} // namespace Powheg

