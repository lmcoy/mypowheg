#include "phasespace.h"

#include "math/fourmomentum.h"
#include <cmath>
#include <iostream>

static Math::FourMomentum boost_from_rest(const Math::FourMomentum &q,
                                          const Math::FourMomentum &ref) {
    double b[] = {ref.At(1) / ref.At(0), ref.At(2) / ref.At(0),
                  ref.At(3) / ref.At(0)};
    double b2 = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
    double gamma = 1 / sqrt(1.0 - b2);
    double bp = b[0] * q.At(1) + b[1] * q.At(2) + b[2] * q.At(3);
    double gamma2 = (gamma - 1) / b2;

    double t = gamma * (q.At(0) + bp);
    double x = q.At(1) + gamma2 * bp * b[0] + gamma * b[0] * q.At(0);
    double y = q.At(2) + gamma2 * bp * b[1] + gamma * b[1] * q.At(0);
    double z = q.At(3) + gamma2 * bp * b[2] + gamma * b[2] * q.At(0);
    return Math::FourMomentum(t, x, y, z);
}

/**
 * @brief Generate 3 particle phase space
 *
 * GenPhasespace generates a phase space with 3 final state particles.
 * Particle 3 and 4 are generated as a decay product of an intermediate
 * particle Q2. All phase space parameters v have to be in [0,1].
 *
 * @param S proton CMS
 * @param v[0] momentum fraction of parton 1
 * @param v[1] momentum fraction of parton 2
 * @param v[2] mass^2 of intermediate particle Q2 -> p1 p2
 * @param v[3] cos theta of particle 5 in partonic CMS
 * @param v[4] phi of particle 5 in partonic CMS
 * @param v[5] cos theta of particle 3 in rest frame of Q2
 * @param v[6] phi of particle 3 in rest frame of Q2
 */
Phasespace::Phasespace GenPhasespace(double S, const std::array<double, 7> &v) {
    assert(v[0] >= 0.0 && v[0] <= 1.0);
    assert(v[1] >= 0.0 && v[1] <= 1.0);
    assert(v[2] >= 0.0 && v[2] <= 1.0);
    assert(v[3] >= 0.0 && v[3] <= 1.0);
    assert(v[4] >= 0.0 && v[4] <= 1.0);
    assert(v[5] >= 0.0 && v[5] <= 1.0);
    assert(v[6] >= 0.0 && v[6] <= 1.0);
    double x1 = v[0];
    double x2 = v[1];
    double s = x1 * x2 * S;
    double p = sqrt(s) / 2.0;

    Math::FourMomentum k1(p, 0.0, 0.0, p);
    Math::FourMomentum k2(p, 0.0, 0.0, -p);

    double Q2 = v[2] * s;

    double E1 = (s - Q2) / sqrt(s) / 2.0;
    double costheta1 = (2 * v[3] - 1.0);
    double sintheta1 = sqrt(1 - costheta1 * costheta1);
    double phi1 = 2. * M_PI * v[4];
    double costheta2 = (2 * v[5] - 1.0);
    double sintheta2 = sqrt(1 - costheta2 * costheta2);
    double phi2 = 2. * M_PI * v[6];

    Math::FourMomentum p1(E1, E1 * sintheta1 * cos(phi1),
                          E1 * sintheta1 * sin(phi1), E1 * costheta1);
    Math::FourMomentum q(2. * p - p1.At(0), -p1.At(1), -p1.At(2), -p1.At(3));

    // p2, p3 in rest frame of q
    double t = q.Dot(q);
    double E2 = 0.5 * sqrt(t);
    Math::FourMomentum p2(E2, E2 * sintheta2 * cos(phi2),
                          E2 * sintheta2 * sin(phi2), E2 * costheta2);
    Math::FourMomentum p3(p2.At(0), -p2.At(1), -p2.At(2), -p2.At(3));

    // boost p2, p3 to frame of p1 and q
    p2 = boost_from_rest(p2, q);
    p3 = boost_from_rest(p3, q);

    double p2p2 = p2.Dot(p2);
    if (p2p2 < 0.0) {
        double x = sqrt(-p2p2 + p2.MomentumMagnitudeSqr());
        p2.SetE(x);
    }
    double p3p3 = p3.Dot(p3);
    if (p3p3 < 0.0) {
        double x = sqrt(-p3p3 + p3.MomentumMagnitudeSqr());
        p3.SetE(x);
    }

    double jac = (1 - v[2]) * s / (8 * M_PI * M_PI * M_PI * 16.0);

    Phasespace::Phasespace ps;
    ps.N = 3;
    ps.S = S;
    ps.X1 = x1;
    ps.X2 = x2;
    ps.Jacobian = jac;
    ps.Momenta[0] = k1;
    ps.Momenta[1] = k2;
    ps.Momenta[2] = p2;
    ps.Momenta[3] = p3;
    ps.Momenta[4] = p1;

    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;

    return ps;
}

Phasespace::Phasespace BornPSGenerator::Gen(int n, double *x,
                                            UserProcess::Data *userdata) {
    double SqrtS = userdata->SqrtS;
    std::array<double, 7> v = {0.0};
    assert(n >= 7);
    for (int i = 0; i < 7; i++) {
        v[i] = x[i];
    }
    return GenPhasespace(SqrtS * SqrtS, v);
}