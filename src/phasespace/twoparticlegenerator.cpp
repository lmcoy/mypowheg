#include "phasespace/twoparticlegenerator.h"

#include "math/math.h"
#include <cmath>

namespace Phasespace {
/**
@brief Generate 2 particle phase space.

A two particle phase space is created from the input parameters
-# x1 - momentum fraction of the 1st inital state parton
-# x2 - momentum fraction of the 2nd inital state parton
-# y  - cos(theta) of the first final state particle
-# phi - azimuthal angle of the 1st final state particle

These parameters have to be in the above order in the array x.
All parameters have to be from 0 to 1, i.e. the value for phi is multiplied
internally with 2pi to obtain the angle. (y is converted internally as well).

@param[out] ps the phase space is written to ps.
@param n number of parameters has to be 4.
@param x n parameters.
*/
void TwoParticleGenerator::operator()(Phasespace *ps, int n,
                                      double const *const x, double S, int m,
                                      double *masses) const {
    assert(n == 4 && "need four parameters");
    const double x1 = x[0];
    const double x2 = x[1];
    assert(x1 > 0.0 && x1 <= 1.0 && "x1 not in range (0,1]");
    assert(x2 > 0.0 && x2 <= 1.0 && "x2 not in range (0,1]");
    const double CosTheta = 2.0 * x[2] - 1.0;
    const double phi = 2.0 * Math::Pi * x[3];
    assert(fabs(CosTheta) <= 1.0 && "cos theta not in range [-1.0, 1.0]");
    // assert(phi >= 0.0 && phi < 2.0 * Math::Pi &&
    //        "azimuth angle not in range [0, 2pi)");

    ps->N = 2;
    ps->S = S;
    const double s = x1 * x2 * S;
    const double sqrts = sqrt(s);

    ps->X1 = x1;
    ps->X2 = x2;

    // initial state particles
    const double p = sqrts / 2.0;
    ps->Momenta[0].Set(p, 0.0, 0.0, p);
    ps->Momenta[1].Set(p, 0.0, 0.0, -p);

    ps->Masses[0] = masses[0];
    ps->Masses[1] = masses[1];
    const double m1 = ps->Masses[0];
    const double m2 = ps->Masses[1];
    // final state particles
    const double m1_2 = m1 * m1;
    const double m2_2 = m2 * m2;
    const double m1_4 = m1_2 * m1_2;
    const double m2_4 = m2_2 * m2_2;
    const double p_final = sqrt(m1_4 - 2.0 * m1_2 * m2_2 + m2_4 -
                                2.0 * m1_2 * s - 2.0 * m2_2 * s + s * s) /
                           (2.0 * sqrts);
    ps->Momenta[2].SetFromPMCosThetaPhi(p_final, m1, CosTheta, phi);
    const double px = -ps->Momenta[2].PX();
    const double py = -ps->Momenta[2].PY();
    const double pz = -ps->Momenta[2].PZ();
    ps->Momenta[3].Set(sqrt(p_final * p_final + m2_2), px, py, pz);

    ps->Jacobian = 1.0/(8.0*Math::Pi);
}

} // end namespace
