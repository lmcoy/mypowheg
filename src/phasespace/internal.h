#ifndef PHASESPACE_INTERNAL_H_
#define PHASESPACE_INTERNAL_H_

#include <cmath>

namespace Phasespace {

namespace {

/// gamma returns gamma for a Lorentz boost for parton CMS -> lab frame.
double gamma(double x1, double x2) {
    return 0.5 * (x1 + x2) / sqrt(x1 * x2);
}

/// betagamma returns beta*gamma for a Lorentz boost for parton CMS -> lab
/// frame.
double betagamma(double x1, double x2) {
    return -0.5 * (x1 - x2) / sqrt(x1 * x2);
}

/// boost_z boosts the 4-momentm \p p in z direction.
void boost_z(double gamma, double betagamma, Math::FourMomentum *p) {
    const double E = p->E();
    const double pz = p->PZ();
    p->SetE(E * gamma - betagamma * pz);
    p->SetPZ(gamma * pz - E * betagamma);
}
} // end namespace

} // end namespace Phasespace

#endif
