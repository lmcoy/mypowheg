#ifndef FKS_COLOR_H_
#define FKS_COLOR_H_

#include "math/math.h"

namespace {

const double CA = 3.0;
const double CF = 4.0 / 3.0;
const double nf = 5.0;
const double TF = 0.5;

double gamma(int flavour) {
    static const double gamma_gluon = (11.0 * CA - 4.0 * TF * nf) / 6.0;
    static const double gamma_quark = (3.0 / 2.0) * CF;
    if (flavour == 0 || flavour == 21) {
        return gamma_gluon;
    }
    if (flavour > -6 && flavour < 6) {
        return gamma_quark;
    }
    return 0.0;
}

double gamma_prime(int flavour) {
    static const double gamma_prime_gluon =
        ((67.0 / 9.0) - (2.0 / 3.0) * Math::Pi * Math::Pi) * CA -
        (23.0 / 9.0) * TF * nf;
    static const double gamma_prime_quark =
        (13.0 / 2.0 - 2 * Math::Pi * Math::Pi / 3.0) * CF;
    if (flavour == 0 || flavour == 21) {
        return gamma_prime_gluon;
    }
    if (flavour > -6 && flavour < 6) {
        return gamma_prime_quark;
    }

    return 0.0;
}

double ColorC(int flavour) {
    if (flavour == 0 || flavour == 21) {
        return CA;
    }
    if (flavour > -6 && flavour < 6) {
        return CF;
    }
    return 0.0;
}

} // end namespace

#endif
