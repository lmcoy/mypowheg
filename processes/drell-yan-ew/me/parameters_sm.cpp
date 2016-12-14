//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph5_aMC@NLO v. 2.0.1, 2014-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "parameters_sm.h"
#include <iostream>
#include <cassert>

using namespace std;

void Parameters_sm::Set(double Gf, double aEWM1, double mz, double wz) {
    WZ = wz;
    MZ = mz;
    double MZ__exp__2 = pow(MZ, 2.);
    double MW__exp__2 = pow(MW, 2.);
    alpha = sqrt(2.0) / M_PI * Gf * MW__exp__2 * (1.0 - MW__exp__2 / MZ__exp__2);
    ee = 2. * sqrt(alpha) * sqrt(M_PI);
    std::complex<double> mu_w_2 = std::complex<double>(MW__exp__2, -MW * WidthW);
    std::complex<double> mu_z_2 = std::complex<double>(MZ__exp__2, -MZ * WZ);

    MuZ = sqrt(mu_z_2);
    std::complex<double> cw2 = mu_w_2 / mu_z_2;
    cw = sqrt(cw2);
    sw = sqrt(1.0 - cw2);

    GC_1 = std::complex<double>(0.0, -ee*(1./ 3.));
    GC_2 = std::complex<double>(0.0, ee*(2./ 3.));
    GC_3 = std::complex<double>(0.0,-ee);
    std::complex<double> cI(0.0, 1.0);
    GC_50 = -cI * (cw * ee) / (2. * sw);
    GC_51 = cI * (cw * ee) / (2. * sw);
    GC_58 = -cI * (ee * sw) / (6. * cw);
    GC_59 = cI * (ee * sw) / (2. * cw);
}

double Parameters_sm::Mass(int pdg) const {
    int apdg = std::abs(pdg);
    switch (apdg) {
    case 1:
        return MDQuark;
    case 2:
        return MUQuark;
    case 3:
        return MSQuark;
    case 4:
        return MCQuark;
    case 5:
        return MBquark;
    case 6:
        return MTQuark;
    case 11:
        return MElectron;
    case 12:
        return 0.0;
    case 13:
        return MMuon;
    case 14:
        return 0.0;
    case 15:
        return MTau;
    case 16:
        return 0.0;
    case 21:
        return 0.0;
    case 22:
        return 0.0;
    case 23:
        return MZ;
    case 24:
        return MW;
    case 25:
        return MH;
    }
    assert(0 && "mass not implemented");

    return 0.0;
}

void Parameters_alphaS::Set(double alphaS) {
  aS = alphaS;
  sqrt__aS = sqrt(aS);
    G = 2. * sqrt__aS * sqrt(M_PI);
    G__exp__2 = pow(G, 2.);
    GC_11 = std::complex<double>(0.0, G);
}


