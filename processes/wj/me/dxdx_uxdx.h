//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.2, 2016-06-10
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//
// 2016 07 04: modified by Lennart Oymanns
//
//==========================================================================

#ifndef DXDX_UXDX_D9AWC8AF_H
#define DXDX_UXDX_D9AWC8AF_H

#include <complex>

#include "parameters_sm.h"
#include "phasespace/phasespace.h"

//==========================================================================
// A class for calculating the matrix elements for
// Process: d~ d~ > mu+ vm u~ d~ WEIGHTED<=6 @1
// Process: s~ s~ > mu+ vm c~ s~ WEIGHTED<=6 @1
//--------------------------------------------------------------------------

class DXDX_UXDX {
  public:
    double Calculate(const Phasespace::Phasespace &ps, int perm[],
                     const Parameters_sm &param);

    // Constants for array limits
    static const int ninitial = 2;
    static const int nexternal = 6;
    static const int nprocesses = 1;

  private:
    double calculate_wavefunctions(const int perm[], const int hel[],
                                   const Parameters_sm &pars);

    // vector with external particle masses
    double mME[nexternal] = { 0.0 };

    // vector with momenta (to be changed each event)
    double momenta[nexternal][4];
};

#endif
