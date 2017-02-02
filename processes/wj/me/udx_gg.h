//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.2, 2016-06-10
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//
// 2016 07 11: modified by Lennart Oymanns
//
//==========================================================================

#ifndef UDX_GG_SJS72SDQ_H
#define UDX_GG_SJS72SDQ_H

#include <complex>

#include "parameters_sm.h"
#include "phasespace/phasespace.h"

//==========================================================================
// A class for calculating the matrix elements for
// Process: u d~ > mu+ vm g g WEIGHTED<=6 @1
// Process: c s~ > mu+ vm g g WEIGHTED<=6 @1
//--------------------------------------------------------------------------

class UDX_GG {
  public:
    double Calculate(const Phasespace::Phasespace &ps, int perm[],
                     const Parameters_sm &param);

    // Constants for array limits
    static const int ninitial = 2;
    static const int nexternal = 6;
    static const int nprocesses = 2;

  private:
    double calculate_wavefunctions(const int perm[], const int hel[],
                                   const Parameters_sm &pars,
                                   double momenta[][4]);

    // vector with external particle masses
    double mME[nexternal] = { 0.0 };
};

#endif
