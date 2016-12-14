//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.0.1, 2014-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef GUX_MUPMUMUX_H
#define GUX_MUPMUMUX_H

#include <complex>

#include "phasespace/phasespace.h"
#include "parameters_sm.h"


//==========================================================================
// A class for calculating the matrix elements for
// Process: g u~ > mu+ mu- u~ WEIGHTED=5
// Process: g c~ > mu+ mu- c~ WEIGHTED=5
//--------------------------------------------------------------------------

class ME_gux_mupmumux {
  public:
    // Constructor.
    ME_gux_mupmumux()  {
    }

    ~ME_gux_mupmumux() {
    }


    // Calculate flavour-independent parts of cross section.
    double Calculate(const Phasespace::Phasespace &ps, int perm[],
                   const Parameters_sm &param,
                   const Parameters_alphaS &param_aS);



    // Constants for array limits
    static const int ninitial = 2;
    static const int nexternal = 5;
    static const int nprocesses = 2;

  private:
    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[],const Parameters_sm &par,
                                 const Parameters_alphaS &param_aS);
    static const int nwavefuncs = 9;
    std::complex<double> w[nwavefuncs][18];
    static const int namplitudes = 4;
    std::complex<double> amp[namplitudes];
    double matrix_gux_mupmumux();

    // Store the matrix element value from sigmaKin
    double matrix_element;

    // Color flows, used when selecting color
    double jamp2;

    // vector with external particle masses
    double mME[nexternal] = { 0.0 };

    // vector with momenta (to be changed each event)
    double momenta[nexternal][4];
};


#endif
