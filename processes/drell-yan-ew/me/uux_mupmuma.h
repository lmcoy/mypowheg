//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.0.1, 2014-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef uux_mupmuma_H
#define uux_mupmuma_H

#include <complex> 

#include "phasespace/phasespace.h"
#include "parameters_sm.h"

//==========================================================================
// A class for calculating the matrix elements for
// Process: u u~ > mu+ mu- a WEIGHTED=6
// Process: c c~ > mu+ mu- a WEIGHTED=6
//--------------------------------------------------------------------------

class ME_uux_mupmuma
{
  public:
// Calculate flavour-independent parts of cross section.
    double Calculate(const Phasespace::Phasespace &ps, int perm[],
                     const Parameters_sm &param,
                     const Parameters_alphaS &param_aS, int flag = DEFAULT);

    // Constants for array limits
    static const int ninitial = 2; 
    static const int nexternal = 5; 
    static const int nprocesses = 2; 

    static const int DEFAULT = 0;
    static const int ONLYISR = 1;
    static const int ONLYFSR = 2;

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
        void calculate_wavefunctions(const int perm[], const int hel[],
                                 const Parameters_sm &par,
                                 const Parameters_alphaS &param_aS, int flag);
    static const int nwavefuncs = 13; 
    std::complex<double> w[nwavefuncs][18]; 
    static const int namplitudes = 8; 
    std::complex<double> amp[namplitudes]; 
    double matrix_uux_mupmuma(); 

    // Store the matrix element value
    double matrix_element; 

    // Color flows, used when selecting color
    double  jamp2; 

// vector with external particle masses
    double mME[5] = { 0.0 };

    // vector with momenta (to be changed each event)
    double momenta[nexternal][4];

}; 


#endif  // uux_mupmuma_H
