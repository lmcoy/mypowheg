//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.0, 2016-05-12
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//
// 2016 06 07: modified by Lennart Oymanns
//
//==========================================================================

#ifndef WA_H_3HJ80JK2K
#define WA_H_3HJ80JK2K

#include <complex>

#include "phasespace/phasespace.h"
#include "parameters_sm.h"

/**
 * @brief W gamma matrix element
 *
 * Wa calculates the matrix element for 
 *     q q'~ -> mu+ nu a 
 * and 
 *     q q'~ -> nu~ mu- a.
 * (a=gamma)
 *
 * Note that the CKM matrix is not included.
 */
class Wa {
  public:
    // Calculate flavour-independent parts of cross section.
    double Calculate(const Phasespace::Phasespace &ps, int perm[],
                     const Parameters_sm &param);

    // Constants for array limits
    static const int ninitial = 2;
    static const int nexternal = 5;
    static const int nprocesses = 2;

  private:
    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[],
                                 const Parameters_sm &par);
    static const int nwavefuncs = 10;
    std::complex<double> w[nwavefuncs][18];
    static const int namplitudes = 4;
    std::complex<double> amp[namplitudes];
    double matrix();

    // vector with external particle masses
    double mME[nexternal] = { 0.0 };

    // vector with momenta (to be changed each event)
    double momenta[nexternal][4];
};

#endif
