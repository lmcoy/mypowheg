//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.0, 2016-05-12
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//
// 2016 06 07: modified by Lennart Oymanns
//
//==========================================================================

#ifndef WQX_H_HS89KJAQ
#define WQX_H_HS89KJAQ

#include <complex>

namespace Phasespace {
class Phasespace;
}
class Parameters_sm;

/**
 * @brief W q~ matrix element
 *
 * Wqx calculates the matrix element for 
 *     g q~ -> mu+ nu q'~
 * and 
 *     g q~ -> nu~ mu- q'~.
 *
 * Note that the CKM matrix is not included.
 */
class Wqx {
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
    static const int nwavefuncs = 8;
    std::complex<double> w[nwavefuncs][18];
    static const int namplitudes = 2;
    std::complex<double> amp[namplitudes];
    double matrix();

    // vector with external particle masses
    double mME[nexternal] = { 0.0 };

    // vector with momenta (to be changed each event)
    double momenta[nexternal][4];
};

#endif 
