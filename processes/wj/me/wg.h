//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.0, 2016-05-12
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//
// 2016 06 07: modified by Lennart Oymanns
//
//==========================================================================

#ifndef WG_H_BSHD23OW
#define WG_H_BSHD23OW

#include <complex>

#include "parameters_sm.h"
#include "phasespace/phasespace.h"

#include "process/matrixelement.h"

/**
 * @brief W g matrix element
 *
 * Wg calculates the matrix element for
 *      q q'~ -> mu+ nu g
 * and
 *      q q'~ -> nu~ mu- g.
 *
 * Note that the CKM matrix is not included.
 */
class Wg {
  public:
    // Calculate flavour-independent parts of cross section.
    UserProcess::SpinCorrelated Calculate(const Phasespace::Phasespace &ps,
                                          int perm[],
                                          const Parameters_sm &param);

    // Constants for array limits
    static const int ninitial = 2;
    static const int nexternal = 5;
    static const int nprocesses = 2;

  private:
    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    std::array<std::complex<double>, 5>
    calculate_wavefunctions(const int perm[], const int hel[],
                            const Parameters_sm &par, int helicity_gluon);

    // vector with external particle masses
    double mME[nexternal] = { 0.0 };

    // vector with momenta (to be changed each event)
    double momenta[nexternal][4];
};

#endif
