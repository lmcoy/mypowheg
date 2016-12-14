//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.0, 2016-05-12
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//
// 2016 06 07: modified by Lennart Oymanns
//
//==========================================================================

#ifndef WQ_H_HSJ4SJK1
#define WQ_H_HSJ4SJK1

#include <complex> 

namespace Phasespace {
class Phasespace;
}
class Parameters_sm;

/**
 * @brief W q matrix element
 *
 * Wg calculates the matrix element for 
 *      g q -> mu+ nu q' 
 * and 
 *      g q -> nu~ mu- q'.
 *
 * Note that the CKM matrix is not included.
 */
class Wq
{
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
    double matrix(const Phasespace::Phasespace &ps,const int perm[], const int hel[],const Parameters_sm &par ); 


    // vector with external particle masses
    double mME[nexternal] = { 0.0 };
    

}; 


#endif  // MG5_Sigma_sm_gu_mupvmd_H
