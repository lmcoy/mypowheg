//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.0, 2016-05-12
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//
// 2016 06 07: modified by Lennart Oymanns
//
//==========================================================================

#include "wq.h"

#include "../external/aloha_cm/aloha.h"

#include "phasespace/phasespace.h"
#include "parameters_sm.h"


//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.
double Wq::Calculate(const Phasespace::Phasespace &ps, int perm[],
                      const Parameters_sm &param) {




  // Local variables and constants
  const int ncomb = 32; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0; 

 

 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1}, {-1,
      -1, -1, -1, 1}, {-1, -1, -1, 1, -1}, {-1, -1, -1, 1, 1}, {-1, -1, 1, -1,
      -1}, {-1, -1, 1, -1, 1}, {-1, -1, 1, 1, -1}, {-1, -1, 1, 1, 1}, {-1, 1,
      -1, -1, -1}, {-1, 1, -1, -1, 1}, {-1, 1, -1, 1, -1}, {-1, 1, -1, 1, 1},
      {-1, 1, 1, -1, -1}, {-1, 1, 1, -1, 1}, {-1, 1, 1, 1, -1}, {-1, 1, 1, 1,
      1}, {1, -1, -1, -1, -1}, {1, -1, -1, -1, 1}, {1, -1, -1, 1, -1}, {1, -1,
      -1, 1, 1}, {1, -1, 1, -1, -1}, {1, -1, 1, -1, 1}, {1, -1, 1, 1, -1}, {1,
      -1, 1, 1, 1}, {1, 1, -1, -1, -1}, {1, 1, -1, -1, 1}, {1, 1, -1, 1, -1},
      {1, 1, -1, 1, 1}, {1, 1, 1, -1, -1}, {1, 1, 1, -1, 1}, {1, 1, 1, 1, -1},
      {1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const double denominators = 96; 

  ntry +=  1; 

    double matrix_element = 0.; 




    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        double t = matrix(ps, perm, helicities[ihel],param); 
	matrix_element += t;
            if (ntry < 2) {
        // Mirror initial state momenta for mirror process 
	int p0 = perm[0];
                int p1 = perm[1];
                perm[0] = p1;
                perm[1] = p0;
        // Calculate wavefunctions
        double tsum = t + matrix(ps, perm, helicities[ihel],param); 
        // Mirror back 
	perm[0] = p0;
                perm[1] = p1;
         
        
        // Store which helicities give non-zero result
        if (tsum != 0. )
        {
          goodhel[ihel] = true; 

        }
	    }
      }
    }
 
  
  


    return matrix_element / denominators; 



}



//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess
double Wq::matrix(const Phasespace::Phasespace &ps, const int perm[], const int hel[], const Parameters_sm &pars)
{
    constexpr std::complex<double> ZERO(0.0, 0.0);
    constexpr int nwavefuncs = 8; 
    std::complex<double> w[nwavefuncs][18]; 
    
        // vector with momenta 
    double momenta[nexternal][4];
    
        for (int i = 0; i < nexternal; i++) {
        momenta[i][0] = ps.Momenta[i].E();
        momenta[i][1] = ps.Momenta[i].PX();
        momenta[i][2] = ps.Momenta[i].PY();
        momenta[i][3] = ps.Momenta[i].PZ();
    }

  // Calculate all wavefunctions
  vxxxxx(momenta[perm[0]], mME[0], hel[0], -1, w[0]); 
  ixxxxx(momenta[perm[1]], mME[1], hel[1], +1, w[1]); 
  ixxxxx(momenta[perm[2]], mME[2], hel[2], -1, w[2]); 
  oxxxxx(momenta[perm[3]], mME[3], hel[3], +1, w[3]); 
  oxxxxx(momenta[perm[4]], mME[4], hel[4], +1, w[4]); 
  FFV1_2(w[1], w[0], pars.GC_11, ZERO, w[5]); 
  FFV2_3(w[2], w[3], pars.GC_100, pars.MuW, w[6]); 
  FFV1_1(w[4], w[0], pars.GC_11, ZERO, w[7]); 

  // Calculate all amplitudes
 constexpr int namplitudes = 2; 
    std::complex<double> amp[namplitudes]; 
  FFV2_0(w[5], w[4], w[6], pars.GC_100, amp[0]); 
  FFV2_0(w[1], w[7], w[6], pars.GC_100, amp[1]); 


  // The color matrix;
  constexpr double cf = 4; 

  // Calculate color flows
  std::complex<double> jamp = amp[0] + amp[1]; 

  // Sum and square the color flows to get the matrix element
   return cf * std::norm(jamp); 
}



