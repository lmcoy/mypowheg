// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// 1
// 
#include "SSSS1_4.h"

void SSSS1_4(std::complex<double> S1[], std::complex<double> S2[], std::complex<double> S3[],
    std::complex<double> COUP, double M4, double W4, std::complex<double> S4[])
{
  std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> denom; 
  double P4[4]; 
  S4[0] = +S1[0] + S2[0] + S3[0]; 
  S4[1] = +S1[1] + S2[1] + S3[1]; 
  P4[0] = -S4[0].real(); 
  P4[1] = -S4[1].real(); 
  P4[2] = -S4[1].imag(); 
  P4[3] = -S4[0].imag(); 
  denom = COUP/(pow(P4[0], 2) - pow(P4[1], 2) - pow(P4[2], 2) - pow(P4[3], 2) -
      M4 * (M4 - cI * W4));
  S4[2] = denom * cI * S3[2] * S2[2] * S1[2]; 
}

