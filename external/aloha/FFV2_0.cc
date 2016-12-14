// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// Gamma(3,2,-1)*ProjM(-1,1)
// 
#include "FFV2_0.h"

void FFV2_0(std::complex<double> F1[], std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, std::complex<double> & vertex)
{
  std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP4; 
  TMP4 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])));
  vertex = COUP * - cI * TMP4; 
}

