// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// Metric(1,4)*Metric(2,3) - Metric(1,2)*Metric(3,4)
// 
#include "VVVV3_0.h"

void VVVV3_0(std::complex<double> V1[], std::complex<double> V2[], std::complex<double> V3[],
    std::complex<double> V4[], std::complex<double> COUP, std::complex<double> & vertex)
{
  std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP21; 
  std::complex<double> TMP27; 
  std::complex<double> TMP19; 
  std::complex<double> TMP13; 
  TMP27 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP19 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP21 = (V1[2] * V4[2] - V1[3] * V4[3] - V1[4] * V4[4] - V1[5] * V4[5]); 
  TMP13 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]); 
  vertex = COUP * (-cI * (TMP19 * TMP21) + cI * (TMP13 * TMP27)); 
}

