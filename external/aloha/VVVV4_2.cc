// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// Metric(1,3)*Metric(2,4) - Metric(1,2)*Metric(3,4)
// 
#include "VVVV4_2.h"

void VVVV4_2(std::complex<double> V1[], std::complex<double> V3[], std::complex<double> V4[],
    std::complex<double> COUP, double M2, double W2, std::complex<double> V2[])
{
  std::complex<double> cI = std::complex<double> (0., 1.); 
  double P2[4]; 
  std::complex<double> TMP17; 
  std::complex<double> denom; 
  std::complex<double> TMP27; 
  std::complex<double> TMP24; 
  std::complex<double> TMP9; 
  double OM2; 
  OM2 = 0.; 
  if (M2 != 0.)
    OM2 = 1./pow(M2, 2); 
  V2[0] = +V1[0] + V3[0] + V4[0]; 
  V2[1] = +V1[1] + V3[1] + V4[1]; 
  P2[0] = -V2[0].real(); 
  P2[1] = -V2[1].real(); 
  P2[2] = -V2[1].imag(); 
  P2[3] = -V2[0].imag(); 
  TMP24 = (P2[0] * V4[2] - P2[1] * V4[3] - P2[2] * V4[4] - P2[3] * V4[5]); 
  TMP17 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP27 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP9 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  V2[2] = denom * (OM2 * P2[0] * (-cI * (TMP9 * TMP27) + cI * (TMP17 * TMP24))
      + (-cI * (TMP17 * V4[2]) + cI * (V1[2] * TMP27)));
  V2[3] = denom * (OM2 * P2[1] * (-cI * (TMP9 * TMP27) + cI * (TMP17 * TMP24))
      + (-cI * (TMP17 * V4[3]) + cI * (V1[3] * TMP27)));
  V2[4] = denom * (OM2 * P2[2] * (-cI * (TMP9 * TMP27) + cI * (TMP17 * TMP24))
      + (-cI * (TMP17 * V4[4]) + cI * (V1[4] * TMP27)));
  V2[5] = denom * (OM2 * P2[3] * (-cI * (TMP9 * TMP27) + cI * (TMP17 * TMP24))
      + (-cI * (TMP17 * V4[5]) + cI * (V1[5] * TMP27)));
}

