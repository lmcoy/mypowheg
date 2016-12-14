// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// Metric(1,3)*Metric(2,4) - Metric(1,2)*Metric(3,4)
// 
#include "VVVV4_1.h"

void VVVV4_1(std::complex<double> V2[], std::complex<double> V3[], std::complex<double> V4[],
    std::complex<double> COUP, double M1, double W1, std::complex<double> V1[])
{
  std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> denom; 
  double P1[4]; 
  std::complex<double> TMP22; 
  std::complex<double> TMP15; 
  std::complex<double> TMP14; 
  std::complex<double> TMP27; 
  double OM1; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./pow(M1, 2); 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP15 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP14 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP27 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  TMP22 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (OM1 * P1[0] * (-cI * (TMP14 * TMP27) + cI * (TMP15 * TMP22))
      + (-cI * (V3[2] * TMP22) + cI * (V2[2] * TMP27)));
  V1[3] = denom * (OM1 * P1[1] * (-cI * (TMP14 * TMP27) + cI * (TMP15 * TMP22))
      + (-cI * (V3[3] * TMP22) + cI * (V2[3] * TMP27)));
  V1[4] = denom * (OM1 * P1[2] * (-cI * (TMP14 * TMP27) + cI * (TMP15 * TMP22))
      + (-cI * (V3[4] * TMP22) + cI * (V2[4] * TMP27)));
  V1[5] = denom * (OM1 * P1[3] * (-cI * (TMP14 * TMP27) + cI * (TMP15 * TMP22))
      + (-cI * (V3[5] * TMP22) + cI * (V2[5] * TMP27)));
}

