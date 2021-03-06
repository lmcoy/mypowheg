// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// P(3,1)*Metric(1,2) - P(3,2)*Metric(1,2) - P(2,1)*Metric(1,3) +
// P(2,3)*Metric(1,3) + P(1,2)*Metric(2,3) - P(1,3)*Metric(2,3)
// 
#include "VVV1_1.h"

void VVV1_1(std::complex<double> V2[], std::complex<double> V3[], std::complex<double> COUP,
    double M1, double W1, std::complex<double> V1[])
{
  std::complex<double> cI = std::complex<double> (0., 1.); 
  std::complex<double> TMP12; 
  std::complex<double> TMP11; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  std::complex<double> TMP16; 
  std::complex<double> TMP15; 
  std::complex<double> TMP14; 
  std::complex<double> denom; 
  double OM1; 
  std::complex<double> TMP19; 
  std::complex<double> TMP18; 
  OM1 = 0.; 
  if (M1 != 0.)
    OM1 = 1./pow(M1, 2); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V1[0] = +V2[0] + V3[0]; 
  V1[1] = +V2[1] + V3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP19 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP18 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]); 
  TMP15 = (P1[0] * V3[2] - P1[1] * V3[3] - P1[2] * V3[4] - P1[3] * V3[5]); 
  TMP14 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP16 = (P2[0] * V3[2] - P2[1] * V3[3] - P2[2] * V3[4] - P2[3] * V3[5]); 
  TMP11 = (P1[0] * P2[0] - P1[1] * P2[1] - P1[2] * P2[2] - P1[3] * P2[3]); 
  TMP12 = (P1[0] * P3[0] - P1[1] * P3[1] - P1[2] * P3[2] - P1[3] * P3[3]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (OM1 * P1[0] * (TMP19 * (-cI * (TMP12) + cI * (TMP11)) + (-cI
      * (TMP14 * TMP16) + cI * (TMP15 * TMP18))) + (TMP19 * (-cI * (P2[0]) + cI
      * (P3[0])) + (V2[2] * (-cI * (TMP15) + cI * (TMP16)) + V3[2] * (-cI *
      (TMP18) + cI * (TMP14)))));
  V1[3] = denom * (OM1 * P1[1] * (TMP19 * (-cI * (TMP12) + cI * (TMP11)) + (-cI
      * (TMP14 * TMP16) + cI * (TMP15 * TMP18))) + (TMP19 * (-cI * (P2[1]) + cI
      * (P3[1])) + (V2[3] * (-cI * (TMP15) + cI * (TMP16)) + V3[3] * (-cI *
      (TMP18) + cI * (TMP14)))));
  V1[4] = denom * (OM1 * P1[2] * (TMP19 * (-cI * (TMP12) + cI * (TMP11)) + (-cI
      * (TMP14 * TMP16) + cI * (TMP15 * TMP18))) + (TMP19 * (-cI * (P2[2]) + cI
      * (P3[2])) + (V2[4] * (-cI * (TMP15) + cI * (TMP16)) + V3[4] * (-cI *
      (TMP18) + cI * (TMP14)))));
  V1[5] = denom * (OM1 * P1[3] * (TMP19 * (-cI * (TMP12) + cI * (TMP11)) + (-cI
      * (TMP14 * TMP16) + cI * (TMP15 * TMP18))) + (TMP19 * (-cI * (P2[3]) + cI
      * (P3[3])) + (V2[5] * (-cI * (TMP15) + cI * (TMP16)) + V3[5] * (-cI *
      (TMP18) + cI * (TMP14)))));
}

