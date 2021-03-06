// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// Metric(1,3)*Metric(2,4) - Metric(1,2)*Metric(3,4)
//
#include "VVVV4_3.h"

void VVVV4_3(std::complex<double> V1[], std::complex<double> V2[],
             std::complex<double> V4[], std::complex<double> COUP,
             std::complex<double> M3, std::complex<double> V3[]) {
    std::complex<double> cI = std::complex<double>(0., 1.);
    std::complex<double> TMP22;
    std::complex<double> TMP10;
    std::complex<double> OM3;
    double P3[4];
    std::complex<double> denom;
    std::complex<double> TMP28;
    std::complex<double> TMP13;
    OM3 = 0.;
    if (M3 != 0.)
        OM3 = 1. / pow(M3, 2);
    V3[0] = +V1[0] + V2[0] + V4[0];
    V3[1] = +V1[1] + V2[1] + V4[1];
    P3[0] = -V3[0].real();
    P3[1] = -V3[1].real();
    P3[2] = -V3[1].imag();
    P3[3] = -V3[0].imag();
    TMP22 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]);
    TMP28 = (P3[0] * V4[2] - P3[1] * V4[3] - P3[2] * V4[4] - P3[3] * V4[5]);
    TMP10 = (P3[0] * V1[2] - P3[1] * V1[3] - P3[2] * V1[4] - P3[3] * V1[5]);
    TMP13 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]);
    denom = COUP / (pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) -
                    pow(P3[3], 2) - pow(M3, 2));
    V3[2] =
        denom * (OM3 * P3[0] * (-cI * (TMP13 * TMP28) + cI * (TMP10 * TMP22)) +
                 (-cI * (V1[2] * TMP22) + cI * (TMP13 * V4[2])));
    V3[3] =
        denom * (OM3 * P3[1] * (-cI * (TMP13 * TMP28) + cI * (TMP10 * TMP22)) +
                 (-cI * (V1[3] * TMP22) + cI * (TMP13 * V4[3])));
    V3[4] =
        denom * (OM3 * P3[2] * (-cI * (TMP13 * TMP28) + cI * (TMP10 * TMP22)) +
                 (-cI * (V1[4] * TMP22) + cI * (TMP13 * V4[4])));
    V3[5] =
        denom * (OM3 * P3[3] * (-cI * (TMP13 * TMP28) + cI * (TMP10 * TMP22)) +
                 (-cI * (V1[5] * TMP22) + cI * (TMP13 * V4[5])));
}
