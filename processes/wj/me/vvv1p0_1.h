#ifndef VVV1P0_1_AJSUHA2U_H
#define VVV1P0_1_AJSUHA2U_H

#include <complex>

static void VVV1P0_1(std::complex<double> V2[], std::complex<double> V3[],
                     std::complex<double> COUP, std::complex<double> M1,
                     std::complex<double> V1[]) {
    std::complex<double> cI = std::complex<double>(0., 1.);
    std::complex<double> denom;
    double P1[4];
    double P2[4];
    std::complex<double> TMP7;
    double P3[4];
    std::complex<double> TMP6;
    std::complex<double> TMP5;
    std::complex<double> TMP4;
    std::complex<double> TMP3;
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
    TMP5 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]);
    TMP4 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]);
    TMP7 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]);
    TMP6 = (P3[0] * V2[2] - P3[1] * V2[3] - P3[2] * V2[4] - P3[3] * V2[5]);
    TMP3 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]);
    denom = COUP / (pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) -
                    pow(P1[3], 2) - M1 * M1);
    V1[2] = denom * (TMP7 * (-cI * (P2[0]) + cI * (P3[0])) +
                     (V2[2] * (-cI * (TMP3) + cI * (TMP4)) +
                      V3[2] * (-cI * (TMP6) + cI * (TMP5))));
    V1[3] = denom * (TMP7 * (-cI * (P2[1]) + cI * (P3[1])) +
                     (V2[3] * (-cI * (TMP3) + cI * (TMP4)) +
                      V3[3] * (-cI * (TMP6) + cI * (TMP5))));
    V1[4] = denom * (TMP7 * (-cI * (P2[2]) + cI * (P3[2])) +
                     (V2[4] * (-cI * (TMP3) + cI * (TMP4)) +
                      V3[4] * (-cI * (TMP6) + cI * (TMP5))));
    V1[5] = denom * (TMP7 * (-cI * (P2[3]) + cI * (P3[3])) +
                     (V2[5] * (-cI * (TMP3) + cI * (TMP4)) +
                      V3[5] * (-cI * (TMP6) + cI * (TMP5))));
}

#endif
