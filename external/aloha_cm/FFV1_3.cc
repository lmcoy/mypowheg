// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// Gamma(3,2,1)
//
#include "FFV1_3.h"

void FFV1_3(std::complex<double> F1[], std::complex<double> F2[],
            std::complex<double> COUP, std::complex<double> M3,
            std::complex<double> V3[]) {
    std::complex<double> cI = std::complex<double>(0., 1.);
    std::complex<double> denom;
    std::complex<double> OM3;
    double P3[4];
    std::complex<double> TMP3;
    OM3 = 0.;
    if (M3 != 0.)
        OM3 = 1. / pow(M3, 2);
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    P3[0] = -V3[0].real();
    P3[1] = -V3[1].real();
    P3[2] = -V3[1].imag();
    P3[3] = -V3[0].imag();
    TMP3 =
        (F1[2] * (F2[4] * (P3[0] + P3[3]) + F2[5] * (P3[1] + cI * (P3[2]))) +
         (F1[3] * (F2[4] * (P3[1] - cI * (P3[2])) + F2[5] * (P3[0] - P3[3])) +
          (F1[4] * (F2[2] * (P3[0] - P3[3]) - F2[3] * (P3[1] + cI * (P3[2]))) +
           F1[5] *
               (F2[2] * (+cI * (P3[2]) - P3[1]) + F2[3] * (P3[0] + P3[3])))));
    denom = COUP / (pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) -
                    pow(P3[3], 2) - pow(M3, 2));
    V3[2] = denom * -cI * (F2[4] * F1[2] + F2[5] * F1[3] + F2[2] * F1[4] +
                           F2[3] * F1[5] - P3[0] * OM3 * TMP3);
    V3[3] = denom * -cI * (F2[3] * F1[4] + F2[2] * F1[5] - F2[5] * F1[2] -
                           F2[4] * F1[3] - P3[1] * OM3 * TMP3);
    V3[4] = denom * -cI *
            (-cI * (F2[5] * F1[2] + F2[2] * F1[5]) +
             cI * (F2[4] * F1[3] + F2[3] * F1[4]) - P3[2] * OM3 * TMP3);
    V3[5] = denom * -cI * (F2[5] * F1[3] + F2[2] * F1[4] - F2[4] * F1[2] -
                           F2[3] * F1[5] - P3[3] * OM3 * TMP3);
}
