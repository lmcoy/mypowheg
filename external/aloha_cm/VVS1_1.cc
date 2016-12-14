// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// Metric(1,2)
//
#include "VVS1_1.h"

void VVS1_1(std::complex<double> V2[], std::complex<double> S3[],
            std::complex<double> COUP, std::complex<double> M1,
            std::complex<double> V1[]) {
    std::complex<double> cI = std::complex<double>(0., 1.);
    double P1[4];
    std::complex<double> OM1;
    std::complex<double> denom;
    std::complex<double> TMP14;
    OM1 = 0.;
    if (M1 != 0.)
        OM1 = 1. / pow(M1, 2);
    V1[0] = +V2[0] + S3[0];
    V1[1] = +V2[1] + S3[1];
    P1[0] = -V1[0].real();
    P1[1] = -V1[1].real();
    P1[2] = -V1[1].imag();
    P1[3] = -V1[0].imag();
    TMP14 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]);
    denom = COUP / (pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) -
                    pow(P1[3], 2) - pow(M1, 2));
    V1[2] = denom * S3[2] * (-cI * (V2[2]) + cI * (P1[0] * OM1 * TMP14));
    V1[3] = denom * S3[2] * (-cI * (V2[3]) + cI * (P1[1] * OM1 * TMP14));
    V1[4] = denom * S3[2] * (-cI * (V2[4]) + cI * (P1[2] * OM1 * TMP14));
    V1[5] = denom * S3[2] * (-cI * (V2[5]) + cI * (P1[3] * OM1 * TMP14));
}