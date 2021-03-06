// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// ProjM(2,1) - ProjP(2,1)
//
#include "FFS2_3.h"

void FFS2_3(std::complex<double> F1[], std::complex<double> F2[],
            std::complex<double> COUP, std::complex<double> M3,
            std::complex<double> S3[]) {
    std::complex<double> cI = std::complex<double>(0., 1.);
    std::complex<double> denom;
    std::complex<double> TMP1;
    std::complex<double> TMP0;
    double P3[4];
    S3[0] = +F1[0] + F2[0];
    S3[1] = +F1[1] + F2[1];
    P3[0] = -S3[0].real();
    P3[1] = -S3[1].real();
    P3[2] = -S3[1].imag();
    P3[3] = -S3[0].imag();
    TMP1 = (F2[4] * F1[4] + F2[5] * F1[5]);
    TMP0 = (F2[2] * F1[2] + F2[3] * F1[3]);
    denom = COUP / (pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) -
                    pow(P3[3], 2) - pow(M3, 2));
    S3[2] = denom * (-cI * (TMP1) + cI * (TMP0));
}
