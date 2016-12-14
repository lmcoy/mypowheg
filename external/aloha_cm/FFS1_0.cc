// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// ProjM(2,1)
//
#include "FFS1_0.h"

void FFS1_0(std::complex<double> F1[], std::complex<double> F2[],
            std::complex<double> S3[], std::complex<double> COUP,
            std::complex<double> &vertex) {
    std::complex<double> cI = std::complex<double>(0., 1.);
    std::complex<double> TMP0;
    TMP0 = (F2[2] * F1[2] + F2[3] * F1[3]);
    vertex = COUP * -cI * TMP0 * S3[2];
}
