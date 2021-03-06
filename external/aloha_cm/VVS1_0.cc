// This File is Automatically generated by ALOHA
// The process calculated in this file is:
// Metric(1,2)
//
#include "VVS1_0.h"

void VVS1_0(std::complex<double> V1[], std::complex<double> V2[],
            std::complex<double> S3[], std::complex<double> COUP,
            std::complex<double> &vertex) {
    std::complex<double> cI = std::complex<double>(0., 1.);
    std::complex<double> TMP13;
    TMP13 = (V1[2] * V2[2] - V1[3] * V2[3] - V1[4] * V2[4] - V1[5] * V2[5]);
    vertex = COUP * -cI * TMP13 * S3[2];
}
