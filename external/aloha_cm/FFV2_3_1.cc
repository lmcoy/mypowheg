// This File is Automatically generated by ALOHA
//
#include "FFV2_3_1.h"
#include "FFV2_1.h"
#include "FFV3_1.h"

void FFV2_3_1(std::complex<double> F2[], std::complex<double> V3[],
              std::complex<double> COUP1, std::complex<double> COUP2,
              std::complex<double> M1, std::complex<double> F1[]) {
    std::complex<double> Ftmp[6];
    FFV2_1(F2, V3, COUP1, M1, F1);
    FFV3_1(F2, V3, COUP2, M1, Ftmp);
    int i = 2;
    while (i < 6) {
        F1[i] = F1[i] + Ftmp[i];
        i++;
    }
}