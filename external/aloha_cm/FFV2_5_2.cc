// This File is Automatically generated by ALOHA
//
#include "FFV2_5_2.h"
#include "FFV2_2.h"
#include "FFV5_2.h"

void FFV2_5_2(std::complex<double> F1[], std::complex<double> V3[],
              std::complex<double> COUP1, std::complex<double> COUP2,
              std::complex<double> M2, std::complex<double> F2[]) {
    std::complex<double> Ftmp[6];
    FFV2_2(F1, V3, COUP1, M2, F2);
    FFV5_2(F1, V3, COUP2, M2, Ftmp);
    int i = 2;
    while (i < 6) {
        F2[i] = F2[i] + Ftmp[i];
        i++;
    }
}
