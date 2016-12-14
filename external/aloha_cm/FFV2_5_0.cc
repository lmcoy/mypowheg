// This File is Automatically generated by ALOHA
//
#include "FFV2_5_0.h"
#include "FFV2_0.h"
#include "FFV5_0.h"

void FFV2_5_0(std::complex<double> F1[], std::complex<double> F2[],
              std::complex<double> V3[], std::complex<double> COUP1,
              std::complex<double> COUP2, std::complex<double> &vertex) {
    std::complex<double> tmp;
    FFV2_0(F1, F2, V3, COUP1, vertex);
    FFV5_0(F1, F2, V3, COUP2, tmp);
    vertex = vertex + tmp;
}
