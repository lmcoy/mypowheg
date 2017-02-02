#ifndef FFV1P0_3_AHBST62G_H_
#define FFV1P0_3_AHBST62G_H_

static void FFV1P0_3(std::complex<double> F1[], std::complex<double> F2[],
                     std::complex<double> COUP, std::complex<double> M3,
                     std::complex<double> V3[]) {
    std::complex<double> cI = std::complex<double>(0., 1.);
    double P3[4];
    std::complex<double> denom;
    V3[0] = +F1[0] + F2[0];
    V3[1] = +F1[1] + F2[1];
    P3[0] = -V3[0].real();
    P3[1] = -V3[1].real();
    P3[2] = -V3[1].imag();
    P3[3] = -V3[0].imag();
    denom = COUP / (pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) -
                    pow(P3[3], 2) - M3 * M3);
    V3[2] = denom * -cI *
            (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5] * F2[3]);
    V3[3] = denom * -cI *
            (F1[4] * F2[3] + F1[5] * F2[2] - F1[2] * F2[5] - F1[3] * F2[4]);
    V3[4] = denom * -cI * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) +
                           cI * (F1[3] * F2[4] + F1[4] * F2[3]));
    V3[5] = denom * -cI *
            (F1[3] * F2[5] + F1[4] * F2[2] - F1[2] * F2[4] - F1[5] * F2[3]);
}

#endif
