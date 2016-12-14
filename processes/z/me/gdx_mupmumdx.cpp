//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.0.1, 2014-01-20
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "gdx_mupmumdx.h"

#include <complex>

#include "../external/aloha_cm/aloha.h"
#include "../external/aloha_cm/aloha_pow2.h"

using namespace std;

namespace {

void FFV1P0_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    complex<double> M3, complex<double> V3[]) {
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag();
  denom =
      COUP / (pow2(P3[0]) - pow2(P3[1]) - pow2(P3[2]) - pow2(P3[3]) - pow2(M3));
  V3[2] = denom * - cI * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5]
      * F2[3]);
  V3[3] = denom * - cI * (F1[4] * F2[3] + F1[5] * F2[2] - F1[2] * F2[5] - F1[3]
      * F2[4]);
  V3[4] = denom * - cI * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) + cI * (F1[3] *
      F2[4] + F1[4] * F2[3]));
  V3[5] = denom * - cI * (F1[3] * F2[5] + F1[4] * F2[2] - F1[2] * F2[4] - F1[5]
      * F2[3]);
}

} // namespace

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g d~ > mu+ mu- d~ WEIGHTED=5

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

double ME_gdx_mupmumdx::Calculate(const Phasespace::Phasespace &ps, int perm[],
                                  const Parameters_sm &param) {
    for (int i = 0; i < nexternal; i++) {
        momenta[i][0] = ps.Momenta[i].E();
        momenta[i][1] = ps.Momenta[i].PX();
        momenta[i][2] = ps.Momenta[i].PY();
        momenta[i][3] = ps.Momenta[i].PZ();
    }

    // Reset color flows
    jamp2 = 0.;

    // Local variables and constants
    const int ncomb = 32;
    static bool goodhel[ncomb] = { ncomb * false };
    static int ntry = 0, sum_hel = 0, ngood = 0;
    static int igood[ncomb];
    static int jhel;
    double t;
    // Helicities for the process
    static const int helicities[ncomb][nexternal] = { { -1, -1, -1, -1, -1 },
                                                      { -1, -1, -1, -1, 1 },
                                                      { -1, -1, -1, 1, -1 },
                                                      { -1, -1, -1, 1, 1 },
                                                      { -1, -1, 1, -1, -1 },
                                                      { -1, -1, 1, -1, 1 },
                                                      { -1, -1, 1, 1, -1 },
                                                      { -1, -1, 1, 1, 1 },
                                                      { -1, 1, -1, -1, -1 },
                                                      { -1, 1, -1, -1, 1 },
                                                      { -1, 1, -1, 1, -1 },
                                                      { -1, 1, -1, 1, 1 },
                                                      { -1, 1, 1, -1, -1 },
                                                      { -1, 1, 1, -1, 1 },
                                                      { -1, 1, 1, 1, -1 },
                                                      { -1, 1, 1, 1, 1 },
                                                      { 1, -1, -1, -1, -1 },
                                                      { 1, -1, -1, -1, 1 },
                                                      { 1, -1, -1, 1, -1 },
                                                      { 1, -1, -1, 1, 1 },
                                                      { 1, -1, 1, -1, -1 },
                                                      { 1, -1, 1, -1, 1 },
                                                      { 1, -1, 1, 1, -1 },
                                                      { 1, -1, 1, 1, 1 },
                                                      { 1, 1, -1, -1, -1 },
                                                      { 1, 1, -1, -1, 1 },
                                                      { 1, 1, -1, 1, -1 },
                                                      { 1, 1, -1, 1, 1 },
                                                      { 1, 1, 1, -1, -1 },
                                                      { 1, 1, 1, -1, 1 },
                                                      { 1, 1, 1, 1, -1 },
                                                      { 1, 1, 1, 1, 1 } };
    // Denominators: spins, colors and identical particles
    const int denominators = 96;

    ntry = ntry + 1;

    // Reset the matrix elements
    matrix_element = 0.;

    if (sum_hel == 0 || ntry < 10) {
        // Calculate the matrix element for all helicities
        for (int ihel = 0; ihel < ncomb; ihel++) {
            if (goodhel[ihel] || ntry < 2) {
                calculate_wavefunctions(perm, helicities[ihel], param);
                t = matrix_gdx_mupmumdx();

                double tsum = 0;

                matrix_element += t;
                tsum += t;

                // Store which helicities give non-zero result
                if (tsum != 0. && !goodhel[ihel]) {
                    goodhel[ihel] = true;
                    ngood++;
                    igood[ngood] = ihel;
                }
            }
        }
        jhel = 0;
        sum_hel = min(sum_hel, ngood);
    } else {
        // Only use the "good" helicities
        for (int j = 0; j < sum_hel; j++) {
            jhel++;
            if (jhel >= ngood)
                jhel = 0;
            double hwgt = double(ngood) / double(sum_hel);
            int ihel = igood[jhel];
            calculate_wavefunctions(perm, helicities[ihel], param);
            t = matrix_gdx_mupmumdx();

            matrix_element += t * hwgt;
        }
    }

    matrix_element /= denominators;

    return matrix_element;
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void ME_gdx_mupmumdx::calculate_wavefunctions(const int perm[], const int hel[],
                                              const Parameters_sm &pars) {

    // Calculate all wavefunctions
    vxxxxx(momenta[perm[0]], mME[0], hel[0], -1, w[0]);
    oxxxxx(momenta[perm[1]], mME[1], hel[1], -1, w[1]);
    ixxxxx(momenta[perm[2]], mME[2], hel[2], -1, w[2]);
    oxxxxx(momenta[perm[3]], mME[3], hel[3], +1, w[3]);
    ixxxxx(momenta[perm[4]], mME[4], hel[4], -1, w[4]);
    FFV1_1(w[1], w[0], pars.GC_11, std::complex<double>(0.0,0.0), w[5]);
    FFV1P0_3(w[2], w[3], pars.GC_3, std::complex<double>(0.0,0.0), w[6]);
    FFV2_4_3(w[2], w[3], pars.GC_50, pars.GC_59, pars.MuZ, w[7]);
    FFV1_2(w[4], w[0], pars.GC_11, std::complex<double>(0.0,0.0), w[8]);

    // Calculate all amplitudes
    // Amplitude(s) for diagram number 0
    FFV1_0(w[4], w[5], w[6], pars.GC_1, amp[0]);
    FFV2_3_0(w[4], w[5], w[7], pars.GC_50, pars.GC_58, amp[1]);
    FFV1_0(w[8], w[1], w[6], pars.GC_1, amp[2]);
    FFV2_3_0(w[8], w[1], w[7], pars.GC_50, pars.GC_58, amp[3]);
}
double ME_gdx_mupmumdx::matrix_gdx_mupmumdx() {

    std::complex<double> ztemp;
    std::complex<double> jamp;
    // The color matrix;

    static const double cf = 4;

    // Calculate color flows
    jamp = -amp[0] - amp[1] - amp[2] - amp[3];

    // Sum and square the color flows to get the matrix element
    double matrix = 0;

    ztemp = cf * jamp;
    matrix = real(ztemp * conj(jamp));

    // Store the leading color flows for choice of color

    jamp2 += real(jamp * conj(jamp));

    return matrix;
}

