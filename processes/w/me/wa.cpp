//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.0, 2016-05-12
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//
// 2016 06 07: modified by Lennart Oymanns
//
//==========================================================================

#include "wa.h"

#include "../external/aloha_cm/aloha.h"

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: q q'~ > mu nu a

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

double Wa::Calculate(const Phasespace::Phasespace &ps, int perm[],
                     const Parameters_sm &param) {

    for (int i = 0; i < nexternal; i++) {
        momenta[i][0] = ps.Momenta[i].E();
        momenta[i][1] = ps.Momenta[i].PX();
        momenta[i][2] = ps.Momenta[i].PY();
        momenta[i][3] = ps.Momenta[i].PZ();
    }

    // Local variables and constants
    const int ncomb = 32;
    static bool goodhel[ncomb] = { ncomb * false };
    static int ntry = 0;

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
    const double denominators = 36.0;

    ntry += 1;

    double matrix_element = 0.;

    // Calculate the matrix element for all helicities
    for (int ihel = 0; ihel < ncomb; ihel++) {
        if (goodhel[ihel] || ntry < 2) {
            calculate_wavefunctions(perm, helicities[ihel], param);
            double t0 = matrix();

            matrix_element += t0;
            if (ntry < 2) {
                double tsum = t0;
                // Mirror initial state momenta for mirror process
                int p0 = perm[0];
                int p1 = perm[1];
                perm[0] = p1;
                perm[1] = p0;
                // Calculate wavefunctions
                calculate_wavefunctions(perm, helicities[ihel], param);
                // Mirror back
                perm[0] = p0;
                perm[1] = p1;
                // Calculate matrix elements
                tsum += matrix();

                // Store which helicities give non-zero result
                if (tsum != 0.) {
                    goodhel[ihel] = true;
                }
            }
        }
    }

    return matrix_element / denominators;
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void Wa::calculate_wavefunctions(const int perm[], const int hel[], const Parameters_sm &pars) {
    // Calculate wavefunctions for all processes
    static const std::complex<double> ZERO(0.0, 0.0);

    // Calculate all wavefunctions
    ixxxxx(momenta[perm[0]], mME[0], hel[0], +1, w[0]);
    oxxxxx(momenta[perm[1]], mME[1], hel[1], -1, w[1]);
    ixxxxx(momenta[perm[2]], mME[2], hel[2], -1, w[2]);
    oxxxxx(momenta[perm[3]], mME[3], hel[3], +1, w[3]);
    vxxxxx(momenta[perm[4]], mME[4], hel[4], +1, w[4]);
    FFV2_3(w[0], w[1], pars.GC_100, pars.MuW, w[5]);
    FFV2_3(w[2], w[3], pars.GC_100, pars.MuW, w[6]);
    FFV1_2(w[2], w[4], pars.GC_3, ZERO, w[7]);
    FFV1_2(w[0], w[4], pars.GC_2, ZERO, w[8]);
    FFV1_1(w[1], w[4], pars.GC_1, ZERO, w[9]);

    // Calculate all amplitudes
    // Amplitude(s) for diagram number 0

    // FSR
    VVV1_0(w[4], w[5], w[6], pars.GC_4, amp[0]);
    FFV2_0(w[7], w[3], w[5], pars.GC_100, amp[1]);

    // ISR
    FFV2_0(w[8], w[1], w[6], pars.GC_100, amp[2]);
    FFV2_0(w[0], w[9], w[6], pars.GC_100, amp[3]);
}

double Wa::matrix() {
    // The color matrix;
    static const double cf = 3;

    // Calculate color flows
    std::complex<double> jamp = +amp[0] + amp[1] + amp[2] + amp[3];

    // Sum and square the color flows to get the matrix element
    return cf * std::norm(jamp);
}

