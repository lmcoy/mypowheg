//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.2, 2016-06-10
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//
// 2016 07 04: modified by Lennart Oymanns
//
//==========================================================================

#include "dxdx_uxdx.h"
#include "../external/aloha_cm/aloha.h"
#include "ffv1p0_3.h"

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: d~ d~ > mu+ vm u~ d~ WEIGHTED<=6 @1
// Process: s~ s~ > mu+ vm c~ s~ WEIGHTED<=6 @1

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

double DXDX_UXDX::Calculate(const Phasespace::Phasespace &ps, int perm[],
                            const Parameters_sm &param) {

    for (int i = 0; i < nexternal; i++) {
        momenta[i][0] = ps.Momenta[i].E();
        momenta[i][1] = ps.Momenta[i].PX();
        momenta[i][2] = ps.Momenta[i].PY();
        momenta[i][3] = ps.Momenta[i].PZ();
    }

    // Local variables and constants
    const int ncomb = 64;
    static bool goodhel[ncomb] = { ncomb * false };
    static int ntry = 0;

    // Helicities for the process
    static const int helicities[ncomb][nexternal] = {
        { -1, -1, -1, -1, -1, -1 }, { -1, -1, -1, -1, -1, 1 },
        { -1, -1, -1, -1, 1, -1 },  { -1, -1, -1, -1, 1, 1 },
        { -1, -1, -1, 1, -1, -1 },  { -1, -1, -1, 1, -1, 1 },
        { -1, -1, -1, 1, 1, -1 },   { -1, -1, -1, 1, 1, 1 },
        { -1, -1, 1, -1, -1, -1 },  { -1, -1, 1, -1, -1, 1 },
        { -1, -1, 1, -1, 1, -1 },   { -1, -1, 1, -1, 1, 1 },
        { -1, -1, 1, 1, -1, -1 },   { -1, -1, 1, 1, -1, 1 },
        { -1, -1, 1, 1, 1, -1 },    { -1, -1, 1, 1, 1, 1 },
        { -1, 1, -1, -1, -1, -1 },  { -1, 1, -1, -1, -1, 1 },
        { -1, 1, -1, -1, 1, -1 },   { -1, 1, -1, -1, 1, 1 },
        { -1, 1, -1, 1, -1, -1 },   { -1, 1, -1, 1, -1, 1 },
        { -1, 1, -1, 1, 1, -1 },    { -1, 1, -1, 1, 1, 1 },
        { -1, 1, 1, -1, -1, -1 },   { -1, 1, 1, -1, -1, 1 },
        { -1, 1, 1, -1, 1, -1 },    { -1, 1, 1, -1, 1, 1 },
        { -1, 1, 1, 1, -1, -1 },    { -1, 1, 1, 1, -1, 1 },
        { -1, 1, 1, 1, 1, -1 },     { -1, 1, 1, 1, 1, 1 },
        { 1, -1, -1, -1, -1, -1 },  { 1, -1, -1, -1, -1, 1 },
        { 1, -1, -1, -1, 1, -1 },   { 1, -1, -1, -1, 1, 1 },
        { 1, -1, -1, 1, -1, -1 },   { 1, -1, -1, 1, -1, 1 },
        { 1, -1, -1, 1, 1, -1 },    { 1, -1, -1, 1, 1, 1 },
        { 1, -1, 1, -1, -1, -1 },   { 1, -1, 1, -1, -1, 1 },
        { 1, -1, 1, -1, 1, -1 },    { 1, -1, 1, -1, 1, 1 },
        { 1, -1, 1, 1, -1, -1 },    { 1, -1, 1, 1, -1, 1 },
        { 1, -1, 1, 1, 1, -1 },     { 1, -1, 1, 1, 1, 1 },
        { 1, 1, -1, -1, -1, -1 },   { 1, 1, -1, -1, -1, 1 },
        { 1, 1, -1, -1, 1, -1 },    { 1, 1, -1, -1, 1, 1 },
        { 1, 1, -1, 1, -1, -1 },    { 1, 1, -1, 1, -1, 1 },
        { 1, 1, -1, 1, 1, -1 },     { 1, 1, -1, 1, 1, 1 },
        { 1, 1, 1, -1, -1, -1 },    { 1, 1, 1, -1, -1, 1 },
        { 1, 1, 1, -1, 1, -1 },     { 1, 1, 1, -1, 1, 1 },
        { 1, 1, 1, 1, -1, -1 },     { 1, 1, 1, 1, -1, 1 },
        { 1, 1, 1, 1, 1, -1 },      { 1, 1, 1, 1, 1, 1 }
    };
    // Denominators: spins, colors and identical particles
    const double denominator = 36;

    ntry += 1;

    double matrix_element = 0.0;

    // Calculate the matrix element for all helicities
    for (int ihel = 0; ihel < ncomb; ihel++) {
        if (goodhel[ihel] || ntry < 2) {
            double t = calculate_wavefunctions(perm, helicities[ihel], param);
            matrix_element += t;

            // Store which helicities give non-zero result
            if (t != 0.) {
                goodhel[ihel] = true;
            }
        }
    }

    return matrix_element / denominator;
}

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

double DXDX_UXDX::calculate_wavefunctions(const int perm[], const int hel[],
                                          const Parameters_sm &pars) {
    constexpr std::complex<double> ZERO(0.0, 0.0);
    constexpr int nwavefuncs = 12;
    std::complex<double> w[nwavefuncs][18];

    // Calculate all wavefunctions
    oxxxxx(momenta[perm[0]], mME[0], hel[0], -1, w[0]);
    oxxxxx(momenta[perm[1]], mME[1], hel[1], -1, w[1]);
    ixxxxx(momenta[perm[2]], mME[2], hel[2], -1, w[2]);
    oxxxxx(momenta[perm[3]], mME[3], hel[3], +1, w[3]);
    ixxxxx(momenta[perm[4]], mME[4], hel[4], -1, w[4]);
    ixxxxx(momenta[perm[5]], mME[5], hel[5], -1, w[5]);
    FFV1P0_3(w[5], w[0], pars.GC_11, ZERO, w[6]);
    FFV2_3(w[2], w[3], pars.GC_100, pars.MuW, w[7]);
    FFV1_2(w[4], w[6], pars.GC_11, ZERO, w[8]);
    FFV2_2(w[4], w[7], pars.GC_100, ZERO, w[9]);
    FFV1P0_3(w[5], w[1], pars.GC_11, ZERO, w[10]);
    FFV1_2(w[4], w[10], pars.GC_11, ZERO, w[11]);

    // Calculate all amplitudes
    // Amplitude(s) for diagram number 0
    constexpr int namplitudes = 4;
    std::complex<double> amp[namplitudes];
    FFV2_0(w[8], w[1], w[7], pars.GC_100, amp[0]);
    FFV1_0(w[9], w[1], w[6], pars.GC_11, amp[1]);
    FFV2_0(w[11], w[0], w[7], pars.GC_100, amp[2]);
    FFV1_0(w[9], w[0], w[10], pars.GC_11, amp[3]);

    // Calculate color flows
    const int ncolor = 2;
    std::complex<double> jamp[ncolor];
    jamp[0] =
        +1. / 2. * (+amp[0] + amp[1] + 1. / 3. * amp[2] + 1. / 3. * amp[3]);
    jamp[1] =
        +1. / 2. * (-1. / 3. * amp[0] - 1. / 3. * amp[1] - amp[2] - amp[3]);

    // Sum and square the color flows to get the matrix element
    constexpr double cf[ncolor][ncolor] = { { 9, 3 }, { 3, 9 } };
    double matrix =
        cf[0][0] * std::norm(jamp[0]) + cf[1][1] * std::norm(jamp[1]) +
        (cf[0][1] + cf[1][0]) * std::real(jamp[0] * std::conj(jamp[1]));

    return matrix;
}
