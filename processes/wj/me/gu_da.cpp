//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.2, 2016-06-10
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//
// 2016 10 13: modified by Lennart Oymanns
//
//==========================================================================

#include "gu_da.h"
#include "../external/aloha_cm/aloha.h"

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g u > mu+ vm d a WEIGHTED<=7 @1
// Process: g c > mu+ vm s a WEIGHTED<=7 @1

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

double GU_DA::Calculate(const Phasespace::Phasespace &ps, int perm[],
                        const Parameters_sm &param) {
    double momenta[nexternal][4];
    for (int i = 0; i < nexternal; i++) {
        momenta[i][0] = ps.Momenta[i].E();
        momenta[i][1] = ps.Momenta[i].PX();
        momenta[i][2] = ps.Momenta[i].PY();
        momenta[i][3] = ps.Momenta[i].PZ();
    }

    // Local variables and constants
    const int ncomb = 64;
    static bool goodhel[ncomb] = {ncomb * false};
    static int ntry = 0;

    // Helicities for the process
    static const int helicities[ncomb][nexternal] = {
        {-1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, 1},
        {-1, -1, -1, -1, 1, -1},  {-1, -1, -1, -1, 1, 1},
        {-1, -1, -1, 1, -1, -1},  {-1, -1, -1, 1, -1, 1},
        {-1, -1, -1, 1, 1, -1},   {-1, -1, -1, 1, 1, 1},
        {-1, -1, 1, -1, -1, -1},  {-1, -1, 1, -1, -1, 1},
        {-1, -1, 1, -1, 1, -1},   {-1, -1, 1, -1, 1, 1},
        {-1, -1, 1, 1, -1, -1},   {-1, -1, 1, 1, -1, 1},
        {-1, -1, 1, 1, 1, -1},    {-1, -1, 1, 1, 1, 1},
        {-1, 1, -1, -1, -1, -1},  {-1, 1, -1, -1, -1, 1},
        {-1, 1, -1, -1, 1, -1},   {-1, 1, -1, -1, 1, 1},
        {-1, 1, -1, 1, -1, -1},   {-1, 1, -1, 1, -1, 1},
        {-1, 1, -1, 1, 1, -1},    {-1, 1, -1, 1, 1, 1},
        {-1, 1, 1, -1, -1, -1},   {-1, 1, 1, -1, -1, 1},
        {-1, 1, 1, -1, 1, -1},    {-1, 1, 1, -1, 1, 1},
        {-1, 1, 1, 1, -1, -1},    {-1, 1, 1, 1, -1, 1},
        {-1, 1, 1, 1, 1, -1},     {-1, 1, 1, 1, 1, 1},
        {1, -1, -1, -1, -1, -1},  {1, -1, -1, -1, -1, 1},
        {1, -1, -1, -1, 1, -1},   {1, -1, -1, -1, 1, 1},
        {1, -1, -1, 1, -1, -1},   {1, -1, -1, 1, -1, 1},
        {1, -1, -1, 1, 1, -1},    {1, -1, -1, 1, 1, 1},
        {1, -1, 1, -1, -1, -1},   {1, -1, 1, -1, -1, 1},
        {1, -1, 1, -1, 1, -1},    {1, -1, 1, -1, 1, 1},
        {1, -1, 1, 1, -1, -1},    {1, -1, 1, 1, -1, 1},
        {1, -1, 1, 1, 1, -1},     {1, -1, 1, 1, 1, 1},
        {1, 1, -1, -1, -1, -1},   {1, 1, -1, -1, -1, 1},
        {1, 1, -1, -1, 1, -1},    {1, 1, -1, -1, 1, 1},
        {1, 1, -1, 1, -1, -1},    {1, 1, -1, 1, -1, 1},
        {1, 1, -1, 1, 1, -1},     {1, 1, -1, 1, 1, 1},
        {1, 1, 1, -1, -1, -1},    {1, 1, 1, -1, -1, 1},
        {1, 1, 1, -1, 1, -1},     {1, 1, 1, -1, 1, 1},
        {1, 1, 1, 1, -1, -1},     {1, 1, 1, 1, -1, 1},
        {1, 1, 1, 1, 1, -1},      {1, 1, 1, 1, 1, 1}};
    // Denominators: spins, colors and identical particles
    const double denominator = 96.0;

    ntry += 1;

    double matrix_element = 0.0;

    // Calculate the matrix element for all helicities
    for (int ihel = 0; ihel < ncomb; ihel++) {
        if (goodhel[ihel] || ntry < 2) {
            double t1 =
                calculate_wavefunctions(perm, helicities[ihel], param, momenta);
            matrix_element += t1;

            // Mirror initial state momenta for mirror process
            int perm0 = perm[0];
            int perm1 = perm[1];
            // Calculate wavefunctions
            double t2 =
                calculate_wavefunctions(perm, helicities[ihel], param, momenta);
            // Mirror back
            perm[0] = perm0;
            perm[1] = perm1;
            // Calculate matrix elements

            double tsum = t1 + t2;

            // Store which helicities give non-zero result
            if (tsum != 0.0) {
                goodhel[ihel] = true;
            }
        }
    }

    return matrix_element / denominator;
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

double GU_DA::calculate_wavefunctions(const int perm[], const int hel[],
                                      const Parameters_sm &pars,
                                      double momenta[][4]) {

    constexpr std::complex<double> ZERO(0.0, 0.0);
    static const int nwavefuncs = 18;
    std::complex<double> w[nwavefuncs][18];

    // Calculate all wavefunctions
    vxxxxx(momenta[perm[0]], mME[0], hel[0], -1, w[0]);
    ixxxxx(momenta[perm[1]], mME[1], hel[1], +1, w[1]);
    ixxxxx(momenta[perm[2]], mME[2], hel[2], -1, w[2]);
    oxxxxx(momenta[perm[3]], mME[3], hel[3], +1, w[3]);
    oxxxxx(momenta[perm[4]], mME[4], hel[4], +1, w[4]);
    vxxxxx(momenta[perm[5]], mME[5], hel[5], +1, w[5]);
    FFV1_2(w[1], w[0], pars.GC_11, ZERO, w[6]);
    FFV2_3(w[2], w[3], pars.GC_100, pars.MuW, w[7]);
    FFV2_3(w[6], w[4], pars.GC_100, pars.MuW, w[8]);
    FFV1_2(w[6], w[5], pars.GC_2, ZERO, w[9]);
    FFV1_1(w[4], w[5], pars.GC_1, ZERO, w[10]);
    FFV1_2(w[2], w[5], pars.GC_3, ZERO, w[11]);
    FFV1_1(w[4], w[0], pars.GC_11, ZERO, w[12]);
    FFV1_2(w[1], w[5], pars.GC_2, ZERO, w[13]);
    FFV2_3(w[1], w[12], pars.GC_100, pars.MuW, w[14]);
    FFV1_1(w[12], w[5], pars.GC_1, ZERO, w[15]);
    FFV1_2(w[13], w[0], pars.GC_11, ZERO, w[16]);
    FFV1_1(w[10], w[0], pars.GC_11, ZERO, w[17]);

    // Calculate all amplitudes
    // Amplitude(s) for diagram number 0
    static const int namplitudes = 10;
    std::complex<double> amp[namplitudes];
    VVV1_0(w[5], w[8], w[7], pars.GC_4, amp[0]);
    FFV2_0(w[9], w[4], w[7], pars.GC_100, amp[1]);
    FFV2_0(w[6], w[10], w[7], pars.GC_100, amp[2]);
    FFV2_0(w[11], w[3], w[8], pars.GC_100, amp[3]);
    FFV2_0(w[13], w[12], w[7], pars.GC_100, amp[4]);
    VVV1_0(w[5], w[14], w[7], pars.GC_4, amp[5]);
    FFV2_0(w[1], w[15], w[7], pars.GC_100, amp[6]);
    FFV2_0(w[11], w[3], w[14], pars.GC_100, amp[7]);
    FFV2_0(w[16], w[4], w[7], pars.GC_100, amp[8]);
    FFV2_0(w[1], w[17], w[7], pars.GC_100, amp[9]);

    // Calculate color flows
    std::complex<double> jamp = +amp[0] + amp[1] + amp[2] + amp[3] + amp[4] +
                                amp[5] + amp[6] + amp[7] + amp[8] + amp[9];

    // Sum and square the color flows to get the matrix element
    constexpr double cf = 4.0;
    double matrix = cf * std::norm(jamp);

    return matrix;
}
