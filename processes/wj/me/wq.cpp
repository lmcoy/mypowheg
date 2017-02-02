//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.4.0, 2016-05-12
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//
// 2016 06 07: modified by Lennart Oymanns
//
//==========================================================================

#include "wq.h"

#include "../external/aloha_cm/aloha.h"

#include "parameters_sm.h"
#include "phasespace/phasespace.h"

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.
UserProcess::SpinCorrelated Wq::Calculate(const Phasespace::Phasespace &ps,
                                          int perm[],
                                          const Parameters_sm &param) {

    // Local variables and constants
    const int ncomb = 16;
    static bool goodhel[ncomb] = { ncomb * false };
    static int ntry = 0;

    // Helicities for the process
    static const int helicities[ncomb][nexternal] = {
        { 0, -1, -1, -1, -1 }, { 0, -1, -1, -1, 1 }, { 0, -1, -1, 1, -1 },
        { 0, -1, -1, 1, 1 },   { 0, -1, 1, -1, -1 }, { 0, -1, 1, -1, 1 },
        { 0, -1, 1, 1, -1 },   { 0, -1, 1, 1, 1 },   { 0, 1, -1, -1, -1 },
        { 0, 1, -1, -1, 1 },   { 0, 1, -1, 1, -1 },  { 0, 1, -1, 1, 1 },
        { 0, 1, 1, -1, -1 },   { 0, 1, 1, -1, 1 },   { 0, 1, 1, 1, -1 },
        { 0, 1, 1, 1, 1 }
    };
    // Denominators: spins, colors and identical particles
    const double denominators = 96;

    ntry += 1;

    double matrix_element = 0.;
    UserProcess::SpinCorrelated spincorr;

    // Calculate the matrix element for all helicities
    for (int ihel = 0; ihel < ncomb; ihel++) {
        if (goodhel[ihel] || ntry < 2) {
            auto m_plus = matrix(ps, perm, helicities[ihel], param, +1);
            auto m_minus = matrix(ps, perm, helicities[ihel], param, -1);
            double t = std::norm(m_plus[0]) + std::norm(m_minus[0]);
            matrix_element += t;

            for (int mu = 0; mu < 4; mu++) {
                for (int nu = mu; nu < 4; nu++) {
                    std::complex<double> ceps_mu_p = std::conj(m_plus[mu + 1]);
                    std::complex<double> ceps_mu_m = std::conj(m_minus[mu + 1]);
                    std::complex<double> eps_nu_p = m_plus[nu + 1];
                    std::complex<double> eps_nu_m = m_minus[nu + 1];
                    spincorr[mu][nu] +=
                        ceps_mu_p * eps_nu_p * std::norm(m_plus[0]) +
                        ceps_mu_p * eps_nu_m * m_plus[0] *
                            std::conj(m_minus[0]) +
                        ceps_mu_m * eps_nu_p * m_minus[0] *
                            std::conj(m_plus[0]) +
                        ceps_mu_m * eps_nu_m * std::norm(m_minus[0]);
                }
            }

            if (ntry < 2) {
                // Mirror initial state momenta for mirror process
                int p0 = perm[0];
                int p1 = perm[1];
                perm[0] = p1;
                perm[1] = p0;
                // Calculate wavefunctions
                auto mp = matrix(ps, perm, helicities[ihel], param, +1);
                auto mm = matrix(ps, perm, helicities[ihel], param, -1);
                double tsum = t + std::norm(mp[0]) + std::norm(mm[0]);
                // Mirror back
                perm[0] = p0;
                perm[1] = p1;

                // Store which helicities give non-zero result
                if (tsum != 0.) {
                    goodhel[ihel] = true;
                }
            }
        }
    }
    spincorr[1][0] = std::conj(spincorr[0][1]);
    spincorr[2][0] = std::conj(spincorr[0][2]);
    spincorr[2][1] = std::conj(spincorr[1][2]);
    spincorr[3][0] = std::conj(spincorr[0][3]);
    spincorr[3][1] = std::conj(spincorr[1][3]);
    spincorr[3][2] = std::conj(spincorr[2][3]);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            spincorr[i][j] /= denominators;
        }
    }

    return spincorr;
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess
std::array<std::complex<double>, 5>
Wq::matrix(const Phasespace::Phasespace &ps, const int perm[], const int hel[],
           const Parameters_sm &pars, int helicity_gluon) {
    constexpr std::complex<double> ZERO(0.0, 0.0);
    constexpr int nwavefuncs = 8;
    std::complex<double> w[nwavefuncs][18];

    // vector with momenta
    double momenta[nexternal][4];

    for (int i = 0; i < nexternal; i++) {
        momenta[i][0] = ps.Momenta[i].E();
        momenta[i][1] = ps.Momenta[i].PX();
        momenta[i][2] = ps.Momenta[i].PY();
        momenta[i][3] = ps.Momenta[i].PZ();
    }

    // Calculate all wavefunctions
    vxxxxx(momenta[perm[0]], mME[0], helicity_gluon, -1, w[0]);
    ixxxxx(momenta[perm[1]], mME[1], hel[1], +1, w[1]);
    ixxxxx(momenta[perm[2]], mME[2], hel[2], -1, w[2]);
    oxxxxx(momenta[perm[3]], mME[3], hel[3], +1, w[3]);
    oxxxxx(momenta[perm[4]], mME[4], hel[4], +1, w[4]);
    FFV1_2(w[1], w[0], pars.GC_11, ZERO, w[5]);
    FFV2_3(w[2], w[3], pars.GC_100, pars.MuW, w[6]);
    FFV1_1(w[4], w[0], pars.GC_11, ZERO, w[7]);

    // Calculate all amplitudes
    constexpr int namplitudes = 2;
    std::complex<double> amp[namplitudes];
    FFV2_0(w[5], w[4], w[6], pars.GC_100, amp[0]);
    FFV2_0(w[1], w[7], w[6], pars.GC_100, amp[1]);

    // The color matrix;
    constexpr double cf = 4;

    // Calculate color flows
    std::complex<double> jamp = sqrt(cf) * (amp[0] + amp[1]);

    return { { jamp, w[0][2], w[0][3], w[0][4], w[0][5] } };
}
