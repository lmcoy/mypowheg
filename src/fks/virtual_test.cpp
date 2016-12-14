#include "gtest/gtest.h"

#include "fks/scales.h"
#include "math/fourmomentum.h"

#include "fks/virtual.cpp"

TEST(QED_Q, UUX_MUMMUP) {
    // test u ux -> mu+ mu- with PS in partonic CMS
    double sqrts = 500.0;
    double sqrtsh = sqrts / 2.0;
    Math::FourMomentum momenta[4];
    momenta[0].Set(sqrtsh, 0.0, 0.0, sqrtsh);
    momenta[1].Set(sqrtsh, 0.0, 0.0, -sqrtsh);
    momenta[2].Set(sqrtsh, sqrtsh / sqrt(2), -sqrtsh / sqrt(2), 0.0);
    momenta[3].Set(sqrtsh, -sqrtsh / sqrt(2), sqrtsh / sqrt(2), 0.0);

    // u u~ -> mu+ mu-
    int pdgs[] = { 2, -2, 13, -13 };

    // use scales such that the logs in Q are 1
    FKS::Scales scales;
    scales.muF = sqrts;
    scales.Q2 = sqrts * sqrts / exp(1.0);

    double Q = Qfin<Type::QED>(4, pdgs, momenta, sqrts, scales);

    EXPECT_DOUBLE_EQ(Q, 2.0 * (5.0 - 2.0 / 3.0 * M_PI * M_PI) - 4.0 / 3.0);
}

TEST(QED_Q, SU_CD) {
    // test s u -> c d with PS in partonic CMS
    double sqrts = 500.0;
    double sqrtsh = sqrts / 2.0;
    Math::FourMomentum momenta[4];
    momenta[0].Set(sqrtsh, 0.0, 0.0, sqrtsh);
    momenta[1].Set(sqrtsh, 0.0, 0.0, -sqrtsh);
    momenta[2].Set(sqrtsh, sqrtsh / sqrt(2), -sqrtsh / sqrt(2), 0.0);
    momenta[3].Set(sqrtsh, -sqrtsh / sqrt(2), sqrtsh / sqrt(2), 0.0);

    // s u -> c d
    int pdgs[] = { 3, 2, 4, 1 };

    FKS::Scales scales;
    scales.muF = sqrts;
    scales.Q2 = 2.0 * sqrts * sqrts;

    double Q = Qfin<Type::QED>(4, pdgs, momenta, sqrts, scales);

    EXPECT_DOUBLE_EQ(Q, 5.0 / 9.0 * (13.0 / 2.0 - 2.0 * M_PI * M_PI / 3.0) +
                            5.0 / 3.0 * log(2.0));
}

TEST(I,TwoParticlePS) {
    double sqrts = 500.0;
    double sqrtsh = sqrts / 2.0;
    Math::FourMomentum momenta[4];
    momenta[0].Set(sqrtsh, 0.0, 0.0, sqrtsh);
    momenta[1].Set(sqrtsh, 0.0, 0.0, -sqrtsh);
    momenta[2].Set(sqrtsh, sqrtsh / sqrt(2), -sqrtsh / sqrt(2), 0.0);
    momenta[3].Set(sqrtsh, -sqrtsh / sqrt(2), sqrtsh / sqrt(2), 0.0);

    double log_s_over_Q2 = log(0.5);

    {
        double I = FiniteI(momenta[0], momenta[1], log_s_over_Q2);
        EXPECT_DOUBLE_EQ(I, 0.5 * log(0.5) * log(0.5) - M_PI * M_PI / 6.0);
    }
    {
        double I = FiniteI(momenta[2], momenta[3], log_s_over_Q2);
        EXPECT_DOUBLE_EQ(I, 0.5 * log(0.5) * log(0.5) - M_PI * M_PI / 6.0);
    }
    {
        double I = FiniteI(momenta[0], momenta[2], log_s_over_Q2);
        EXPECT_DOUBLE_EQ(I, 1.5 * log(2.0) * log(2.0) - M_PI * M_PI / 12.0);
    }
    {
        double I = FiniteI(momenta[0], momenta[3], log_s_over_Q2);
        EXPECT_DOUBLE_EQ(I, 1.5 * log(2.0) * log(2.0) - M_PI * M_PI / 12.0);
    }
    {
        double I = FiniteI(momenta[1], momenta[2], log_s_over_Q2);
        EXPECT_DOUBLE_EQ(I, 1.5 * log(2.0) * log(2.0) - M_PI * M_PI / 12.0);
    }
    {
        double I = FiniteI(momenta[1], momenta[3], log_s_over_Q2);
        EXPECT_DOUBLE_EQ(I, 1.5 * log(2.0) * log(2.0) - M_PI * M_PI / 12.0);
    }
}

    
TEST(I,TwoParticlePS_2) {
    double sqrts = 500.0;
    double sqrtsh = sqrts / 2.0;
    Math::FourMomentum momenta[2];
    momenta[0].Set(sqrtsh, 0.0, 0.0, sqrtsh);
    momenta[1].Set(sqrtsh, sqrtsh / sqrt(2), 0.0, -sqrtsh / sqrt(2));

    double log_s_over_Q2 = log(0.5);

    double I = FiniteI(momenta[0], momenta[1], log_s_over_Q2);

    double log_arg = 0.5*(1.0 + 1.0/sqrt(2.0));
    double Log = log(log_arg);

    double dilog = 1.1885451988366879863;
    double pre = 0.5 * pow(log_s_over_Q2, 2) + log_s_over_Q2 * Log +
                 0.5 * Log * Log - Log * log(1.0 - log_arg) - dilog;
    EXPECT_DOUBLE_EQ(I, pre);
}

TEST( Eps2Pole, QED ) {
    int pdgs_dy[4] = { -2, 2, -13, 13 };
    double eps2 = FKS::QED::Eps2Pole(4, pdgs_dy, 1.0);
    // sum of squared charges
    EXPECT_DOUBLE_EQ((2.0 + 8.0 / 9.0), eps2);

    int pdgs_dy2[4] = { 2, -2, -13, 13 };
    eps2 = FKS::QED::Eps2Pole(4, pdgs_dy2, 1.0);
    // sum of squared charges
    EXPECT_DOUBLE_EQ((2.0 + 8.0 / 9.0), eps2);

    int pdgs_dy3[4] = { -1, 1, -13, 13 };
    eps2 = FKS::QED::Eps2Pole(4, pdgs_dy3, 1.0);
    // sum of squared charges
    EXPECT_DOUBLE_EQ((2.0 + 2.0 / 9.0), eps2);

    int pdgs_dy4[4] = { 1, -1, -13, 13 };
    eps2 = FKS::QED::Eps2Pole(4, pdgs_dy4, 1.0);
    // sum of squared charges
    EXPECT_DOUBLE_EQ((2.0 + 2.0 / 9.0), eps2);
}

TEST(Eps2Pole, QCD) {
    int pdgs_dy[4] = { -2, 2, -13, 13 };
    double eps2 = FKS::QCD::Eps2Pole(4, pdgs_dy, 1.0);
    EXPECT_DOUBLE_EQ(2.0 * (4.0 / 3.0), eps2);

    int pdgs_dy2[4] = { -1, 1, -13, 13 };
    eps2 = FKS::QCD::Eps2Pole(4, pdgs_dy2, 1.0);
    EXPECT_DOUBLE_EQ(2.0 * (4.0 / 3.0), eps2);
}

TEST(EpsPole, QCD_DY) {
    int pdgs[] = {-2, 2, -13, 13};
    double sqrts = 500.0;
    double sqrtsh = sqrts / 2.0;
    Math::FourMomentum momenta[4];
    momenta[0].Set(sqrtsh, 0.0, 0.0, sqrtsh);
    momenta[1].Set(sqrtsh, 0.0, 0.0, -sqrtsh);
    momenta[2].Set(sqrtsh, sqrtsh / sqrt(2), -sqrtsh / sqrt(2), 0.0);
    momenta[3].Set(sqrtsh, -sqrtsh / sqrt(2), sqrtsh / sqrt(2), 0.0);

    double born = 1.23;
    Util::Matrix2 colcorr(4, 0.0);
    colcorr.Set(0, 1, (4.0 / 3.0) * born);
    colcorr.Set(1, 0, (4.0 / 3.0) * born);

    double Q2 = 91.188 * 91.188;
    double epspole = FKS::QCD::EpsPole(4, pdgs, momenta, born, colcorr, Q2);
    // analytical result for the virtual 1/eps Pole of DY QCD
    // CF * (-3 + 2 * log(s/Q^2)) * born (c.f. 0709.2092 (7.190))
    // The real pole which is returned by EpsPole cancels this pole.
    double V_pole = born * (-4.0 + 8.0 / 3.0 * log(sqrts * sqrts / Q2));
    EXPECT_NEAR(epspole + V_pole, 0.0, 1e-15);
}
