#include "gtest/gtest.h"

#include "fks/limits.h"
#include "math/fourmomentum.h"
#include "process/matrixelement.h"

double diff_in_perc(double a, double b) { return 100.0 * (a - b) / b; }

TEST(QEDCollinearLimitFSR, UUX_MUPMUMA) {
    // matrix elment u u~ -> mu+ mu-
    //
    // phase space:
    // (20.493902,  0.0000000,  0.0000000,  20.493902)
    // (20.493902,  0.0000000,  0.0000000, -20.493902)
    // (20.493902, -16.579915, -12.046013,  0.0000000)
    // (20.493902,  16.579915,  12.046013, -0.0000000)
    double bornme = 0.001345504968737048;

    double s = 1680.0000767424158;
    Phasespace::Phasespace fake_ps;
    fake_ps.X1 = 1.0;
    fake_ps.X2 = 1.0;
    fake_ps.S = s;
    fake_ps.N = 2;
    fake_ps.Momenta[3].SetE(20.493902);
    double xi = 0.9;
    double limit =
        FKS::QED::CollinearLimitFSR(13, 13, fake_ps, 3, xi, 0.0, 1.0 / 132.507,
                                    bornme, UserProcess::SpinCorrelated());

    double g_collinear = 3.068399029586553e-06;
    EXPECT_NEAR(limit / g_collinear, 1.0, 1e-4);
}

TEST(QEDCollinearLimitISR, UUX_MUPMUMA_1) {
    // matrix elment u u~ -> mu+ mu-
    //
    // phase space:
    // ( 15.6525,        0,        0,  15.6525)
    // ( 15.6525,        0,        0, -15.6525)
    // ( 15.6525,   3.8695,  11.9091,  9.39149)
    // ( 15.6525,  -3.8695, -11.9091, -9.39149)
    double bornme = 0.00204141958196957;

    double s = 979.9999999999999;
    double xi = 0.2;
    int bornpdgs[] = {2, -2, -13, 13};
    int realpdgs[] = {2, -2, -13, 13, 22};
    double limit = FKS::QED::CollinearLimitISR(4, realpdgs, bornpdgs, xi, 1,
                                               0.0, s, 1.0 / 132.507, bornme,
                                               UserProcess::SpinCorrelated());

    double g_collinear = 1.151999135097969e-06;
    EXPECT_NEAR(limit / g_collinear, 1.0, 1e-4);

    // second ps point
    // phase space:
    // ( 28.9828,        0,        0,  28.9828)
    // ( 28.9828,        0,        0, -28.9828)
    // ( 28.9828,   -21.49, -15.6134, -11.5931)
    // ( 28.9828,    21.49,  15.6134,  11.5931)
    bornme = 0.0009841120664355969;
    s = 3360.0;
    xi = 0.6;
    limit = FKS::QED::CollinearLimitISR(4, realpdgs, bornpdgs, xi, 1, 0.0, s,
                                        1.0 / 132.507, bornme,
                                        UserProcess::SpinCorrelated());
    g_collinear = 1.145562470886872e-07;
    EXPECT_NEAR(limit / g_collinear, 1.0, 1e-4);
}

TEST(QEDCollinearLimitISR, UUX_MUPMUMA_2) {
    // matrix elment u u~ -> mu+ mu-
    //
    // y -> -1
    //
    // phase space:
    // ( 15.6525,        0,        0,  15.6525)
    // ( 15.6525,        0,        0, -15.6525)
    // ( 15.6525,   3.8695,  11.9091,  9.39149)
    // ( 15.6525,  -3.8695, -11.9091, -9.39149)
    double bornme = 0.00204141958196957;

    double s = 979.9999999999999;
    double xi = 0.2;
    int bornpdgs[] = {2, -2, -13, 13};
    int realpdgs[] = {2, -2, -13, 13, 22};
    double limit = FKS::QED::CollinearLimitISR(4, realpdgs, bornpdgs, xi, -1,
                                               0.0, s, 1.0 / 132.507, bornme,
                                               UserProcess::SpinCorrelated());

    double g_collinear = 1.151999135097969e-06;
    EXPECT_NEAR(limit / g_collinear, 1.0, 1e-4);

    // second ps point
    // phase space:
    // ( 28.9828,        0,        0,  28.9828)
    // ( 28.9828,        0,        0, -28.9828)
    // ( 28.9828,   -21.49, -15.6134, -11.5931)
    // ( 28.9828,    21.49,  15.6134,  11.5931)
    bornme = 0.0009841120664355969;
    s = 3360.0;
    xi = 0.6;
    limit = FKS::QED::CollinearLimitISR(4, realpdgs, bornpdgs, xi, -1, 0.0, s,
                                        1.0 / 132.507, bornme,
                                        UserProcess::SpinCorrelated());
    g_collinear = 1.145562470886872e-07;
    EXPECT_NEAR(limit / g_collinear, 1.0, 1e-4);
}

TEST(QEDSoftLimit, UUX_MUPMUMA_FSR) {
    // matrix element u u~ > mu+ mu-

    // phase space:
    Math::FourMomentum born_ps[4];
    born_ps[0].Set(20.493902, 0.000000, 0.000000, 20.493902);
    born_ps[1].Set(20.493902, 0.000000, 0.000000, -20.493902);
    born_ps[2].Set(20.493902, -15.195743, -11.040353, -8.197561);
    born_ps[3].Set(20.493902, 15.195743, 11.040353, 8.197561);
    Util::Matrix2 bornme(1, 0.001271991648686304);

    double s = 1680.0000767424158;
    double y = 0.8;
    double phi = 0.6;

    int pdg[] = {2, -2, -13, 13};
    double limit = FKS::QED::SoftLimit(4, born_ps, pdg, 22, s, 2, 1.0 / 132.507,
                                       bornme, y, phi);

    double g_soft = 4.280439223725195e-07;
    EXPECT_NEAR(limit / g_soft, 1.0, 1e-5);
}

TEST(QEDSoftLimit, UUX_MUPMUMA_ISR) {
    // matrix element u u~ > mu+ mu-

    // phase space:
    Math::FourMomentum born_ps[4];
    born_ps[0].Set(15.65247584, 0, 0, 15.65247584);
    born_ps[1].Set(15.65247584, 0, 0, -15.65247584);
    born_ps[2].Set(15.65247584, 3.869504832, 11.90911132, 9.391485505);
    born_ps[3].Set(15.65247584, -3.869504832, -11.90911132, -9.391485505);
    Util::Matrix2 bornme(1, 0.00204141958196957);

    double s = 979.9999999999999;
    double y = 0.2;
    double phi = 0.2;

    int pdg[] = {2, -2, -13, 13};
    double limit = FKS::QED::SoftLimit(4, born_ps, pdg, 22, s, 0, 1.0 / 132.507,
                                       bornme, y, phi);

    double g_soft = 2.656416321608588e-06;
    EXPECT_NEAR(limit / g_soft, 1.0, 1e-6);
}

TEST(QEDSoftCollinearLimitFSR, UUX_MUPMUMA) {
    // matrix elment u u~ -> mu+ mu-
    //
    // phase space:
    // (15.652476,  0.000000,  0.000000,  15.652476)
    // (15.652476,  0.000000,  0.000000, -15.652476)
    // (15.652476,  3.869505,  11.909111,  9.391486)
    // (15.652476, -3.869505, -11.909111, -9.391486)
    double bornme = 0.00204141958196957;

    double s = 979.9999999999999;
    Phasespace::Phasespace fake_ps;
    fake_ps.X1 = 1.0;
    fake_ps.X2 = 1.0;
    fake_ps.S = s;
    fake_ps.N = 2;

    double limit =
        FKS::QED::SoftCollinearLimitFSR(13, 13, fake_ps, 3, 0.0, 1.0 / 132.507,
                                        bornme, UserProcess::SpinCorrelated());

    double g_collinear = 1.579386188279321e-06;
    EXPECT_NEAR(limit / g_collinear, 1.0, 1e-3);
}

TEST(QEDSoftCollinearLimitISR, UUX_MUPMUMA) {
    // matrix element u u~ > mu+ mu-
    //
    // y -> 1 (y = 0.999999), xi -> 0 (xi = 1e-7)
    //
    // phase space:
    // ( 28.9828,        0,        0,  28.9828)
    // ( 28.9828,        0,        0, -28.9828)
    // ( 28.9828,   -21.49, -15.6134, -11.5931)
    // ( 28.9828,    21.49,  15.6134,  11.5931)
    double bornme = 0.0009841120664355969;
    int bornpdgs[] = {2, -2, -13, 13};
    int realpdgs[] = {2, -2, -13, 13, 22};
    double s = 3360.0;
    double limit = FKS::QED::SoftCollinearLimitISR(
        4, realpdgs, bornpdgs, 1, 0.0, s, 1.0 / 132.507, bornme,
        UserProcess::SpinCorrelated());
    double g_collinear = 1.97399703591074e-07;

    EXPECT_NEAR(limit / g_collinear, 1.0, 1e-3);
}

TEST(Consistency, SoftAndSoftColl_ISR) {
    // matrix element u u~ > mu+ mu-

    // phase space:
    Math::FourMomentum born_ps[4];
    born_ps[0].Set(15.65247584, 0, 0, 15.65247584);
    born_ps[1].Set(15.65247584, 0, 0, -15.65247584);
    born_ps[2].Set(15.65247584, 3.869504832, 11.90911132, 9.391485505);
    born_ps[3].Set(15.65247584, -3.869504832, -11.90911132, -9.391485505);
    double bornme = 0.00204141958196957;

    double s = 979.9999999999999;
    double y = 0.999999999999;
    double phi = 0.6;

    int pdg[] = {2, -2, 13, -13};
    int realpdgs[] = {2, -2, -13, 13, 22};
    // ISR
    double limit_s =
        FKS::QED::SoftLimit(4, born_ps, pdg, 22, s, 0, 1.0 / 132.507,
                            Util::Matrix2(1, bornme), y, phi);
    double limit_sc = FKS::QED::SoftCollinearLimitISR(
        4, realpdgs, pdg, 1, 0.0, s, 1.0 / 132.507, bornme,
        UserProcess::SpinCorrelated());
    EXPECT_NEAR(limit_sc / limit_s, 1.0, 2e-4);
}

TEST(Consistency, SoftAndSoftColl_FSR) {
    // matrix element u u~ > mu+ mu-

    // phase space:
    Math::FourMomentum born_ps[4];
    born_ps[0].Set(15.65247584, 0, 0, 15.65247584);
    born_ps[1].Set(15.65247584, 0, 0, -15.65247584);
    born_ps[2].Set(15.65247584, 3.869504832, 11.90911132, 9.391485505);
    born_ps[3].Set(15.65247584, -3.869504832, -11.90911132, -9.391485505);
    double bornme = 0.00204141958196957;

    double s = 979.9999999999999;
    double y = 0.999999;
    double phi = 0.6;

    Phasespace::Phasespace fake_ps;
    fake_ps.X1 = 1.0;
    fake_ps.X2 = 1.0;
    fake_ps.S = s;
    fake_ps.N = 2;

    int pdg[] = {2, -2, 13, -13};
    double limit_s_2 =
        FKS::QED::SoftLimit(4, born_ps, pdg, 22, s, 2, 1.0 / 132.507,
                            Util::Matrix2(1, bornme), y, phi);
    double limit_s_3 =
        FKS::QED::SoftLimit(4, born_ps, pdg, 22, s, 3, 1.0 / 132.507,
                            Util::Matrix2(1, bornme), y, phi);
    double limit_sc =
        FKS::QED::SoftCollinearLimitFSR(13, 13, fake_ps, 3, phi, 1.0 / 132.507,
                                        bornme, UserProcess::SpinCorrelated());
    EXPECT_NEAR(diff_in_perc(limit_sc, limit_s_2), 0.0, 0.3)
        << "more than 0.3 % difference";
    EXPECT_NEAR(diff_in_perc(limit_sc, limit_s_3), 0.0, 0.3)
        << "more than 0.3 % difference";
}

TEST(QCDSoft, UXU_MUPMUMG) {
    // phase space:
    Math::FourMomentum born_ps[4];
    born_ps[0].Set(57.43015617, 0, 0, 57.43015617);
    born_ps[1].Set(57.43015617, 0, 0, -57.43015617);
    born_ps[2].Set(57.43015617, -13.14151116, -41.83464177, 37.08620054);
    born_ps[3].Set(57.43015617, 13.14151116, 41.83464177, -37.08620054);

    Util::Matrix2 colorcor(4);
    colorcor.Set(0, 1, 0.0167499194703149);
    colorcor.Set(1, 0, 0.0167499194703149);

    double y = -0.654882011935115;
    double phi = 5.6839306132448;
    double s = 13192.8913521231;
    double alpha_s = 0.129783;
    int pdgs[] = {-2, 2, -13, 13};
    double limit = FKS::QCD::SoftLimit(4, born_ps, pdgs, 21, s, 0, alpha_s,
                                       colorcor, y, phi);

    EXPECT_NEAR(limit, 3.31299204483349e-05, 1e-7);
}

TEST(QCDSoft, GU_MUPMUMU) {
    // phase space:
    Math::FourMomentum born_ps[4];
    born_ps[0].Set(57.43015617, 0, 0, 57.43015617);
    born_ps[1].Set(57.43015617, 0, 0, -57.43015617);
    born_ps[2].Set(57.43015617, -13.14151116, -41.83464177, 37.08620054);
    born_ps[3].Set(57.43015617, 13.14151116, 41.83464177, -37.08620054);

    Util::Matrix2 colorcor(4);
    colorcor.Set(0, 1, 0.0167499194703149);
    colorcor.Set(1, 0, 0.0167499194703149);

    double y = -0.654882011935115;
    double phi = 5.6839306132448;
    double s = 13192.8913521231;
    double alpha_s = 0.129783;
    int pdgs[] = {21, 2, -13, 13};
    double limit = FKS::QCD::SoftLimit(4, born_ps, pdgs, 2, s, 0, alpha_s,
                                       colorcor, y, phi);

    EXPECT_NEAR(limit, 0.0, 1e-7);
}

TEST(QCDSoft, UXG_MUPMUMUX) {
    // phase space:
    Math::FourMomentum born_ps[4];
    born_ps[0].Set(57.43015617, 0, 0, 57.43015617);
    born_ps[1].Set(57.43015617, 0, 0, -57.43015617);
    born_ps[2].Set(57.43015617, -13.14151116, -41.83464177, 37.08620054);
    born_ps[3].Set(57.43015617, 13.14151116, 41.83464177, -37.08620054);

    Util::Matrix2 colorcor(4);
    colorcor.Set(0, 1, 0.0167499194703149);
    colorcor.Set(1, 0, 0.0167499194703149);

    double y = -0.654882011935115;
    double phi = 5.6839306132448;
    double s = 13192.8913521231;
    double alpha_s = 0.129783;
    int pdgs[] = {-2, 21, -13, 13};
    double limit = FKS::QCD::SoftLimit(4, born_ps, pdgs, -2, s, 0, alpha_s,
                                       colorcor, y, phi);

    EXPECT_NEAR(limit, 0.0, 1e-7);
}

TEST(QCDCollinearISR1, UXU_MUPMUMG) {
    // y == 1
    // phase space:
    // (    46.30332265,               0,               0,     46.30332265)
    // (    46.30332265,               0,               0,    -46.30332265)
    // (    46.30332265,    -42.46700328,     12.05383173,     13.97341981)
    // (    46.30332265,     42.46700328,    -12.05383173,    -13.97341981)
    double bornme = 0.35916345530503124284;

    double s = 8575.99075378023;
    double xi = 0.285088492732579;
    double alpha_s = 0.129783;

    int bornpdgs[] = {-2, 2, -13, 13};
    int realpdgs[] = {-2, 2, -13, 13, 21};
    double limit = FKS::QCD::CollinearLimitISR(4, realpdgs, bornpdgs, xi, 1,
                                               0.0, s, alpha_s, bornme,
                                               UserProcess::SpinCorrelated());

    double pre = 0.00110092313279678;
    EXPECT_NEAR(limit / pre, 1.0, 1e-5);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    xi = 0.134784312759763;
    limit = FKS::QCD::CollinearLimitISR(4, realpdgs, bornpdgs, xi, 1, 0.0, s,
                                        alpha_s, bornme,
                                        UserProcess::SpinCorrelated());
    pre = 0.000342406169151627;
    EXPECT_NEAR(limit / pre, 1.0, 1e-5);
}

TEST(QCDCollinearISR2, UXU_MUPMUMG) {
    // y == -1
    // phase space:
    // (    46.30332265,               0,               0,     46.30332265)
    // (    46.30332265,               0,               0,    -46.30332265)
    // (    46.30332265,    -42.46700328,     12.05383173,     13.97341981)
    // (    46.30332265,     42.46700328,    -12.05383173,    -13.97341981)
    double bornme = 0.35916345530503124284;

    double s = 8575.99075378023;
    double xi = 0.626275281886555;
    double alpha_s = 0.129783;

    int bornpdgs[] = {-2, 2, -13, 13};
    int realpdgs[] = {-2, 2, -13, 13, 21};
    double limit = FKS::QCD::CollinearLimitISR(4, realpdgs, bornpdgs, xi, -1,
                                               0.0, s, alpha_s, bornme,
                                               UserProcess::SpinCorrelated());

    double pre = 0.000830316011077274;
    EXPECT_NEAR(limit / pre, 1.0, 1e-5);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    xi = 0.896328949986896;
    limit = FKS::QCD::CollinearLimitISR(4, realpdgs, bornpdgs, xi, -1, 0.0, s,
                                        alpha_s, bornme,
                                        UserProcess::SpinCorrelated());
    pre = 0.00019792211054073;
    EXPECT_NEAR(limit / pre, 1.0, 1e-5);
}

TEST(QCDCollinearISR1, UXG_MUPMUMUX) {
    // y == 1
    // phase space:
    // (    46.30332265,               0,               0,     46.30332265)
    // (    46.30332265,               0,               0,    -46.30332265)
    // (    46.30332265,    -42.46700328,     12.05383173,     13.97341981)
    // (    46.30332265,     42.46700328,    -12.05383173,    -13.97341981)
    double bornme = 0.35916345530503124284;

    double s = 8575.99075378023;
    double xi = 0.285088492732579;
    double alpha_s = 0.129783;

    int bornpdgs[] = {-2, 2, -13, 13};
    int realpdgs[] = {-2, 21, -13, 13, -2};
    double limit = FKS::QCD::CollinearLimitISR(4, realpdgs, bornpdgs, xi, 1,
                                               0.0, s, alpha_s, bornme,
                                               UserProcess::SpinCorrelated());

    EXPECT_NEAR(limit, 0.0, 1e-15);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    xi = 0.134784312759763;
    limit = FKS::QCD::CollinearLimitISR(4, realpdgs, bornpdgs, xi, 1, 0.0, s,
                                        alpha_s, bornme,
                                        UserProcess::SpinCorrelated());
    EXPECT_NEAR(limit, 0.0, 1e-15);
}

TEST(QCDCollinearISR2, UXG_MUPMUMUX) {
    // y == -1
    // phase space:
    // (    46.30332265,               0,               0,     46.30332265)
    // (    46.30332265,               0,               0,    -46.30332265)
    // (    46.30332265,    -42.46700328,     12.05383173,     13.97341981)
    // (    46.30332265,     42.46700328,    -12.05383173,    -13.97341981)
    double bornme = 0.35916345530503124284;

    double s = 8575.99075378023;
    double xi = 0.626275281886555;
    double alpha_s = 0.129783;

    int bornpdgs[] = {-2, 2, -13, 13};
    int realpdgs[] = {-2, 21, -13, 13, -2};
    double limit = FKS::QCD::CollinearLimitISR(4, realpdgs, bornpdgs, xi, -1,
                                               0.0, s, alpha_s, bornme,
                                               UserProcess::SpinCorrelated());

    double pre = 9.10087872524386e-05;
    EXPECT_NEAR(limit / pre, 1.0, 1e-5);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    xi = 0.896328949986896;
    limit = FKS::QCD::CollinearLimitISR(4, realpdgs, bornpdgs, xi, -1, 0.0, s,
                                        alpha_s, bornme,
                                        UserProcess::SpinCorrelated());
    pre = 5.35866270885596e-05;
    EXPECT_NEAR(limit / pre, 1.0, 1e-5);
}

TEST(QCDCollinearISR1, GU_MUPMUMU) {
    // y == 1
    // phase space:
    // (    46.30332265,               0,               0,     46.30332265)
    // (    46.30332265,               0,               0,    -46.30332265)
    // (    46.30332265,    -42.46700328,     12.05383173,     13.97341981)
    // (    46.30332265,     42.46700328,    -12.05383173,    -13.97341981)
    double bornme = 0.35916345530503124284;

    double s = 8575.99075378023;
    double xi = 0.285088492732579;
    double alpha_s = 0.129783;

    int bornpdgs[] = {-2, 2, -13, 13};
    int realpdgs[] = {21, 2, -13, 13, -2};
    double limit = FKS::QCD::CollinearLimitISR(4, realpdgs, bornpdgs, xi, 1,
                                               0.0, s, alpha_s, bornme,
                                               UserProcess::SpinCorrelated());

    double pre = 4.61393118576831e-05;
    EXPECT_NEAR(limit / pre, 1.0, 1e-5);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    xi = 0.134784312759763;
    limit = FKS::QCD::CollinearLimitISR(4, realpdgs, bornpdgs, xi, 1, 0.0, s,
                                        alpha_s, bornme,
                                        UserProcess::SpinCorrelated());
    pre = 7.58899824526624e-06;
    EXPECT_NEAR(limit / pre, 1.0, 1e-5);
}

TEST(QCDCollinearISR2, GU_MUPMUMU) {
    // y == -1
    // phase space:
    // (    46.30332265,               0,               0,     46.30332265)
    // (    46.30332265,               0,               0,    -46.30332265)
    // (    46.30332265,    -42.46700328,     12.05383173,     13.97341981)
    // (    46.30332265,     42.46700328,    -12.05383173,    -13.97341981)
    double bornme = 0.35916345530503124284;

    double s = 8575.99075378023;
    double xi = 0.626275281886555;
    double alpha_s = 0.129783;

    int bornpdgs[] = {-2, 2, -13, 13};
    int realpdgs[] = {21, 2, -13, 13, -2};
    double limit = FKS::QCD::CollinearLimitISR(4, realpdgs, bornpdgs, xi, -1,
                                               0.0, s, alpha_s, bornme,
                                               UserProcess::SpinCorrelated());

    EXPECT_NEAR(limit, 0.0, 1e-15);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    xi = 0.896328949986896;
    limit = FKS::QCD::CollinearLimitISR(4, realpdgs, bornpdgs, xi, -1, 0.0, s,
                                        alpha_s, bornme,
                                        UserProcess::SpinCorrelated());
    EXPECT_NEAR(limit, 0.0, 1e-15);
}

TEST(QCDSoftCollinearISR1, UXU_MUPMUMG) {
    // y == 1
    // phase space:
    // (    46.30332265,               0,               0,     46.30332265)
    // (    46.30332265,               0,               0,    -46.30332265)
    // (    46.30332265,    -42.46700328,     12.05383173,     13.97341981)
    // (    46.30332265,     42.46700328,    -12.05383173,    -13.97341981)
    double bornme = 0.35916345530503124284;

    double s = 8575.99075378023;
    double alpha_s = 0.129783;

    int bornpdgs[] = {-2, 2, -13, 13};
    int realpdgs[] = {-2, 2, -13, 13, 21};
    double limit = FKS::QCD::SoftCollinearLimitISR(
        4, realpdgs, bornpdgs, 1, 0.0, s, alpha_s, bornme,
        UserProcess::SpinCorrelated());

    double pre = 0.00145711634230421;
    EXPECT_NEAR(limit / pre, 1.0, 1e-5);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    limit = FKS::QCD::SoftCollinearLimitISR(4, realpdgs, bornpdgs, 1, 0.0, s,
                                            alpha_s, bornme,
                                            UserProcess::SpinCorrelated());
    pre = 0.000391635050295075;
    EXPECT_NEAR(limit / pre, 1.0, 1e-5);
}

TEST(QCDSoftCollinearISR2, UXU_MUPMUMG) {
    // y == -1
    // phase space:
    // (    46.30332265,               0,               0,     46.30332265)
    // (    46.30332265,               0,               0,    -46.30332265)
    // (    46.30332265,    -42.46700328,     12.05383173,     13.97341981)
    // (    46.30332265,     42.46700328,    -12.05383173,    -13.97341981)
    double bornme = 0.35916345530503124284;

    double s = 8575.99075378023;
    double alpha_s = 0.129783;

    int bornpdgs[] = {-2, 2, -13, 13};
    int realpdgs[] = {-2, 2, -13, 13, 21};
    double limit = FKS::QCD::SoftCollinearLimitISR(
        4, realpdgs, bornpdgs, -1, 0.0, s, alpha_s, bornme,
        UserProcess::SpinCorrelated());

    double pre = 0.00145711634230421;
    EXPECT_NEAR(limit / pre, 1.0, 1e-5);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    limit = FKS::QCD::SoftCollinearLimitISR(4, realpdgs, bornpdgs, -1, 0.0, s,
                                            alpha_s, bornme,
                                            UserProcess::SpinCorrelated());
    pre = 0.000391635050295075;
    EXPECT_NEAR(limit / pre, 1.0, 1e-5);
}

TEST(QCDSoftCollinearISR1, UXG_MUPMUMUX) {
    // y == 1
    // phase space:
    // (    46.30332265,               0,               0,     46.30332265)
    // (    46.30332265,               0,               0,    -46.30332265)
    // (    46.30332265,    -42.46700328,     12.05383173,     13.97341981)
    // (    46.30332265,     42.46700328,    -12.05383173,    -13.97341981)
    double bornme = 0.35916345530503124284;

    double s = 8575.99075378023;
    double alpha_s = 0.129783;

    int bornpdgs[] = {-2, 2, -13, 13};
    int realpdgs[] = {-2, 21, -13, 13, -2};
    double limit = FKS::QCD::SoftCollinearLimitISR(
        4, realpdgs, bornpdgs, 1, 0.0, s, alpha_s, bornme,
        UserProcess::SpinCorrelated());

    EXPECT_NEAR(limit, 0.0, 1e-15);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    limit = FKS::QCD::SoftCollinearLimitISR(4, realpdgs, bornpdgs, 1, 0.0, s,
                                            alpha_s, bornme,
                                            UserProcess::SpinCorrelated());
    EXPECT_NEAR(limit, 0.0, 1e-15);
}

TEST(QCDSoftCollinearISR2, UXG_MUPMUMUX) {
    // y == -1
    // phase space:
    // (    46.30332265,               0,               0,     46.30332265)
    // (    46.30332265,               0,               0,    -46.30332265)
    // (    46.30332265,    -42.46700328,     12.05383173,     13.97341981)
    // (    46.30332265,     42.46700328,    -12.05383173,    -13.97341981)
    double bornme = 0.35916345530503124284;

    double s = 8575.99075378023;
    double alpha_s = 0.129783;

    int bornpdgs[] = {-2, 2, -13, 13};
    int realpdgs[] = {-2, 21, -13, 13, -2};
    double limit = FKS::QCD::SoftCollinearLimitISR(
        4, realpdgs, bornpdgs, -1, 0.0, s, alpha_s, bornme,
        UserProcess::SpinCorrelated());

    EXPECT_NEAR(limit, 0.0, 1e-15);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    limit = FKS::QCD::SoftCollinearLimitISR(4, realpdgs, bornpdgs, -1, 0.0, s,
                                            alpha_s, bornme,
                                            UserProcess::SpinCorrelated());
    EXPECT_NEAR(limit, 0.0, 1e-15);
}

TEST(QCDSoftCollinearISR1, GU_MUPMUMU) {
    // y == 1
    // phase space:
    // (    46.30332265,               0,               0,     46.30332265)
    // (    46.30332265,               0,               0,    -46.30332265)
    // (    46.30332265,    -42.46700328,     12.05383173,     13.97341981)
    // (    46.30332265,     42.46700328,    -12.05383173,    -13.97341981)
    double bornme = 0.35916345530503124284;

    double s = 8575.99075378023;
    double alpha_s = 0.129783;

    int bornpdgs[] = {-2, 2, -13, 13};
    int realpdgs[] = {21, 2, -13, 13, -2};
    double limit = FKS::QCD::SoftCollinearLimitISR(
        4, realpdgs, bornpdgs, 1, 0.0, s, alpha_s, bornme,
        UserProcess::SpinCorrelated());

    EXPECT_NEAR(limit, 0.0, 1e-15);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    limit = FKS::QCD::SoftCollinearLimitISR(4, realpdgs, bornpdgs, 1, 0.0, s,
                                            alpha_s, bornme,
                                            UserProcess::SpinCorrelated());

    EXPECT_NEAR(limit, 0.0, 1e-15);
}

TEST(QCDSoftCollinearISR2, GU_MUPMUMU) {
    // y == -1
    // phase space:
    // (    46.30332265,               0,               0,     46.30332265)
    // (    46.30332265,               0,               0,    -46.30332265)
    // (    46.30332265,    -42.46700328,     12.05383173,     13.97341981)
    // (    46.30332265,     42.46700328,    -12.05383173,    -13.97341981)
    double bornme = 0.35916345530503124284;

    double s = 8575.99075378023;
    double alpha_s = 0.129783;

    int bornpdgs[] = {-2, 2, -13, 13};
    int realpdgs[] = {21, 2, -13, 13, -2};
    double limit = FKS::QCD::SoftCollinearLimitISR(
        4, realpdgs, bornpdgs, -1, 0.0, s, alpha_s, bornme,
        UserProcess::SpinCorrelated());

    EXPECT_NEAR(limit, 0.0, 1e-15);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    limit = FKS::QCD::SoftCollinearLimitISR(4, realpdgs, bornpdgs, -1, 0.0, s,
                                            alpha_s, bornme,
                                            UserProcess::SpinCorrelated());
    EXPECT_NEAR(limit, 0.0, 1e-15);
}

// Test for quark to gluon ISR splitting
//                    q
//                   /
//                  /
//                 /
//   q -----------X/\/\/\/\/\/\(BORN)
//
// born process = g u >  mu+ nu d
// real process = u u >  mu+ nu d u
TEST(QCDCollinearISR1, QuarkToGluonSplitting) {
    double alphas = 0.11799999999999999;

    Phasespace::Phasespace ps;

    ps.Momenta[0].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, 7.500000000000000000000000e+02);
    ps.Momenta[1].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, -7.500000000000000000000000e+02);
    ps.Momenta[2].Set(
        6.878681818281603455034201e+02, 2.541798304645197390527755e+02,
        5.694804931172980104747694e+02, -2.902537119753787351328356e+02);
    ps.Momenta[3].Set(
        5.460999311052265738908318e+02, -2.749480393978777570396232e+01,
        -5.215564519790506210483727e+02, 1.595244116380622187989502e+02);
    ps.Momenta[4].Set(
        2.660318870666130806057481e+02, -2.266850265247319953232363e+02,
        -4.792404113824733968840519e+01, 1.307293003373165447555948e+02);

    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;

    double x1 = 0.5;
    double x2 = 0.2;
    ps.X1 = x1;
    ps.X2 = x2;
    auto p0 = ps.Momenta[0].Plus(ps.Momenta[1]);
    double s = p0.Dot(p0);
    ps.S = s / x1 / x2;
    ps.N = 3;

    double xi = 0.3;
    double y = 0.9999999999;
    double phi = 3.1415 / 3.0;

    UserProcess::SpinCorrelated spincorr;
    spincorr[0][0] = std::complex<double>(0, 0);
    spincorr[0][1] = std::complex<double>(0, 0);
    spincorr[0][2] = std::complex<double>(0, 0);
    spincorr[0][3] = std::complex<double>(0, 0);
    spincorr[1][0] = std::complex<double>(0, -0);
    spincorr[1][1] = std::complex<double>(6.83389304101289358e-09, 0);
    spincorr[1][2] =
        std::complex<double>(3.80254365238820558e-09, 5.49374912883507661e-09);
    spincorr[1][3] = std::complex<double>(0, 0);
    spincorr[2][0] = std::complex<double>(0, -0);
    spincorr[2][1] =
        std::complex<double>(3.80254365238820558e-09, -5.49374912883507661e-09);
    spincorr[2][2] = std::complex<double>(6.53223827926307841e-09, 0);
    spincorr[2][3] = std::complex<double>(0, 0);
    spincorr[3][0] = std::complex<double>(0, -0);
    spincorr[3][1] = std::complex<double>(0, -0);
    spincorr[3][2] = std::complex<double>(0, -0);
    spincorr[3][3] = std::complex<double>(0, 0);
    double born = 1.33661313202759728e-08;
    double me_real = 0.00308482655576152376;

    int bornpdgs[] = {21, 2, -13, 14, 1};
    int realpdgs[] = {2, 2, -13, 14, 1, 2};
    double result = FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1,
                                                phi, ps.X1 * ps.X2 * ps.S,
                                                alphas, born, spincorr);

    double pre = xi * xi * (1 - y * y) * me_real;
    ASSERT_NEAR(result / pre, 1.0, 5e-5);
}

// Test for quark to gluon ISR splitting
//                    q
//                   /
//                  /
//                 /
//   q -----------X/\/\/\/\/\/\(BORN)
//
// born process = g u >  mu+ nu d
// real process = u u >  mu+ nu d u
//
// anti collinear (y=-1). No collinear singularity.
TEST(QCDCollinearISR2, QuarkToGluonSplitting) {
    double alphas = 0.11799999999999999;

    double xi = 0.3;
    double phi = 3.1415 / 3.0;
    double born = 1.0;
    UserProcess::SpinCorrelated spincorr;

    double s = 1.0;

    int bornpdgs[] = {21, 2, -13, 14, 1};
    int realpdgs[] = {2, 2, -13, 14, 1, 2};
    double result = FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1,
                                                phi, s, alphas, born, spincorr);

    ASSERT_NEAR(result, 0.0, 5e-5);
}

// Test for gluon to gluon ISR splitting
//                     _ g
//                   _|
//                 _|
//                |
//   g /\/\/\/\/\X/\/\/\/\/\/\(BORN)
//
// born process: g u > mu+ nu d
// real process: g u > mu+ nu d g
TEST(QCDCollinearISR1, GluonToGluonSplitting) {
    double alphas = 0.11799999999999999;

    Phasespace::Phasespace ps;

    ps.Momenta[0].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, 7.500000000000000000000000e+02);
    ps.Momenta[1].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, -7.500000000000000000000000e+02);
    ps.Momenta[2].Set(
        6.878681818281603455034201e+02, 2.541798304645197390527755e+02,
        5.694804931172980104747694e+02, -2.902537119753787351328356e+02);
    ps.Momenta[3].Set(
        5.460999311052265738908318e+02, -2.749480393978777570396232e+01,
        -5.215564519790506210483727e+02, 1.595244116380622187989502e+02);
    ps.Momenta[4].Set(
        2.660318870666130806057481e+02, -2.266850265247319953232363e+02,
        -4.792404113824733968840519e+01, 1.307293003373165447555948e+02);

    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;

    double x1 = 0.5;
    double x2 = 0.2;
    ps.X1 = x1;
    ps.X2 = x2;
    auto p0 = ps.Momenta[0].Plus(ps.Momenta[1]);
    double s = p0.Dot(p0);
    ps.S = s / x1 / x2;
    ps.N = 3;

    double xi = 0.3;
    double y = 0.9999999999;
    double phi = 3.1415 / 3.0;

    UserProcess::SpinCorrelated spincorr;
    spincorr[0][0] = std::complex<double>(0, 0);
    spincorr[0][1] = std::complex<double>(0, 0);
    spincorr[0][2] = std::complex<double>(0, 0);
    spincorr[0][3] = std::complex<double>(0, 0);
    spincorr[1][0] = std::complex<double>(0, -0);
    spincorr[1][1] = std::complex<double>(6.83389304101289358e-09, 0);
    spincorr[1][2] =
        std::complex<double>(3.80254365238820558e-09, 5.49374912883507661e-09);
    spincorr[1][3] = std::complex<double>(0, 0);
    spincorr[2][0] = std::complex<double>(0, -0);
    spincorr[2][1] =
        std::complex<double>(3.80254365238820558e-09, -5.49374912883507661e-09);
    spincorr[2][2] = std::complex<double>(6.53223827926307841e-09, 0);
    spincorr[2][3] = std::complex<double>(0, 0);
    spincorr[3][0] = std::complex<double>(0, -0);
    spincorr[3][1] = std::complex<double>(0, -0);
    spincorr[3][2] = std::complex<double>(0, -0);
    spincorr[3][3] = std::complex<double>(0, 0);
    double born = 1.33661313202759728e-08;
    double me_real = 0.0223972592924080449;

    int bornpdgs[] = {21, 2, -13, 14, 1};
    int realpdgs[] = {21, 2, -13, 14, 1, 21};
    double result = FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1,
                                                phi, ps.X1 * ps.X2 * ps.S,
                                                alphas, born, spincorr);

    double pre = xi * xi * (1 - y * y) * me_real;
    ASSERT_NEAR(result / pre, 1.0, 5e-5);
}

// Test for quark to gluon ISR splitting y = -1
//                     _ g
//                   _|
//                 _|
//                |
//   g -----------X----------(BORN)
//
// born process: g u > mu+ nu d
// real process: g u > mu+ nu d g
TEST(QCDCollinearISR2, QuarkToGluonSplitting2) {
    double alphas = 0.11799999999999999;

    Phasespace::Phasespace ps;

    ps.Momenta[0].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, 7.500000000000000000000000e+02);
    ps.Momenta[1].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, -7.500000000000000000000000e+02);
    ps.Momenta[2].Set(
        6.878681818281603455034201e+02, 2.541798304645197390527755e+02,
        5.694804931172980104747694e+02, -2.902537119753787351328356e+02);
    ps.Momenta[3].Set(
        5.460999311052265738908318e+02, -2.749480393978777570396232e+01,
        -5.215564519790506210483727e+02, 1.595244116380622187989502e+02);
    ps.Momenta[4].Set(
        2.660318870666130806057481e+02, -2.266850265247319953232363e+02,
        -4.792404113824733968840519e+01, 1.307293003373165447555948e+02);

    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;

    double x1 = 0.5;
    double x2 = 0.2;
    ps.X1 = x1;
    ps.X2 = x2;
    auto p0 = ps.Momenta[0].Plus(ps.Momenta[1]);
    double s = p0.Dot(p0);
    ps.S = s / x1 / x2;
    ps.N = 3;

    double xi = 0.3;
    double y = -0.9999999999;
    double phi = 3.1415 / 3.0;

    UserProcess::SpinCorrelated spincorr;
    spincorr[0][0] = std::complex<double>(0, 0);
    spincorr[0][1] = std::complex<double>(0, 0);
    spincorr[0][2] = std::complex<double>(0, 0);
    spincorr[0][3] = std::complex<double>(0, 0);
    spincorr[1][0] = std::complex<double>(0, -0);
    spincorr[1][1] = std::complex<double>(6.83389304101289358e-09, 0);
    spincorr[1][2] =
        std::complex<double>(3.80254365238820558e-09, 5.49374912883507661e-09);
    spincorr[1][3] = std::complex<double>(0, 0);
    spincorr[2][0] = std::complex<double>(0, -0);
    spincorr[2][1] =
        std::complex<double>(3.80254365238820558e-09, -5.49374912883507661e-09);
    spincorr[2][2] = std::complex<double>(6.53223827926307841e-09, 0);
    spincorr[2][3] = std::complex<double>(0, 0);
    spincorr[3][0] = std::complex<double>(0, -0);
    spincorr[3][1] = std::complex<double>(0, -0);
    spincorr[3][2] = std::complex<double>(0, -0);
    spincorr[3][3] = std::complex<double>(0, 0);
    double born = 1.33661313202759728e-08;
    double me_real = 0.00777780227985162224;

    int bornpdgs[] = {21, 2, -13, 14, 1};
    int realpdgs[] = {21, 2, -13, 14, 1, 21};
    double result = FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1,
                                                phi, ps.X1 * ps.X2 * ps.S,
                                                alphas, born, spincorr);

    double pre = xi * xi * (1 - y * y) * me_real;
    ASSERT_NEAR(result / pre, 1.0, 5e-5);
}

TEST(QCDSoftCollinearISR1, QuarkToGluonSplitting) {
    double alphas = 0.11799999999999999;

    Phasespace::Phasespace ps;

    ps.Momenta[0].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, 7.500000000000000000000000e+02);
    ps.Momenta[1].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, -7.500000000000000000000000e+02);
    ps.Momenta[2].Set(
        6.878681818281603455034201e+02, 2.541798304645197390527755e+02,
        5.694804931172980104747694e+02, -2.902537119753787351328356e+02);
    ps.Momenta[3].Set(
        5.460999311052265738908318e+02, -2.749480393978777570396232e+01,
        -5.215564519790506210483727e+02, 1.595244116380622187989502e+02);
    ps.Momenta[4].Set(
        2.660318870666130806057481e+02, -2.266850265247319953232363e+02,
        -4.792404113824733968840519e+01, 1.307293003373165447555948e+02);

    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;

    double x1 = 0.5;
    double x2 = 0.2;
    ps.X1 = x1;
    ps.X2 = x2;
    auto p0 = ps.Momenta[0].Plus(ps.Momenta[1]);
    double s = p0.Dot(p0);
    ps.S = s / x1 / x2;
    ps.N = 3;

    double phi = 3.1415 / 3.0;

    UserProcess::SpinCorrelated spincorr;
    spincorr[0][0] = std::complex<double>(0, 0);
    spincorr[0][1] = std::complex<double>(0, 0);
    spincorr[0][2] = std::complex<double>(0, 0);
    spincorr[0][3] = std::complex<double>(0, 0);
    spincorr[1][0] = std::complex<double>(0, -0);
    spincorr[1][1] = std::complex<double>(6.83389304101289358e-09, 0);
    spincorr[1][2] =
        std::complex<double>(3.80254365238820558e-09, 5.49374912883507661e-09);
    spincorr[1][3] = std::complex<double>(0, 0);
    spincorr[2][0] = std::complex<double>(0, -0);
    spincorr[2][1] =
        std::complex<double>(3.80254365238820558e-09, -5.49374912883507661e-09);
    spincorr[2][2] = std::complex<double>(6.53223827926307841e-09, 0);
    spincorr[2][3] = std::complex<double>(0, 0);
    spincorr[3][0] = std::complex<double>(0, -0);
    spincorr[3][1] = std::complex<double>(0, -0);
    spincorr[3][2] = std::complex<double>(0, -0);
    spincorr[3][3] = std::complex<double>(0, 0);
    double born = 1.33661313202759728e-08;

    int bornpdgs[] = {21, 2, -13, 14, 1};
    int realpdgs[] = {2, 2, -13, 14, 1, 2};
    double result = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                    phi, ps.X1 * ps.X2 * ps.S,
                                                    alphas, born, spincorr);

    ASSERT_NEAR(result, 0.0, 1e-15);
}

// Test for gluon to gluon ISR splitting
//                     _ g
//                   _|
//                 _|
//                |
//   g /\/\/\/\/\X/\/\/\/\/\/\(BORN)
//
// real process: g u > mu+ nu d g
// born process: g u > mu+ nu d
TEST(QCDSoftCollinearISR1, GluonToGluonSplitting) {
    double alphas = 0.11799999999999999;

    Phasespace::Phasespace ps;

    ps.Momenta[0].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, 7.500000000000000000000000e+02);
    ps.Momenta[1].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, -7.500000000000000000000000e+02);
    ps.Momenta[2].Set(
        6.878681818281603455034201e+02, 2.541798304645197390527755e+02,
        5.694804931172980104747694e+02, -2.902537119753787351328356e+02);
    ps.Momenta[3].Set(
        5.460999311052265738908318e+02, -2.749480393978777570396232e+01,
        -5.215564519790506210483727e+02, 1.595244116380622187989502e+02);
    ps.Momenta[4].Set(
        2.660318870666130806057481e+02, -2.266850265247319953232363e+02,
        -4.792404113824733968840519e+01, 1.307293003373165447555948e+02);

    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;

    double x1 = 0.5;
    double x2 = 0.2;
    ps.X1 = x1;
    ps.X2 = x2;
    auto p0 = ps.Momenta[0].Plus(ps.Momenta[1]);
    double s = p0.Dot(p0);
    ps.S = s / x1 / x2;
    ps.N = 3;

    double xi = 1e-3;
    double y = 0.9999999999;
    double phi = 3.1415 / 3.0;

    UserProcess::SpinCorrelated spincorr;

    spincorr[0][0] = std::complex<double>(0, 0);
    spincorr[0][1] = std::complex<double>(0, 0);
    spincorr[0][2] = std::complex<double>(0, 0);
    spincorr[0][3] = std::complex<double>(0, 0);
    spincorr[1][0] = std::complex<double>(0, -0);
    spincorr[1][1] = std::complex<double>(6.83389304101289358e-09, 0);
    spincorr[1][2] =
        std::complex<double>(3.80254365238820558e-09, 5.49374912883507661e-09);
    spincorr[1][3] = std::complex<double>(0, 0);
    spincorr[2][0] = std::complex<double>(0, -0);
    spincorr[2][1] =
        std::complex<double>(3.80254365238820558e-09, -5.49374912883507661e-09);
    spincorr[2][2] = std::complex<double>(6.53223827926307841e-09, 0);
    spincorr[2][3] = std::complex<double>(0, 0);
    spincorr[3][0] = std::complex<double>(0, -0);
    spincorr[3][1] = std::complex<double>(0, -0);
    spincorr[3][2] = std::complex<double>(0, -0);
    spincorr[3][3] = std::complex<double>(0, 0);
    double born = 1.33661313202759728e-08;
    double me_real = 2113.44156045684531;

    int bornpdgs[] = {21, 2, -13, 14, 1};
    int realpdgs[] = {21, 2, -13, 14, 1, 21};
    double result = FKS::QCD::SoftCollinearLimitISR(5, realpdgs, bornpdgs, 1,
                                                    phi, ps.X1 * ps.X2 * ps.S,
                                                    alphas, born, spincorr);

    double pre = xi * xi * (1 - y * y) * me_real;
    ASSERT_NEAR(result / pre, 1.0, 5e-4);
}

//
//                         q
//                         /
//                        /
//                       /
//                      /
//   (BORN)/\/\/\/\/\/\X---------- q
//
// real process: u d~ > mu+ nu c c~
// born process: u d~ > mu+ nu g
TEST(QCDCollinearFSR, GluonSplittingToQuarks) {

    double alphas = 0.11799999999999999;

    Phasespace::Phasespace ps;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.S = 1.96e+08;
    ps.X1 = 0.191096;
    ps.X2 = 0.068618;
    ps.Jacobian = 614.835;
    ps.N = 3;
    ps.Momenta[0].Set(801.5728381, 0, 0, 801.5728381);
    ps.Momenta[1].Set(801.5728381, 0, 0, -801.5728381);
    ps.Momenta[2].Set(105.7409648, 66.71985503, 75.49816789, -32.08799204);
    ps.Momenta[3].Set(736.3506575, 364.2504358, -372.2172065, -520.5653292);
    ps.Momenta[4].Set(761.0540539, -430.9702909, 296.7190386, 552.6533212);

    UserProcess::SpinCorrelated spincorr;
    spincorr[0][0] = std::complex<double>(0, 0);
    spincorr[0][1] = std::complex<double>(0, 0);
    spincorr[0][2] = std::complex<double>(0, 0);
    spincorr[0][3] = std::complex<double>(0, 0);
    spincorr[1][0] = std::complex<double>(0, -0);
    spincorr[1][1] = std::complex<double>(7.4992605567097649e-08, 0);
    spincorr[1][2] =
        std::complex<double>(4.31494876871413944e-08, -1.86946995660878943e-08);
    spincorr[1][3] =
        std::complex<double>(3.53138392285661161e-08, 1.00371662837759417e-08);
    spincorr[2][0] = std::complex<double>(0, -0);
    spincorr[2][1] =
        std::complex<double>(4.31494876871413944e-08, 1.86946995660878943e-08);
    spincorr[2][2] =
        std::complex<double>(2.94878416719430758e-08, 2.9410866224107649e-24);
    spincorr[2][3] =
        std::complex<double>(1.78168534444796286e-08, 1.45785066355718212e-08);
    spincorr[3][0] = std::complex<double>(0, -0);
    spincorr[3][1] =
        std::complex<double>(3.53138392285661161e-08, -1.00371662837759417e-08);
    spincorr[3][2] =
        std::complex<double>(1.78168534444796286e-08, -1.45785066355718212e-08);
    spincorr[3][3] =
        std::complex<double>(1.79725979365160808e-08, -1.47054331120538245e-24);

    double bornme = 1.22453045175556812e-07;

    double xi = 0.1;
    double y = 0.9999999;
    double phi = M_PI / 3.0;

    double me_real = 1.06658e-05;

    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, 4, -4};
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, bornme, spincorr);
    double pre = me_real * (1 - y) * xi * xi;

    ASSERT_NEAR(pre / limit, 1.0, 5e-4);
}

//
//                         q
//                         /
//                        /
//                       /
//                      /
//   (BORN)/\/\/\/\/\/\X---------- q
//
// real process: u d~ > mu+ nu c c~
// born process: u d~ > mu+ nu g
TEST(QCDSoftCollinearFSR, GluonSplittingToQuarks) {

    double alphas = 0.11799999999999999;

    Phasespace::Phasespace ps;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.S = 1.96e+08;
    ps.X1 = 0.191096;
    ps.X2 = 0.068618;
    ps.Jacobian = 614.835;
    ps.N = 3;
    ps.Momenta[0].Set(801.5728381, 0, 0, 801.5728381);
    ps.Momenta[1].Set(801.5728381, 0, 0, -801.5728381);
    ps.Momenta[2].Set(105.7409648, 66.71985503, 75.49816789, -32.08799204);
    ps.Momenta[3].Set(736.3506575, 364.2504358, -372.2172065, -520.5653292);
    ps.Momenta[4].Set(761.0540539, -430.9702909, 296.7190386, 552.6533212);

    UserProcess::SpinCorrelated spincorr;

    spincorr[0][0] = std::complex<double>(0, 0);
    spincorr[0][1] = std::complex<double>(0, 0);
    spincorr[0][2] = std::complex<double>(0, 0);
    spincorr[0][3] = std::complex<double>(0, 0);
    spincorr[1][0] = std::complex<double>(0, -0);
    spincorr[1][1] = std::complex<double>(7.4992605567097649e-08, 0);
    spincorr[1][2] =
        std::complex<double>(4.31494876871413944e-08, -1.86946995660878943e-08);
    spincorr[1][3] =
        std::complex<double>(3.53138392285661161e-08, 1.00371662837759417e-08);
    spincorr[2][0] = std::complex<double>(0, -0);
    spincorr[2][1] =
        std::complex<double>(4.31494876871413944e-08, 1.86946995660878943e-08);
    spincorr[2][2] =
        std::complex<double>(2.94878416719430758e-08, 2.9410866224107649e-24);
    spincorr[2][3] =
        std::complex<double>(1.78168534444796286e-08, 1.45785066355718212e-08);
    spincorr[3][0] = std::complex<double>(0, -0);
    spincorr[3][1] =
        std::complex<double>(3.53138392285661161e-08, -1.00371662837759417e-08);
    spincorr[3][2] =
        std::complex<double>(1.78168534444796286e-08, -1.45785066355718212e-08);
    spincorr[3][3] =
        std::complex<double>(1.79725979365160808e-08, -1.47054331120538245e-24);
    double bornme = 1.22453045175556812e-07;

    double phi = M_PI / 3.0;

    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, 4, -4};
    double limit = FKS::QCD::SoftCollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, phi, alphas, bornme, spincorr);

    ASSERT_NEAR(limit, 0.0, 1e-14);
}

// Test for gluon to gluon FSR splitting
//                        _ g
//                      _|
//                    _|
//                   |
//   (BORN)/\/\/\/\/\X/\/\/\/\/\/\ g
//
// real process: u d~ > mu+ nu g g
// born process: u d~ > mu+ nu g
TEST(QCDCollinearFSR, GluonSplittingToGluons) {

    double alphas = 0.11799999999999999;

    Phasespace::Phasespace ps;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.S = 1.96e+08;
    ps.X1 = 0.191096;
    ps.X2 = 0.068618;
    ps.Jacobian = 614.835;
    ps.N = 3;
    ps.Momenta[0].Set(801.5728381, 0, 0, 801.5728381);
    ps.Momenta[1].Set(801.5728381, 0, 0, -801.5728381);
    ps.Momenta[2].Set(105.7409648, 66.71985503, 75.49816789, -32.08799204);
    ps.Momenta[3].Set(736.3506575, 364.2504358, -372.2172065, -520.5653292);
    ps.Momenta[4].Set(761.0540539, -430.9702909, 296.7190386, 552.6533212);

    UserProcess::SpinCorrelated spincorr;
    spincorr[0][0] = std::complex<double>(0, 0);
    spincorr[0][1] = std::complex<double>(0, 0);
    spincorr[0][2] = std::complex<double>(0, 0);
    spincorr[0][3] = std::complex<double>(0, 0);
    spincorr[1][0] = std::complex<double>(0, -0);
    spincorr[1][1] = std::complex<double>(7.4992605567097649e-08, 0);
    spincorr[1][2] =
        std::complex<double>(4.31494876871413944e-08, -1.86946995660878943e-08);
    spincorr[1][3] =
        std::complex<double>(3.53138392285661161e-08, 1.00371662837759417e-08);
    spincorr[2][0] = std::complex<double>(0, -0);
    spincorr[2][1] =
        std::complex<double>(4.31494876871413944e-08, 1.86946995660878943e-08);
    spincorr[2][2] =
        std::complex<double>(2.94878416719430758e-08, 2.9410866224107649e-24);
    spincorr[2][3] =
        std::complex<double>(1.78168534444796286e-08, 1.45785066355718212e-08);
    spincorr[3][0] = std::complex<double>(0, -0);
    spincorr[3][1] =
        std::complex<double>(3.53138392285661161e-08, -1.00371662837759417e-08);
    spincorr[3][2] =
        std::complex<double>(1.78168534444796286e-08, -1.45785066355718212e-08);
    spincorr[3][3] =
        std::complex<double>(1.79725979365160808e-08, -1.47054331120538245e-24);

    double bornme = 1.22453045175556812e-07;

    double xi = 0.1;
    double y = 0.9999999;
    double phi = M_PI / 3.0;

    double me_real = 0.000877256509901880028;

    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, 0, 0};

    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, bornme, spincorr);
    double pre = me_real * (1 - y) * xi * xi;

    ASSERT_NEAR(pre / limit, 1.0, 5e-4);
}

// Test for gluon to gluon FSR splitting
//                        _ g
//                      _|
//                    _|
//                   |
//   (BORN)/\/\/\/\/\X/\/\/\/\/\/\ g
//
// real process: u d~ > mu+ nu g g
// born process: u d~ > mu+ nu g
TEST(QCDSoftCollinearFSR, GluonSplittingToGluons) {

    double alphas = 0.11799999999999999;

    Phasespace::Phasespace ps;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.S = 1.96e+08;
    ps.X1 = 0.191096;
    ps.X2 = 0.068618;
    ps.Jacobian = 614.835;
    ps.N = 3;
    ps.Momenta[0].Set(801.5728381, 0, 0, 801.5728381);
    ps.Momenta[1].Set(801.5728381, 0, 0, -801.5728381);
    ps.Momenta[2].Set(105.7409648, 66.71985503, 75.49816789, -32.08799204);
    ps.Momenta[3].Set(736.3506575, 364.2504358, -372.2172065, -520.5653292);
    ps.Momenta[4].Set(761.0540539, -430.9702909, 296.7190386, 552.6533212);

    double xi = 1e-5;
    double y = 0.9999999;
    double phi = M_PI / 3.0;

    UserProcess::SpinCorrelated spincorr;

    spincorr[0][0] = std::complex<double>(0, 0);
    spincorr[0][1] = std::complex<double>(0, 0);
    spincorr[0][2] = std::complex<double>(0, 0);
    spincorr[0][3] = std::complex<double>(0, 0);
    spincorr[1][0] = std::complex<double>(0, -0);
    spincorr[1][1] = std::complex<double>(7.4992605567097649e-08, 0);
    spincorr[1][2] =
        std::complex<double>(4.31494876871413944e-08, -1.86946995660878943e-08);
    spincorr[1][3] =
        std::complex<double>(3.53138392285661161e-08, 1.00371662837759417e-08);
    spincorr[2][0] = std::complex<double>(0, -0);
    spincorr[2][1] =
        std::complex<double>(4.31494876871413944e-08, 1.86946995660878943e-08);
    spincorr[2][2] =
        std::complex<double>(2.94878416719430758e-08, 2.9410866224107649e-24);
    spincorr[2][3] =
        std::complex<double>(1.78168534444796286e-08, 1.45785066355718212e-08);
    spincorr[3][0] = std::complex<double>(0, -0);
    spincorr[3][1] =
        std::complex<double>(3.53138392285661161e-08, -1.00371662837759417e-08);
    spincorr[3][2] =
        std::complex<double>(1.78168534444796286e-08, -1.45785066355718212e-08);
    spincorr[3][3] =
        std::complex<double>(1.79725979365160808e-08, -1.47054331120538245e-24);
    double bornme = 1.22453045175556812e-07;
    double me_real = 84749.3272381366114;

    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, 0, 0};

    double limit = FKS::QCD::SoftCollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, phi, alphas, bornme, spincorr);
    double pre = me_real * (1 - y) * xi * xi;

    ASSERT_NEAR(pre / limit, 1.0, 5e-4);
}

TEST(QCDCollinearFSR, QuarkSplitting) {

    double alphas = 0.11799999999999999;

    Phasespace::Phasespace ps;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.S = 1.96e+08;
    ps.X1 = 0.191096;
    ps.X2 = 0.068618;
    ps.Jacobian = 614.835;
    ps.N = 3;
    ps.Momenta[0].Set(801.5728381, 0, 0, 801.5728381);
    ps.Momenta[1].Set(801.5728381, 0, 0, -801.5728381);
    ps.Momenta[2].Set(105.7409648, 66.71985503, 75.49816789, -32.08799204);
    ps.Momenta[3].Set(736.3506575, 364.2504358, -372.2172065, -520.5653292);
    ps.Momenta[4].Set(761.0540539, -430.9702909, 296.7190386, 552.6533212);

    UserProcess::SpinCorrelated spincorr;
    double bornme = 1.11473e-06;

    double xi = 0.1;
    double y = 0.99999999999;
    double phi = M_PI / 3.0;

    double me_real = 69.0295;

    int bornpdgs[] = {0, 2, -13, 14, 1};
    int realpdgs[] = {0, 2, -13, 14, 1, 0};
    double limit = FKS::QCD::CollinearLimitFSR(
        realpdgs[4], bornpdgs[4], ps, 4, xi, phi, alphas, bornme, spincorr);
    double pre = me_real * (1 - y) * xi * xi;

    ASSERT_NEAR(pre / limit, 1.0, 5e-5);
}

// Test for soft gluon radiation
//
// multiplies the matrix element with the ISR FKS prefactor xi^2*(1-y^2).
TEST(QCDSoftLimit, GluonRadiationWj) {

    double alphas = 0.11799999999999999;

    Phasespace::Phasespace ps;

    ps.Momenta[0].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, 7.500000000000000000000000e+02);
    ps.Momenta[1].Set(
        7.500000000000000000000000e+02, 0.000000000000000000000000e+00,
        0.000000000000000000000000e+00, -7.500000000000000000000000e+02);
    ps.Momenta[2].Set(
        6.878681818281603455034201e+02, 2.541798304645197390527755e+02,
        5.694804931172980104747694e+02, -2.902537119753787351328356e+02);
    ps.Momenta[3].Set(
        5.460999311052265738908318e+02, -2.749480393978777570396232e+01,
        -5.215564519790506210483727e+02, 1.595244116380622187989502e+02);
    ps.Momenta[4].Set(
        2.660318870666130806057481e+02, -2.266850265247319953232363e+02,
        -4.792404113824733968840519e+01, 1.307293003373165447555948e+02);

    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;

    double x1 = 0.5;
    double x2 = 0.2;
    ps.X1 = x1;
    ps.X2 = x2;
    auto p0 = ps.Momenta[0].Plus(ps.Momenta[1]);
    double s = p0.Dot(p0);
    ps.S = s / x1 / x2;
    ps.N = 3;

    double xi = 1e-8;
    double y = 0.1;
    double phi = 3.1415 / 3.0;

    double born = 1.33661313202759728e-08;
    double me_real = 2421.79770727493587;

    int bornpdgs[] = {21, 2, -13, 14, 1};
    int realpdgs[] = {21, 2, -13, 14, 1, 21};

    Util::Matrix2 ColorCorr;
    ColorCorr.SetLen(5);
    double CA = 3.0;
    // born color factor: N / 2
    ColorCorr.Set(4, 1, -1. / 6. * born); // real color factor: -1/4
    ColorCorr.Set(1, 4, -1. / 6. * born);
    ColorCorr.Set(4, 0, CA / 2 * born); // real color factor: N^2 / 4
    ColorCorr.Set(0, 4, CA / 2 * born);
    ColorCorr.Set(0, 1, CA / 2 * born);
    ColorCorr.Set(1, 0, CA / 2 * born);

    double limit =
        FKS::QCD::SoftLimit(5, ps.Momenta.data(), bornpdgs, realpdgs[5],
                            ps.X1 * ps.X2 * ps.S, 0, alphas, ColorCorr, y, phi);
    double pre = xi * xi * (1 - y * y) * me_real;
    ASSERT_NEAR(limit / pre, 1.0, 1e-8);
}

TEST(QCDCollinearFSR, TwoGluonFinalState) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.Momenta[0].Set(59.735574318254578, 0, 0, 59.735574318254578);
    ps.Momenta[1].Set(59.735574318254578, 0, 0, -59.735574318254578);
    ps.Momenta[2].Set(58.816829706527031, -49.254881082388366,
                      -30.510807882710672, 10.122586064376122);
    ps.Momenta[3].Set(18.886606419725116, 11.502312673385969,
                      14.914157238664279, -1.403074865309847);
    ps.Momenta[4].Set(41.767712510257013, 37.752568409002393,
                      15.596650644046395, -8.7195111990662753);
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.N = 3;
    ps.X1 = 0.041562390198298459;
    ps.X2 = 0.0020320709083419168;
    ps.S = 169000000;
    ps.Jacobian = 2.5146302639035372;

    double alphas = 0.11799999999999999;

    double y = -0.99999998387895999996;
    double phi = 2.9061114346473165;

    double xi = 0.0497308378401237067;

    int bornpdgs[] = {2, -1, -13, 14, 0};
    int realpdgs[] = {2, -1, -13, 14, 21, 21};

    double born = 8.69230960775116906e-05;
    double real = 1149.16750749380208;

    UserProcess::SpinCorrelated spin; // unused
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, -1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);

    double me = xi * xi * (1 - y * y) * real;
    EXPECT_NEAR(limit / me, 1.0, 2e-4);
}

TEST(QCDCollinearISR1, NAME) {
    // born phase space
    Phasespace::Phasespace ps;
    ps.S = 169000000;
    ps.X1 = 3.907985046680551e-14;
    ps.X2 = 0.0009853946746503084;
    ps.Jacobian = 1.571531091663275e-12;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.Momenta[0].Set(4.033619322210627e-05, 0, 0, 4.033619322210627e-05);
    ps.Momenta[1].Set(4.033619322210627e-05, 0, 0, -4.033619322210627e-05);
    ps.Momenta[2].Set(1.197745634288066e-05, 1.068563839979435e-05,
                      -5.06149328043127e-06, 1.912558026989396e-06);
    ps.Momenta[3].Set(3.003797300357577e-05, 8.75773984736797e-06,
                      -1.710411577839625e-05, 2.30874649619352e-05);
    ps.Momenta[4].Set(3.865695709775614e-05, -1.944337824716231e-05,
                      2.216560905882752e-05, -2.500002298892459e-05);

    // parameters
    double alphas = 0.11799999999999999;

    // Born matrix element
    double born = 9.205883526290226e-20;

    // spin correlated born
    UserProcess::SpinCorrelated spin;
    spin[0][0] = std::complex<double>(0, 0);
    spin[0][1] = std::complex<double>(0, 0);
    spin[0][2] = std::complex<double>(0, 0);
    spin[0][3] = std::complex<double>(0, 0);
    spin[1][0] = std::complex<double>(0, -0);
    spin[1][1] =
        std::complex<double>(1.709225873475582e-20, 1.003088512701685e-36);
    spin[1][2] =
        std::complex<double>(1.119475294245184e-20, -3.400037689549069e-20);
    spin[1][3] = std::complex<double>(0, 0);
    spin[2][0] = std::complex<double>(0, -0);
    spin[2][1] =
        std::complex<double>(1.119475294245184e-20, 3.400037689549069e-20);
    spin[2][2] =
        std::complex<double>(7.496657652814644e-20, -1.003088512701685e-36);
    spin[2][3] = std::complex<double>(0, 0);
    spin[3][0] = std::complex<double>(0, -0);
    spin[3][1] = std::complex<double>(0, -0);
    spin[3][2] = std::complex<double>(0, -0);
    spin[3][3] = std::complex<double>(0, 0);

    // radiation variables
    double phi = 0.4544334237382444;
    double xi = 0.487217223946826;

    int bornpdgs[] = {0, -1, -13, 14, -2};
    int realpdgs[] = {2, -1, -13, 14, -2, 2};

    // FKS g = real matrix element * (1.0 - y * y) * xi * xi
    double fks_g = 2.227407152119047e-10;
    double limit =
        FKS::QCD::CollinearLimitISR(5, realpdgs, bornpdgs, xi, 1, phi,
                                    ps.X1 * ps.X2 * ps.S, alphas, born, spin);
    // = 2.22742e-10
    EXPECT_NEAR(limit / fks_g, 1.0, 1e-5) << "limit = " << limit
                                          << ", fks_g = " << fks_g;
}

#include "wjlimits.cpp"
