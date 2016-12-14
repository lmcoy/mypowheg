#include "gtest/gtest.h"

#include "math/fourmomentum.h"
#include "fks/limits.h"

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
    double xi = 0.9;
    double limit =
        FKS::QED::CollinearLimitFSR(13, 13, xi, s, 1.0 / 132.507, bornme);

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
    int bornpdgs[] = { 2, -2, -13, 13 };
    int realpdgs[] = { 2, -2, -13, 13, 22 };
    double limit = FKS::QED::CollinearLimitISR(realpdgs, bornpdgs, xi, 1, s,
                                               1.0 / 132.507, bornme);

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
    limit = FKS::QED::CollinearLimitISR(realpdgs, bornpdgs, xi, 1, s,
                                        1.0 / 132.507, bornme);
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
    int bornpdgs[] = { 2, -2, -13, 13 };
    int realpdgs[] = { 2, -2, -13, 13, 22 };
    double limit = FKS::QED::CollinearLimitISR(realpdgs, bornpdgs, xi, -1, s,
                                               1.0 / 132.507, bornme);

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
    limit = FKS::QED::CollinearLimitISR(realpdgs, bornpdgs, xi, -1, s,
                                        1.0 / 132.507, bornme);
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

    int pdg[] = { 2, -2, -13, 13 };
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

    int pdg[] = { 2, -2, -13, 13 };
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
    double limit =
        FKS::QED::SoftCollinearLimitFSR(13, 13, s, 1.0 / 132.507, bornme);

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
    int bornpdgs[] = { 2, -2, -13, 13 };
    int realpdgs[] = { 2, -2, -13, 13, 22 };
    double s = 3360.0;
    double limit = FKS::QED::SoftCollinearLimitISR(realpdgs, bornpdgs, 1, s,
                                                   1.0 / 132.507, bornme);
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

    int pdg[] = { 2, -2, 13, -13 };
    int realpdgs[] = { 2, -2, -13, 13, 22 };
    // ISR
    double limit_s =
        FKS::QED::SoftLimit(4, born_ps, pdg, 22, s, 0, 1.0 / 132.507,
                            Util::Matrix2(1, bornme), y, phi);
    double limit_sc = FKS::QED::SoftCollinearLimitISR(realpdgs, pdg, 1, s,
                                                      1.0 / 132.507, bornme);
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

    int pdg[] = { 2, -2, 13, -13 };
    double limit_s_2 =
        FKS::QED::SoftLimit(4, born_ps, pdg, 22, s, 2, 1.0 / 132.507,
                            Util::Matrix2(1, bornme), y, phi);
    double limit_s_3 =
        FKS::QED::SoftLimit(4, born_ps, pdg, 22, s, 3, 1.0 / 132.507,
                            Util::Matrix2(1, bornme), y, phi);
    double limit_sc =
        FKS::QED::SoftCollinearLimitFSR(13, 13, s, 1.0 / 132.507, bornme);
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
    int pdgs[] = { -2, 2, -13, 13 };
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
    int pdgs[] = { 21, 2, -13, 13 };
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
    int pdgs[] = { -2, 21, -13, 13 };
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

    int bornpdgs[] = { -2, 2, -13, 13 };
    int realpdgs[] = { -2, 2, -13, 13, 21 };
    double limit = FKS::QCD::CollinearLimitISR(realpdgs, bornpdgs, xi, 1, s,
                                               alpha_s, bornme);

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
    limit = FKS::QCD::CollinearLimitISR(realpdgs, bornpdgs, xi, 1, s, alpha_s,
                                        bornme);
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

    int bornpdgs[] = { -2, 2, -13, 13 };
    int realpdgs[] = { -2, 2, -13, 13, 21 };
    double limit = FKS::QCD::CollinearLimitISR(realpdgs, bornpdgs, xi, -1, s,
                                               alpha_s, bornme);

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
    limit = FKS::QCD::CollinearLimitISR(realpdgs, bornpdgs, xi, -1, s, alpha_s,
                                        bornme);
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

    int bornpdgs[] = { -2, 2, -13, 13 };
    int realpdgs[] = { -2, 21, -13, 13, -2 };
    double limit = FKS::QCD::CollinearLimitISR(realpdgs, bornpdgs, xi, 1, s,
                                               alpha_s, bornme);

    EXPECT_NEAR(limit, 0.0, 1e-15);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    xi = 0.134784312759763;
    limit = FKS::QCD::CollinearLimitISR(realpdgs, bornpdgs, xi, 1, s, alpha_s,
                                        bornme);
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

    int bornpdgs[] = { -2, 2, -13, 13 };
    int realpdgs[] = { -2, 21, -13, 13, -2 };
    double limit = FKS::QCD::CollinearLimitISR(realpdgs, bornpdgs, xi, -1, s,
                                               alpha_s, bornme);

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
    limit = FKS::QCD::CollinearLimitISR(realpdgs, bornpdgs, xi, -1, s, alpha_s,
                                        bornme);
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

    int bornpdgs[] = { -2, 2, -13, 13 };
    int realpdgs[] = { 21, 2, -13, 13, -2 };
    double limit = FKS::QCD::CollinearLimitISR(realpdgs, bornpdgs, xi, 1, s,
                                               alpha_s, bornme);

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
    limit = FKS::QCD::CollinearLimitISR(realpdgs, bornpdgs, xi, 1, s, alpha_s,
                                        bornme);
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

    int bornpdgs[] = { -2, 2, -13, 13 };
    int realpdgs[] = { 21, 2, -13, 13, -2 };
    double limit = FKS::QCD::CollinearLimitISR(realpdgs, bornpdgs, xi, -1, s,
                                               alpha_s, bornme);

    EXPECT_NEAR(limit, 0.0, 1e-15);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    xi = 0.896328949986896;
    limit = FKS::QCD::CollinearLimitISR(realpdgs, bornpdgs, xi, -1, s, alpha_s,
                                        bornme);
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

    int bornpdgs[] = { -2, 2, -13, 13 };
    int realpdgs[] = { -2, 2, -13, 13, 21 };
    double limit = FKS::QCD::SoftCollinearLimitISR(realpdgs, bornpdgs, 1, s,
                                                   alpha_s, bornme);

    double pre = 0.00145711634230421;
    EXPECT_NEAR(limit / pre, 1.0, 1e-5);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    limit = FKS::QCD::SoftCollinearLimitISR(realpdgs, bornpdgs, 1, s, alpha_s,
                                            bornme);
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

    int bornpdgs[] = { -2, 2, -13, 13 };
    int realpdgs[] = { -2, 2, -13, 13, 21 };
    double limit = FKS::QCD::SoftCollinearLimitISR(realpdgs, bornpdgs, -1, s,
                                                   alpha_s, bornme);

    double pre = 0.00145711634230421;
    EXPECT_NEAR(limit / pre, 1.0, 1e-5);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    limit = FKS::QCD::SoftCollinearLimitISR(realpdgs, bornpdgs, -1, s, alpha_s,
                                            bornme);
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

    int bornpdgs[] = { -2, 2, -13, 13 };
    int realpdgs[] = { -2, 21, -13, 13, -2 };
    double limit = FKS::QCD::SoftCollinearLimitISR(realpdgs, bornpdgs, 1, s,
                                                   alpha_s, bornme);

    EXPECT_NEAR(limit, 0.0, 1e-15);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    limit = FKS::QCD::SoftCollinearLimitISR(realpdgs, bornpdgs, 1, s, alpha_s,
                                            bornme);
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

    int bornpdgs[] = { -2, 2, -13, 13 };
    int realpdgs[] = { -2, 21, -13, 13, -2 };
    double limit = FKS::QCD::SoftCollinearLimitISR(realpdgs, bornpdgs, -1, s,
                                                   alpha_s, bornme);

    EXPECT_NEAR(limit, 0.0, 1e-15);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    limit = FKS::QCD::SoftCollinearLimitISR(realpdgs, bornpdgs, -1, s, alpha_s,
                                            bornme);
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

    int bornpdgs[] = { -2, 2, -13, 13 };
    int realpdgs[] = { 21, 2, -13, 13, -2 };
    double limit = FKS::QCD::SoftCollinearLimitISR(realpdgs, bornpdgs, 1, s,
                                                   alpha_s, bornme);

    EXPECT_NEAR(limit, 0.0, 1e-15);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    limit = FKS::QCD::SoftCollinearLimitISR(realpdgs, bornpdgs, 1, s, alpha_s,
                                            bornme);

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

    int bornpdgs[] = { -2, 2, -13, 13 };
    int realpdgs[] = { 21, 2, -13, 13, -2 };
    double limit = FKS::QCD::SoftCollinearLimitISR(realpdgs, bornpdgs, -1, s,
                                                   alpha_s, bornme);

    EXPECT_NEAR(limit, 0.0, 1e-15);

    // second point
    // (     44.1287843,               0,               0,      44.1287843)
    // (     44.1287843,               0,               0,     -44.1287843)
    // (     44.1287843,    -33.72487942,     26.89162591,      9.31786283)
    // (     44.1287843,     33.72487942,    -26.89162591,     -9.31786283)
    bornme = 0.087679704752541803847;
    s = 7789.39841421111;
    limit = FKS::QCD::SoftCollinearLimitISR(realpdgs, bornpdgs, -1, s, alpha_s,
                                            bornme);
    EXPECT_NEAR(limit, 0.0, 1e-15);
}
