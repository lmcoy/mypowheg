#include "gtest/gtest.h"

#include <array>
#include <string>
#include <map>

#include "phasespace/realphasespace.h"
#include "phasespace/twoparticlegenerator.h"
#include "math/math.h"

class ISRTest : public ::testing::Test {
  protected:
    void addPS(const std::string &name, double S, double x1, double x2,
               double y, double phi, double xi_fks, double y_fks,
               double phi_fks) {
        /*ps_real[name] = Phasespace::Phasespace<3>(14000.0 * 14000.0,
                                                    { { 0.0, 0.0, 0.0 } });*/
        ps_born.insert(std::make_pair(name, new Phasespace::Phasespace()));
        ps_real.insert(std::make_pair(name, new Phasespace::Phasespace()));
        double yn = (y + 1.0) / 2.0;
        double params[] = { x1, x2, yn, phi / (2.0 * Math::Pi) };
        Phasespace::TwoParticleGenerator ps_gen;
        double masses[2] = { 0.0 };
        ps_gen(ps_born[name], 4, params, 14000.0 * 14000.0, 2, masses);
        Phasespace::GenRealPhasespaceISR(ps_real[name], ps_born[name],
                                            xi_fks, y_fks, phi_fks, false);
    }

    virtual void SetUp() {
        addPS("ps1", 14000.0 * 14000.0, 0.5, 0.5, 0.3, 0.0, 0.2, 0.3, 0.2);
        addPS("y=1,xi=0.1", 14000.0 * 14000.0, 0.4, 0.3, -0.2, 0.2, 0.1, 1.0,
              1.3);
        addPS("y=-1,xi=0.8", 14000.0 * 14000.0, 0.4, 0.3, -0.2, 0.2, 0.8, -1.0,
              1.3);
    }

    virtual void TearDown() {
        for (auto it = ps_born.begin(); it != ps_born.end(); ++it) {
            delete it->second;
        }
        for (auto it = ps_real.begin(); it != ps_real.end(); ++it) {
            delete it->second;
        }
    }

    std::map<std::string, Phasespace::Phasespace *> ps_born;
    std::map<std::string, Phasespace::Phasespace *> ps_real;
};

TEST_F(ISRTest, Rapidity_of_ktot) {
    for (auto &x : ps_real) {
        std::string name = x.first;
        Math::FourMomentum ktot =
            ps_real[name]->Momenta[0].Plus(ps_real[name]->Momenta[1]).Minus(
                ps_real[name]->Momenta[4]);

        double y = (ktot.E() + ktot.PZ()) / (ktot.E() - ktot.PZ());
        EXPECT_NEAR(y, 1.0, 1e-9) << "rapidity of ktot should be zero! ps = "
                                  << name; // <=> rapidity == 0
    }
}

TEST_F(ISRTest, InvMass_of_ktot) {
    for (auto &x : ps_real) {
        std::string name = x.first;
        Math::FourMomentum ktot =
            ps_real[name]->Momenta[0].Plus(ps_real[name]->Momenta[1]).Minus(
                ps_real[name]->Momenta[4]);

        double m2 = ktot.E() * ktot.E() - ktot.PX() * ktot.PX() -
                    ktot.PY() * ktot.PY() - ktot.PZ() * ktot.PZ();

        EXPECT_DOUBLE_EQ(m2, ps_born[name]->X1 * ps_born[name]->X2 *
                                 ps_born[name]->S)
            << "inv mass squared of ktot should be x1*x2*S, where x=born! ps = "
            << name;
    }
}

TEST_F(ISRTest, MomentumConservation) {
    for (auto &x : ps_real) {
        std::string name = x.first;
        Math::FourMomentum k_in =
            ps_real[name]->Momenta[0].Plus(ps_real[name]->Momenta[1]);
        Math::FourMomentum k_out =
            ps_real[name]->Momenta[2].Plus(ps_real[name]->Momenta[3]).Plus(
                ps_real[name]->Momenta[4]);
        EXPECT_NEAR(k_in.E(), k_out.E(), 1e-9);
        EXPECT_NEAR(k_in.PX(), k_out.PX(), 1e-9);
        EXPECT_NEAR(k_in.PY(), k_out.PY(), 1e-9);
        EXPECT_NEAR(k_in.PZ(), k_out.PZ(), 1e-9);
    }
}

TEST_F(ISRTest, CollinearLimit1) {
    double S = ps_born["y=1,xi=0.1"]->S;
    double x1 = ps_born["y=1,xi=0.1"]->X1;
    double x2 = ps_born["y=1,xi=0.1"]->X2;
    double xi = 0.1;
    // check x1, x2 of the n+1 phase space
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->X1, x1 / (1.0 - xi));
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->X2, x2);

    // test first initial state momentum
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[0].E(),
                     sqrt(S * x1 * x2) / (1.0 - xi) * 0.5);
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[0].PX(), 0.0);
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[0].PY(), 0.0);
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[0].PZ(),
                     sqrt(S * x1 * x2) / (1.0 - xi) * 0.5);

    // test second inital state momentum
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[1].E(),
                     sqrt(S * x1 * x2) * 0.5);
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[1].PX(), 0.0);
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[1].PY(), 0.0);
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[1].PZ(),
                     -sqrt(S * x1 * x2) * 0.5);

    // test if the fks parton is correct
    double k = 0.5 * sqrt(S * x2 * x1 / (1.0 - xi)) * xi;
    double ch = 1.0 / sqrt(1.0 - xi * xi / ((2.0 - xi) * (2.0 - xi)));
    double sh = 0.5 * xi / sqrt(1.0 - xi);
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[4].E(), k * (ch + sh));
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[4].PX(), 0.0);
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[4].PY(), 0.0);
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[4].PZ(), k * (ch + sh));

    // check if final state momenta are equal
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[2].E(),
                     ps_born["y=1,xi=0.1"]->Momenta[2].E());
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[2].PX(),
                     ps_born["y=1,xi=0.1"]->Momenta[2].PX());
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[2].PY(),
                     ps_born["y=1,xi=0.1"]->Momenta[2].PY());
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[2].PZ(),
                     ps_born["y=1,xi=0.1"]->Momenta[2].PZ());

    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[3].E(),
                     ps_born["y=1,xi=0.1"]->Momenta[3].E());
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[3].PX(),
                     ps_born["y=1,xi=0.1"]->Momenta[3].PX());
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[3].PY(),
                     ps_born["y=1,xi=0.1"]->Momenta[3].PY());
    EXPECT_DOUBLE_EQ(ps_real["y=1,xi=0.1"]->Momenta[3].PZ(),
                     ps_born["y=1,xi=0.1"]->Momenta[3].PZ());
}

TEST_F(ISRTest, CollinearLimit2) {
    double S = ps_born["y=-1,xi=0.8"]->S;
    double x1 = ps_born["y=-1,xi=0.8"]->X1;
    double x2 = ps_born["y=-1,xi=0.8"]->X2;
    double xi = 0.8;
    // check x1, x2 of the n+1 phase space
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->X1, x1);
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->X2, x2 / (1.0 - xi));

    // test first initial state momentum
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[0].E(),
                     sqrt(S * x1 * x2) * 0.5);
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[0].PX(), 0.0);
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[0].PY(), 0.0);
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[0].PZ(),
                     sqrt(S * x1 * x2) * 0.5);

    // test second inital state momentum
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[1].E(),
                     sqrt(S * x1 * x2) / (1.0 - xi) * 0.5);
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[1].PX(), 0.0);
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[1].PY(), 0.0);
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[1].PZ(),
                     -sqrt(S * x1 * x2) / (1.0 - xi) * 0.5);

    // test if the fks parton is correct
    double k = 0.5 * sqrt(S * x2 * x1 / (1.0 - xi)) * xi;
    double ch = 1.0 / sqrt(1.0 - xi * xi / ((2.0 - xi) * (2.0 - xi)));
    double sh = 0.5 * xi / sqrt(1.0 - xi);
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[4].E(), k * (ch + sh));
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[4].PX(), 0.0);
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[4].PY(), 0.0);
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[4].PZ(), -k * (ch + sh));

    // check if final state momenta are equal
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[2].E(),
                     ps_born["y=-1,xi=0.8"]->Momenta[2].E());
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[2].PX(),
                     ps_born["y=-1,xi=0.8"]->Momenta[2].PX());
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[2].PY(),
                     ps_born["y=-1,xi=0.8"]->Momenta[2].PY());
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[2].PZ(),
                     ps_born["y=-1,xi=0.8"]->Momenta[2].PZ());

    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[3].E(),
                     ps_born["y=-1,xi=0.8"]->Momenta[3].E());
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[3].PX(),
                     ps_born["y=-1,xi=0.8"]->Momenta[3].PX());
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[3].PY(),
                     ps_born["y=-1,xi=0.8"]->Momenta[3].PY());
    EXPECT_DOUBLE_EQ(ps_real["y=-1,xi=0.8"]->Momenta[3].PZ(),
                     ps_born["y=-1,xi=0.8"]->Momenta[3].PZ());
}

class FSRTest : public ::testing::Test {
  protected:
    void addPS(const std::string &name, double S, double x1, double x2,
               double y, double phi, double xi_fks, double y_fks,
               double phi_fks) {
        /*ps_real[name] = Phasespace::Phasespace<3>(14000.0 * 14000.0,
                                                    { { 0.0, 0.0, 0.0 } });*/
        ps_born.insert(std::make_pair(name, new Phasespace::Phasespace()));
        ps_real.insert(std::make_pair(name, new Phasespace::Phasespace()));
        double yn = (y + 1.0) / 2.0;
        double params[] = { x1, x2, yn, phi / (2.0 * Math::Pi) };
        Phasespace::TwoParticleGenerator ps_gen;
        double masses[2] = { 0.0 };
        ps_gen(ps_born[name], 4, params, 14000.0 * 14000.0, 2, masses);
        Phasespace::GenRealPhasespaceFSR(ps_real[name], ps_born[name], 2,
                                            xi_fks, y_fks, phi_fks);
    }

    virtual void SetUp() {
        addPS("ps1", 14000.0 * 14000.0, 0.5, 0.5, 0.3, 0.0, 0.2, 0.3, 0.2);
        addPS("y=1,xi=0.1", 14000.0 * 14000.0, 0.4, 0.3, -0.2, 0.2, 0.1, 1.0,
              1.3);
        addPS("y=-1,xi=0.8", 14000.0 * 14000.0, 0.4, 0.3, -0.2, 0.2, 0.8, -1.0,
              1.3);
    }

    virtual void TearDown() {
        for (auto it = ps_born.begin(); it != ps_born.end(); ++it) {
            delete it->second;
        }
        for (auto it = ps_real.begin(); it != ps_real.end(); ++it) {
            delete it->second;
        }
    }

    std::map<std::string, Phasespace::Phasespace *> ps_born;
    std::map<std::string, Phasespace::Phasespace *> ps_real;
};

TEST_F(FSRTest, MomentumConservation) {
    for (auto &x : ps_real) {
        std::string name = x.first;
        Math::FourMomentum k_in =
            ps_real[name]->Momenta[0].Plus(ps_real[name]->Momenta[1]);
        Math::FourMomentum k_out =
            ps_real[name]->Momenta[2].Plus(ps_real[name]->Momenta[3]).Plus(
                ps_real[name]->Momenta[4]);

        EXPECT_NEAR(k_in.E(), k_out.E(), 1e-11);
        EXPECT_NEAR(k_in.PX(), k_out.PX(), 1e-11);
        EXPECT_NEAR(k_in.PY(), k_out.PY(), 1e-11);
        EXPECT_NEAR(k_in.PZ(), k_out.PZ(), 1e-11);
    }
}

TEST(FSR, MG5) {
    Phasespace::Phasespace ps2;
    Phasespace::Phasespace ps3;

    ps2.X1 = 1.0;
    ps2.X2 = 1.0;
    ps2.Momenta[0].Set(7000.0, 0.0, 0.0, 7000.0);
    ps2.Momenta[1].Set(7000.0, 0.0, 0.0, -7000.0);
    ps2.Momenta[2].Set(7000, -6951.9736739841346, -557.16607125299265,
                       599.68992427423916);
    ps2.Momenta[3]
        .Set(7000, 6951.9736739841346, 557.16607125299265, -599.68992427423916);
    ps2.Masses[0] = 0.0;
    ps2.Masses[1] = 0.0;
    ps2.N = 2;
    ps2.S = 14000.0 * 14000.0;

    double xi_fks = 0.10245306005308354;
    double y_fks = 0.54946672918388273;
    double phi_fks = 2.4399270425265396;
    Phasespace::GenRealPhasespaceFSR(&ps3, &ps2, 3, xi_fks, y_fks, phi_fks);

    EXPECT_NEAR(ps3.Momenta[2].E(), 6851.571368, 1e-6);
    EXPECT_NEAR(ps3.Momenta[2].PX(), -6804.563396, 1e-6);
    EXPECT_NEAR(ps3.Momenta[2].PY(), -545.351872, 1e-6);
    EXPECT_NEAR(ps3.Momenta[2].PZ(), 586.974045, 1e-6);

    EXPECT_NEAR(ps3.Momenta[3].E(), 6431.257212, 1e-6);
    EXPECT_NEAR(ps3.Momenta[3].PX(), 6354.979745, 1e-6);
    EXPECT_NEAR(ps3.Momenta[3].PY(), 145.100174, 1e-6);
    EXPECT_NEAR(ps3.Momenta[3].PZ(), -976.856034, 1e-6);

    EXPECT_NEAR(ps3.Momenta[4].E(), 717.171420, 1e-6);
    EXPECT_NEAR(ps3.Momenta[4].PX(), 449.583651, 1e-6);
    EXPECT_NEAR(ps3.Momenta[4].PY(), 400.251697, 1e-6);
    EXPECT_NEAR(ps3.Momenta[4].PZ(), 389.881989, 1e-6);
}

TEST(FSR, MG5_2) {
    Phasespace::Phasespace ps2;
    Phasespace::Phasespace ps3;

    ps2.X1 = 1.0;
    ps2.X2 = 1.0;
    ps2.Momenta[0].Set(7000.0, 0.0, 0.0, 7000.0);
    ps2.Momenta[1].Set(7000.0, 0.0, 0.0, -7000.0);
    ps2.Momenta[2].Set(7000, -4782.8060804346505, -2561.3445346656608,
                       4423.1527411672869);
    ps2.Momenta[3]
        .Set(7000, 4782.8060804346505, 2561.3445346656608, -4423.1527411672869);
    ps2.Masses[0] = 0.0;
    ps2.Masses[1] = 0.0;
    ps2.N = 2;
    ps2.S = 14000.0 * 14000.0;

    double xi_fks = 0.042751806635932948;
    double y_fks = 0.7832812878648463;
    double phi_fks = 3.5384559880513438;
    Phasespace::GenRealPhasespaceFSR(&ps3, &ps2, 3, xi_fks, y_fks, phi_fks);

    EXPECT_NEAR(ps3.Momenta[2].E(), 6968.8139728867091, 1e-6);
    EXPECT_NEAR(ps3.Momenta[2].PX(), -4761.4979775629299, 1e-6);
    EXPECT_NEAR(ps3.Momenta[2].PY(), -2549.933368936438, 1e-6);
    EXPECT_NEAR(ps3.Momenta[2].PZ(), 4403.4469466941055, 1e-6);

    EXPECT_NEAR(ps3.Momenta[3].E(), 6731.9233806617603, 1e-6);
    EXPECT_NEAR(ps3.Momenta[3].PX(), 4472.8789400922742, 1e-6);
    EXPECT_NEAR(ps3.Momenta[3].PY(), 2474.1686239631535, 1e-6);
    EXPECT_NEAR(ps3.Momenta[3].PZ(), -4380.7118155130629, 1e-6);

    EXPECT_NEAR(ps3.Momenta[4].E(), 299.26264645153066, 1e-6);
    EXPECT_NEAR(ps3.Momenta[4].PX(), 288.61903747065594, 1e-6);
    EXPECT_NEAR(ps3.Momenta[4].PY(), 75.764744973283655, 1e-6);
    EXPECT_NEAR(ps3.Momenta[4].PZ(), -22.735131181043158, 1e-6);
}

TEST(FSR, MG5_3) {
    Phasespace::Phasespace ps2;
    Phasespace::Phasespace ps3;

    ps2.X1 = 1.0;
    ps2.X2 = 1.0;
    ps2.Momenta[0].Set(7000.0, 0.0, 0.0, 7000.0);
    ps2.Momenta[1].Set(7000.0, 0.0, 0.0, -7000.0);
    ps2.Momenta[2].Set(7000, -2459.2041467673243, -6526.8283439648185,
                       -593.99219939311274);
    ps2.Momenta[3]
        .Set(7000, 2459.2041467673243, 6526.8283439648185, 593.99219939311274);
    ps2.Masses[0] = 0.0;
    ps2.Masses[1] = 0.0;
    ps2.N = 2;
    ps2.S = 14000.0 * 14000.0;

    double xi_fks = 0.17290220634641654;
    double y_fks = 0.89945170467629232;
    double phi_fks = 1.4765429233333769;
    Phasespace::GenRealPhasespaceFSR(&ps3, &ps2, 2, xi_fks, y_fks, phi_fks);

    EXPECT_NEAR(ps3.Momenta[2].E(), 5840.4527557081346, 1e-6);
    EXPECT_NEAR(ps3.Momenta[2].PX(), -2461.2814457686359, 1e-6);
    EXPECT_NEAR(ps3.Momenta[2].PY(), -5277.1440862179461, 1e-6);
    EXPECT_NEAR(ps3.Momenta[2].PZ(), -452.47356792276128, 1e-6);

    EXPECT_NEAR(ps3.Momenta[3].E(), 6949.2317998669487, 1e-6);
    EXPECT_NEAR(ps3.Momenta[3].PX(), 2441.368522725737, 1e-6);
    EXPECT_NEAR(ps3.Momenta[3].PY(), 6479.4918685933226, 1e-6);
    EXPECT_NEAR(ps3.Momenta[3].PZ(), 589.68421155650412, 1e-6);

    EXPECT_NEAR(ps3.Momenta[4].E(), 1210.3154444249158, 1e-6);
    EXPECT_NEAR(ps3.Momenta[4].PX(), 19.912923042898399, 1e-6);
    EXPECT_NEAR(ps3.Momenta[4].PY(), -1202.3477823753756, 1e-6);
    EXPECT_NEAR(ps3.Momenta[4].PZ(), -137.21064363374273, 1e-6);
}

TEST(FSR, SoftLimit) {
    Phasespace::Phasespace ps;
    Phasespace::Phasespace ps_real;

    double params[] = { 0.3, 0.4, 0.5, 0.6 };
    Phasespace::TwoParticleGenerator ps_gen;
    double masses[2] = { 0.0 };
    ps_gen(&ps, 4, params, 14000.*14000., 2 , masses);

    double y = 0.6;
    double xi = 1e-25;
    Phasespace::GenRealPhasespaceFSR(&ps_real, &ps, 2, xi, y, 0.6);

    EXPECT_NEAR(ps_real.Momenta[4].E(), 0.0, 1e-20);
    EXPECT_NEAR(ps_real.Momenta[4].PX(), 0.0, 1e-20);
    EXPECT_NEAR(ps_real.Momenta[4].PY(), 0.0, 1e-20);
    EXPECT_NEAR(ps_real.Momenta[4].PZ(), 0.0, 1e-20);
}

TEST(FSR, CollinearLimit) {
    Phasespace::Phasespace ps;
    Phasespace::Phasespace ps_real;

    double params[] = { 0.3, 0.4, 0.5, 0.6 };
    Phasespace::TwoParticleGenerator ps_gen;
    double masses[2] = { 0.0 };
    ps_gen(&ps, 4, params, 14000.*14000., 2 , masses);

    double y = 1.0;
    double xi = 0.8;
    Phasespace::GenRealPhasespaceFSR(&ps_real, &ps, 2, xi, y, 0.6);

    EXPECT_NEAR(ps_real.Momenta[2].PX() / ps_real.Momenta[2].E(),
                ps_real.Momenta[2].PX() / ps_real.Momenta[2].E(), 1e-20);
    EXPECT_NEAR(ps_real.Momenta[2].PY() / ps_real.Momenta[2].E(),
                ps_real.Momenta[2].PY() / ps_real.Momenta[2].E(), 1e-20);
    EXPECT_NEAR(ps_real.Momenta[2].PZ() / ps_real.Momenta[2].E(),
                ps_real.Momenta[2].PZ() / ps_real.Momenta[2].E(), 1e-20);
}
