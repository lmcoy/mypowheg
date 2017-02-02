#include "gtest/gtest.h"

#include <array>
#include <map>
#include <string>

#include "math/math.h"
#include "phasespace/realphasespace.h"
#include "phasespace/twoparticlegenerator.h"

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
        double params[] = {x1, x2, yn, phi / (2.0 * Math::Pi)};
        Phasespace::TwoParticleGenerator ps_gen;
        double masses[2] = {0.0};
        ps_gen(ps_born[name], 4, params, 14000.0 * 14000.0, 2, masses);
        Phasespace::GenRealPhasespaceISR(ps_real[name], ps_born[name], xi_fks,
                                         y_fks, phi_fks, false);
    }

    void addPS3(const std::string &name, double xi_fks, double y_fks,
                double phi_fks) {
        ps_born.insert(std::make_pair(name, new Phasespace::Phasespace()));
        ps_real.insert(std::make_pair(name, new Phasespace::Phasespace()));
        ps_born[name]->S = 196000000;
        ps_born[name]->X1 = 0.9011708913176371;
        ps_born[name]->X2 = 0.3714679254657938;
        ps_born[name]->Jacobian = 2193.672575057868;
        ps_born[name]->N = 3;
        ps_born[name]->Momenta[0].Set(4050.067652880409, 0, 0,
                                      4050.067652880409);
        ps_born[name]->Momenta[1].Set(4050.067652880409, 0, 0,
                                      -4050.067652880409);
        ps_born[name]->Momenta[2].Set(3967.564009595429, 2243.913315465149,
                                      -3208.644423221378, -641.106986582977);
        ps_born[name]->Momenta[3].Set(3595.157112400982, -1913.998123283016,
                                      2897.798247481882, 929.8016798260696);
        ps_born[name]->Momenta[4].Set(537.4141837644091, -329.915192182132,
                                      310.8461757394962, -288.6946932430926);
        Phasespace::GenRealPhasespaceISR(ps_real[name], ps_born[name], xi_fks,
                                         y_fks, phi_fks, false);
    }

    virtual void SetUp() {
        addPS("ps1", 14000.0 * 14000.0, 0.5, 0.5, 0.3, 0.0, 0.2, 0.3, 0.2);
        addPS("y=1,xi=0.1", 14000.0 * 14000.0, 0.4, 0.3, -0.2, 0.2, 0.1, 1.0,
              1.3);
        addPS("y=-1,xi=0.8", 14000.0 * 14000.0, 0.4, 0.3, -0.2, 0.2, 0.8, -1.0,
              1.3);
        addPS3("3 particle", 0.2, 0.4, 3.0);
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
        int N = ps_real[name]->N + 1;
        Math::FourMomentum ktot = ps_real[name]
                                      ->Momenta[0]
                                      .Plus(ps_real[name]->Momenta[1])
                                      .Minus(ps_real[name]->Momenta[N]);

        double y = (ktot.E() + ktot.PZ()) / (ktot.E() - ktot.PZ());
        EXPECT_NEAR(y, 1.0, 1e-9) << "rapidity of ktot should be zero! ps = "
                                  << name; // <=> rapidity == 0
    }
}

TEST_F(ISRTest, InvMass_of_ktot) {
    for (auto &x : ps_real) {
        std::string name = x.first;
        int N = ps_real[name]->N + 1;
        Math::FourMomentum ktot = ps_real[name]
                                      ->Momenta[0]
                                      .Plus(ps_real[name]->Momenta[1])
                                      .Minus(ps_real[name]->Momenta[N]);

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
        Math::FourMomentum k_out(0.0, 0.0, 0.0, 0.0);
        for (int f = 2; f < ps_real[name]->N + 2; f++) {
            k_out = k_out.Plus(ps_real[name]->Momenta[f]);
        }
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
                                                    { { 0.0, 0.0, 0.0 }
           });*/
        ps_born.insert(std::make_pair(name, new Phasespace::Phasespace()));
        ps_real.insert(std::make_pair(name, new Phasespace::Phasespace()));
        double yn = (y + 1.0) / 2.0;
        double params[] = {x1, x2, yn, phi / (2.0 * Math::Pi)};
        Phasespace::TwoParticleGenerator ps_gen;
        double masses[2] = {0.0};
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
        Math::FourMomentum k_out = ps_real[name]
                                       ->Momenta[2]
                                       .Plus(ps_real[name]->Momenta[3])
                                       .Plus(ps_real[name]->Momenta[4]);

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
    ps2.Momenta[3].Set(7000, 6951.9736739841346, 557.16607125299265,
                       -599.68992427423916);
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
    ps2.Momenta[3].Set(7000, 4782.8060804346505, 2561.3445346656608,
                       -4423.1527411672869);
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
    ps2.Momenta[3].Set(7000, 2459.2041467673243, 6526.8283439648185,
                       593.99219939311274);
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

    double params[] = {0.3, 0.4, 0.5, 0.6};
    Phasespace::TwoParticleGenerator ps_gen;
    double masses[2] = {0.0};
    ps_gen(&ps, 4, params, 14000. * 14000., 2, masses);

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

    double params[] = {0.3, 0.4, 0.5, 0.6};
    Phasespace::TwoParticleGenerator ps_gen;
    double masses[2] = {0.0};
    ps_gen(&ps, 4, params, 14000. * 14000., 2, masses);

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

// Compare to MadFKS real phase spaces
TEST(FSR, ThreeParticleBorn) {
    Phasespace::Phasespace ps;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;

    ps.S = 169000000.00000000;
    ps.X1 = 1.9973163356326647E-004;
    ps.X2 = 0.53547277258683934;
    ps.Jacobian = 1.7289366081869456E-006;
    ps.Momenta[0].Set(67.221079134383388, 0.0000000000000000,
                      0.0000000000000000, 67.221079134383388);
    ps.Momenta[1].Set(67.221079134383388, -0.0000000000000000,
                      -0.0000000000000000, -67.221079134383388);
    ps.Momenta[2].Set(29.073144891019879, 2.9341793871702926,
                      -10.707444535334993, 26.869852561948775);
    ps.Momenta[3].Set(61.770698373824040, -39.711434164250534,
                      -0.71631537177801174, -47.308646847024661);
    ps.Momenta[4].Set(43.598315003922856, 36.777254777080245,
                      11.423759907113004, 20.438794285075883);
    double xi_fks = 0.13223109638959257;
    double y_fks = 0.12475850527718457;
    double phi_fks = 2.4220908743555327;

    Phasespace::Phasespace ps_real;
    Phasespace::GenRealPhasespaceFSR(&ps_real, &ps, 4, xi_fks, y_fks, phi_fks);

    EXPECT_NEAR(ps_real.Momenta[0].E(), 67.221079134383388, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[0].PX(), 0.0, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[0].PY(), 0.0, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[0].PZ(), 67.221079134383388, 1e-12);

    EXPECT_NEAR(ps_real.Momenta[1].E(), 67.221079134383388, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[1].PX(), 0.0, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[1].PY(), 0.0, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[1].PZ(), -67.221079134383388, 1e-12);

    EXPECT_NEAR(ps_real.Momenta[2].E(), 29.745846866169902, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[2].PX(), 4.2154066445473806, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[2].PY(), -10.309469410126729, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[2].PZ(), 27.581888860054473, 1e-12);

    EXPECT_NEAR(ps_real.Momenta[3].E(), 58.966086220276239, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[3].PX(), -37.081479879610512, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[3].PY(), 0.10060167373062345, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[3].PZ(), -45.847061556639950, 1e-12);

    EXPECT_NEAR(ps_real.Momenta[4].E(), 36.841508187889566, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[4].PX(), 34.709263141478004, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[4].PY(), 5.0266663797619033, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[4].PZ(), 11.282570754953923, 1e-12);

    EXPECT_NEAR(ps_real.Momenta[5].E(), 8.8887169944310802, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[5].PX(), -1.8431899064148523, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[5].PY(), 5.1822013566342084, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[5].PZ(), 6.9826019416315628, 1e-12);
}

// compare results from MadFKS
//
// You have to modify genps_fks.f in the SubProcesses directory to print the
// phase spaces.
// The easiest way to compile the new genps_fks.f is to run the normal cross
// section calculation for fixed order NLO. Then change to
// SubProcesses/P*/all_G1 and run ../madevent_mintFO.
// Note that the real jacobian is  xpswgt*xi_i_fks.
TEST(ThreeParticleBorn, MadFKS1) {
    Phasespace::Phasespace ps;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;

    ps.S = 169000000.00000000;
    ps.X1 = 7.7511082039846799E-004;
    ps.X2 = 0.60384494721413340;
    ps.Jacobian = 3.9507002589821781E-007;
    ps.Momenta[0].Set(140.62352324596711, 0.0000000000000000,
                      0.0000000000000000, 140.62352324596711);
    ps.Momenta[1].Set(140.62352324596711, -0.0000000000000000,
                      -0.0000000000000000, -140.62352324596711);
    ps.Momenta[2].Set(109.90677488326632, 81.164604262401951,
                      50.794232920195824, 53.960653092948746);
    ps.Momenta[3].Set(41.620030490572447, 21.156407719843184,
                      -27.953174686826962, 22.433309505190600);
    ps.Momenta[4].Set(129.72024111809543, -102.32101198224514,
                      -22.841058233368866, -76.393962598139353);

    double xi_fks = 0.26597711351648368;
    double y_fks = 0.78232389650413814;
    double phi_fks = 1.2670037155846734;
    int j_fks = 4;
    int i_fks = 5;

    Phasespace::Phasespace ps_real;
    Phasespace::GenRealPhasespace(&ps_real, &ps, i_fks, j_fks, xi_fks, y_fks,
                                  phi_fks);

    EXPECT_NEAR(ps_real.X1, 7.7511082039846799E-004, 1e-12);
    EXPECT_NEAR(ps_real.X2, 0.60384494721413340, 1e-12);
    EXPECT_NEAR(ps_real.Jacobian, 3.1612357056018420E-006, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[0].At(0), 140.62352324596711, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[0].At(1), 0.0000000000000000, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[0].At(2), 0.0000000000000000, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[0].At(3), 140.62352324596711, 1e-12);

    EXPECT_NEAR(ps_real.Momenta[1].At(0), 140.62352324596711, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[1].At(1), -0.0000000000000000, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[1].At(2), -0.0000000000000000, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[1].At(3), -140.62352324596711, 1e-12);

    EXPECT_NEAR(ps_real.Momenta[2].At(0), 107.68169279590836, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[2].At(1), 79.321049780039189, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[2].At(2), 50.382697357047341, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[2].At(3), 52.584235586467372, 1e-12);

    EXPECT_NEAR(ps_real.Momenta[3].At(0), 41.092993317552370, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[3].At(1), 20.455607868556580, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[3].At(2), -28.109613817979721, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[3].At(3), 21.910084836172860, 1e-12);

    EXPECT_NEAR(ps_real.Momenta[4].At(0), 95.069721572992975, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[4].At(1), -80.357048341669184, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[4].At(2), -0.81920968266448213, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[4].At(3), -50.796905784475463, 1e-12);

    EXPECT_NEAR(ps_real.Momenta[5].At(0), 37.402638805480471, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[5].At(1), -19.419609306926571, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[5].At(2), -21.453873856403163, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[5].At(3), -23.697414638164762, 1e-12);
}

TEST(ThreeParticleBorn, MadFKS2) {
    Phasespace::Phasespace ps;
    ps.N = 3;
    ps.Masses[0] = 0.0;
    ps.Masses[1] = 0.0;
    ps.Masses[2] = 0.0;
    ps.S = 169000000.00000000;
    ps.X1 = 6.3029924593260302E-004;
    ps.X2 = 8.2646357127841015E-004;
    ps.Jacobian = 3.5497095089005390E-004;
    ps.Momenta[0].Set(4.6913583537903927, 0.0000000000000000,
                      0.0000000000000000, 4.6913583537903927);
    ps.Momenta[1].Set(4.6913583537903927, -0.0000000000000000,
                      -0.0000000000000000, -4.6913583537903927);
    ps.Momenta[2].Set(4.6214879548064696, 1.8654124695158638,
                      2.5403532117984744, -3.3800876903271595);
    ps.Momenta[3].Set(4.2562682411694439, -1.9813438347292822,
                      -2.1332742194222538, 3.1047120726592787);
    ps.Momenta[4].Set(0.50496051160487099, 0.11593136521341833,
                      -0.40707899237622081, 0.27537561766788071);

    double xi_fks = 0.86069984658528009;
    double y_fks = 0.26260143699680361;
    double phi_fks = 4.4613542977711003;
    int j_fks = 0; //            0
    int i_fks = 5;

    Phasespace::Phasespace ps_real;
    Phasespace::GenRealPhasespace(&ps_real, &ps, i_fks, j_fks, xi_fks, y_fks,
                                  phi_fks, false);

    EXPECT_NEAR(ps_real.X1, 2.0648410898136565E-003, 1e-12);
    EXPECT_NEAR(ps_real.X2, 1.8110576472135846E-003, 1e-12);
    EXPECT_NEAR(ps_real.Jacobian, 6.9850410979795552E-004, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[0].At(0), 15.368746763474254, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[0].At(1), 0.0000000000000000, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[0].At(2), 0.0000000000000000, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[0].At(3), 15.368746763474254, 1e-12);

    EXPECT_NEAR(ps_real.Momenta[1].At(0), 10.280332633789223, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[1].At(1), 0.0000000000000000, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[1].At(2), 0.0000000000000000, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[1].At(3), -10.280332633789223, 1e-12);

    EXPECT_NEAR(ps_real.Momenta[2].At(0), 10.166753011120161, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[2].At(1), 3.5028915001434799, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[2].At(2), 8.9256722510431796, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[2].At(3), -3.3800876903271595, 1e-12);

    EXPECT_NEAR(ps_real.Momenta[3].At(0), 3.5204843293515244, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[3].At(1), -1.1202360606574138, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[3].At(2), 1.2245995376024537, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[3].At(3), 3.1047120726592787, 1e-12);

    EXPECT_NEAR(ps_real.Momenta[4].At(0), 0.34871843857827256, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[4].At(1), 0.21045791391072000, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[4].At(2), -3.8474472976891494E-002, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[4].At(3), 0.27537561766788071, 1e-12);

    EXPECT_NEAR(ps_real.Momenta[5].At(0), 11.613123618213518, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[5].At(1), -2.5931133533967867, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[5].At(2), -10.111797315668744, 1e-12);
    EXPECT_NEAR(ps_real.Momenta[5].At(3), 5.0884141296850309, 1e-12);
}
