#include "gtest/gtest.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include <mpi.h>

#include "LHAPDF/LHAPDF.h"

#include "math/fourmomentum.h"
#include "me/parameters_sm.h"
#include "me/wg.h"
#include "phasespace/phasespace.h"
#include "random/rnd.h"
#include "vegas/vegas.h"

#include "phasespace.h"

class PSTest : public ::testing::TestWithParam<std::array<double, 7>> {};

struct Data {
    double SqrtS = 13000.0;
    LHAPDF::PDF *pdf = 0;
};

const double CrossSectionToPb = 3.8937944929e8;

bool applycuts(const Phasespace::Phasespace &ps, const Data &data) {
    // cut on inv. mass of final state particles
    const auto &momenta = ps.Momenta;
    double pT = momenta[4].PT();
    double eta = momenta[4].Eta();
    if (fabs(eta) > 1e8) {
        return false;
    }
    if (pT < 10.0) {
        return false;
    }

    return true;
}

double ME2(const Phasespace::Phasespace &ps, const Parameters_sm &p) {
    Wg me;

    int perm[] = {0, 1, 2, 3, 4};

    auto m = me.Calculate(ps, perm, p);

    return m.Born();
}

// x1, x2, x3 -- for NLO real radiation, can be ignored here
// wgt -- VEGAS weight
int XSecBorn(const Phasespace::Phasespace &ps, double x1, double x2, double x3,
             double wgt, double *out, Data *params) {

    Phasespace::Phasespace ps_lab;
    ps_lab.SetToLabFromCMS(&ps);
    bool cuts = applycuts(ps_lab, *params);
    if (!cuts) {
        *out = 0.0;
        return 0;
    }

    double muF = 91.188;
    double f1 = params->pdf->xfxQ(2, ps.X1, muF) / ps.X1;  // LHAPDF
    double f2 = params->pdf->xfxQ(-1, ps.X2, muF) / ps.X2; // LHAPDF
    double lumi = f1 * f2;

    f1 = params->pdf->xfxQ(4, ps.X1, muF) / ps.X1;  // LHAPDF
    f2 = params->pdf->xfxQ(-3, ps.X2, muF) / ps.X2; // LHAPDF
    lumi += f1 * f2;

    Parameters_sm param;
    double aew = 1.325070e+02;
    double mz = 9.1188e+01;
    param.MW = 8.04190024457561634108e+01;
    param.WidthW = 2.0476;

    double Gf = M_PI / sqrt(2.0) / aew / param.MW / param.MW *
                (1.0 / (1.0 - param.MW * param.MW / mz / mz));

    param.Set(Gf, aew, mz, 2.441404);

    param.SetAlphaS(params->pdf->alphasQ(muF));

    double me2 = ME2(ps, param); // squared matrix element

    double s = ps.X1 * ps.X2 * ps.S;
    double dsigma = 1.0 / (2.0 * s) * me2;

    double result_born = dsigma * CrossSectionToPb;

    *out = lumi * result_born;
    return 0;
}

int One(const Phasespace::Phasespace &ps, double x1, double x2, double x3,
        double wgt, double *out, Data *params) {

    *out = 1.0;
    return 0;
}

typedef int (*function)(const Phasespace::Phasespace &, double, double, double,
                        double, double *, Data *);

double clip(double x) {
    if (x >= 1.0 && x < 1.0 + 1e-9) {
        return 1.0 - 1e-12;
    }
    if (x <= 0.0 && x > -1e-9) {
        return 1e-12;
    }
    return x;
}

template <function func>
static int integrand_simple(int n, const double *x, const double wgt,
                            void *params, int threadid, double *out) {

    Data *userdata = (Data *)params;
    double xx[n];

    for (int i = 0; i < n; i++) {
        xx[i] = clip(x[i]);
    }

    std::array<double, 7> v;
    for (int i = 0; i < 7; i++) {
        v[i] = xx[i];
    }

    double S = userdata->SqrtS * userdata->SqrtS;
    Phasespace::Phasespace ps = GenPhasespace(S, v);

    double integrand = 0.0;
    int ret = func(ps, 0.0, 0.0, 0.0, wgt, &integrand, userdata);

    *out = integrand * ps.Jacobian;
    return ret;
}

template <function func>
static int integrand_breitwigner(int n, const double *x, const double wgt,
                                 void *params, int threadid, double *out) {

    Data *userdata = (Data *)params;
    double xx[n];

    for (int i = 0; i < n; i++) {
        xx[i] = clip(x[i]);
    }

    std::array<double, 7> v;
    for (int i = 0; i < 7; i++) {
        v[i] = xx[i];
    }

    double S = userdata->SqrtS * userdata->SqrtS;
    double M = 80.385;
    double G = 2.085;

    double a = v[1];
    double b = v[2];

    v[1] = b * (1.0 - a) + a;
    v[2] = a / v[1];

    double jac = (1.0 - a) / v[1];

    assert(v[1] >= 0.0 && v[1] <= 1.0);
    assert(v[2] >= 0.0 && v[2] <= 1.0);

    double Mll = 1.0;
    double Delta =
        atan((S - M * M) / (M * G)) - atan((Mll * Mll - M * M) / (M * G));
    double s = v[0];
    double z = v[1];
    double ts = Delta * s + atan((Mll * Mll - M * M) / (M * G));
    double tau = clip((G * tan(ts) + M) * M / S);
    LIB_ASSERT(tau >= 0.0 && tau <= 1.0, "tau = %.16g\n", tau);
    double x1 = (1.0 - tau) * z + tau;
    v[0] = clip(x1);
    v[1] = clip(tau / x1);

    double tmp = tau * S - M * M;
    double inv_breit_wigner = (M * M * G * G + tmp * tmp) / (M * G * S);
    jac *= Delta * (1.0 - tau) / x1 * inv_breit_wigner;

    Phasespace::Phasespace ps = GenPhasespace(S, v);

    double integrand = 0.0;
    int ret = func(ps, 0.0, 0.0, 0.0, wgt * jac, &integrand, userdata);

    *out = integrand * ps.Jacobian * jac;
    return ret;
}

TEST(CrossSection, SimpleTrafo) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // init pdfs
    if (rank != 0) {
        LHAPDF::setVerbosity(0);
    }
    auto pdf = LHAPDF::mkPDF(244600);

    const int NDIM = 10; // 4d integral x_1, x_2, and 4 angles and Q^2 of
                         // intermediate particle
    VegasState *state = vegas_new(NDIM, 50);
    vegas_set_random(state, Random::init_rnd, NULL, Random::get_rnd,
                     Random::free_rnd);

    vegas_integrand integrand_ptr = 0;

    integrand_ptr = integrand_simple<XSecBorn>;

    Data userdata;
    userdata.pdf = pdf;
    userdata.SqrtS = 13000.0;

    // -----------------------------------------
    // VEGAS grid setup
    // -----------------------------------------
    int verbosity = 1;

    int seed1 = 1; // random seed
    vegas_seed_random(state, seed1);

    // integration result
    double integral = 0.0, error = 0.0, prob = 0.0;
    // number of vegas iterations and integrand evaluations per interation
    int N_setup_points = 100000;
    int N_setup_iterations = 5;
    if (rank == 0) {
        printf("***************************************************************"
               "****\n");
        printf("*                                                              "
               "   *\n");
        printf("* setup VEGAS grid                                             "
               "   *\n");
        printf("*                                                              "
               "   *\n");
        printf("***************************************************************"
               "****\n");
    }
    vegas_integrate(state, seed1, integrand_ptr, (void *)&userdata, 1,
                    N_setup_points, N_setup_iterations, verbosity, &integral,
                    &error, &prob);
    // values in integral, error etc. should be discarded.

    // -----------------------------------------
    // VEGAS integration
    // -----------------------------------------
    vegas_reset_int(state); // only keep grid from previous integration
    vegas_disable_update_grid(state); // don't update grid during integration

    int seed2 = 2; // seed for VEGAS integration
    vegas_seed_random(state, seed2);
    int N_points = 1000000;
    int N_iterations = 10;
    if (rank == 0) {
        printf("***************************************************************"
               "****\n");
        printf("*                                                              "
               "   *\n");
        printf("* VEGAS integration                                            "
               "   *\n");
        printf("*                                                              "
               "   *\n");
        printf("***************************************************************"
               "****\n");
    }
    vegas_integrate(state, seed2, integrand_ptr, (void *)&userdata, 1, N_points,
                    N_iterations, verbosity, &integral, &error, &prob);

    if (rank == 0) {
        printf("***************************************************************"
               "****\n");
        printf("*                                                              "
               "   *\n");
        printf("* result                                                       "
               "   *\n");
        printf("*                                                              "
               "   *\n");
        printf("***************************************************************"
               "****\n");
        printf("\n");
        printf("I = %g +- %g (chi^2/ndf = %g)\n", integral, error, prob);
        double pre = 839.5; // MadGraph result
        double e = sqrt(0.65 * 0.65 + error * error);
        ASSERT_NEAR(integral, pre, 5 * e)
            << "result doesn't agree with MadGraph";
    }

    delete pdf;
}

TEST(CrossSection, BreitWignerTrafo) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // init pdfs
    if (rank != 0) {
        LHAPDF::setVerbosity(0);
    }
    auto pdf = LHAPDF::mkPDF(244600);

    const int NDIM = 10; // 4d integral x_1, x_2, and 4 angles and Q^2 of
                         // intermediate particle
    VegasState *state = vegas_new(NDIM, 50);
    vegas_set_random(state, Random::init_rnd, NULL, Random::get_rnd,
                     Random::free_rnd);

    vegas_integrand integrand_ptr = 0;

    integrand_ptr = integrand_breitwigner<XSecBorn>;

    Data userdata;
    userdata.pdf = pdf;
    userdata.SqrtS = 13000.0;

    // -----------------------------------------
    // VEGAS grid setup
    // -----------------------------------------
    int verbosity = 1;

    int seed1 = 1; // random seed
    vegas_seed_random(state, seed1);

    // integration result
    double integral = 0.0, error = 0.0, prob = 0.0;
    // number of vegas iterations and integrand evaluations per interation
    int N_setup_points = 100000;
    int N_setup_iterations = 5;
    if (rank == 0) {
        printf("***************************************************************"
               "****\n");
        printf("*                                                              "
               "   *\n");
        printf("* setup VEGAS grid                                             "
               "   *\n");
        printf("*                                                              "
               "   *\n");
        printf("***************************************************************"
               "****\n");
    }
    vegas_integrate(state, seed1, integrand_ptr, (void *)&userdata, 1,
                    N_setup_points, N_setup_iterations, verbosity, &integral,
                    &error, &prob);
    // values in integral, error etc. should be discarded.

    // -----------------------------------------
    // VEGAS integration
    // -----------------------------------------
    vegas_reset_int(state); // only keep grid from previous integration
    vegas_disable_update_grid(state); // don't update grid during integration

    int seed2 = 2; // seed for VEGAS integration
    vegas_seed_random(state, seed2);
    int N_points = 1000000;
    int N_iterations = 10;
    if (rank == 0) {
        printf("***************************************************************"
               "****\n");
        printf("*                                                              "
               "   *\n");
        printf("* VEGAS integration                                            "
               "   *\n");
        printf("*                                                              "
               "   *\n");
        printf("***************************************************************"
               "****\n");
    }
    vegas_integrate(state, seed2, integrand_ptr, (void *)&userdata, 1, N_points,
                    N_iterations, verbosity, &integral, &error, &prob);

    if (rank == 0) {
        printf("***************************************************************"
               "****\n");
        printf("*                                                              "
               "   *\n");
        printf("* result                                                       "
               "   *\n");
        printf("*                                                              "
               "   *\n");
        printf("***************************************************************"
               "****\n");
        printf("\n");
        printf("I = %g +- %g (chi^2/ndf = %g)\n", integral, error, prob);
        double pre = 839.5; // MadGraph result
        double e = sqrt(0.65 * 0.65 + error * error);
        ASSERT_NEAR(integral, pre, 5 * e)
            << "result doesn't agree with MadGraph";
    }

    delete pdf;
}

TEST(PS, Volume) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int NDIM = 7; // 4d integral x_1, x_2, and 4 angles and Q^2 of
                        // intermediate particle
    VegasState *state = vegas_new(NDIM, 50);
    vegas_set_random(state, Random::init_rnd, NULL, Random::get_rnd,
                     Random::free_rnd);

    vegas_integrand integrand_ptr = 0;

    integrand_ptr = integrand_simple<One>;

    Data userdata;
    userdata.SqrtS = 13000.0;

    // -----------------------------------------
    // VEGAS grid setup
    // -----------------------------------------
    int verbosity = 1;

    int seed1 = 1; // random seed
    vegas_seed_random(state, seed1);

    // integration result
    double integral = 0.0, error = 0.0, prob = 0.0;
    // number of vegas iterations and integrand evaluations per interation
    int N_setup_points = 10000;
    int N_setup_iterations = 5;
    if (rank == 0) {
        printf("***************************************************************"
               "****\n");
        printf("*                                                              "
               "   *\n");
        printf("* setup VEGAS grid                                             "
               "   *\n");
        printf("*                                                              "
               "   *\n");
        printf("***************************************************************"
               "****\n");
    }
    vegas_integrate(state, seed1, integrand_ptr, (void *)&userdata, 1,
                    N_setup_points, N_setup_iterations, verbosity, &integral,
                    &error, &prob);
    // values in integral, error etc. should be discarded.

    // -----------------------------------------
    // VEGAS integration
    // -----------------------------------------
    vegas_reset_int(state); // only keep grid from previous integration
    vegas_disable_update_grid(state); // don't update grid during integration

    int seed2 = 2; // seed for VEGAS integration
    vegas_seed_random(state, seed2);
    int N_points = 500000;
    int N_iterations = 5;
    if (rank == 0) {
        printf("***************************************************************"
               "****\n");
        printf("*                                                              "
               "   *\n");
        printf("* VEGAS integration                                            "
               "   *\n");
        printf("*                                                              "
               "   *\n");
        printf("***************************************************************"
               "****\n");
    }
    vegas_integrate(state, seed2, integrand_ptr, (void *)&userdata, 1, N_points,
                    N_iterations, verbosity, &integral, &error, &prob);

    if (rank == 0) {
        printf("***************************************************************"
               "****\n");
        printf("*                                                              "
               "   *\n");
        printf("* result                                                       "
               "   *\n");
        printf("*                                                              "
               "   *\n");
        printf("***************************************************************"
               "****\n");
        printf("\n");
        printf("I = %g +- %g (chi^2/ndf = %g)\n", integral, error, prob);

        ASSERT_NEAR(integral, userdata.SqrtS * userdata.SqrtS /
                                  (1024.0 * M_PI * M_PI * M_PI),
                    error * 5.0)
            << "phase space volume not correct";
    }
}

TEST_P(PSTest, PS) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) {
        SUCCEED();
        return;
    }
    double S = 14000.0 * 14000.0;
    std::array<double, 7> v = GetParam();
    Phasespace::Phasespace ps = GenPhasespace(S, v);

    auto initial = ps.Momenta[0].Plus(ps.Momenta[1]);
    auto final = ps.Momenta[2].Plus(ps.Momenta[3].Plus(ps.Momenta[4]));
    ASSERT_NEAR(initial.PX(), 0.0, 1e-13) << "not in CMS frame";
    ASSERT_NEAR(initial.PY(), 0.0, 1e-13) << "not in CMS frame";
    ASSERT_NEAR(initial.PZ(), 0.0, 1e-13) << "not in CMS frame";

    ASSERT_NEAR(S / ps.S, 1.0, 1e-16) << "S not set correctly";

    double s = initial.Dot(initial);
    ASSERT_NEAR(ps.X1 * ps.X2 * ps.S / s, 1.0, 1e-14)
        << "partonic s not consistent";

    ASSERT_NEAR(final.PX(), initial.PX(), 1e-12)
        << "momentum in x direction not conserved";
    ASSERT_NEAR(final.PY(), initial.PY(), 1e-12)
        << "momentum in y direction not conserved";
    ASSERT_NEAR(final.PZ(), initial.PZ(), 1e-12)
        << "momentum in z direction not conserved";

    ASSERT_NEAR(final.E() / initial.E(), 1.0, 1e-12) << "energy not conserved";

    ASSERT_NEAR(ps.Momenta[2].Dot(ps.Momenta[2]), 0.0, 1e-8)
        << "particle 3 not massless";
    ASSERT_NEAR(ps.Momenta[3].Dot(ps.Momenta[3]), 0.0, 1e-8)
        << "particle 4 not massless";
    ASSERT_NEAR(ps.Momenta[4].Dot(ps.Momenta[4]), 0.0, 1e-8)
        << "particle 5 not massless";

    auto res = ps.Momenta[2].Plus(ps.Momenta[3]);
    double Q2 = res.Dot(res);
    ASSERT_NEAR(Q2 / (ps.X1 * ps.X2 * ps.S), v[2], 1e-15)
        << "particle 3 and 4 are not daughters of the resonance Q";
}

TEST_P(PSTest, LAB) {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank != 0) {
        SUCCEED();
        return;
    }
    double S = 14000.0 * 14000.0;
    std::array<double, 7> v = GetParam();
    Phasespace::Phasespace ps = GenPhasespace(S, v);

    Phasespace::Phasespace ps_lab;
    ps_lab.SetToLabFromCMS(&ps);

    auto initial = ps_lab.Momenta[0].Plus(ps_lab.Momenta[1]);
    auto final =
        ps_lab.Momenta[2].Plus(ps_lab.Momenta[3].Plus(ps_lab.Momenta[4]));

    ASSERT_NEAR(final.PX(), initial.PX(), 1e-12)
        << "momentum in x direction not conserved";
    ASSERT_NEAR(final.PY(), initial.PY(), 1e-12)
        << "momentum in y direction not conserved";
    ASSERT_NEAR(final.PZ(), initial.PZ(), 1e-12)
        << "momentum in z direction not conserved";

    ASSERT_NEAR(final.E() / initial.E(), 1.0, 1e-12) << "energy not conserved";

    ASSERT_NEAR(ps_lab.Momenta[2].Dot(ps_lab.Momenta[2]), 0.0, 1e-7)
        << "particle 3 not massless";
    ASSERT_NEAR(ps_lab.Momenta[3].Dot(ps_lab.Momenta[3]), 0.0, 1e-7)
        << "particle 4 not massless";
    ASSERT_NEAR(ps_lab.Momenta[4].Dot(ps_lab.Momenta[4]), 0.0, 1e-7)
        << "particle 5 not massless";

    auto P1 = ps_lab.Momenta[0];
    P1.Scale(1.0 / ps.X1);

    auto P2 = ps_lab.Momenta[1];
    P2.Scale(1.0 / ps.X2);

    ASSERT_NEAR(P1.E(), sqrt(ps.S) / 2.0, 1e-10)
        << "proton energy not consistent with S";
    ASSERT_NEAR(P1.PX(), 0.0, 1e-10) << "proton not in z direction?";
    ASSERT_NEAR(P1.PY(), 0.0, 1e-10) << "proton not in z direction?";
    ASSERT_NEAR(P1.PZ(), sqrt(ps.S) / 2.0, 1e-10)
        << "proton momentum not consistent with S";

    ASSERT_NEAR(P2.E(), sqrt(ps.S) / 2.0, 1e-10)
        << "proton energy not consistent with S";
    ASSERT_NEAR(P2.PX(), 0.0, 1e-10) << "proton not in z direction?";
    ASSERT_NEAR(P2.PY(), 0.0, 1e-10) << "proton not in z direction?";
    ASSERT_NEAR(P2.PZ(), -sqrt(ps.S) / 2.0, 1e-10)
        << "proton momentum not consistent with S";
}

const std::vector<std::array<double, 7>> values = {
    {0.561402, 0.029818, 0.640930, 0.264198, 0.543440, 0.841203, 0.340566},
    {0.154373, 0.398442, 0.206282, 0.378498, 0.072933, 0.840558, 0.132368},
    {0.266698, 0.972760, 0.663167, 0.018735, 0.368691, 0.161293, 0.559502},
    {0.373244, 0.615696, 0.142081, 0.335989, 0.658287, 0.013480, 0.521004},
    {0.392427, 0.051733, 0.481191, 0.080760, 0.336827, 0.904680, 0.101807},
    {0.009186, 0.251938, 0.380905, 0.441995, 0.851722, 0.753045, 0.120762},
    {0.748755, 0.914471, 0.090683, 0.818798, 0.646400, 0.051276, 0.397170},
    {0.035944, 0.430185, 0.943689, 0.607959, 0.991922, 0.895173, 0.600650},
    {0.210407, 0.083549, 0.275740, 0.337554, 0.289098, 0.060247, 0.777669},
    {0.812021, 0.683626, 0.748137, 0.932160, 0.294042, 0.345454, 0.250985},
    {0.633557, 0.874955, 0.365076, 0.128658, 0.341022, 0.275305, 0.066921},
    {0.787235, 0.364737, 0.631376, 0.633001, 0.648149, 0.464041, 0.881205},
    {0.008057, 0.565911, 0.255062, 0.582305, 0.913922, 0.812485, 0.088187},
    {0.878265, 0.586692, 0.328029, 0.849482, 0.091584, 0.491634, 0.997333},
    {0.183296, 0.758704, 0.986985, 0.221967, 0.807672, 0.038221, 0.598484},
    {0.341324, 0.056130, 0.461773, 0.090374, 0.778656, 0.179248, 0.569609},
    {0.680662, 0.828892, 0.353992, 0.486139, 0.926868, 0.073614, 0.133767},
    {0.995050, 0.639159, 0.586587, 0.098900, 0.599138, 0.145336, 0.948993},
    {0.654994, 0.792519, 0.121640, 0.272302, 0.828204, 0.913109, 0.311823},
    {0.958739, 0.566272, 0.361485, 0.501750, 0.135560, 0.825924, 0.241276},
};
INSTANTIATE_TEST_CASE_P(DifferentPSValues, PSTest, ::testing::ValuesIn(values));

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank != 0) { // print only messages from rank 0
        ::testing::TestEventListeners &listeners =
            ::testing::UnitTest::GetInstance()->listeners();
        delete listeners.Release(listeners.default_result_printer());
        // Adds a listener to the end.  Google Test takes the ownership.
        listeners.Append(new ::testing::EmptyTestEventListener);
    }
    auto ret = RUN_ALL_TESTS();
    MPI_Finalize();
    return ret;
}
