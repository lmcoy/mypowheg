#include <cmath>
#include <iostream>
#include <random>

#include <mpi.h>

#include "gtest/gtest.h"

#include "vegas.h"

// random seeds
int seed1 = 232;
int seed2 = 62731;

class MPIEnvironment : public ::testing::Environment {
  public:
    virtual ~MPIEnvironment() {}
    virtual void SetUp() { MPI_Init(NULL, NULL); }
    virtual void TearDown() { MPI_Finalize(); }
};

::testing::Environment *const mpi_env =
    ::testing::AddGlobalTestEnvironment(new MPIEnvironment);
/*********************************************************************
 *                     random number generation                      *
 *********************************************************************/
struct RNG {
    RNG(int seed) {
        r = new std::mt19937(seed);
        dis = new std::uniform_real_distribution<>(0.0, 1.0);
    }
    ~RNG() {
        delete r;
        delete dis;
    }
    std::mt19937 *r;
    std::uniform_real_distribution<> *dis;
};

void *init_rnd(void *args, int seed) { return new RNG(seed); }

void free_rnd(void *r) {
    RNG *rng = (RNG *)r;
    delete rng;
}

double get_rnd(void *r) {
    RNG *rng = (RNG *)r;
    return (*rng->dis)(*rng->r);
}

/*********************************************************************
 *                          test functions                           *
 *********************************************************************/

// integral ~ 0.73389261543701282506 (mathematica)
int f1(int n, const double *x, const double wgt, void *userdata, int threadid,
       double *out) {
    double cx =
        0.5 * x[0] + 0.3 * x[1] + 0.6 * x[2] + 0.9 * x[3] + 0.125 * x[4];
    double w1 = 0.7;
    *out = cos(cx + 2.0 * M_PI * w1);
    return 0;
}

// integral ~ 1.0543373707555682624e-6 (mathematica)
int f2(int n, const double *x, const double wgt, void *userdata, int threadid,
       double *out) {
    double w[8] = {0.3, 0.4, 0.9, 0.1, 0.4, 0.5, 0.2, 0.8};
    double c[8] = {0.1, 0.9, 0.4, 0.5, 0.7, 0.4, 0.8, 0.3};

    double ret = 1.0;
    for (int i = 0; i < 8; i++) {
        ret /= ((x[i] - w[i]) * (x[i] - w[i]) + 1. / c[i] / c[i]);
    }
    *out = ret;
    return 0;
}

// integral = 15/34
int f3(int n, const double *x, const double wgt, void *userdata, int threadid,
       double *out) {
    *out = 1.0 / pow(1 + 0.2 * x[0] + 0.5 * x[1], 3.0);
    return 0;
}

// integral = 144/25 * pi * erf(5/24)^2
int f4(int n, const double *x, const double wgt, void *userdata, int threadid,
       double *out) {
    *out = exp(-25.0 / 144.0 * (pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)));
    return 0;
}

// integral ~ 0.675283390170236356
int f5(int n, const double *x, const double wgt, void *userdata, int threadid,
       double *out) {
    *out = exp(-sqrt(pow(x[0] - 0.5, 2.0) + pow(x[1] - 1.0 / 3.0, 2.0)));
    return 0;
}

// integral = 144/25 * pi * erf(5/24)^2
int f6(int n, const double *x, const double wgt, void *userdata, int threadid,
       double *out) {
    static int counter = 0;
    counter += 1;
    if (counter % 4 == 0) {
        return -1;
    }
    *out = exp(-25.0 / 144.0 * (pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)));
    return 0;
}

// integral = 1
int f7(int n, const double *x, const double wgt, void *userdata, int threadid,
       double *out) {
    *out = 1.0;
    return 0;
}

/*********************************************************************
 *                               tests                               *
 *********************************************************************/

TEST(VEGAS, F1) {
    double tgral, sd, chi2a;
    int ncall = 100000; // number of function calls
    int itmx = 10;      // number of iterations
    int verbosity = 0;

    RecordProperty("dimensions", 5);

    struct VegasState *state = vegas_new(5, 50);
    vegas_set_random(state, init_rnd, NULL, get_rnd, free_rnd);
    vegas_seed_random(state, seed1);

    vegas_integrate(state, seed1, f1, NULL, 6, ncall, itmx, verbosity, &tgral,
                    &sd, &chi2a);

    // forget integral value
    vegas_reset_int(state);
    // do the integration with the previously obtained grid
    vegas_integrate(state, seed2, f1, NULL, 6, ncall, itmx, verbosity, &tgral,
                    &sd, &chi2a);
    vegas_free(state);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        EXPECT_NEAR(0.73389261543701282506, tgral, 5.0 * sd);
    }
}

TEST(VEGAS, F2) {
    double tgral, sd, chi2a;
    int ncall = 50000; // number of function calls
    int itmx = 10;     // number of iterations
    int verbosity = 0;

    RecordProperty("dimensions", 8);

    struct VegasState *state = vegas_new(8, 50);
    vegas_set_random(state, init_rnd, NULL, get_rnd, free_rnd);
    vegas_seed_random(state, seed1);

    // setup grid
    vegas_integrate(state, seed1, f2, NULL, 6, 10000, 5, verbosity, &tgral, &sd,
                    &chi2a);

    // forget integral value
    vegas_reset_int(state);
    // do the integration with the previously obtained grid
    vegas_integrate(state, seed2, f2, NULL, 6, ncall, itmx, verbosity, &tgral,
                    &sd, &chi2a);
    vegas_free(state);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        EXPECT_NEAR(1.0543373707555682624e-6, tgral, 5.0 * sd);
    }
}

TEST(VEGAS, F3) {
    double tgral, sd, chi2a;
    int ncall = 50000; // number of function calls
    int itmx = 10;     // number of iterations
    int verbosity = 0;

    RecordProperty("dimensions", 2);

    struct VegasState *state = vegas_new(2, 50);
    vegas_set_random(state, init_rnd, NULL, get_rnd, free_rnd);
    vegas_seed_random(state, seed1);

    // setup grid
    vegas_integrate(state, seed1, f3, NULL, 6, 10000, 5, verbosity, &tgral, &sd,
                    &chi2a);

    // forget integral value
    vegas_reset_int(state);
    // do the integration with the previously obtained grid
    vegas_integrate(state, seed2, f3, NULL, 6, ncall, itmx, verbosity, &tgral,
                    &sd, &chi2a);
    vegas_free(state);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        EXPECT_NEAR(15.0 / 34.0, tgral, 5.0 * sd);
    }
}

TEST(VEGAS, F4) {
    double tgral, sd, chi2a;
    int ncall = 50000; // number of function calls
    int itmx = 10;     // number of iterations
    int verbosity = 0;

    RecordProperty("dimensions", 2);

    struct VegasState *state = vegas_new(2, 50);
    vegas_set_random(state, init_rnd, NULL, get_rnd, free_rnd);
    vegas_seed_random(state, seed1);

    // setup grid
    vegas_integrate(state, seed1, f4, NULL, 6, 10000, 5, verbosity, &tgral, &sd,
                    &chi2a);

    // forget integral value
    vegas_reset_int(state);
    // do the integration with the previously obtained grid
    vegas_integrate(state, seed2, f4, NULL, 6, ncall, itmx, verbosity, &tgral,
                    &sd, &chi2a);
    vegas_free(state);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        EXPECT_NEAR(144.0 / 25.0 * M_PI * pow(erf(5.0 / 24.0), 2.0), tgral,
                    5.0 * sd);
    }
}

struct UserData {
    int *data;
};

void create_msg(void *userdata, int n, char *msg) {
    UserData *ud = (UserData *)userdata;
    printf("%d\n", *ud->data);
    *(ud->data) = 5;
    snprintf(msg, n, "%d", 1);
}

void process_msg(int num, int n, char *msg, int ncall, double integral,
                 double sigma, void *args) {
    int *i = (int *)args;
    *i += atoi(msg);
}

TEST(VEGAS, F5) {
    double tgral, sd, chi2a;
    int ncall = 100000; // number of function calls
    int itmx = 10;      // number of iterations
    int verbosity = 0;

    RecordProperty("dimensions", 2);

    UserData ud;
    ud.data = new int;
    *(ud.data) = 4;

    struct VegasState *state = vegas_new(2, 50);
    vegas_set_random(state, init_rnd, NULL, get_rnd, free_rnd);
    vegas_seed_random(state, seed1);

    // setup grid
    vegas_integrate(state, seed1, f5, (void *)&ud, 6, 10000, 5, verbosity,
                    &tgral, &sd, &chi2a);

    // forget integral value
    vegas_reset_int(state);
    int number = 0;
    vegas_register(state, AFTER_ITERATION, create_msg, process_msg, 100,
                   &number);
    // do the integration with the previously obtained grid
    vegas_integrate(state, seed2, f5, (void *)&ud, 6, ncall, itmx, verbosity,
                    &tgral, &sd, &chi2a);
    vegas_free(state);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    if (rank == 0) {
        EXPECT_EQ(itmx * num_procs, number);
        EXPECT_NEAR(0.675283390170236356, tgral, 5.0 * sd);
    }
}

TEST(VEGAS, F6) {
    double tgral, sd, chi2a;
    int ncall = 50000; // number of function calls
    int itmx = 10;     // number of iterations
    int verbosity = 0;

    RecordProperty("dimensions", 2);

    struct VegasState *state = vegas_new(2, 50);
    vegas_set_random(state, init_rnd, NULL, get_rnd, free_rnd);
    vegas_seed_random(state, seed1);

    // setup grid
    vegas_integrate(state, seed1, f6, NULL, 6, 10000, 5, verbosity, &tgral, &sd,
                    &chi2a);

    // forget integral value
    vegas_reset_int(state);
    // do the integration with the previously obtained grid
    vegas_integrate(state, seed2, f6, NULL, 6, ncall, itmx, verbosity, &tgral,
                    &sd, &chi2a);
    vegas_free(state);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        EXPECT_NEAR(144.0 / 25.0 * M_PI * pow(erf(5.0 / 24.0), 2.0), tgral,
                    5.0 * sd);
    }
}

TEST(VEGAS, F7) {
    double tgral, sd, chi2a;
    int ncall = 100000; // number of function calls
    int itmx = 10;      // number of iterations
    int verbosity = 0;

    RecordProperty("dimensions", 10);

    struct VegasState *state = vegas_new(10, 50);
    vegas_set_random(state, init_rnd, NULL, get_rnd, free_rnd);
    vegas_seed_random(state, seed1);

    // setup grid
    vegas_integrate(state, seed1, f7, NULL, 6, 10000, 5, verbosity, &tgral, &sd,
                    &chi2a);

    // forget integral value
    vegas_reset_int(state);
    // do the integration with the previously obtained grid
    vegas_integrate(state, seed2, f7, NULL, 6, ncall, itmx, verbosity, &tgral,
                    &sd, &chi2a);
    vegas_free(state);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        EXPECT_NEAR(1.0, tgral, 5.0 * sd);
    }
}
