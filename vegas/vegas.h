#ifndef LO_VEGAS_H_
#define LO_VEGAS_H_

#ifdef __cplusplus
#include <cstdio>
extern "C" {
#else
#include <stdio.h>
#endif

#define VEGAS_MPI_TAG 45239

/**
 * rnd_gen is a pointer to some random number generator state.
 */
typedef void *rnd_gen;

/**
 * vegas_getrandom should returns uniformly distributed random number in [0,1).
 */
typedef double (*vegas_getrandom)(rnd_gen);

/**
 * vegas_initrandom should be a function which returns a random number
 * generator.
 *
 * Parameters:
 * args -- arguments which are used to init the random number generator.
 * seed -- seed of the random number generator
 */
typedef rnd_gen (*vegas_initrandom)(void *args, int seed);

/**
 * vegas_freerandom should be a function which frees the random number generator
 * initialized by vegas_initrandom
 */
typedef void (*vegas_freerandom)(rnd_gen);

/**
 * vegas_create_msg creates a message which is sent by the worker
 * over MPI.
 * vegas_integrate calls this function with the pointer to userdata.
 *
 * See vegas_register.
 */
typedef void (*vegas_create_msg)(void *userdata, int n, char *msg);

/**
 * vegas_process_msg is called in the MPI master and processes the
 * received message msg. The size of the message is given by n.
 *
 * See vegas_register.
 */
typedef void (*vegas_process_msg)(int num, int n, char *msg, int ncall,
                                  double integral, double sigma, void *args);

/**
 * VegasState represents the state of a vegas integration.
 */
struct VegasState;

/**
 * vegas_integrand is a function pointer to the integrand.
 *
 * Parameters:
 * n         -- number of dimensions of the integrand
 * x         -- n function arguments from 0 to 1
 * wgt       -- VEGAS weight
 * userdata  -- pointer to the user data
 * thread_id -- id of the thread which evaluates the integral.
 * out       -- result of the integrand evaluation
 *
 * Returns:
 * The integrand should return 0. If the function returns -1 the VEGAS
 * integrator uses the integrand value but increases the total number of points
 * by 1. If the function returns -2, the VEGAS integrator omits the integrand
 * value and draws a new point instead. If the function returns a value greater
 * than 0, the VEGAS integrator decreases the number of points by this value.
 */
typedef int (*vegas_integrand)(int n, const double *x, double wgt,
                               void *userdata, int thread_id, double *out);

/**
 * vegas_new initialises a new VegasState. Note that you have to call
 * vegas_set_random, vegas_init_grid before you can use the state.
 *
 * Parameters:
 * ndim -- dimension of the integral you want to integrate
 * nbin -- number of bins which VEGAS should use in every dimension.
 *
 * This function allocates new memory. Therefore, you should use
 * vegas_free to free the memory when it is noth needed anymore.
 */
struct VegasState *vegas_new(int ndim, int nbin);

/**
 * vegas_free frees the memory which was allocated with vegas_new.
 */
void vegas_free(struct VegasState *state);

/**
 * vegas_set_random sets the random number generator which will be used
 * for the Monte Carlo integration.
 *
 * Parameters:
 * state        -- pointer to the VegasState
 * initrandom   -- function which initializes the random number generator
 * init_rnd_arg -- argument of the initrandom
 * getrandom    -- function to obtain a random number
 * freerandom   -- function to free the random number generator
 */
void vegas_set_random(struct VegasState *state, vegas_initrandom initrandom,
                      void *init_rnd_arg, vegas_getrandom getrandom,
                      vegas_freerandom freerandom);

/**
 * vegas_seed_random seeds the random number generator.
 */
void vegas_seed_random(struct VegasState *state, int seed);

/**
 * vegas_init_grid initializes the VEGAS grid. This function is called by
 * vegas_new.
 */
void vegas_init_grid(struct VegasState *state);

/**
 * vegas_reset_int sets all so far obtained integration results to 0.
 * This can be useful if you use a first run to obtain the VEGAS grid
 * and don't want to use the integration results from the grid setup
 * stage for your final result.
 */
void vegas_reset_int(struct VegasState *state);

/**
 * vegas_get_integral writes the result of the last integration to `integral`,
 * `error` and `chi2ndf`.
 */
void vegas_get_integral(struct VegasState *state, double *integral,
                        double *error, double *chi2ndf);

/**
 * vegas_get_integralabs writes the integral of the absolute value of the
 * integrand to `integralabs`, `error` and `chi2ndf`.
 */
void vegas_get_integralabs(struct VegasState *state, double *integralabs,
                        double *error, double *chi2ndf);

/**
 * vegas_reset_grid resets the grid to an aequidistant grid in all dimensions.
 */
void vegas_reset_grid(struct VegasState *state);

enum When { AFTER_ITERATION };

/**
 * vegas_register registers callback functions which are used to compose/process
 * messages which vegas_integrate sends during integration.
 *
 * Parameters:
 * state     -- pointer to the VegasState
 * when      -- when to send the messages
 * f         -- creates the message from  userdata of vegas_integrate
 * g         -- processes in the MPI master the received message.
 * bufferlen -- buffer length for message
 */
int vegas_register(struct VegasState *state, enum When when, vegas_create_msg f,
                   vegas_process_msg g, int bufferlen, void *args);

/**
 * vegas_integrate integrates the integrand fxn and stores the integral value
 * in tgral, the error in sd and the chi^2/ndf in chi2a.
 * It uses Monte Carlo integration with the VEGAS algorithm.
 *
 * vegas_integrate uses MPI to parallelize the integration.
 * Every MPI worker uses its own random number generator which uses the
 * same random seed. The used random numbers do not depend
 * on the number of MPI workers. Therefore, the integration result is
 * independent of num_threads.
 *
 * Parameters:
 * state       -- VegasState which contains the grid etc.
 * seed        -- random number generator seed
 * fxn         -- function pointer to the integrand
 * userdata    -- pointer to data structure which is passed to fxn
 * num_threads -- number of pthreads (unused! only in pthread version relevant)
 * ncall       -- number of integrand calls per iteration
 * itmx        -- number of iterations
 * verbosity   -- controls the intermediate output to stdout.
 * tgral       -- output for the integral value
 * sd          -- output for the error of the integral value
 * chi2_ndf    -- chi^2/ndf
 *
 * Verbosity parameter:
 * <= 0 -- no output to stdout (only output to stderr if there is a pthread
 *         error)
 * >  0 -- print intermediate results for every iteration
 *
 * Return values:
 *  0 -- no error
 * -1 -- error in pthread_create
 * -2 -- error in pthread_join
 */
int vegas_integrate(struct VegasState *state, int seed, vegas_integrand fxn,
                    void *userdata, int num_threads, const int ncall,
                    const int itmx, const int verbosity, double *tgral,
                    double *sd, double *chi2_ndf);

/**
 * vegas_write_grid_1d writes the grid in dimension dim to stream.
 *
 * The output is a gnuplot script. The plot can be used to check if the grid
 * looks sensible.
 */
void vegas_write_grid_1d(struct VegasState *state, int dim, FILE *stream);

/**
 * vegas_write_grid_2d writes the grid in dimension dim1 and dimension dim2
 * to stream.
 *
 * The output is a gnuplot script. The plot can be used to check if the grid
 * looks sensible.
 */
void vegas_write_grid_2d(struct VegasState *state, int dim1, int dim2,
                         FILE *stream);

/**
 * vegas_read_grid_from_file reads the grid written with
 * vegas_write_grid_to_file from filename.
 */
int vegas_read_grid_from_file(struct VegasState *state, const char *filename);

/**
 * vegas_write_grid_to_file writes the grid of state to filename.
 */
int vegas_write_grid_to_file(const struct VegasState *state,
                             const char *filename);

int vegas_write_max_to_file(const struct VegasState * state, const char* filename);

/**
 * vegas_enable_update_grid enables the updating of the vegas grid.
 */
void vegas_enable_update_grid(struct VegasState *state);

/**
 * vegas_disable_update_grid disables the updating of the vegas grid.
 */
void vegas_disable_update_grid(struct VegasState *state);

/**
 * vegas_set_relative_accuracy sets the relative accuracy goal to `relerorr`.
 * If the integral error divided by the absolute value of the integral is
 * less than zero, the integration is stopped. The default value is 0.0.
 */
void vegas_set_relative_accuracy(struct VegasState *state, double relerror);

/**
 * vegas_set_absolute_accuracy sets the absolute accuracy goal to `abserorr`.
 * If the integral error is less than `abserror`, the integration is stopped.
 * The default value is 0.0.
 */
void vegas_set_absolute_accuracy(struct VegasState *state, double abserror);

#ifdef __cplusplus
}
#endif

#endif
