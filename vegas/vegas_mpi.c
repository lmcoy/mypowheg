#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>

#include <mpi.h>

#include "vegas.h"

struct Matrix {
    double *m;
    int ncol;
    int nrow;
};

struct Matrix *matrix_new(int rows, int cols) {
    struct Matrix *m = malloc(sizeof(struct Matrix));
    m->ncol = cols;
    m->nrow = rows;
    m->m = malloc(sizeof(double) * rows * cols);
    int i;
    for (i = 0; i < rows * cols; i++) {
        m->m[i] = 0.0;
    }
    return m;
}

static void matrix_free(struct Matrix *m) {
    free(m->m);
    m->m = NULL;
    free(m);
    m = NULL;
}

static void matrix_set(struct Matrix *m, int row, int col, double v) {
    m->m[col + row * m->ncol] = v;
}

static double matrix_get(struct Matrix *m, int row, int col) {
    return m->m[col + row * m->ncol];
}

static void matrix_add_to_elem(struct Matrix *m, int row, int col, double v) {
    m->m[col + row * m->ncol] += v;
}

static void matrix_set_zero(struct Matrix *m) {
    int i;
    for (i = 0; i < m->nrow * m->ncol; i++) {
        m->m[i] = 0.0;
    }
}

static void matrix_add(struct Matrix *m, struct Matrix *p) {
    int i;
    for (i = 0; i < m->nrow * m->ncol; i++) {
        m->m[i] += p->m[i];
    }
}

static void matrix_mul_factor(struct Matrix *m, double v) {
    int i;
    for (i = 0; i < m->nrow * m->ncol; i++) {
        m->m[i] *= v;
    }
}

static void rebin(const double rc, const int nd, const double *r, double *xin,
                  struct Matrix *xi, const int j) {
    int k = 0;
    double dr = 0.0, xn = 0.0, xo = 0.0;

    int i;
    for (i = 0; i < nd - 1; i++) {
        while (rc > dr) {
            dr += r[k];
            k += 1;
        }
        if (k > 1) {
            xo = matrix_get(xi, j, k - 2);
        }
        xn = matrix_get(xi, j, k - 1);
        dr -= rc;
        xin[i] = xn - (xn - xo) * dr / r[k - 1];
    }
    for (i = 0; i < nd - 1; i++) {
        matrix_set(xi, j, i, xin[i]);
    }
    matrix_set(xi, j, nd - 1, 1.0);
}

struct VegasRandom {
    vegas_initrandom init;
    vegas_freerandom free;
    vegas_getrandom get;
    void * init_arg;
};

struct VegasState {
  struct Matrix *xi;
  int nbin;
  int ndim;
  double si;
  double schi;
  double swgt;
  int has_grid;
  int refine_grid; ///< If zero the grid is not updated
  int ncalls_last_it; ///< Number of points which were used in the last
                      ///iteration. Equals ncall if no event was rejected.

  struct VegasRandom random;
  rnd_gen *generators;
  vegas_create_msg msg_after_it;
  vegas_process_msg proc_after_it;
  int after_it_bufferlen;
  void * proc_args;
};

void vegas_enable_update_grid(struct VegasState * state) {
    state->refine_grid = 1;
}

void vegas_disable_update_grid(struct VegasState * state) {
    state->refine_grid = 0;
}

void vegas_reset_grid(struct VegasState *state) {
    int ndim = state->ndim;
    int j;
    for (j = 0; j < ndim; j++) {
        matrix_set(state->xi, j, 0, 1.0);
    }
    state->has_grid = 0;
}

void vegas_reset_int(struct VegasState *state) {
    state->si = 0.0;
    state->schi = 0.0;
    state->swgt = 0.0;
}

void vegas_init_grid(struct VegasState *state) {
    int ndim = state->ndim;
    int nbin = state->nbin;
    double *r = malloc(sizeof(double) * nbin);
    int i, j;
    for (i = 0; i < nbin; i++) {
        r[i] = 1.0;
    }
    double *xin = malloc(sizeof(double) * nbin);
    for (j = 0; j < ndim; j++) {
        rebin(1.0 / (double)nbin, nbin, r, xin, state->xi, j);
    }
    free(xin);
    xin = NULL;
    free(r);
    r = NULL;
    state->has_grid = 1;
}

void vegas_set_random(struct VegasState *state, vegas_initrandom initrandom,
                      void *init_rnd_arg, vegas_getrandom getrandom,
                      vegas_freerandom freerandom) {
    state->random.init_arg = init_rnd_arg;
    state->random.init = initrandom;
    state->random.get = getrandom;
    state->random.free = freerandom;
}

static void vegas_free_random(struct VegasState * state) {
    if (state->generators == NULL) {
        return;
    }
    int dim;
    int ndim = state->ndim;
    struct VegasRandom * rnd = &state->random;
    for (dim = 0; dim < ndim; dim++) {
        rnd->free(state->generators[dim]);
    }
    free(state->generators);
    state->generators = NULL;
}

void vegas_seed_random(struct VegasState *state, int seed) {
    if (state->generators != NULL) {
        vegas_free_random(state);
    }
    int ndim = state->ndim;
    struct VegasRandom * rnd = &state->random;

    state->generators = malloc(sizeof(rnd_gen) * ndim);
    int dim;
    for (dim = 0; dim < ndim; dim++) {
        state->generators[dim] = rnd->init(rnd->init_arg, seed + dim * 127);
    }
}


struct VegasState *vegas_new(int ndim, int nbin) {
    struct VegasState *state = malloc(sizeof(struct VegasState));
    state->xi = matrix_new(ndim, nbin);
    state->nbin = nbin;
    state->ndim = ndim;
    state->si = 0.0;
    state->schi = 0.0;
    state->swgt = 0.0;
    vegas_reset_grid(state);
    vegas_init_grid(state);
    state->random.get = NULL;
    state->random.init = NULL;
    state->random.free = NULL;
    state->random.init_arg = NULL;
    state->msg_after_it = NULL;
    state->proc_after_it = NULL;
    state->proc_args =NULL;
    state->after_it_bufferlen = 0;
    state->refine_grid = 1;
    state->ncalls_last_it = 0;
    state->generators = NULL;
    return state;
}

int vegas_read_grid_from_file(struct VegasState *state, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        return -1;
    }
    int ndim, nbin;
    int ret = fscanf(file, "%d %d", &ndim, &nbin);
    if (ret != 2) {
        goto fail_pre;
    }
    struct Matrix *m = matrix_new(ndim, nbin);
    int i;
    for (i = 0; i < ndim * nbin; i++) {
        double value;
        int ret = fscanf(file, "%lf ", &value);
        if (ret != 1) {
            goto fail;
        }
        m->m[i] = value;
    }
    fclose(file);
    matrix_free(state->xi);
    state->xi = m;
    state->ndim = ndim;
    state->nbin = nbin;

    return 0;

fail:
    matrix_free(m);
fail_pre:
    fclose(file);
    return -1;
}

int vegas_write_grid_to_file(const struct VegasState *state,
                             const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        return -1;
    }
    int ret = fprintf(file, "%d %d ", state->ndim, state->nbin);
    if (ret < 4) {
        goto fail;
    }
    int i;
    for (i = 0; i < state->ndim * state->nbin; i++) {
        int ret = fprintf(file, "%.12g ", state->xi->m[i]);
        if (ret < 2) {
            goto fail;
        }
    }
    fclose(file);
    return 0;
fail:
    fclose(file);
    return -1;
}

void vegas_free(struct VegasState *state) {
    matrix_free(state->xi);
    vegas_free_random(state);
    state->generators = NULL;
    free(state);
}

int vegas_register(struct VegasState *state, enum When when, vegas_create_msg f,
                   vegas_process_msg g, int bufferlen, void *pargs) {
    switch (when) {
    case AFTER_ITERATION:
        state->msg_after_it = f;
        state->proc_after_it = g;
        state->proc_args = pargs;
        state->after_it_bufferlen = bufferlen;
        return 0;
    }
    return -1;
}


void vegas_write_grid_1d(struct VegasState *state, int dim, FILE *stream) {
    fprintf(stream,
            "set yrange [0:1]\nset xrange [0:1]\nplot \"-\" with impulses\n");
    int i;
    for (i = 0; i < state->nbin; i++) {
        fprintf(stream, "%20g    1\n", matrix_get(state->xi, dim, i));
    }
}

void vegas_write_grid_2d(struct VegasState *state, int dim1, int dim2,
                         FILE *stream) {
    fprintf(stream, "set yrange[0:1]\nset xrange[0:1]\n");
    fprintf(stream, "set title \"vegas grid\"\n");
    fprintf(stream, "set xlabel \"axis %d\"\n", dim1);
    fprintf(stream, "set ylabel \"axis %d\"\n", dim2);
    fprintf(stream,
            "plot \"-\" w impulse lc 1 ti \"\", \"-\" w lines lc 1 ti\"\"\n");
    int i;
    for (i = 0; i < state->nbin; i++) {
        fprintf(stream, "%20g    1\n", matrix_get(state->xi, dim1, i));
    }
    fprintf(stream, "EOF\n");
    for (i = 0; i < state->nbin; i++) {
        fprintf(stream, "0 %20g\n", matrix_get(state->xi, dim2, i));
        fprintf(stream, "1 %20g\n\n", matrix_get(state->xi, dim2, i));
    }
    fprintf(stream, "EOF\n");
}

static int min(int a, int b) {
    if (a < b) {
        return a;
    }
    return b;
}

struct int_state {
    int *ia;
    double *x;
    struct Matrix * xi;
    int ndim;
    rnd_gen * generators;
    struct Matrix * d;
    struct Matrix * di;
    double fb;
    double f2b;
    int ncalls;
};

static void int_state_init(struct int_state *istate, int ndim, int nbin,
                           rnd_gen* rnd_gen,
                           int seed) {
    istate->di = matrix_new(nbin, ndim);
    istate->d = matrix_new(nbin, ndim);
    istate->ia = malloc(sizeof(int) * ndim);
    istate->x = malloc(sizeof(double) * ndim);

    istate->xi = matrix_new(ndim, nbin);

    istate->generators = rnd_gen;
    istate->ndim = ndim;
    istate->ncalls = 0;
}

static void int_state_free(struct int_state *istate) {
    free(istate->x);
    istate->x = NULL;
    free(istate->ia);
    istate->ia = NULL;
    matrix_free(istate->d);
    matrix_free(istate->di);
    matrix_free(istate->xi);
}

static void integration(struct int_state *istate, int ncall, int additional,
                        int ndim, int nbin, int rank, int num_procs,
                        struct VegasRandom *rnd, vegas_integrand fxn,
                        void *userdata) {
    matrix_set_zero(istate->di);
    matrix_set_zero(istate->d);
    istate->fb = 0.0;
    istate->f2b = 0.0;
    istate->ncalls = 0;
    int k, j, i;
    int n_max = ncall + additional;
    int rej = 0;
    for (k = 0; k < n_max; k++) {
        double wgt = 1.0;
        int reject = 0;
        double ff;
        for (j = 0; j < ndim; j++) {
            double u = 0.0;
            if (k < ncall || rej > 0) {
                for (i = 0; i < num_procs; i++) {
                    double tmp = rnd->get(istate->generators[j]);
                    if (i == rank) {
                        u = tmp;
                    }
                }
            } else {
                // only worker 0 is left but we have to draw the same random
                // numbers in every worker to keep the generators in sync.
                u = rnd->get(istate->generators[j]);
            }
            double xn = (1 - u) * nbin + 1.0;
            istate->ia[j] = min((int)(xn), nbin);
            double x_1 = matrix_get(istate->xi, j, istate->ia[j] - 1);
            double x_2 = 0.0;
            if (istate->ia[j] > 1) {
                x_2 = matrix_get(istate->xi, j, istate->ia[j] - 2);
            }
            double delta_x = x_1 - x_2;
            istate->x[j] = x_2 + (xn - istate->ia[j]) * delta_x;
            wgt *= delta_x * (double)nbin;
        }
        if (rej == 0 && k >= ncall && rank != 0) {
            // only worker 0 is responsible for the random numbers which
            // we couldn't distribute equally to all workers.
            continue;
        }
        reject = fxn(ndim, istate->x, wgt, userdata, rank, &ff);
        if (reject == -1) {
            // use this point but increase number of points
            rej += 1;
            if (rank != 0 && rej == 1) {
                n_max -= additional;
            }
            n_max += 1;
        } else if(reject == -2) {
            // redo point and don't use integrand value
            k -= 1;
            continue;
        } else if (reject > 0) {
            n_max -= reject;
        }

        double f = wgt * ff;
        double f2 = f * f;
        istate->fb += f;
        istate->f2b += f2;
        for (j = 0; j < ndim; j++) {
            int iaj = istate->ia[j];
            matrix_add_to_elem(istate->di, iaj - 1, j, f);
            matrix_add_to_elem(istate->d, iaj - 1, j, f2);
        }
    }
    if (rej == 0 && rank != 0) {
        istate->ncalls = ncall;
    } else {
        istate->ncalls = n_max;
    }
}

static void worker(struct VegasState *state, int ndim, int nbin, int rank,
                   int num_procs, int ncall, int additional, int itmx,
                   vegas_integrand fxn, void *userdata, struct VegasRandom *rnd,
                   int seed, vegas_create_msg cmsg, int bufferlen) {
    struct int_state istate;
    int_state_init(&istate, ndim, nbin, state->generators, seed);

    char *buffer = malloc(sizeof(char) * bufferlen);

    int it;
    for (it = 0; it < itmx; it++) {
        MPI_Bcast(istate.xi->m, ndim * nbin, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        integration(&istate, ncall, additional, ndim, nbin, rank, num_procs, 
                    rnd, fxn, userdata);

        MPI_Send(&istate.fb, 1, MPI_DOUBLE, 0, VEGAS_MPI_TAG, MPI_COMM_WORLD);
        MPI_Send(&istate.f2b, 1, MPI_DOUBLE, 0, VEGAS_MPI_TAG, MPI_COMM_WORLD);
        int di_size = istate.di->nrow * istate.di->ncol;
        int d_size = istate.d->nrow * istate.d->ncol;
        MPI_Send(istate.di->m, di_size, MPI_DOUBLE, 0, VEGAS_MPI_TAG,
                 MPI_COMM_WORLD);
        MPI_Send(istate.d->m, d_size, MPI_DOUBLE, 0, VEGAS_MPI_TAG,
                 MPI_COMM_WORLD);
        MPI_Send(&istate.ncalls, 1, MPI_INT, 0, VEGAS_MPI_TAG, MPI_COMM_WORLD);

        if(cmsg != NULL) {
            cmsg(userdata, bufferlen, buffer);
            MPI_Send(buffer, bufferlen, MPI_CHAR, 0, VEGAS_MPI_TAG,
                     MPI_COMM_WORLD);
        }
    }
    int_state_free(&istate);
    free(buffer);
    buffer = NULL;
}

static void master(struct VegasState *state, int rank, int num_procs, int ncall,
                   int worker_ncall, int additional, int itmx,
                   vegas_integrand fxn, void *userdata, struct VegasRandom *rnd,
                   int seed, int verbosity, double *tgral, double *sd,
                   double *chi2_ndf) {
    static const double ALPH = 1.5, TINY = 1.0e-30;
    int ndim = state->ndim;
    int nbin = state->nbin;

    double *r = malloc(sizeof(double) * nbin);
    double *dt = malloc(sizeof(double) * ndim);
    double *xin = malloc(sizeof(double) * nbin);
    struct Matrix *di_add = matrix_new(nbin, ndim);
    struct Matrix *d_add = matrix_new(nbin, ndim);

    struct int_state istate;
    int_state_init(&istate, ndim, nbin, state->generators, seed);

    char *buffer = malloc(sizeof(char) * state->after_it_bufferlen);
    int it;
    for (it = 0; it < itmx; it++) {
        int ncalls_total = 0;
        clock_t t_start = clock();
        int i;
        for (i = 0; i < ndim * nbin; i++) {
            istate.xi->m[i] = state->xi->m[i];
        }
        MPI_Bcast(state->xi->m, ndim * nbin, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        istate.fb = 0.0;
        istate.f2b = 0.0;
        struct Matrix *d = istate.d;
        struct Matrix *di = istate.di;
        matrix_set_zero(istate.d);
        matrix_set_zero(istate.di);
        integration(&istate, worker_ncall, additional, ndim, nbin, rank,
                    num_procs, rnd, fxn, userdata);
        double fb = istate.fb;
        double f2b = istate.f2b;
        ncalls_total += istate.ncalls;
        state->ncalls_last_it += istate.ncalls;
        for (i = 1; i < num_procs; i++) {
            MPI_Status status;
            double fb_add;
            double f2b_add;
            int ncalls_worker;
            MPI_Recv(&fb_add, 1, MPI_DOUBLE, i, VEGAS_MPI_TAG, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(&f2b_add, 1, MPI_DOUBLE, i, VEGAS_MPI_TAG, MPI_COMM_WORLD,
                     &status);
            MPI_Recv(di_add->m, nbin * ndim, MPI_DOUBLE, i, VEGAS_MPI_TAG,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(d_add->m, nbin * ndim, MPI_DOUBLE, i, VEGAS_MPI_TAG,
                     MPI_COMM_WORLD, &status);
            MPI_Recv(&ncalls_worker, 1, MPI_INT, i, VEGAS_MPI_TAG, MPI_COMM_WORLD,
                     &status);

            fb += fb_add;
            f2b += f2b_add;
            matrix_add(di, di_add);
            matrix_add(d, d_add);
            state->ncalls_last_it += ncalls_worker;
            ncalls_total += ncalls_worker;
        }
        double jac = 1.0/((double)ncalls_total);
        double dv2g = 1.0 / (ncalls_total - 1.0);
        fb *= jac;
        f2b *= jac * jac;
        matrix_mul_factor(di, jac);
        matrix_mul_factor(d, jac * jac);
        f2b = sqrt(f2b * ncalls_total);
        f2b = (f2b - fb) * (f2b + fb);
        if (f2b <= 0.0) {
            f2b = TINY;
        }
        double integral_it = fb;
        double sigma_it = f2b * dv2g;
        double wgt = 1.0 / sigma_it;
        state->si += wgt * integral_it;
        state->schi += wgt * integral_it * integral_it;
        state->swgt += wgt;
        *tgral = state->si / state->swgt;
        *chi2_ndf = (state->schi - state->si * *tgral) / (it + 0.0001);
        if (*chi2_ndf < 0.0) {
            *chi2_ndf = 0.0;
        }
        *sd = sqrt(1.0 / state->swgt);
        sigma_it = sqrt(sigma_it);
        if (state->msg_after_it != NULL) {
            state->msg_after_it(userdata, state->after_it_bufferlen, buffer);
            state->proc_after_it(num_procs, state->after_it_bufferlen, buffer,
                                 ncalls_total, integral_it, sigma_it,
                                 state->proc_args);
        }
        for (i = 1; i < num_procs; i++) {
            MPI_Status status;
            if (state->msg_after_it != NULL) {
                MPI_Recv(buffer, state->after_it_bufferlen, MPI_CHAR, i,
                         VEGAS_MPI_TAG, MPI_COMM_WORLD, &status);
                state->proc_after_it(num_procs, state->after_it_bufferlen,
                                     buffer, ncalls_total, integral_it,
                                     sigma_it, state->proc_args);
            }
        }
        clock_t t_end = clock();
        if (verbosity > 0) {
            printf("integral: %g +- %g (chi^2/ndf = %g)\n", *tgral, *sd,
                   *chi2_ndf);
            if(verbosity > 1) {
                printf("    last it = %d: %g +- %g\n", it + 1, integral_it, sigma_it);
                double p_rej =
                    ((double)ncall) / ((double)(ncalls_total)) * 100.0;
                printf("    requested points: %d, sampled points: %d used %2.f "
                       "%%\n",
                       ncall, ncalls_total, p_rej);
            }
            double elapsed_secs = (double)(t_end - t_start) / CLOCKS_PER_SEC;
            time_t now = time(0);
            int it_left = itmx - it - 1;
            now += it_left * elapsed_secs;
            char time_buffer[30];
            struct tm now_tm = *localtime(&now);
            strftime(time_buffer, 30, "%c", &now_tm);
            printf("    iteration %d, progress = %.1f%% (%s)\n", it + 1,
                   100.0 * (it + 1) / itmx, time_buffer);
        }
        /* refine grid */
        if (state->refine_grid) {
            int j;
            for (j = 0; j < ndim; j++) {
                double xo = matrix_get(d, 0, j);
                double xn = matrix_get(d, 1, j);
                matrix_set(d, 0, j, (xo + xn) / 2.0);
                dt[j] = matrix_get(d, 0, j);
                for (i = 2; i < nbin; i++) {
                    double rc = xo + xn;
                    xo = xn;
                    xn = matrix_get(d, i, j);
                    matrix_set(d, i - 1, j, (rc + xn) / 3.0);
                    dt[j] += matrix_get(d, i - 1, j);
                }
                matrix_set(d, nbin - 1, j, (xo + xn) / 2.0);
                dt[j] += matrix_get(d, nbin - 1, j);
            }
            for (j = 0; j < ndim; j++) {
                double rc = 0.0;
                for (i = 0; i < nbin; i++) {
                    double dij = matrix_get(d, i, j);
                    if (dij < TINY) {
                        matrix_set(d, i, j, TINY);
                        dij = TINY;
                    }
                    r[i] = pow((1.0 - dij / dt[j]) / (log(dt[j]) - log(dij)),
                               ALPH);
                    rc += r[i];
                }
                rebin(rc / (double)nbin, nbin, r, xin, state->xi, j);
            }
        }
    }
    int_state_free(&istate);
    free(xin);
    xin = NULL;
    free(r);
    r = NULL;
    free(dt);
    dt = NULL;
    matrix_free(d_add);
    matrix_free(di_add);
    free(buffer);
    buffer = NULL;
    free(buffer);
    buffer = NULL;
}

int vegas_integrate(struct VegasState *state, int seed, vegas_integrand fxn,
                    void *userdata, int num_threads, int ncall, int itmx,
                    int verbosity, double *tgral, double *sd,
                    double *chi2_ndf) {

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    assert(state->generators != NULL &&
           "You must use vegas_seed_random before calling vegas_integrate");

    int worker_ncall = ncall / num_procs;
    int additional = ncall % num_procs;
    state->ncalls_last_it = 0;
    if (rank == 0) {
        master(state, rank, num_procs, ncall, worker_ncall, additional, itmx,
               fxn, userdata, &state->random, seed, verbosity, tgral, sd,
               chi2_ndf);
    } else {
        worker(state, state->ndim, state->nbin, rank, num_procs, worker_ncall,
               additional, itmx, fxn, userdata, &state->random, seed,
               state->msg_after_it, state->after_it_bufferlen);
    }
    return 0;
}

