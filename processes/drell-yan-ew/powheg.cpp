#include <cmath>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <cstring>
#include <ctime>
#include <cassert>
#include <cfloat>
#include <cstdio>
#include <cfenv>
#include <memory>

#include "pdf/lhapdf.h"  
#include "physics/pdgcode.h"
#include "vegas.h"
#include "util/databuffer.h"
#include "util/stringutil.h"
#include "random/rnd.h"
#include "math/math.h"
#include "powheg/generateevents.h"
#include "powheg/btilde.h"
#include "powheg/findnorm.h"
#include "lhe/lhefile.h"
#include "util/stringutil.h"

#include "zdata.h"
#include "config.h"
#include "processlist.h"
#include "process/data.h"
#include "genbornps.h"
#include "run.h"
#include "cuts.h"

double clip(double x) {
    if (x >= 1.0 && x < 1.0 + 1e-9) {
        return 1.0 - 1e-12;
    }
    if (x <= 0.0 && x > -1e-9) {
        return 1e-12;
    }
    return x;
}

typedef int (*function)(const Phasespace::Phasespace &, double, double, double,
                        double, double *, UserProcess::Data *);

template <function func>
static int integrand3(int n, const double *x, const double wgt, void *params,
                      int threadid, double *out) {
    constexpr double M = 91.188;
    constexpr double G = 2.45;
    UserProcess::Data *userdata = (UserProcess::Data *)params;
    // implement mll cut in Breit-Wigner trafo.
    // can be set to 0.0 for small cut values (mll < ~100).
    auto cuts = std::static_pointer_cast<DrellYanCuts>(userdata->cuts);
    double Mll = cuts->mllmin;
    double S = userdata->SqrtS * userdata->SqrtS;
    double Delta =
        atan((S - M * M) / (M * G)) - atan((Mll * Mll - M * M) / (M * G));
    double s = x[0];
    double z = x[1];
    double ts = Delta * s + atan((Mll * Mll - M * M) / (M * G));
    double tau = clip((G * tan(ts) + M) * M / S);
    LIB_ASSERT(tau >= 0.0 && tau <= 1.0, "tau = %.16g\n", tau);
    double xx[n];
    double x1 = (1.0 - tau) * z + tau;
    xx[0] = clip(x1);
    xx[1] = clip(tau / x1);

    for (int i = 2; i < n; i++) {
        xx[i] = clip(x[i]);
    }

    double tmp = tau * S - M * M;
    double inv_breit_wigner = (M * M * G * G + tmp * tmp) / (M * G * S);
    double jac = Delta * (1.0 - tau) / x1 * inv_breit_wigner;

    Phasespace::Phasespace ps;
    GenBornPhasespace(&ps, n, xx, userdata);
    double integrand = 0.0;
    int ret = func(ps, xx[3], xx[4], xx[5], jac * wgt, &integrand, userdata);

    *out = integrand * jac;
    return ret;
}

#define NDIM 7

void fill_histograms(void *params, int n, char *msg) {
    UserProcess::Data *ud = (UserProcess::Data *)params;
    Util::DataBuffer *buffer = new Util::DataBuffer(n);
    buffer->AddUInt64((uint64_t)ud->hists->N());
    ud->hists->WriteBinaryToBuffer(buffer);
    memcpy(msg, buffer->Data(), sizeof(char) * n);
    delete buffer;
    ud->hists->Reset();
}

void update_histograms(int num, int n, char *msg, int ncall, double integral,
                       double sigma, void *args) {
    UserProcess::Data *userdata = (UserProcess::Data *)args;
    Util::DataBuffer *buffer = new Util::DataBuffer(n);
    memcpy(buffer->Data(), msg, sizeof(char) * n);
    uint64_t nn = buffer->GetUInt64();
    Histograms *h = new Histograms((int)nn);
    h->InitFrom(*userdata->hists);
    h->ReadBinaryFromBuffer(buffer);
    delete buffer;
    double wgt = 1.0 / ((double)ncall * sigma * sigma);
    userdata->hists->Add(*h, wgt);
    userdata->hists->AddWgt(wgt/(double)num);
    delete h;
}

struct Hists {
    Hists() : hist_lo(0), hist_nlo_qcd(0), hist_nlo_qed(0) {}
    ~Hists() {
        if (hist_lo) {
            delete hist_lo;
        }
        if (hist_nlo_qcd) {
            delete hist_nlo_qcd;
        }
        if (hist_nlo_qed) {
            delete hist_nlo_qed;
        }
    }
    void Init(const Histograms *h) {
        int n = h->N();
        hist_lo = new Histograms(n);
        hist_lo->InitFrom(*h);
        hist_nlo_qcd = new Histograms(n);
        hist_nlo_qcd->InitFrom(*h);
        hist_nlo_qed = new Histograms(n);
        hist_nlo_qed->InitFrom(*h);
    }
    Histograms *hist_lo;
    Histograms *hist_nlo_qcd;
    Histograms *hist_nlo_qed;
};

struct IntegrationResult {
    double I;
    double sigma;
    double chipndf;
};

void fail(int i) {
    MPI_Finalize();
    exit(i);
}

int IntegrateProcess(int rank, int process, bool print_proc, bool load_previous,
                     UserProcess::Data *userdata, IntegrationResult *result,
                     Hists *hists, Run & run) {

    userdata->Reset();

    vegas_integrand integrand_ptr = 0;
    /*********************************************************************
     *                            integration                            *
     *********************************************************************/
    integrand_ptr = integrand3<Powheg::Btilde>;

    VegasState *state = vegas_new(NDIM, 50);
    vegas_set_random(state, Random::init_rnd, NULL, Random::get_rnd,
                     Random::free_rnd);

    // setup grid
    int verbosity = 0;
    int seed1 = userdata->Seed1;
    if (seed1 == 0) {
        if (rank == 0) {
            seed1 = Random::GetSeed();
        }
        MPI_Bcast(&seed1, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    double integral, error, prob;
    if (!load_previous) {
        if (rank == 0) {
            printf(" --- setup grid "
                   "---------------------------------------------------"
                   "---\n");
            printf("    random seed = %d\n", seed1);
        }
        vegas_seed_random(state, seed1);
        vegas_integrate(state, seed1, integrand_ptr, (void *)userdata, 1,
                        userdata->IntParams.nevents_setup,
                        userdata->IntParams.iterations_setup, verbosity,
                        &integral, &error, &prob);
        if (rank == 0) {
            printf("setup result: %g +- %g (chi^2/ndf = %.2f)\n", integral,
                   error, prob);
        }
        vegas_write_grid_to_file(state, "vegas_grid.dat");
    } else {
        const char *gridfile = "vegas_grid.dat";
        int ret = vegas_read_grid_from_file(state, gridfile);
        if (ret != 0) {
            fprintf(stderr, "process %d: error while loading grid from %s\n",
                    rank, gridfile);
            return -1;
        }
        if (rank == 0) {
            printf("using grid from %s\n", gridfile);
        }
        vegas_disable_update_grid(state);
    }

    int seed2 = userdata->Seed2;
    if (!load_previous) {
        vegas_reset_int(state);
        userdata->hists->Reset();
        vegas_disable_update_grid(state);

        vegas_register(state, AFTER_ITERATION, fill_histograms,
                       update_histograms, 1024 * 1024, (void *)userdata);

        if (seed2 == 0) {
            if (rank == 0) {
                seed2 = Random::GetSeed();
            }
            MPI_Bcast(&seed2, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }
        if (rank == 0) {
            printf(" --- integrate "
                   "----------------------------------------------------"
                   "---\n");
            printf("    random seed = %d\n", seed2);
        }
        verbosity = 1;
        userdata->BtildeState.Max = 0.0;
        vegas_seed_random(state, seed2);
        vegas_integrate(state, seed2, integrand_ptr, (void *)userdata, 1,
                        userdata->IntParams.nevents,
                        userdata->IntParams.iterations, verbosity, &integral,
                        &error, &prob);

        double max = -1.0;
        MPI_Reduce(&userdata->BtildeState.Max, &max, 1, MPI_DOUBLE, MPI_MAX,
                   0, MPI_COMM_WORLD);
        if (rank == 0) {
            std::cout << "MAX = " << max << "\n";
            userdata->BtildeState.Max = max;
        }
        MPI_Bcast(&userdata->BtildeState.Max, 1, MPI_DOUBLE, 0,
                  MPI_COMM_WORLD);

        MPI_Reduce(&userdata->BtildeState.MaxVoverB, &max, 1, MPI_DOUBLE, MPI_MAX,
                   0, MPI_COMM_WORLD);
        if (rank == 0) {
            std::cout << "MAX V/B = " << max << "\n";
            userdata->BtildeState.MaxVoverB = 1.1 * max;
        }
        MPI_Bcast(&userdata->BtildeState.MaxVoverB, 1, MPI_DOUBLE, 0,
                  MPI_COMM_WORLD);

        run.Maximum = userdata->BtildeState.Max;
        run.MaximumVoverB = userdata->BtildeState.MaxVoverB;
        run.XSec = integral;
        run.XSecErr = error;
    } else {
        userdata->BtildeState.Max = run.Maximum;
        userdata->BtildeState.MaxVoverB = run.MaximumVoverB;
        if (rank == 0) {
            std::cout << "MAX = " << userdata->BtildeState.Max << "\n";
            std::cout << "MAX V/B = " << userdata->BtildeState.MaxVoverB
                      << "\n";
        }
    }
    if (rank == 0) {
        // always write the run to file because the last run number has to be
        // updated.
        run.Write("integration.dat");
    }

    if (rank == 0) {
        std::cout
            << " --- generate events --------------------------------------\n";
    }
    int seed3 = userdata->GenEvent.seed3;
    if (seed3 == 0) {
        if (rank == 0) {
            seed3 = Random::GetSeed();
        }
        MPI_Bcast(&seed3, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    if (rank == 0) {
        std::cout << "random seed: " << seed3 << "\n";
    }
    double integral_pow;
    double error_pow;
    double prob_pow;
    vegas_reset_int(state);
    userdata->InitEventGeneration();
    LHE::File eventfile;
    if (rank == 0) { // init file
        std::string filename = "events";
        if (run.Number < 10) {
            filename += "0";
        }
        filename += std::to_string(run.Number);
        filename += ".lhe";
        eventfile.Init(filename.c_str());
        std::ifstream config_in("config.txt");
        eventfile.WriteHeader(config_in, false);
        std::stringstream ss;
        ss << "# actually used seeds:\n#  seed1 = " << seed1
           << "\n#  seed2 = " << seed2 << "\n#  seed3 = " << seed3 << "\n";
        eventfile.WriteHeader(ss, false);
        config_in.close();
        LHE::Process proc;
        proc.XSec = run.XSec;
        proc.XSecErr = run.XSecErr;
        eventfile.WriteProcessInfo(userdata->SqrtS, userdata->Lhaid,
                                   { { proc } });
    }
    int num_gen = userdata->GenEvent.N;
    int num_per_it = userdata->GenEvent.NperIt;
    int num_it = num_gen / num_per_it;
    int add = num_gen % num_per_it;
    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    // userdata->EventBuffer.Allocate(num_per_it / 1000 / num_procs * 1.3);
    userdata->EventBuffer.Allocate(10);
    vegas_seed_random(state, seed3);
    for (int i = 0; i < num_it; i++) {
        clock_t t_start = clock();
        int num_events = num_per_it;
        if (i == 0) {
            num_events += add;
        }
        userdata->EventBuffer.SetNMax(num_events);

        vegas_integrand genevents_ptr = 0;
        genevents_ptr = integrand3<Powheg::GenerateEvents>;
        vegas_disable_update_grid(state);
        vegas_integrate(state, seed3, genevents_ptr, (void *)userdata, 1,
                        num_events, 1, 0, &integral_pow, &error_pow, &prob_pow);

        if (rank == 0) {
            eventfile.WriteEvents(userdata->EventBuffer.Length(),
                                  userdata->EventBuffer.Data());
            int num_procs;
            MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
            for (int i = 1; i < num_procs; i++) {
                MPI_Status status;
                int len;
                MPI_Recv(&len, 1, MPI_INT, i, 233, MPI_COMM_WORLD, &status);
                char *buffer = (char *)malloc(sizeof(char) * len);
                MPI_Recv(buffer, len, MPI_CHAR, i, 233, MPI_COMM_WORLD,
                         &status);
                eventfile.WriteEvents(len, buffer);
                free(buffer);
                buffer = NULL;
            }
        } else {
            int len = userdata->EventBuffer.Length();
            MPI_Send(&len, 1, MPI_INT, 0, 233, MPI_COMM_WORLD);
            MPI_Send(userdata->EventBuffer.Data(), len, MPI_CHAR, 0, 233,
                     MPI_COMM_WORLD);
        }
        clock_t t_end = clock();
        double elapsed_secs = (double)(t_end - t_start) / CLOCKS_PER_SEC;
        if (rank == 0) {
            time_t now = time(0);
            int left = num_it - i - 1;
            now += left * elapsed_secs;
            struct tm now_tm = *localtime(&now);
            char time_buffer[40];
            strftime(time_buffer, 40, "%c", &now_tm);
            printf("   events per second: %f\n", num_events / elapsed_secs);
            printf("   end time: %s\n", time_buffer);
        }
        userdata->EventBuffer.Clear();
    }
    int statN = 0;
    int statNenorm = 0;
    int statNeradvar = 0;
    int statNneg = 0;
    int statNmax = 0;
    int statNraderror = 0;
    int statNreject = 0;
    int statNborn = 0;
    int statNreal = 0;
    int statNrejectwovirtual = 0;
    int statN2 = 0;
    int statN3 = 0;
    int statN4 = 0;
    int statNN = 0;
    int statWrongV = 0;

    MPI_Reduce(&userdata->GenEventStatistics.N, &statN, 1, MPI_INT, MPI_SUM, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&userdata->GenEventStatistics.N_ENORM, &statNenorm, 1, MPI_INT,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&userdata->GenEventStatistics.N_NEG, &statNneg, 1, MPI_INT,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&userdata->GenEventStatistics.N_MAX, &statNmax, 1, MPI_INT,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&userdata->GenEventStatistics.N_RADERROR, &statNraderror, 1,
               MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&userdata->GenEventStatistics.N_REJECT, &statNreject, 1, MPI_INT,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&userdata->GenEventStatistics.N_BORN, &statNborn, 1, MPI_INT,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&userdata->GenEventStatistics.N_REAL, &statNreal, 1, MPI_INT,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&userdata->GenEventStatistics.N_ERADVAR, &statNeradvar, 1,
               MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&userdata->GenEventStatistics.N_REJECTWOVIRTUAL,
               &statNrejectwovirtual, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&userdata->GenEventStatistics.N_2TIMES, &statN2, 1, MPI_INT,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&userdata->GenEventStatistics.N_3TIMES, &statN3, 1, MPI_INT,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&userdata->GenEventStatistics.N_4TIMES, &statN4, 1, MPI_INT,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&userdata->GenEventStatistics.N_NTIMES, &statNN, 1, MPI_INT,
               MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&userdata->GenEventStatistics.N_WRONGV, &statWrongV, 1, MPI_INT,
               MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::string efile = eventfile.Filename();
        eventfile.Close();
        if (result) {
            result->I = integral;
            result->sigma = error;
            result->chipndf = prob;
        }
        double swgt = userdata->hists->Swgt();
        hists->hist_nlo_qcd->Add(*userdata->hists, 1.0 / swgt);
        hists->hist_nlo_qcd->SetSwgt(1.0);

        printf(" powheg cross section %f +- %f pb\n", integral_pow, error_pow);
        FILE *out = fopen(efile.c_str(), "a");
        fprintf(out, "# *******************************************\n");
        fprintf(out, "# *  statistics                             *\n");
        fprintf(out, "# *******************************************\n#\n");
        fprintf(out, "# --- Events --------------------------------\n");
        fprintf(out, "# total number of events: %12d\n", statN);
        fprintf(out, "# --------------------------\n");
        fprintf(out, "# rejected events w/o virtual:    %12d\n",
                statNrejectwovirtual);
        fprintf(out, "# rejected events:                %12d\n", statNreject);
        fprintf(out, "# unweighted events:              %12d\n",
                statNreal + statNborn);
        fprintf(out, "#  - born:                        %12d\n", statNborn);
        fprintf(out, "#  - real:                        %12d\n", statNreal);
        fprintf(out, "#  - event written 2 times:       %12d\n", statN2);
        fprintf(out, "#  - event written 3 times:       %12d\n", statN3);
        fprintf(out, "#  - event written 4 times:       %12d\n", statN4);
        fprintf(out, "#  - event written >4 times:      %12d\n", statNN);
        fprintf(out, "# events with errors:\n");
        fprintf(out, "#  - wrong norm:                  %12d\n", statNenorm);
        fprintf(out, "#  - rad. var. out of bounds:     %12d\n", statNeradvar);
        fprintf(out, "#  - neg. weight:                 %12d\n", statNneg);
        fprintf(out, "#  - wrong max:                   %12d\n", statNmax);
        fprintf(out, "#  - rad. error:                  %12d\n", statNraderror);
        fprintf(out, "#  - guessed V too small:         %12d\n", statWrongV);
        fprintf(out, "#\n");
        fprintf(out, "# --- cross section -------------------------\n");
        fprintf(out, "# cross section:        %g +- %g\n", run.XSec,
                run.XSecErr);
        fprintf(out, "# powheg cross section: %g +- %g\n", integral_pow,
                error_pow);
        double diff = integral_pow - run.XSec;
        double diff_err =
            sqrt(run.XSecErr * run.XSecErr + error_pow * error_pow);
        fprintf(out, "#   diff =              %g +- %g (~ %.1f %%)\n", diff,
                diff_err, diff / run.XSec * 100.);
        fclose(out);
    }
    return 0;
}

void print_usage(FILE *stream, char *argv0) {
    fprintf(stream, "usage: %s [options]\n", argv0);
    fprintf(stream, "\noptions:\n");
    fprintf(stream, "  --load   load grid from previous run\n");
}

int main(int argc, char *argv[]) {
    // feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
    MPI_Init(NULL, NULL);

    bool load_previous = false;

    int rank;
    int e = 0;
    Hists hists;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int num_procs;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    if (argc >= 3) {
        if (rank == 0) {
            std::cerr << "error: too many command line arguments\n\n";
                print_usage(stderr, argv[0]);
        }
        fail(1);
    }
    if (argc == 2) {
        if (strcmp(argv[1], "--load") == 0) {
            load_previous = true;
        } else {
            if (rank == 0) {
                std::cerr << "error: unknown command line argument \""
                          << argv[1] << "\"\n\n";
                print_usage(stderr, argv[0]);
            }
            fail(1);
        }
    } else {
        load_previous = false;
    }
    double integral, error, prob;
    VegasState * state;

    std::fstream res_out;
    std::fstream norm_out;

    /*********************************************************************
     *                          set up process                           *
     *********************************************************************/
    auto userdata = std::shared_ptr<ZData>(new ZData);
    if (userdata->Init("config.txt") != 0) {
        fail(1);
    }

    if (rank == 0) {
        std::cout << "Radiation: ";
        if (userdata->RadiationType.QCD) {
            std::cout << "QCD";
            if (userdata->RadiationType.EW) {
                std::cout << ", ";
            }
        }
        if (userdata->RadiationType.EW) {
            std::cout << "EW";
        }
        std::cout << "\n";
    }

    userdata->RadiationRegions = FKS::FindRadiationRegions(userdata->Process);
    if (rank == 0) {
        for (const auto &r : userdata->RadiationRegions) {
            std::cout << "{ fb = " << Physics::PDG::CodesToName(
                                          r.FlavourConfig->Born.Flavours)
                      << ", region = (" << r.Region.I << ", " << r.Region.J
                      << "), real = (";
            for (auto &rr : r.RealFlavour) {
                std::cout << Physics::PDG::CodesToName(rr->Flavours) << " ";
            }
            std::cout << ")\n";
        }
    }

    // add Z resonance
    userdata->ResonanceISR.pdg = 23;
    userdata->ResonanceISR.ID[0] = 2;
    userdata->ResonanceISR.ID[1] = 3;

    userdata->ResonanceFSR.pdg = 23;
    userdata->ResonanceFSR.ID[0] = 2;
    userdata->ResonanceFSR.ID[1] = 3;
    userdata->ResonanceFSR.ID[2] = 4;

    // init pdfs
    std::shared_ptr<PDF::Lhapdf> lhapdf(new PDF::Lhapdf);
    lhapdf->InitByLHAID(userdata->Lhaid);
    userdata->pdf = lhapdf;

    double as_at_mz = lhapdf->AlphaS(91.1876);
    int as_o = lhapdf->GetOrderAlphaS();
    Physics::AlphaSRunning::AlphaSOrder as_order;
    switch (as_o) {
    case 0:
        as_order = Physics::AlphaSRunning::AlphaSOrder::LO;
        break;
    case 1:
        as_order = Physics::AlphaSRunning::AlphaSOrder::NLO;
        break;
    default:
        as_order = Physics::AlphaSRunning::AlphaSOrder::LO;
        fail(1);
    }
    userdata->AlphaS = std::shared_ptr<Physics::IAlphaS>(
        new Physics::AlphaSRunning(as_at_mz, as_order, 5));

    lhapdf->SetThresholdC(sqrt(userdata->AlphaS->ThresholdC2()));
    lhapdf->SetThresholdB(sqrt(userdata->AlphaS->ThresholdB2()));

    if (rank == 0) {
        userdata->Print();
    }

    hists.Init(userdata->hists);

    Run run;
    if (load_previous) {
        if (rank == 0) {
            std::cout << "loading \"integration.dat\"\n";
        }
        if (!run.Read("integration.dat")) {
            if (rank == 0) {
                std::cerr << "    error while loading file!\n";
            }
            fail(1);
        }
        run.Number += 1;
        int i = 0;
        for (auto &r : userdata->RadiationRegions) {
            if (rank == 0) {
                std::cout << Physics::PDG::CodesToName(
                                 r.FlavourConfig->Born.Flavours) << ": ";
                std::cout << "i = " << r.Region.I << ", j = " << r.Region.J;
                std::cout << "{ ";
            }
            for (unsigned j = 0; j < FKS::RadiationRegion::NPDF; j++) {
                r.Norm[j] = run.Norms[i][j];
                if (rank == 0) {
                    std::cout << r.Norm[j];
                    if (j != FKS::RadiationRegion::NPDF - 1) {
                        std::cout << ", ";
                    }
                }
            }
            if (rank == 0) {
                std::cout << " }\n";
            }
            i += 1;
        }
        if (rank == 0) {
            std::cout << "xsec = " << run.XSec << " +- " << run.XSecErr << "\n";
        }
    } else {
        if (rank == 0) {
            std::cout << "--------------------------------------------------\n";
            std::cout << "----- searching norm for upper bounding fct ------\n";
            std::cout << "--------------------------------------------------\n";
        }

        state = vegas_new(NDIM, 4);
        vegas_disable_update_grid(state);
        vegas_set_random(state, Random::init_rnd, NULL, Random::get_rnd,
                         Random::free_rnd);
        vegas_seed_random(state, 0);
        if (rank == 0) {
            std::cout << "  ---- init ----\n";
        }
        vegas_integrand initfindnorm = integrand3<Powheg::InitFindNorm>;
        vegas_integrand findnorm = integrand3<Powheg::FindNorm>;

        vegas_integrate(state, 0, initfindnorm, (void *)&userdata, 1, 50000, 1,
                        1, &integral, &error, &prob);
        for (auto &r : userdata->RadiationRegions) {
            if (rank == 0) {
                for (int i = 1; i < num_procs; i++) {
                    MPI_Status status;
                    int size;
                    MPI_Recv(&size, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
                    Util::DataBuffer buffer(size);
                    MPI_Recv(buffer.Data(), size, MPI_CHAR, i, 1,
                             MPI_COMM_WORLD, &status);
                    HistList histlist;
                    histlist.ReadFromBuffer(&buffer);
                    r.InitHist.Add(histlist);
                }
            } else {
                Util::DataBuffer buffer(Util::DataBuffer::GETSIZE);
                r.InitHist.WriteToBuffer(&buffer);
                int size = buffer.GetDataSize()+100;
                Util::DataBuffer sendbuffer(size);
                r.InitHist.WriteToBuffer(&sendbuffer);
                MPI_Send(&size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
                MPI_Send(sendbuffer.Data(), size, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
            }
            if (rank == 0) {
                Util::DataBuffer buffer(Util::DataBuffer::GETSIZE);
                r.InitHist.WriteToBuffer(&buffer);
                int size = buffer.GetDataSize()+100;
                Util::DataBuffer sendbuffer(size);
                r.InitHist.WriteToBuffer(&sendbuffer);
                for (int i = 1; i < num_procs; i++) {
                    MPI_Send(&size, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
                    MPI_Send(sendbuffer.Data(), size, MPI_CHAR, i, 1,
                             MPI_COMM_WORLD);
                }
            } else {
                MPI_Status status;
                int size;
                MPI_Recv(&size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
                Util::DataBuffer buffer(size);
                MPI_Recv(buffer.Data(), size, MPI_CHAR, 0, 1, MPI_COMM_WORLD,
                         &status);
                r.InitHist.ReadFromBuffer(&buffer);
            }
            r.CreateHistograms(100);
        }
        if (rank == 0) {
            std::cout << "  ---- generating histograms ----\n";
        }
        vegas_integrate(state, 0, findnorm, (void *)userdata.get(), 1, 100000, 1, 1,
                        &integral, &error, &prob);
        vegas_free(state);

        if (rank == 0) {
            // init output file
            norm_out.open("norm.txt", std::ios::out);
            norm_out.close();
            std::cout << "  ---- setting norms ----\n";
        }
        for (auto &r : userdata->RadiationRegions) {
            if (rank == 0) {
                for (int i = 1; i < num_procs; i++) {
                    MPI_Status status;
                    int size;
                    MPI_Recv(&size, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
                    Util::DataBuffer buffer(size);
                    MPI_Recv(buffer.Data(), size, MPI_CHAR, i, 1,
                             MPI_COMM_WORLD, &status);
                    r.MergeNormHistBinaryFromBuffer(&buffer);
                }
                r.ComputeNorm(userdata->UpperBoundingParams.Ratio);
                fprintf(stderr, "new born = %s, i=%d, j=%d, N=%g, %g\n",
                        Physics::PDG::CodesToName(
                            r.FlavourConfig->Born.Flavours).c_str(),
                        r.Region.I, r.Region.J, r.Norm[0], r.Norm[1]);
                std::array<double, FKS::RadiationRegion::NPDF> array;
                for (size_t j = 0; j < array.size(); j++) {
                    array[j] = r.Norm[j];
                }
                run.Norms.push_back(array);
                // write to file
                r.AppendNormHistToFile("norm.txt");
            } else {
                int size = r.GetBinarySize() + 100;
                MPI_Send(&size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
                Util::DataBuffer buffer(size);
                r.WriteNormHistBinaryToBuffer(&buffer);
                MPI_Send(buffer.Data(), size, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
            }
            MPI_Bcast(r.Norm, FKS::RadiationRegion::NPDF, MPI_DOUBLE, 0,
                      MPI_COMM_WORLD);
        }
    }

    IntegrationResult result;
    e = IntegrateProcess(rank, 0, 0, load_previous, userdata.get(), &result, &hists,
                         run);

    if (e != 0) {
        fail(1);
    }
    // write final result
    if (rank == 0) {
        std::string filename = "histograms_pow_qcd.txt";
        res_out.open(filename, std::ios::out);
        if (res_out) {
            hists.hist_nlo_qcd->SetSwgt(1.0);
            res_out << "m_{ll}\n";
            hists.hist_nlo_qcd->WriteToStream(0, res_out, false);
            res_out << "\n\np_{T}\n";
            hists.hist_nlo_qcd->WriteToStream(1, res_out, false);
            res_out << "\n\ny\n";
            hists.hist_nlo_qcd->WriteToStream(2, res_out, false);
            res_out.close();
        }
        filename = "histograms_pow_qed.txt";
        res_out.open(filename, std::ios::out);
        if (res_out) {
            hists.hist_nlo_qed->SetSwgt(1.0);
            res_out << "m_{ll}\n";
            hists.hist_nlo_qed->WriteToStream(0, res_out, false);
            res_out << "\n\np_{T}\n";
            hists.hist_nlo_qed->WriteToStream(1, res_out, false);
            res_out << "\n\ny\n";
            hists.hist_nlo_qed->WriteToStream(2, res_out, false);
            res_out.close();
        }
    }

    MPI_Finalize();
    return 0;
}

