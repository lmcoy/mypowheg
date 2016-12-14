#include "eventgenerator.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <mpi.h>

#include "vegas.h"

#include "physics/pdgcode.h"
#include "util/databuffer.h"
#include "util/stringutil.h"
#include "random/rnd.h"
#include "math/math.h"
#include "powheg/btilde.h"
#include "powheg/findnorm.h"
#include "lhe/lhefile.h"
#include "pdf/lhapdf.h"

#define FAIL -1
#define SUCCESS 0

static int integrand(int n, const double *x, const double wgt, void *params,
                     int threadid, double *out) {
    auto *p = (IntegrandParams *)params;
    auto psgen = p->psgen;
    auto *userdata = p->userdata;
    auto trafo = p->trafo;
    auto func = p->func;

    double xx[n];
    double jac = trafo->Transform(n, x, xx, userdata);

    auto ps = psgen->Gen(n, xx, userdata);
    double integrand = 0.0;
    int ret = func(ps, xx[n - 3], xx[n - 2], xx[n - 1], jac * wgt, &integrand,
                   userdata);

    *out = integrand * jac;
    return ret;
}

EventGenerator::EventGenerator(IntegralTransformationPtr trafo,
                               PhasespaceGeneratorPtr psgen,
                               std::shared_ptr<UserProcess::Data> data)
    : userdata(data), trafo_(trafo), psgen_(psgen) {
    // born + 3 radiation variables
    NDIM = psgen_->Dim() + 3;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
}

EventGenerator::~EventGenerator() {
    if (vegasstate_) {
        vegas_free(vegasstate_);
    }
}

int EventGenerator::Init() {
    if (userdata->Init("config.txt") != 0) {
        return FAIL;
    }

    userdata->RadiationRegions = FKS::FindRadiationRegions(userdata->Process);
    if (verbose && rank == 0) {
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

    // init pdfs
    std::shared_ptr<PDF::Lhapdf> lhapdf(new PDF::Lhapdf);
    lhapdf->InitByLHAID(userdata->Lhaid);
    userdata->pdf = lhapdf;

    // init alpha_s
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
        return FAIL;
    }
    userdata->AlphaS = std::shared_ptr<Physics::IAlphaS>(
        new Physics::AlphaSRunning(as_at_mz, as_order, 5));

    lhapdf->SetThresholdC(sqrt(userdata->AlphaS->ThresholdC2()));
    lhapdf->SetThresholdB(sqrt(userdata->AlphaS->ThresholdB2()));

    if (verbose && rank == 0) {
        userdata->Print();
    }
    return SUCCESS;
}

void EventGenerator::Setup() {
    SearchNormForUpperBounding();
    Integrate();
}

void EventGenerator::SearchNormForUpperBounding() {
    // --------------------------------------------------------------------
    // --------------------------------------------------------------------
    // init vegas etc.
    // --------------------------------------------------------------------
    // --------------------------------------------------------------------
    if (rank == 0) {
        std::cout << "--------------------------------------------------\n";
        std::cout << "----- searching norm for upper bounding fct ------\n";
        std::cout << "--------------------------------------------------\n";
    }
    double integral = 0.0, error = 0.0, prob = 0.0;
    VegasState *state = vegas_new(NDIM, 4);
    vegas_disable_update_grid(state);
    vegas_set_random(state, Random::init_rnd, NULL, Random::get_rnd,
                     Random::free_rnd);
    vegas_seed_random(state, 0);

    // --------------------------------------------------------------------
    // --------------------------------------------------------------------
    // find a rough estimate for the norm
    // --------------------------------------------------------------------
    // --------------------------------------------------------------------
    if (rank == 0) {
        std::cout << "  ---- init ----\n";
    }
    auto init = GenIntegrandParams(Powheg::InitFindNorm);
    vegas_integrate(state, 0, integrand, (void *)&init, 1,
                    userdata->UpperBoundingParams.Ninit, 1, 1, &integral,
                    &error, &prob);
    for (auto &r : userdata->RadiationRegions) {
        if (rank == 0) {
            for (int i = 1; i < num_procs; i++) {
                MPI_Status status;
                int size;
                MPI_Recv(&size, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
                Util::DataBuffer buffer(size);
                MPI_Recv(buffer.Data(), size, MPI_CHAR, i, 1, MPI_COMM_WORLD,
                         &status);
                HistList histlist;
                histlist.ReadFromBuffer(&buffer);
                r.InitHist.Add(histlist);
            }
        } else {
            Util::DataBuffer buffer(Util::DataBuffer::GETSIZE);
            r.InitHist.WriteToBuffer(&buffer);
            int size = buffer.GetDataSize() + 100;
            Util::DataBuffer sendbuffer(size);
            r.InitHist.WriteToBuffer(&sendbuffer);
            MPI_Send(&size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(sendbuffer.Data(), size, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
        }
        if (rank == 0) {
            Util::DataBuffer buffer(Util::DataBuffer::GETSIZE);
            r.InitHist.WriteToBuffer(&buffer);
            int size = buffer.GetDataSize() + 100;
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

    // --------------------------------------------------------------------
    // --------------------------------------------------------------------
    // search better value for the norm
    // --------------------------------------------------------------------
    // --------------------------------------------------------------------
    if (rank == 0) {
        std::cout << "  ---- generating histograms ----\n";
    }
    auto findnorm = GenIntegrandParams(Powheg::FindNorm);
    vegas_integrate(state, 0, integrand, (void *)&findnorm, 1,
                    userdata->UpperBoundingParams.N, 1, 1, &integral, &error,
                    &prob);
    vegas_free(state);

    if (rank == 0) {
        // init output file
        std::fstream norm_out;
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
                MPI_Recv(buffer.Data(), size, MPI_CHAR, i, 1, MPI_COMM_WORLD,
                         &status);
                r.MergeNormHistBinaryFromBuffer(&buffer);
            }
            r.ComputeNorm(userdata->UpperBoundingParams.Ratio);
            fprintf(stderr, "new born = %s, i=%d, j=%d, ",
                    Physics::PDG::CodesToName(r.FlavourConfig->Born.Flavours)
                        .c_str(),
                    r.Region.I, r.Region.J );
            for (size_t i = 0; i < FKS::RadiationRegion::NPDF; i++) {
                fprintf(stderr, ", %g", r.Norm[i]);
            }
            fprintf(stderr, "\n");
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

void EventGenerator::Integrate() {
    userdata->Reset();

    /*********************************************************************
     *                            integration                            *
     *********************************************************************/

    vegasstate_ = vegas_new(NDIM, 50);
    vegas_set_random(vegasstate_, Random::init_rnd, NULL, Random::get_rnd,
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
    if (rank == 0) {
        printf(" --- setup grid "
               "---------------------------------------------------"
               "---\n");
        printf("    random seed = %d\n", seed1);
    }

    auto btilde = GenIntegrandParams(Powheg::Btilde);
    vegas_seed_random(vegasstate_, seed1);
    vegas_integrate(vegasstate_, seed1, integrand, (void *)&btilde, 1,
                    userdata->IntParams.nevents_setup,
                    userdata->IntParams.iterations_setup, verbosity, &integral,
                    &error, &prob);
    if (rank == 0) {
        printf("setup result: %g +- %g (chi^2/ndf = %.2f)\n", integral, error,
               prob);
    }
    vegas_write_grid_to_file(vegasstate_, "vegas_grid.dat");

    int seed2 = userdata->Seed2;
    vegas_reset_int(vegasstate_);
    userdata->hists->Reset();
    vegas_disable_update_grid(vegasstate_);

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
    vegas_seed_random(vegasstate_, seed2);
    vegas_integrate(vegasstate_, seed2, integrand, (void *)&btilde, 1,
                    userdata->IntParams.nevents, userdata->IntParams.iterations,
                    verbosity, &integral, &error, &prob);

    double max = -1.0;
    MPI_Reduce(&userdata->BtildeState.Max, &max, 1, MPI_DOUBLE, MPI_MAX, 0,
               MPI_COMM_WORLD);
    if (rank == 0) {
        std::cout << "MAX = " << max << "\n";
        userdata->BtildeState.Max = max;
    }
    MPI_Bcast(&userdata->BtildeState.Max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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

    if (rank == 0) {
        // always write the run to file because the last run number has to be
        // updated.
        run.Write("integration.dat");
    }
}

int EventGenerator::GenerateEvents() {
    if (rank == 0) {
        std::cout
            << " --- generate events --------------------------------------\n";
    }
    if (vegasstate_ == 0) {
        if (rank == 0) {
            std::cerr << "error: no vegas state. call Integrate() first.\n";
        }
        return FAIL;
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
    vegas_reset_int(vegasstate_);
    userdata->InitEventGeneration();
    LHE::File eventfile;
    if (rank == 0) { // init event file
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
        int seed1 = userdata->Seed1;
        int seed2 = userdata->Seed2;
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

    userdata->EventBuffer.Allocate(10);
    vegas_seed_random(vegasstate_, seed3);
    for (int i = 0; i < num_it; i++) {
        clock_t t_start = clock();
        int num_events = num_per_it;
        if (i == 0) {
            num_events += add;
        }
        userdata->EventBuffer.SetNMax(num_events);

        // generate events
        vegas_disable_update_grid(vegasstate_);
        auto integrandParams = GenIntegrandParams(Powheg::GenerateEvents);
        vegas_integrate(vegasstate_, seed3, integrand, (void *)&integrandParams,
                        1, num_events, 1, 0, &integral_pow, &error_pow,
                        &prob_pow);

        // write events to file
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
        userdata->EventBuffer.Clear();

        // print runtime for this iteration
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
    return SUCCESS;
}
