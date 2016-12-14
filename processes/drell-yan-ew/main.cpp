#include <cmath>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <ctime>
#include <cassert>
#include <cstring>

#include "pdf/pdfinterface.h"
#include "pdf/lhapdf.h"

#include "physics/pdgcode.h"
#include "vegas.h"
#include "util/databuffer.h"
#include "util/stringutil.h"
#include "random/rnd.h"
#include "math/math.h"
#include "fks/xsec.h"

#include "processes.h"
#include "config.h"
#include "processlist.h"
#include "process/data.h"
#include "genbornps.h"

typedef int (*function)(const Phasespace::Phasespace &, double, double, double,
                        double, double *, UserProcess::Data *);

template <function func>
static int integrand(int n, const double *x, const double wgt, void *params,
                     int threadid, double *out) {
    UserProcess::Data *userdata = (UserProcess::Data *)params;
    double SqrtS = userdata->SqrtS;
    // Z breit wigner propagtor
    // transfom uniformly distributed random numbers to breit wigner distributed
    // random numbers
    double m = 91.188;
    double Gamma = 2.45;
    double E = (m + 0.5 * Gamma * tan((2.0 * x[0] - 1.0) * Math::Pi / 2.0));
    double tau = E * E / SqrtS / SqrtS;

    double x1 = x[1];
    if (tau > x1) {
        *out = 0.0;
        return 0;
    }
    double xx[n];
    xx[0] = x1;
    xx[1] = tau / x1;
    for (int i = 2; i < n; i++) {
        xx[i] = x[i];
    }
    double integrand = 0.0;
    double Eprime =
        0.5 * Gamma * Math::Pi * (1.0 + 4.0 * pow(((E - m) / Gamma), 2));
    double jac = 0.5 * SqrtS / sqrt(tau) / Eprime;

    double xjac = 1.0 / (x1 * jac * 2.0);

    Phasespace::Phasespace ps;
    GenBornPhasespace(&ps, n, xx, userdata);
    int ret = func(ps, xx[3], xx[4], xx[5], xjac * wgt, &integrand, userdata);

    *out = integrand * xjac;
    return ret;
}

#if 1

#define NDIM 7

void fill_histograms(void *params, int n, char *msg) {
    UserProcess::Data *ud = (UserProcess::Data *)(params);
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

int IntegrateProcess(int rank, int process, const std::string & type, bool print_proc,
                     UserProcess::Data *userdata, IntegrationResult *result,
                     Hists *hists) {

    userdata->Reset();

    userdata->ProcessID = process;

    if (rank == 0 && print_proc) {
        userdata->Process[process].Print();
    }

    vegas_integrand integrand_ptr = 0;
    /*********************************************************************
     *                            integration                            *
     *********************************************************************/
    if (type == "b") {
        integrand_ptr = integrand<FKS::XSecBorn>;
        if (rank == 0) {
            printf("integrate born part\n");
            printf("===================\n");
        }
    } else if (type == "r") {
        integrand_ptr = integrand<FKS::XSecReal>;
        if (rank == 0) {
            printf("integrate real part\n");
            printf("===================\n");
        }
    } else if (type == "v") {
        integrand_ptr = integrand<FKS::XSecVirtualPlusRemnant>;
        if (rank == 0) {
            printf("integrate virtual part\n");
            printf("======================\n");
        }
    }

    VegasState *state = vegas_new(NDIM, 50);
    vegas_set_random(state, Random::init_rnd, NULL, Random::get_rnd, Random::free_rnd);

    // setup grid
    int verbosity = 0;
    int seed1 = userdata->Seed1;
    if(seed1 == 0) {
        seed1 = time(NULL);
    }
    double integral, error, prob;
    if (rank == 0) {
        printf(" --- setup grid "
               "---------------------------------------------------"
               "---\n");
        printf("    random seed = %d\n", seed1);
    }
    vegas_seed_random(state, seed1);
    vegas_integrate(state, seed1, integrand_ptr, (void *)userdata, 1,
                    userdata->IntParams.nevents_setup,
                    userdata->IntParams.iterations_setup, verbosity, &integral,
                    &error, &prob);
    if (rank == 0) {
        printf("setup result: %g +- %g (chi^2/ndf = %.2f)\n", integral, error,
               prob);
    }

    vegas_reset_int(state);
    userdata->hists->Reset();

    vegas_register(state, AFTER_ITERATION, fill_histograms, update_histograms,
                   1024 * 1024, (void *)userdata);

    int seed2 = userdata->Seed2;
    if(seed2 == 0) {
        seed2 = time(NULL);
    }
    if (rank == 0) {
        printf(" --- integrate "
               "----------------------------------------------------"
               "---\n");
        printf("    random seed = %d\n", seed2);
    }
    verbosity = 1;
    vegas_seed_random(state, seed2);
    vegas_integrate(state, seed2, integrand_ptr, (void *)userdata, 1,
                    userdata->IntParams.nevents, userdata->IntParams.iterations,
                    verbosity, &integral, &error, &prob);

    if (rank == 0) {
        std::fstream res_out;
        res_out.open("results.txt", std::ios::out | std::ios::app);
        if (res_out) {
            std::string process = Physics::PDG::CodesToName(
                userdata->Process.at(userdata->ProcessID).Born.Flavours );
            std::string res =
                Strings::Format("%g +- %g (chi^2/ndf = %.2f)\n",
                                integral, error, prob);
            
            res_out << process << ", (" << type << "): " << res;
            res_out.close();
        }
        printf("final result: %g +- %g (chi^2/ndf = %.2f)\n", integral, error,
               prob);
        if(result) {
            result->I = integral;
            result->sigma = error;
            result->chipndf = prob;
        }
        std::string filename =
            Strings::Format("histograms_%d_%s.txt", process, type.c_str());
        std::fstream ostr;
        ostr.open(filename, std::ios::out);
        if (ostr) {
            ostr << "m_{ll}\n";
            userdata->hists->WriteToStream(0, ostr, false);
            ostr << "\n\np_{T}\n";
            userdata->hists->WriteToStream(1, ostr, false);
            ostr << "\n\ny\n";
            userdata->hists->WriteToStream(2, ostr, false);
            ostr.close();
        }
        double swgt = userdata->hists->Swgt();
        if(type == "b") {
            hists->hist_lo->Add(*(userdata->hists), 1.0 / swgt);
            hists->hist_lo->SetSwgt(1.0);
        } else if (type == "r") {
            hists->hist_nlo_qcd->Add(*userdata->hists, 1.0 / swgt);
            hists->hist_nlo_qcd->SetSwgt(1.0);
        } 
    }
    return 0;
}

int fail(int c) {
    MPI_Finalize();
    exit(c);
}

int main(int argc, char *argv[]) {
    MPI_Init(NULL, NULL);
    int rank;
    int i = 1;
    double I_lo = 0.0;
    double sigma_2_lo = 0.0;
    double I_nlo_qcd = 0.0;
    double sigma_2_nlo_qcd = 0.0;
    double I_nlo_qed = 0.0;
    double sigma_2_nlo_qed = 0.0;
    int previous_proc_id;
    Hists hists;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::fstream res_out;

    /*********************************************************************
     *                          set up process                           *
     *********************************************************************/
    UserProcess::Data userdata;
    ProcessList plist;
    int e = ReadProcessFile(&plist, "processes.txt");
    if (e != 0) {
        fail(1);
    }
    if (userdata.Init("config.txt") != 0) {
        fail(1);
    }

    if (rank == 0) {
        std::cout << "Radiation: ";
        if (userdata.RadiationType.QCD) {
            std::cout << "QCD";
            if (userdata.RadiationType.EW) {
                std::cout << ", ";
            }
        }
        if (userdata.RadiationType.EW) {
            std::cout << "EW";
        }
        std::cout << "\n";
    }
    userdata.Process = GenerateProcesses(userdata.RadiationType.QCD,
                                         userdata.RadiationType.EW);

    std::shared_ptr<PDF::Lhapdf> lhapdf(new PDF::Lhapdf);
    lhapdf->InitByLHAID(userdata.Lhaid);
    userdata.pdf = lhapdf;

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
    userdata.AlphaS = std::shared_ptr<Physics::IAlphaS>(
        new Physics::AlphaSRunning(as_at_mz, as_order, 5));

    if (rank == 0) {
        userdata.Print();
        // init output file
        res_out.open("results.txt", std::ios::out);
        res_out.close();
    }

    i = 1;
    I_lo = 0.0;
    sigma_2_lo = 0.0;
    I_nlo_qcd = 0.0;
    sigma_2_nlo_qcd = 0.0;
    I_nlo_qed = 0.0;
    sigma_2_nlo_qed = 0.0;
    hists.Init(userdata.hists);
    previous_proc_id = -1;
    for (auto &proc : plist) {
        if (rank == 0) {
            printf("/**********************************************************"
                   "***********\n");
            printf("* Process %d / %lu\n", i, plist.size());
            printf(" **********************************************************"
                   "***********\n");
        }
        IntegrationResult result;
        bool print_proc = true;
        if (previous_proc_id == proc.Id) {
            print_proc = false;
        }
        e = IntegrateProcess(rank, proc.Id, proc.Type, print_proc, &userdata, &result, &hists);
        if (proc.Type == "b") {
            I_lo += result.I;
            sigma_2_lo += result.sigma * result.sigma;
        } else if(proc.Type == "r" || proc.Type == "v") {
            I_nlo_qcd += result.I;
            sigma_2_nlo_qcd += result.sigma * result.sigma;
        } 

        if (e != 0) {
            fail(1);
        }
        i += 1;
        previous_proc_id = proc.Id;
    }
    // write final result
    if (rank == 0) {
        res_out.open("results.txt", std::ios::out | std::ios::app);
        if (res_out) {
            res_out << "======================================================="
                       "=========\n";
            res_out << Strings::Format("sigma(LO)           : %10g +- %10g\n",
                                       I_lo, sqrt(sigma_2_lo));
            res_out << Strings::Format("NLO correction (QCD): %10g +- %10g\n",
                                       I_nlo_qcd, sqrt(sigma_2_nlo_qcd));
            res_out << Strings::Format("NLO correction (QED): %10g +- %10g\n",
                                       I_nlo_qed, sqrt(sigma_2_nlo_qed));
            res_out << Strings::Format("sigma(NLO, QCD)     : %10g +- %10g\n",
                                       I_lo + I_nlo_qcd,
                                       sqrt(sigma_2_lo + sigma_2_nlo_qcd));
            res_out << Strings::Format("sigma(NLO, QED)     : %10g +- %10g\n",
                                       I_lo + I_nlo_qed,
                                       sqrt(sigma_2_lo + sigma_2_nlo_qed));
            res_out.close();
        }
        std::string filename = "histograms_lo.txt";
        res_out.open(filename, std::ios::out);
        if (res_out) {
            hists.hist_lo->SetSwgt(1.0);
            res_out << "m_{ll}\n";
            hists.hist_lo->WriteToStream(0, res_out, false);
            res_out << "\n\np_{T}\n";
            hists.hist_lo->WriteToStream(1, res_out, false);
            res_out << "\n\ny\n";
            hists.hist_lo->WriteToStream(2, res_out, false);
            res_out.close();
        }
        filename = "histograms_correction_qcd.txt";
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
        filename = "histograms_correction_qed.txt";
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

#else

int main(int argc, char *argv[]) {
    LHAPDF::setVerbosity(LHAPDF::SILENT);
    LHAPDF::initPDFSet("cteq6ll", LHAPDF::LHPDF, 0);
    UserProcess::Data userdata;
    userdata.Process = GenerateProcesses();
    userdata.FlavourConfig = userdata.Process.Connect(
        "d~ d -> mu+ mu-", 2, { { 2 } }, { { 2 } }, { { 2 } });
    if (userdata.Init("config.txt") != 0) {
        exit(1);
    }

    int ndim = 7;
    int ncomp = 1;
    double result;
    double x[7] = { 0.0250607, 0.00167993, 0.890943, 1.41139e-06,
                    0.999946,  0.937026,   0.40381 };

    mc_real(&ndim, x, 1.0, &result, (void *)&userdata);

    printf("result = %g\n", result);
    return 0;
}
#endif
