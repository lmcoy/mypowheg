#include "process/data.h"

#include <iostream>
#include <complex>
#include "config/file.h"

#include "cuts.h"
#include "me/parameters_sm.h"
#include "myhistograms.h"
#include "scales.h"

using namespace UserProcess;

static void complex_to_str(int n, char * buffer, std::complex<double> c) {
    char sign = '+';
    if (c.imag() < 0.0) {
        sign = '-';
    }
    if (c.real() == 0.0) {
        snprintf(buffer, n, "%g i", c.imag());
    }
    if (c.imag() == 0.0 || c.imag() == -0.0) {
        snprintf(buffer, n, "%g", c.real());
    }
    snprintf(buffer, n, "%g %c %g i", c.real(), sign, fabs(c.imag()));
}

int Data::UserInit(Config::File & cfile) {
    Cuts *icuts = new Cuts();
    if (cfile.GetDoubleInterval("cuts", "mll", &icuts->mllmin, &icuts->mllmax) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }

    if (cfile.GetDoubleInterval("cuts", "pT", &icuts->pTmin, &icuts->pTmax) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    if (cfile.GetDoubleInterval("cuts", "eta", &icuts->EtaMin, &icuts->EtaMax) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    if (cfile.GetDoubleInterval("cuts", "y", &icuts->Ymin, &icuts->Ymax) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }

    Strings::IntList hist_mll;
    if (cfile.GetIntList("histograms", "mll", &hist_mll) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }

    Strings::IntList hist_pt;
    if (cfile.GetIntList("histograms", "pt(mu+)", &hist_pt) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }

    Strings::IntList hist_y;
    if (cfile.GetIntList("histograms", "y(mu+)", &hist_y) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }

    std::string mu_string = "";
    if (cfile.GetString("scales", "mu", &mu_string) != Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    std::string muF_string = "";
    if (cfile.GetString("scales", "muF", &muF_string) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }

    double mu = 4.0;
    int mu_parse = Strings::ParseDouble(mu_string, &mu);
    double muF = 0.0;
    int muF_parse = Strings::ParseDouble(muF_string, &muF);

    if (mu_parse == 0 && muF_parse == 0) {
        Scales = new FixedScales(mu, muF);
    } else if (mu_string == "mll" && muF_string == "mll") {
        Scales = new InvMassScales;
    } else {
        std::cerr << "unknown scales\n";
        return 1;
    }

    cuts = icuts;

    // set params
    auto params = new Parameters_sm;
    params->Set(1.16637E-005, 132.507, 91.1876, 2.4952);
    Params = params;
    params = 0;

    Params_as = new Parameters_alphaS;

    hists = new MyHistograms(3);
    hists->Init1D(0, hist_mll[0], hist_mll[1], hist_mll[2]);
    hists->Init1D(1, hist_pt[0], hist_pt[1], hist_pt[2]);
    hists->Init1D(2, hist_y[0], hist_y[1], hist_y[2]);
    return 0;
}

void Data::UserFree() {
    if (hists) {
        delete hists;
    }
    if (cuts) {
        delete cuts;
    }
    if (Params) {
        delete Params;
    }
    if (Params_as) {
        delete Params_as;
    }
    if (Scales) {
        delete Scales;
    }
}

void Data::UserPrint() const {
    printf("/******************************************************************"
           "***\n"
           " *                               cuts                              "
           "  *\n"
           " ******************************************************************"
           "***/\n");
    printf(" %-6s = (%g, %g)\n", "mll", cuts->mllmin,
           cuts->mllmax);
    printf(" %-6s = (%g, %g)\n", "pT", cuts->pTmin,
           cuts->pTmax);
    printf(" %-6s = (%g, %g)\n", "eta", cuts->EtaMin,
           cuts->EtaMax);
    printf(" %-6s = (%g, %g)\n", "y", cuts->Ymin, cuts->Ymax);
    printf("\n");

    printf("/******************************************************************"
           "***\n"
           " *                           recombination                         "
           "  *\n"
           " ******************************************************************"
           "***/\n");
    printf(" %-6s = %g\n", "dR", Recomb.dR);
    printf("\n");

    printf("/******************************************************************"
           "***\n"
           " *                             Parameter                           "
           "  *\n"
           " ******************************************************************"
           "***/\n");

    Parameters_sm * params = dynamic_cast<Parameters_sm*>(Params);
    printf(" %-6s = %g\n", "MW", params->MW);
    printf(" %-6s = %g\n", "MZ", params->MZ);
    printf(" %-6s = %g\n", "WidthW", params->WidthW);
    printf(" %-6s = %g\n", "WidthZ", params->WZ);
    printf(" %-6s = %g\n", "MH", params->MZ);
    printf(" %-6s = %g\n", "WidthH", params->WidthH);
    printf(" %-6s = %g\n", "alpha", params->alpha);
    printf("\n");

    char buffer[100];
    printf(" --- Complex Masses -----------------------------------------------"
           "---\n");
    complex_to_str(100, buffer, params->MuZ);
    printf(" %-6s = %s\n", "mu_Z", buffer);
    printf("\n");

    printf(" --- Weinberg Angle -----------------------------------------------"
           "---\n");
    complex_to_str(100, buffer, params->cw);
    printf(" %-6s = %s\n", "cos th_w", buffer);
    complex_to_str(100, buffer, params->sw);
    printf(" %-6s = %s\n", "sin th_w", buffer);
    printf("\n");

    printf(" --- Quark Masses -------------------------------------------------"
           "---\n");
    printf(" %-6s = %g\n", "u", params->MUQuark);
    printf(" %-6s = %g\n", "d", params->MDQuark);
    printf(" %-6s = %g\n", "c", params->MCQuark);
    printf(" %-6s = %g\n", "s", params->MSQuark);
    printf(" %-6s = %g\n", "t", params->MTQuark);
    printf(" %-6s = %g\n", "b", params->MBquark);
    printf("\n");

    printf(" --- Lepton Masses ------------------------------------------------"
           "---\n");
    printf(" %-6s = %g\n", "e", params->MElectron);
    printf(" %-6s = %g\n", "mu", params->MMuon);
    printf(" %-6s = %g\n", "tau", params->MTau);
}
