#include "wdata.h"

#include "config/file.h"

#include "cuts.h"
#include "me/parameters_sm.h"
#include "myhistograms.h"
#include "scales.h"
#include "matrixelement.h"

using namespace UserProcess;

static FKS::ProcessList generateProcesses(bool QCD, bool EW);

WData::~WData() {
    if (Params) {
        delete Params;
        Params = 0;
    }

    if (hists) {
        delete hists;
        hists = 0;
    }

    if (Scales) {
        delete Scales;
        Scales = 0;
    }
}

int WData::ProcessInit(Config::File &cfile) {
    auto icuts = std::make_shared<WCuts>();

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
    if (cfile.GetString("scales", "mu", &mu_string) !=
        Config::File::Error::NoError) {
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

    auto me = std::make_shared<MatrixElementW>();
    me->Init(Params, AlphaScheme::Gmu, .0);

    MatrixElement = me;

    hists = new MyHistograms(3);
    hists->Init1D(0, hist_mll[0], hist_mll[1], hist_mll[2]);
    hists->Init1D(1, hist_pt[0], hist_pt[1], hist_pt[2]);
    hists->Init1D(2, hist_y[0], hist_y[1], hist_y[2]);

    Process = generateProcesses(RadiationType.QCD, RadiationType.EW);

    // add W+ resonance
    ResonanceISR.pdg = 24;
    ResonanceISR.ID[0] = 2;
    ResonanceISR.ID[1] = 3;

    ResonanceFSR.pdg = 24;
    ResonanceFSR.ID[0] = 2;
    ResonanceFSR.ID[1] = 3;
    ResonanceFSR.ID[2] = 4;
    return 0;
}

static void complex_to_str(int n, char *buffer, std::complex<double> c) {
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

void WData::ProcessPrint() const {

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

    Parameters_sm *params = dynamic_cast<Parameters_sm *>(Params);
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

FKS::ProcessList generateProcesses(bool QCD, bool EW) {
    using R = MatrixElementW::SubProcesses;
    FKS::ProcessList list;
    double Vud = 0.9748;
    double Vus = 0.2225;
    double Vub = 0.0036;
    double Vcd = 0.2225;
    double Vcs = 0.9740;
    double Vcb = 0.041;
    double Vud2 = Vud * Vud;
    double Vcd2 = Vcd * Vcd;
    double Vus2 = Vus * Vus;
    double Vcs2 = Vcs * Vcs;
    double Vub2 = Vub * Vub;
    double Vcb2 = Vcb * Vcb;
    {
        std::vector<int> pdfs = {
            { -1, 2, -1, 4, -3, 2, -3, 4, -5, 2, -5, 4 }
        };
        std::vector<double> ckm = { { Vud2, Vcd2, Vus2, Vcs2, Vub2, Vcb2 } };
        FKS::FlavourConfig fl1(R::QXQ_MUPNU, { { -1, 2, -13, 14 } }, pdfs, ckm);
        if (EW) {
            FKS::ResonanceList resonances;
            size_t isr_res = resonances.Add(
                FKS::Resonance({ { 2, 3 } }, 80.385, 2.085, 2.0 / 3.0));
            size_t fsr_res = resonances.Add(
                FKS::Resonance({ { 2, 3, 4 } }, 80.385, 2.085, 1.0));
            FKS::RegionList rlist;
            // add all real flavour structures
            rlist.push_back(FKS::Region(4, 2, fsr_res));
            rlist.push_back(FKS::Region(4, 0, isr_res));
            fl1.AddReal(R::QXQ_MUPNU_A, FKS::Type_t::EW,
                        { { -1, 2, -13, 14, 22 } }, pdfs, rlist, resonances);
        }
        if (QCD) {
            FKS::RegionList rlist;
            // add all real flavour structures
            rlist.push_back(FKS::Region(4, 0));
            fl1.AddReal(R::QXQ_MUPNU_G, FKS::Type_t::QCD,
                        { { -1, 2, -13, 14, 21 } }, pdfs, rlist);
            std::vector<int> pdfs_g2 = { { -1, 0, -1, 0, -3, 0, -3, 0,-5, 0,
            -5, 0 } };
            fl1.AddReal(R::QXG_MUPNU_QX, FKS::Type_t::QCD,
                        { { -1, 21, -13, 14, -2 } }, pdfs_g2, rlist);
            std::vector<int> pdfs_g1 = { { 0, 2, 0, 4, 0, 2, 0, 4, 0, 2, 0, 4 }
            };
            fl1.AddReal(R::GQ_MUPNU_Q, FKS::Type_t::QCD, { { 21, 2, -13, 14, 1 }
            },
                        pdfs_g1, rlist);
        }
        list.push_back(fl1);
    }

    {
        std::vector<int> pdfs = {
	  { 2, -1, 2,-3, 2,-5, 4, -1, 4,-3, 4,-5 }
        };
        std::vector<double> ckm = { { Vud2, Vus2, Vub2, Vcd2, Vcs2, Vcb2 } };
        FKS::FlavourConfig fl1(R::QQX_MUPNU, { { 2, -1, -13, 14 } }, pdfs, ckm);
        if (EW) {
            FKS::ResonanceList resonances;
            size_t isr_res = resonances.Add(
                FKS::Resonance({ { 2, 3 } }, 80.385, 2.085, 2.0 / 3.0));
            size_t fsr_res = resonances.Add(
                FKS::Resonance({ { 2, 3, 4 } }, 80.385, 2.085, 1.0));
            FKS::RegionList rlist;
            // add all real flavour structures
            rlist.push_back(FKS::Region(4, 2, fsr_res));
            rlist.push_back(FKS::Region(4, 0, isr_res));
            fl1.AddReal(R::QQX_MUPNU_A, FKS::Type_t::EW,
                        { { 2, -1, -13, 14, 22 } }, pdfs, rlist, resonances);
        }
        if (QCD) {
            FKS::RegionList rlist;
            // add all real flavour structures
            rlist.push_back(FKS::Region(4, 0));
            fl1.AddReal(R::QQX_MUPNU_G, FKS::Type_t::QCD,
                        { { 2, -1, -13, 14, 21 } }, pdfs, rlist);
            std::vector<int> pdfs_g2 = {{ 2, 0, 2,0, 2,0, 4, 0, 4,0, 4,0 }};
            fl1.AddReal(R::QG_MUPNU_Q, FKS::Type_t::QCD,
                        { { 2, 21, -13, 14, 1 } }, pdfs_g2, rlist);
            std::vector<int> pdfs_g1 ={{ 0, -1, 0, -3, 0, -5, 0, -1, 0, -3, 0, -5 }};
            fl1.AddReal(R::GQX_MUPNU_QX, FKS::Type_t::QCD, { { 21, -1, -13, 14,
            -2 } },
                        pdfs_g1, rlist);
        }
        list.push_back(fl1);
    }

    return list;
}
