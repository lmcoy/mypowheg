#include "wjdata.h"

#include "config/file.h"

#include "cuts.h"
#include "matrixelement.h"
#include "me/parameters_sm.h"
#include "myhistograms.h"
#include "scales.h"

using namespace UserProcess;

static FKS::ProcessList generateProcesses(bool QCD, bool EW, bool CutOnReal,
                                          FKS::Param *params);

WJData::~WJData() {
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

int WJData::ProcessInit(Config::File &cfile) {
    auto icuts = std::make_shared<WjCuts>();

    {
        double PTCUT = 0.0;
        auto status = cfile.GetDouble("cuts", "pT", &PTCUT); 
        switch (status) {
            case Config::File::Error::CategoryNotFound:
            case Config::File::Error::KeyNotFound:
                break;
            case Config::File::Error::NoError:
                icuts->PTCUT = PTCUT;
                break;
            default:
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (icuts->PTCUT <= 0.0) {
            std::cerr << "error: PTCUT = " << icuts->PTCUT << "\n";
            return 1;
        }
    }

    {
        double jetDR = 0.0;
        auto status = cfile.GetDouble("cuts", "jetDR", &jetDR); 
        switch (status) {
            case Config::File::Error::CategoryNotFound:
            case Config::File::Error::KeyNotFound:
                break;
            case Config::File::Error::NoError:
                icuts->jetDR = jetDR;
                break;
            default:
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (icuts->jetDR <= 0.0) {
            std::cerr << "error: jetDR <= 0.0\n";
            return 1;
        }
    }


    {
        double photonDR = 0.0;
        auto status = cfile.GetDouble("cuts", "photonDR", &photonDR); 
        switch (status) {
            case Config::File::Error::CategoryNotFound:
            case Config::File::Error::KeyNotFound:
                break;
            case Config::File::Error::NoError:
                icuts->photonDR = photonDR;
                break;
            default:
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (icuts->photonDR <= 0.0) {
            std::cerr << "error: photonDR <= 0.0\n";
            return 1;
        }
    }

    {
        double z_thr = 0.0;
        auto status = cfile.GetDouble("cuts", "z_thr", &z_thr); 
        switch (status) {
            case Config::File::Error::CategoryNotFound:
            case Config::File::Error::KeyNotFound:
                break;
            case Config::File::Error::NoError:
                icuts->z_thr = z_thr;
                break;
            default:
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (icuts->z_thr < 0.5 || icuts->z_thr > 1.0) {
            std::cerr << "error: z_thr < 0.5 || z_thr > 1.0\n";
            return 1;
        }
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

    if (mu_string == "PTW" && muF_string == "PTW" ) {
        Scales = new PTWScales();
    } else {
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
    }

    Strings::IntList hist_ptj;
    if (cfile.GetIntList("histograms", "ptj", &hist_ptj) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    if (hist_ptj.size() != 3) {
        std::cerr << "error: histograms->ptj: expected exactly 3 integers "
                     "(nbins, xlow, xup).\n";
        return 1;
    }
    if (hist_ptj[0] <= 0) {
        std::cerr << "error: histograms->ptj: nbins has to be > 0.\n";
        return 1;
    }

    Strings::IntList hist_ptl;
    if (cfile.GetIntList("histograms", "ptl", &hist_ptl) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    if (hist_ptl.size() != 3) {
        std::cerr << "error: histograms->ptl: expected exactly 3 integers "
                     "(nbins, xlow, xup).\n";
        return 1;
    }
    if (hist_ptl[0] <= 0) {
        std::cerr << "error: histograms->ptl: nbins has to be > 0.\n";
        return 1;
    }

    cuts = icuts;

    // set params
    auto params = new Parameters_sm;
    params->Set(1.16637E-005, 132.507, 91.1876, 2.4952);
    Params = params;
    params = 0;

    bool MadGraphParams = false;
    {
        long usemad = 0;
        bool in = false;
        auto status = cfile.GetInt("Parameter", "UseMadGraphDefault", &usemad); 
        switch (status) {
            case Config::File::Error::CategoryNotFound:
            case Config::File::Error::KeyNotFound:
                break;
            case Config::File::Error::NoError:
                in = true;
                break;
            default:
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (in && usemad > 0) {
            MadGraphParams = true;
        } 
    }
    if (MadGraphParams) {
        params = new Parameters_sm;
        double aew = 1.325070e+02;
        double mz = 9.1188e+01;
        params->MW = 8.04190024457561634108e+01;
        params->WidthW = 2.0476;

        double Gf = M_PI / sqrt(2.0) / aew / params->MW / params->MW *
                    (1.0 / (1.0 - params->MW * params->MW / mz / mz));

        params->Set(Gf, aew, mz, 2.441404);

        params->Vud = 0.97461995;
        params->Vus = 0.2253;
        params->Vub = 0.0033788488;
        params->Vcd = 0.2253;
        params->Vcs = 0.97461995;
        params->Vcb = 0.04101415;
        delete Params;
        Params = params;
        params = 0;
    }

    auto me = std::make_shared<MatrixElementWj>();
    me->Init(Params, AlphaScheme::Gmu, .0);

    MatrixElement = me;

    hists = new MyHistograms(2);
    hists->Init1D(0, hist_ptj[0], hist_ptj[1], hist_ptj[2]);
    hists->SetName(0, "p_T of the leading jet");
    hists->Init1D(1, hist_ptl[0], hist_ptl[1], hist_ptl[2]);
    hists->SetName(1, "p_T of muon");

    Process = generateProcesses(RadiationType.QCD, RadiationType.EW, CutOnReal,
                                Params);

    // add W+ resonance
    Resonance.pdg = 24;
    Resonance.ID[0] = 2;
    Resonance.ID[1] = 3;

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

void WJData::ProcessPrint() const {

    auto wjcuts = std::static_pointer_cast<WjCuts>(cuts);
    printf("/******************************************************************"
           "***\n"
           " *                           cuts                                  "
           "  *\n"
           " ******************************************************************"
           "***/\n");
    printf(" %-6s = %g\n", "pT(l)", wjcuts->PTCUT);
    printf(" NLO only cuts:\n");
    printf(" %-6s = %g\n", "jetDR", wjcuts->jetDR);
    printf(" QCD jet definition:\n");
    printf(" %-6s = %g\n", "photonDR", wjcuts->photonDR);
    printf(" %-6s = %g\n", "z_thr", wjcuts->z_thr);
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

    Parameters_sm *params = dynamic_cast<Parameters_sm *>(Params);
    printf(" %-6s = %g\n", "MW", params->MW);
    printf(" %-6s = %g\n", "MZ", params->MZ);
    printf(" %-6s = %g\n", "WidthW", params->WidthW);
    printf(" %-6s = %g\n", "WidthZ", params->WZ);
    printf(" %-6s = %g\n", "MH", params->MH);
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

FKS::ProcessList generateProcesses(bool QCD, bool EW, bool CutOnReal,
                                   FKS::Param *params) {

    FKS::ProcessList list;
    Parameters_sm * p = (Parameters_sm*)params;
    double Vud = p->Vud;
    double Vus = p->Vus;
    double Vub = p->Vub;
    double Vcd = p->Vcd;
    double Vcs = p->Vcs;
    double Vcb = p->Vcb;
    double Vud2 = Vud * Vud;
    double Vcd2 = Vcd * Vcd;
    double Vus2 = Vus * Vus;
    double Vcs2 = Vcs * Vcs;
    double Vub2 = Vub * Vub;
    double Vcb2 = Vcb * Vcb;
    {
        std::vector<FKS::PDGList> born = {
            {2, -1, -13, 14, 21}, {4, -1, -13, 14, 21}, {4, -5, -13, 14, 21},
            {2, -5, -13, 14, 21}, {4, -3, -13, 14, 21}, {2, -3, -13, 14, 21}};
        std::vector<double> ckm = {{Vud2, Vcd2, Vcb2, Vub2, Vcs2, Vus2}};
        FKS::ColorFlow color1 = {501, 0, 0, 0, 501};
        FKS::ColorFlow color2 = {0, 502, 0, 0, 502};
        FKS::FlavourConfig fl1(3, born, 0, color1, color2, ckm);
        if (QCD) {
            FKS::RegionList rlist1 = {FKS::Region(5, 4), FKS::Region(4, 5),
                                      FKS::Region(4, 0), FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real1 = {
                {2, -1, -13, 14, 21, 21}, {4, -1, -13, 14, 21, 21},
                {4, -5, -13, 14, 21, 21}, {2, -5, -13, 14, 21, 21},
                {4, -3, -13, 14, 21, 21}, {2, -3, -13, 14, 21, 21}};
            fl1.AddReal(32, FKS::Type_t::QCD, real1, rlist1, rlist1);

            FKS::RegionList rlist2 = {FKS::Region(5, 4)};
            FKS::RegionList rlist2all = {FKS::Region(5, 4), FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real2 = {
                {2, -1, -13, 14, -2, 2}, {4, -1, -13, 14, -4, 4},
                {4, -5, -13, 14, -4, 4}, {2, -5, -13, 14, -2, 2},
                {4, -3, -13, 14, -4, 4}, {2, -3, -13, 14, -2, 2}};
            fl1.AddReal(25, FKS::Type_t::QCD, real2, rlist2, rlist2all);

            FKS::RegionList rlist3 = {FKS::Region(5, 4)};
            FKS::RegionList rlist3all = {FKS::Region(5, 4), FKS::Region(4, 0)};
            std::vector<FKS::PDGList> real3 = {
                {2, -1, -13, 14, -1, 1}, {4, -1, -13, 14, -1, 1},
                {4, -5, -13, 14, -5, 5}, {2, -5, -13, 14, -5, 5},
                {4, -3, -13, 14, -3, 3}, {2, -3, -13, 14, -3, 3}};
            fl1.AddReal(4, FKS::Type_t::QCD, real3, rlist3, rlist3all);

            FKS::RegionList rlist4 = {FKS::Region(5, 4)};
            std::vector<FKS::PDGList> real4 = {
                {2, -1, -13, 14, -5, 5}, {4, -1, -13, 14, -5, 5},
                {4, -5, -13, 14, -3, 3}, {2, -5, -13, 14, -4, 4},
                {4, -3, -13, 14, -5, 5}, {2, -3, -13, 14, -5, 5}};
            fl1.AddReal(13, FKS::Type_t::QCD, real4, rlist4, rlist4);

            FKS::RegionList rlist5 = {FKS::Region(5, 4)};
            std::vector<FKS::PDGList> real5 = {
                {2, -1, -13, 14, -4, 4}, {4, -1, -13, 14, -3, 3},
                {4, -5, -13, 14, -2, 2}, {2, -5, -13, 14, -3, 3},
                {4, -3, -13, 14, -2, 2}, {2, -3, -13, 14, -4, 4}};
            fl1.AddReal(13, FKS::Type_t::QCD, real5, rlist5, rlist5);

            FKS::RegionList rlist6 = {FKS::Region(5, 4)};
            std::vector<FKS::PDGList> real6 = {
                {2, -1, -13, 14, -3, 3}, {4, -1, -13, 14, -2, 2},
                {4, -5, -13, 14, -1, 1}, {2, -5, -13, 14, -1, 1},
                {4, -3, -13, 14, -1, 1}, {2, -3, -13, 14, -1, 1}};
            fl1.AddReal(13, FKS::Type_t::QCD, real6, rlist6, rlist6);

            FKS::RegionList rlist7 = {FKS::Region(5, 0)};
            FKS::RegionList rlist7all = {FKS::Region(5, 0), FKS::Region(4, 0),
                                         FKS::Region(5, 4)};
            std::vector<FKS::PDGList> real7 = {
                {2, 21, -13, 14, 21, 1}, {4, 21, -13, 14, 21, 1},
                {4, 21, -13, 14, 21, 5}, {2, 21, -13, 14, 21, 5},
                {4, 21, -13, 14, 21, 3}, {2, 21, -13, 14, 21, 3}};
            fl1.AddReal(36, FKS::Type_t::QCD, real7, rlist7, rlist7all);

            FKS::RegionList rlist8 = {FKS::Region(5, -1)};
            FKS::RegionList rlist8all = {FKS::Region(5, -1), FKS::Region(4, 0),
                                         FKS::Region(5, 4)};
            std::vector<FKS::PDGList> real8 = {
                {21, -1, -13, 14, 21, -2}, {21, -1, -13, 14, 21, -4},
                {21, -5, -13, 14, 21, -4}, {21, -5, -13, 14, 21, -2},
                {21, -3, -13, 14, 21, -4}, {21, -3, -13, 14, 21, -2}};
            fl1.AddReal(11, FKS::Type_t::QCD, real8, rlist8, rlist8all);
        }
        if (EW) {
            FKS::RegionList rlist = {FKS::Region(5, 2), FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real = {{2, -1, -13, 14, 21, 22},
                                              {4, -1, -13, 14, 21, 22},
                                              {4, -5, -13, 14, 21, 22},
                                              {2, -5, -13, 14, 21, 22},
                                              {4, -3, -13, 14, 21, 22},
                                              {2, -3, -13, 14, 21, 22}};
            fl1.AddReal(40, FKS::Type_t::EW, real, rlist, rlist);
        }

        list.push_back(fl1);
    }

    {
        std::vector<FKS::PDGList> born = {
            {21, -1, -13, 14, -2}, {21, -1, -13, 14, -4}, {21, -5, -13, 14, -4},
            {21, -5, -13, 14, -2}, {21, -3, -13, 14, -4}, {21, -3, -13, 14, -2}};
        std::vector<double> ckm = {{Vud2, Vcd2, Vcb2, Vub2, Vcs2, Vus2}};
        FKS::ColorFlow color1 = {501, 0, 0, 0, 0};
        FKS::ColorFlow color2 = {502, 501, 0, 0, 502};
        FKS::FlavourConfig fl1(1, born, 0, color1, color2, ckm);

        if (QCD) {
            FKS::RegionList rlist1 = {FKS::Region(5, 0)};
            FKS::RegionList rlist1all = {FKS::Region(5, 4), FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real1 = {
                {2, -1, -13, 14, -2, 2}, {4, -1, -13, 14, -4, 4},
                {4, -5, -13, 14, -4, 4}, {2, -5, -13, 14, -2, 2},
                {4, -3, -13, 14, -4, 4}, {2, -3, -13, 14, -2, 2}};
            fl1.AddReal(25, FKS::Type_t::QCD, real1, rlist1, rlist1all);

            FKS::RegionList rlist2 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real2 = {
                {1, -1, -13, 14, -2, 1}, {1, -1, -13, 14, -4, 1},
                {5, -5, -13, 14, -4, 5}, {5, -5, -13, 14, -2, 5},
                {3, -3, -13, 14, -4, 3}, {3, -3, -13, 14, -2, 3}};
            fl1.AddReal(7, FKS::Type_t::QCD, real2, rlist2, rlist2);

            FKS::RegionList rlist3 = {FKS::Region(4, 0), FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real3 = {
                {-2, -1, -13, 14, -2, -2}, {-4, -1, -13, 14, -4, -4},
                {-4, -5, -13, 14, -4, -4}, {-2, -5, -13, 14, -2, -2},
                {-4, -3, -13, 14, -4, -4}, {-2, -3, -13, 14, -2, -2}};
            fl1.AddReal(34, FKS::Type_t::QCD, real3, rlist3, rlist3);

            FKS::RegionList rlist4 = {FKS::Region(5, -1)};
	    FKS::RegionList rlist4all = {FKS::Region(5, -1), FKS::Region(5, -2)};
            std::vector<FKS::PDGList> real4 = {
                {-1, -1, -13, 14, -2, -1}, {-1, -1, -13, 14, -4, -1},
                {-5, -5, -13, 14, -4, -5}, {-5, -5, -13, 14, -2, -5},
                {-3, -3, -13, 14, -4, -3}, {-3, -3, -13, 14, -2, -3}};
            fl1.AddReal(2, FKS::Type_t::QCD, real4, rlist4, rlist4all);

            FKS::RegionList rlist5 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real5 = {
                {3, -1, -13, 14, -2, 3}, {2, -1, -13, 14, -4, 2},
                {1, -5, -13, 14, -4, 1}, {1, -5, -13, 14, -2, 1},
                {1, -3, -13, 14, -4, 1}, {1, -3, -13, 14, -2, 1}};
            fl1.AddReal(17, FKS::Type_t::QCD, real5, rlist5, rlist5);

            FKS::RegionList rlist6 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real6 = {
                {4, -1, -13, 14, -2, 4}, {3, -1, -13, 14, -4, 3},
                {2, -5, -13, 14, -4, 2}, {3, -5, -13, 14, -2, 3},
                {2, -3, -13, 14, -4, 2}, {4, -3, -13, 14, -2, 4}};
            fl1.AddReal(17, FKS::Type_t::QCD, real6, rlist6, rlist6);

            FKS::RegionList rlist7 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real7 = {
                {5, -1, -13, 14, -2, 5}, {5, -1, -13, 14, -4, 5},
                {3, -5, -13, 14, -4, 3}, {4, -5, -13, 14, -2, 4},
                {5, -3, -13, 14, -4, 5}, {5, -3, -13, 14, -2, 5}};
            fl1.AddReal(17, FKS::Type_t::QCD, real7, rlist7, rlist7);

            FKS::RegionList rlist8 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real8 = {
                {-5, -1, -13, 14, -2, -5}, {-5, -1, -13, 14, -4, -5},
                {-3, -5, -13, 14, -4, -3}, {-4, -5, -13, 14, -2, -4},
                {-5, -3, -13, 14, -4, -5}, {-5, -3, -13, 14, -2, -5}};
            fl1.AddReal(12, FKS::Type_t::QCD, real8, rlist8, rlist8);

            FKS::RegionList rlist9 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real9 = {
                {-4, -1, -13, 14, -2, -4}, {-3, -1, -13, 14, -4, -3},
                {-2, -5, -13, 14, -4, -2}, {-3, -5, -13, 14, -2, -3},
                {-2, -3, -13, 14, -4, -2}, {-4, -3, -13, 14, -2, -4}};
            fl1.AddReal(12, FKS::Type_t::QCD, real9, rlist9, rlist9);

            FKS::RegionList rlist10 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real10 = {
                {-3, -1, -13, 14, -2, -3}, {-2, -1, -13, 14, -4, -2},
                {-1, -5, -13, 14, -4, -1}, {-1, -5, -13, 14, -2, -1},
                {-1, -3, -13, 14, -4, -1}, {-1, -3, -13, 14, -2, -1}};
            fl1.AddReal(12, FKS::Type_t::QCD, real10, rlist10, rlist10);

            FKS::RegionList rlist11 = {FKS::Region(5, 4), FKS::Region(5, 0)};
            FKS::RegionList rlist11all = {FKS::Region(5, 0), FKS::Region(4, -1),
                                          FKS::Region(5, 4)};
            std::vector<FKS::PDGList> real11 = {
                {21, -1, -13, 14, -2, 21}, {21, -1, -13, 14, -4, 21},
                {21, -5, -13, 14, -4, 21}, {21, -5, -13, 14, -2, 21},
                {21, -3, -13, 14, -4, 21}, {21, -3, -13, 14, -2, 21}};
            fl1.AddReal(30, FKS::Type_t::QCD, real11, rlist11, rlist11all);

            FKS::RegionList rlist12 = {FKS::Region(5, -2)};
	    FKS::RegionList rlist12all = {FKS::Region(5, -1), FKS::Region(5, -2), 
	      FKS::Region(4,-1), FKS::Region(4,-2)};
            std::vector<FKS::PDGList> real12 = {
                {21, 21, -13, 14, -2, 1}, {21, 21, -13, 14, -4, 1},
                {21, 21, -13, 14, -4, 5}, {21, 21, -13, 14, -2, 5},
                {21, 21, -13, 14, -4, 3}, {21, 21, -13, 14, -2, 3}};
            fl1.AddReal(35, FKS::Type_t::QCD, real12, rlist12, rlist12all);
        }
        if (EW) {
            FKS::RegionList rlist = {FKS::Region(5, 4), FKS::Region(5, 2),
                                     FKS::Region(5, -2)};
            FKS::RegionList rlistall;
            if(CutOnReal) {
                rlistall = rlist;
            } else {
                rlistall = {FKS::Region(5, 4),  FKS::Region(5, 2),
                            FKS::Region(5, -2), FKS::Region(4, -1)};
            }

            std::vector<FKS::PDGList> real = {{21, -1, -13, 14, -2, 22},
                                              {21, -1, -13, 14, -4, 22},
                                              {21, -5, -13, 14, -4, 22},
                                              {21, -5, -13, 14, -2, 22},
                                              {21, -3, -13, 14, -4, 22},
                                              {21, -3, -13, 14, -2, 22}};
            fl1.AddReal(44, FKS::Type_t::EW, real, rlist, rlistall);
        }

        list.push_back(fl1);
    }
    {
        std::vector<FKS::PDGList> born = {
            {2, 21, -13, 14, 1}, {4, 21, -13, 14, 1}, {4, 21, -13, 14, 5},
            {2, 21, -13, 14, 5}, {4, 21, -13, 14, 3}, {2, 21, -13, 14, 3}};
        FKS::ColorFlow color1 = {501, 502, 0, 0, 502};
        FKS::ColorFlow color2 = {0, 501, 0, 0, 0};
        std::vector<double> ckm = {{Vud2, Vcd2, Vcb2, Vub2, Vcs2, Vus2}};
        FKS::FlavourConfig fl1(2, born, 0, color1, color2, ckm);

        if (QCD) {
            FKS::RegionList rlist1 = {FKS::Region(4, 0), FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real1 = {
                {2, 1, -13, 14, 1, 1}, {4, 1, -13, 14, 1, 1},
                {4, 5, -13, 14, 5, 5}, {2, 5, -13, 14, 5, 5},
                {4, 3, -13, 14, 3, 3}, {2, 3, -13, 14, 3, 3}};
            fl1.AddReal(33, FKS::Type_t::QCD, real1, rlist1, rlist1);

            FKS::RegionList rlist2 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real2 = {
                {2, -2, -13, 14, 1, -2}, {4, -4, -13, 14, 1, -4},
                {4, -4, -13, 14, 5, -4}, {2, -2, -13, 14, 5, -2},
                {4, -4, -13, 14, 3, -4}, {2, -2, -13, 14, 3, -2}};
            fl1.AddReal(10, FKS::Type_t::QCD, real2, rlist2, rlist2);

            FKS::RegionList rlist3 = {FKS::Region(5, 0)};
            FKS::RegionList rlist3all = {FKS::Region(5, 4), FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real3 = {
                {2, -1, -13, 14, 1, -1}, {4, -1, -13, 14, 1, -1},
                {4, -5, -13, 14, 5, -5}, {2, -5, -13, 14, 5, -5},
                {4, -3, -13, 14, 3, -3}, {2, -3, -13, 14, 3, -3}};
            fl1.AddReal(20, FKS::Type_t::QCD, real3, rlist3, rlist3all);

            FKS::RegionList rlist4 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real4 = {
                {2, 3, -13, 14, 1, 3}, {4, 2, -13, 14, 1, 2},
                {4, 1, -13, 14, 5, 1}, {2, 1, -13, 14, 5, 1},
                {4, 1, -13, 14, 3, 1}, {2, 1, -13, 14, 3, 1}};
            fl1.AddReal(27, FKS::Type_t::QCD, real4, rlist4, rlist4);

            FKS::RegionList rlist5 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real5 = {
                {2, 4, -13, 14, 1, 4}, {4, 3, -13, 14, 1, 3},
                {4, 2, -13, 14, 5, 2}, {2, 3, -13, 14, 5, 3},
                {4, 2, -13, 14, 3, 2}, {2, 4, -13, 14, 3, 4}};
            fl1.AddReal(27, FKS::Type_t::QCD, real5, rlist5, rlist5);

            FKS::RegionList rlist6 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real6 = {
                {2, 5, -13, 14, 1, 5}, {4, 5, -13, 14, 1, 5},
                {4, 3, -13, 14, 5, 3}, {2, 4, -13, 14, 5, 4},
                {4, 5, -13, 14, 3, 5}, {2, 5, -13, 14, 3, 5}};
            fl1.AddReal(27, FKS::Type_t::QCD, real6, rlist6, rlist6);

            FKS::RegionList rlist7 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real7 = {
                {2, -5, -13, 14, 1, -5}, {4, -5, -13, 14, 1, -5},
                {4, -3, -13, 14, 5, -3}, {2, -4, -13, 14, 5, -4},
                {4, -5, -13, 14, 3, -5}, {2, -5, -13, 14, 3, -5}};
            fl1.AddReal(38, FKS::Type_t::QCD, real7, rlist7, rlist7);

            FKS::RegionList rlist8 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real8 = {
                {2, -4, -13, 14, 1, -4}, {4, -3, -13, 14, 1, -3},
                {4, -2, -13, 14, 5, -2}, {2, -3, -13, 14, 5, -3},
                {4, -2, -13, 14, 3, -2}, {2, -4, -13, 14, 3, -4}};
            fl1.AddReal(38, FKS::Type_t::QCD, real8, rlist8, rlist8);

            FKS::RegionList rlist9 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real9 = {
                {2, -3, -13, 14, 1, -3}, {4, -2, -13, 14, 1, -2},
                {4, -1, -13, 14, 5, -1}, {2, -1, -13, 14, 5, -1},
                {4, -1, -13, 14, 3, -1}, {2, -1, -13, 14, 3, -1}};
            fl1.AddReal(38, FKS::Type_t::QCD, real9, rlist9, rlist9);
            

            FKS::RegionList rlist10 = {FKS::Region(5,4), FKS::Region(5, 0)};
            FKS::RegionList rlist10all = {FKS::Region(5, 4), FKS::Region(5, 0),
                                          FKS::Region(4, 0)};
            std::vector<FKS::PDGList> real10 = {
                {2, 21, -13, 14, 1, 21}, {4, 21, -13, 14, 1, 21},
                {4, 21, -13, 14, 5, 21}, {2, 21, -13, 14, 5, 21},
                {4, 21, -13, 14, 3, 21}, {2, 21, -13, 14, 3, 21}};
            fl1.AddReal(31, FKS::Type_t::QCD, real10, rlist10, rlist10all);

            
            FKS::RegionList rlist11 = {FKS::Region(5, -1)};
	    FKS::RegionList rlist11all = {FKS::Region(5, -1), FKS::Region(5, -2), 
	      FKS::Region(4,-1), FKS::Region(4,-2)};
            std::vector<FKS::PDGList> real11 = {
                {21, 21, -13, 14, 1, -2}, {21, 21, -13, 14, 1, -4},
                {21, 21, -13, 14, 5, -4}, {21, 21, -13, 14, 5, -2},
                {21, 21, -13, 14, 3, -4}, {21, 21, -13, 14, 3, -2}};
            fl1.AddReal(26, FKS::Type_t::QCD, real11, rlist11, rlist11all);

            FKS::RegionList rlist12 = {FKS::Region(5, -2)};
	    FKS::RegionList rlist12all = {FKS::Region(5, -2), FKS::Region(5, -1)};
            std::vector<FKS::PDGList> real12 = {
                {2, 2, -13, 14, 1, 2}, {4, 4, -13, 14, 1, 4},
                {4, 4, -13, 14, 5, 4}, {2, 2, -13, 14, 5, 2},
                {4, 4, -13, 14, 3, 4}, {2, 2, -13, 14, 3, 2}};
            fl1.AddReal(18, FKS::Type_t::QCD, real12, rlist12, rlist12all);
            
        }
        if (EW) {
            FKS::RegionList rlist = {FKS::Region(5, 4), FKS::Region(5, 2),
                                     FKS::Region(5, -1)};
            FKS::RegionList rlistall;
            if(CutOnReal) {
                rlistall = rlist;
            } else {
                rlistall = {FKS::Region(5, 4),  FKS::Region(5, 2),
                            FKS::Region(5, -1), FKS::Region(4, -2)};
            }
            std::vector<FKS::PDGList> real = {{2, 21, -13, 14, 1, 22},
                                              {4, 21, -13, 14, 1, 22},
                                              {4, 21, -13, 14, 5, 22},
                                              {2, 21, -13, 14, 5, 22},
                                              {4, 21, -13, 14, 3, 22},
                                              {2, 21, -13, 14, 3, 22}};
            fl1.AddReal(43, FKS::Type_t::EW, real, rlist, rlistall);
        }

        list.push_back(fl1);
    }
    {
        std::vector<FKS::PDGList> born = {
            {-1, 2, -13, 14, 21}, {-1, 4, -13, 14, 21}, {-5, 4, -13, 14, 21},
            {-5, 2, -13, 14, 21}, {-3, 4, -13, 14, 21}, {-3, 2, -13, 14, 21}};
        FKS::ColorFlow color1 = {0, 502, 0, 0, 502};
        FKS::ColorFlow color2 = {501, 0, 0, 0, 501};
        std::vector<double> ckm = {{Vud2, Vcd2, Vcb2, Vub2, Vcs2, Vus2}};
        FKS::FlavourConfig fl1(5, born, 0, color1, color2, ckm);

        if (QCD) {
            FKS::RegionList rlist1 = {FKS::Region(5, 4), FKS::Region(4, 5),
                                      FKS::Region(4, 0), FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real1 = {
                {-1, 2, -13, 14, 21, 21}, {-1, 4, -13, 14, 21, 21},
                {-5, 4, -13, 14, 21, 21}, {-5, 2, -13, 14, 21, 21},
                {-3, 4, -13, 14, 21, 21}, {-3, 2, -13, 14, 21, 21}};
            fl1.AddReal(21, FKS::Type_t::QCD, real1, rlist1, rlist1);

            FKS::RegionList rlist2 = {FKS::Region(5, 4)};
	    FKS::RegionList rlist2all = {FKS::Region(5,0), FKS::Region(5, 4)};
            std::vector<FKS::PDGList> real2 = {
                {-1, 2, -13, 14, -2, 2}, {-1, 4, -13, 14, -4, 4},
                {-5, 4, -13, 14, -4, 4}, {-5, 2, -13, 14, -2, 2},
                {-3, 4, -13, 14, -4, 4}, {-3, 2, -13, 14, -2, 2}};
            fl1.AddReal(6, FKS::Type_t::QCD, real2, rlist2, rlist2all);

            FKS::RegionList rlist3 = {FKS::Region(5, 4)};
            FKS::RegionList rlist3all = {FKS::Region(4, 0), FKS::Region(5, 4)};
            std::vector<FKS::PDGList> real3 = {
                {-1, 2, -13, 14, -1, 1}, {-1, 4, -13, 14, -1, 1},
                {-5, 4, -13, 14, -5, 5}, {-5, 2, -13, 14, -5, 5},
                {-3, 4, -13, 14, -3, 3}, {-3, 2, -13, 14, -3, 3}};
            fl1.AddReal(5, FKS::Type_t::QCD, real3, rlist3, rlist3all);

            FKS::RegionList rlist4 = {FKS::Region(5, 4)};
            std::vector<FKS::PDGList> real4 = {
                {-1, 2, -13, 14, -5, 5}, {-1, 4, -13, 14, -5, 5},
                {-5, 4, -13, 14, -3, 3}, {-5, 2, -13, 14, -4, 4},
                {-3, 4, -13, 14, -5, 5}, {-3, 2, -13, 14, -5, 5}};
            fl1.AddReal(37, FKS::Type_t::QCD, real4, rlist4, rlist4);

            FKS::RegionList rlist5 = {FKS::Region(5, 4)};
            std::vector<FKS::PDGList> real5 = {
                {-1, 2, -13, 14, -4, 4}, {-1, 4, -13, 14, -3, 3},
                {-5, 4, -13, 14, -2, 2}, {-5, 2, -13, 14, -3, 3},
                {-3, 4, -13, 14, -2, 2}, {-3, 2, -13, 14, -4, 4}};
            fl1.AddReal(37, FKS::Type_t::QCD, real5, rlist5, rlist5);

            FKS::RegionList rlist6 = {FKS::Region(5, 4)};
            std::vector<FKS::PDGList> real6 = {
                {-1, 2, -13, 14, -3, 3}, {-1, 4, -13, 14, -2, 2},
                {-5, 4, -13, 14, -1, 1}, {-5, 2, -13, 14, -1, 1},
                {-3, 4, -13, 14, -1, 1}, {-3, 2, -13, 14, -1, 1}};
            fl1.AddReal(37, FKS::Type_t::QCD, real6, rlist6, rlist6);

            FKS::RegionList rlist7 = {FKS::Region(5, 0)};
            FKS::RegionList rlist7all = {FKS::Region(5, 0), FKS::Region(4, 0),
                                         FKS::Region(5, 4)};
            std::vector<FKS::PDGList> real7 = {
                {21, 2, -13, 14, 21, 1}, {21, 4, -13, 14, 21, 1},
                {21, 4, -13, 14, 21, 5}, {21, 2, -13, 14, 21, 5},
                {21, 4, -13, 14, 21, 3}, {21, 2, -13, 14, 21, 3}};
            fl1.AddReal(39, FKS::Type_t::QCD, real7, rlist7, rlist7all);

            FKS::RegionList rlist8 = {FKS::Region(5, 0)};
            FKS::RegionList rlist8all = {FKS::Region(5, 4), FKS::Region(5, 0),
                                         FKS::Region(4, 0)};
            std::vector<FKS::PDGList> real8 = {
                {-1, 21, -13, 14, 21, -2}, {-1, 21, -13, 14, 21, -4},
                {-5, 21, -13, 14, 21, -4}, {-5, 21, -13, 14, 21, -2},
                {-3, 21, -13, 14, 21, -4}, {-3, 21, -13, 14, 21, -2}};
            fl1.AddReal(14, FKS::Type_t::QCD, real8, rlist8, rlist8all);
        }
        if (EW) {
            FKS::RegionList rlist = {FKS::Region(5, 2), FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real = {{-1, 2, -13, 14, 21, 22},
                                              {-1, 4, -13, 14, 21, 22},
                                              {-5, 4, -13, 14, 21, 22},
                                              {-5, 2, -13, 14, 21, 22},
                                              {-3, 4, -13, 14, 21, 22},
                                              {-3, 2, -13, 14, 21, 22}};
            fl1.AddReal(41, FKS::Type_t::EW, real, rlist, rlist);
        }

        list.push_back(fl1);
    }

    {
        std::vector<FKS::PDGList> born = {
            {21, 2, -13, 14, 1}, {21, 4, -13, 14, 1}, {21, 4, -13, 14, 5},
            {21, 2, -13, 14, 5}, {21, 4, -13, 14, 3}, {21, 2, -13, 14, 3}};
        FKS::ColorFlow color1 = {501, 502, 0, 0, 501};
        FKS::ColorFlow color2 = {502, 0, 0, 0, 0};
        std::vector<double> ckm = {{Vud2, Vcd2, Vcb2, Vub2, Vcs2, Vus2}};
        FKS::FlavourConfig fl1(0, born, 0, color1, color2, ckm);

        if (QCD) {
            FKS::RegionList rlist1 = {FKS::Region(5, -1)};
	    FKS::RegionList rlist1all = {FKS::Region(5, -1), FKS::Region(5, -2)};
            std::vector<FKS::PDGList> real1 = {
                {2, 2, -13, 14, 1, 2}, {4, 4, -13, 14, 1, 4},
                {4, 4, -13, 14, 5, 4}, {2, 2, -13, 14, 5, 2},
                {4, 4, -13, 14, 3, 4}, {2, 2, -13, 14, 3, 2}};
            fl1.AddReal(18, FKS::Type_t::QCD, real1, rlist1, rlist1all);

            FKS::RegionList rlist2 = {FKS::Region(4, 0), FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real2 = {
                {1, 2, -13, 14, 1, 1}, {1, 4, -13, 14, 1, 1},
                {5, 4, -13, 14, 5, 5}, {5, 2, -13, 14, 5, 5},
                {3, 4, -13, 14, 3, 3}, {3, 2, -13, 14, 3, 3}};
            fl1.AddReal(28, FKS::Type_t::QCD, real2, rlist2, rlist2);

            FKS::RegionList rlist3 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real3 = {
                {-2, 2, -13, 14, 1, -2}, {-4, 4, -13, 14, 1, -4},
                {-4, 4, -13, 14, 5, -4}, {-2, 2, -13, 14, 5, -2},
                {-4, 4, -13, 14, 3, -4}, {-2, 2, -13, 14, 3, -2}};
            fl1.AddReal(1, FKS::Type_t::QCD, real3, rlist3, rlist3);

            FKS::RegionList rlist4 = {FKS::Region(5, 0)};
            FKS::RegionList rlist4all = {FKS::Region(5, 0), FKS::Region(5, 4)};
            std::vector<FKS::PDGList> real4 = {
                {-1, 2, -13, 14, 1, -1}, {-1, 4, -13, 14, 1, -1},
                {-5, 4, -13, 14, 5, -5}, {-5, 2, -13, 14, 5, -5},
                {-3, 4, -13, 14, 3, -3}, {-3, 2, -13, 14, 3, -3}};
            fl1.AddReal(15, FKS::Type_t::QCD, real4, rlist4, rlist4all);

            FKS::RegionList rlist5 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real5 = {
                {3, 2, -13, 14, 1, 3}, {2, 4, -13, 14, 1, 2},
                {1, 4, -13, 14, 5, 1}, {1, 2, -13, 14, 5, 1},
                {1, 4, -13, 14, 3, 1}, {1, 2, -13, 14, 3, 1}};
            fl1.AddReal(3, FKS::Type_t::QCD, real5, rlist5, rlist5);

            FKS::RegionList rlist6 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real6 = {
                {4, 2, -13, 14, 1, 4}, {3, 4, -13, 14, 1, 3},
                {2, 4, -13, 14, 5, 2}, {3, 2, -13, 14, 5, 3},
                {2, 4, -13, 14, 3, 2}, {4, 2, -13, 14, 3, 4}};
            fl1.AddReal(3, FKS::Type_t::QCD, real6, rlist6, rlist6);

            FKS::RegionList rlist7 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real7 = {
                {5, 2, -13, 14, 1, 5}, {5, 4, -13, 14, 1, 5},
                {3, 4, -13, 14, 5, 3}, {4, 2, -13, 14, 5, 4},
                {5, 4, -13, 14, 3, 5}, {5, 2, -13, 14, 3, 5}};
            fl1.AddReal(3, FKS::Type_t::QCD, real7, rlist7, rlist7);

            FKS::RegionList rlist8 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real8 = {
                {-5, 2, -13, 14, 1, -5}, {-5, 4, -13, 14, 1, -5},
                {-3, 4, -13, 14, 5, -3}, {-4, 2, -13, 14, 5, -4},
                {-5, 4, -13, 14, 3, -5}, {-5, 2, -13, 14, 3, -5}};
            fl1.AddReal(19, FKS::Type_t::QCD, real8, rlist8, rlist8);

            FKS::RegionList rlist9 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real9 = {
                {-4, 2, -13, 14, 1, -4}, {-3, 4, -13, 14, 1, -3},
                {-2, 4, -13, 14, 5, -2}, {-3, 2, -13, 14, 5, -3},
                {-2, 4, -13, 14, 3, -2}, {-4, 2, -13, 14, 3, -4}};
            fl1.AddReal(19, FKS::Type_t::QCD, real9, rlist9, rlist9);

            FKS::RegionList rlist10 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real10 = {
                {-3, 2, -13, 14, 1, -3}, {-2, 4, -13, 14, 1, -2},
                {-1, 4, -13, 14, 5, -1}, {-1, 2, -13, 14, 5, -1},
                {-1, 4, -13, 14, 3, -1}, {-1, 2, -13, 14, 3, -1}};
            fl1.AddReal(19, FKS::Type_t::QCD, real10, rlist10, rlist10);

            FKS::RegionList rlist11 = {FKS::Region(5, 4), FKS::Region(5, 0)};
            FKS::RegionList rlist11all = {FKS::Region(5, 0), FKS::Region(4, 0),
                                          FKS::Region(5, 4)};
            std::vector<FKS::PDGList> real11 = {
                {21, 2, -13, 14, 1, 21}, {21, 4, -13, 14, 1, 21},
                {21, 4, -13, 14, 5, 21}, {21, 2, -13, 14, 5, 21},
                {21, 4, -13, 14, 3, 21}, {21, 2, -13, 14, 3, 21}};
            fl1.AddReal(23, FKS::Type_t::QCD, real11, rlist11, rlist11all);

            FKS::RegionList rlist12 = {FKS::Region(5, -2)};
	    FKS::RegionList rlist12all = {FKS::Region(5, -1), FKS::Region(5, -2), 
	      FKS::Region(4,-1), FKS::Region(4,-2)};
            std::vector<FKS::PDGList> real12 = {
                {21, 21, -13, 14, 1, -2}, {21, 21, -13, 14, 1, -4},
                {21, 21, -13, 14, 5, -4}, {21, 21, -13, 14, 5, -2},
                {21, 21, -13, 14, 3, -4}, {21, 21, -13, 14, 3, -2}};
            fl1.AddReal(26, FKS::Type_t::QCD, real12, rlist12, rlist12all);
        }
        if (EW) {
            FKS::RegionList rlist = {FKS::Region(5, 4), FKS::Region(5, 2),
                                     FKS::Region(5, -2)};
            FKS::RegionList rlistall;
            if(CutOnReal) {
                rlistall = rlist;
            } else {
                rlistall = {FKS::Region(5, 4),  FKS::Region(5, 2),
                            FKS::Region(5, -2), FKS::Region(4, -1)};
            }
            std::vector<FKS::PDGList> real = {{21, 2, -13, 14, 1, 22},
                                              {21, 4, -13, 14, 1, 22},
                                              {21, 4, -13, 14, 5, 22},
                                              {21, 2, -13, 14, 5, 22},
                                              {21, 4, -13, 14, 3, 22},
                                              {21, 2, -13, 14, 3, 22}};
            fl1.AddReal(42, FKS::Type_t::EW, real, rlist, rlistall);
        }

        list.push_back(fl1);
    }
    {
        std::vector<FKS::PDGList> born = {
            {-1, 21, -13, 14, -2}, {-1, 21, -13, 14, -4}, {-5, 21, -13, 14, -4},
            {-5, 21, -13, 14, -2}, {-3, 21, -13, 14, -4}, {-3, 21, -13, 14, -2}};
        FKS::ColorFlow color1 = {0, 501, 0, 0, 0};
        FKS::ColorFlow color2 = {501, 502, 0, 0, 502};
        std::vector<double> ckm = {{Vud2, Vcd2, Vcb2, Vub2, Vcs2, Vus2}};
        FKS::FlavourConfig fl1(4, born, 0, color1, color2, ckm);

        if (QCD) {
            FKS::RegionList rlist1 = {FKS::Region(5, 0)};
	    FKS::RegionList rlist1all = {FKS::Region(5,0), FKS::Region(5, 4)};
            std::vector<FKS::PDGList> real1 = {
                {-1, 2, -13, 14, -2, 2}, {-1, 4, -13, 14, -4, 4},
                {-5, 4, -13, 14, -4, 4}, {-5, 2, -13, 14, -2, 2},
                {-3, 4, -13, 14, -4, 4}, {-3, 2, -13, 14, -2, 2}};
            fl1.AddReal(6, FKS::Type_t::QCD, real1, rlist1, rlist1all);

            FKS::RegionList rlist2 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real2 = {
                {-1, 1, -13, 14, -2, 1}, {-1, 1, -13, 14, -4, 1},
                {-5, 5, -13, 14, -4, 5}, {-5, 5, -13, 14, -2, 5},
                {-3, 3, -13, 14, -4, 3}, {-3, 3, -13, 14, -2, 3}};
            fl1.AddReal(24, FKS::Type_t::QCD, real2, rlist2, rlist2);

            FKS::RegionList rlist3 = {FKS::Region(4, 0), FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real3 = {
                {-1, -2, -13, 14, -2, -2}, {-1, -4, -13, 14, -4, -4},
                {-5, -4, -13, 14, -4, -4}, {-5, -2, -13, 14, -2, -2},
                {-3, -4, -13, 14, -4, -4}, {-3, -2, -13, 14, -2, -2}};
            fl1.AddReal(9, FKS::Type_t::QCD, real3, rlist3, rlist3);

            FKS::RegionList rlist4 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real4 = {
                {-1, 3, -13, 14, -2, 3}, {-1, 2, -13, 14, -4, 2},
                {-5, 1, -13, 14, -4, 1}, {-5, 1, -13, 14, -2, 1},
                {-3, 1, -13, 14, -4, 1}, {-3, 1, -13, 14, -2, 1}};
            fl1.AddReal(16, FKS::Type_t::QCD, real4, rlist4, rlist4);

            FKS::RegionList rlist5 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real5 = {
                {-1, 4, -13, 14, -2, 4}, {-1, 3, -13, 14, -4, 3},
                {-5, 2, -13, 14, -4, 2}, {-5, 3, -13, 14, -2, 3},
                {-3, 2, -13, 14, -4, 2}, {-3, 4, -13, 14, -2, 4}};
            fl1.AddReal(16, FKS::Type_t::QCD, real5, rlist5, rlist5);

            FKS::RegionList rlist6 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real6 = {
                {-1, 5, -13, 14, -2, 5}, {-1, 5, -13, 14, -4, 5},
                {-5, 3, -13, 14, -4, 3}, {-5, 4, -13, 14, -2, 4},
                {-3, 5, -13, 14, -4, 5}, {-3, 5, -13, 14, -2, 5}};
            fl1.AddReal(16, FKS::Type_t::QCD, real6, rlist6, rlist6);

            FKS::RegionList rlist7 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real7 = {
                {-1, -5, -13, 14, -2, -5}, {-1, -5, -13, 14, -4, -5},
                {-5, -3, -13, 14, -4, -3}, {-5, -4, -13, 14, -2, -4},
                {-3, -5, -13, 14, -4, -5}, {-3, -5, -13, 14, -2, -5}};
            fl1.AddReal(8, FKS::Type_t::QCD, real7, rlist7, rlist7);

            FKS::RegionList rlist8 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real8 = {
                {-1, -4, -13, 14, -2, -4}, {-1, -3, -13, 14, -4, -3},
                {-5, -2, -13, 14, -4, -2}, {-5, -3, -13, 14, -2, -3},
                {-3, -2, -13, 14, -4, -2}, {-3, -4, -13, 14, -2, -4}};
            fl1.AddReal(8, FKS::Type_t::QCD, real8, rlist8, rlist8);

            FKS::RegionList rlist9 = {FKS::Region(5, 0)};
            std::vector<FKS::PDGList> real9 = {
                {-1, -3, -13, 14, -2, -3}, {-1, -2, -13, 14, -4, -2},
                {-5, -1, -13, 14, -4, -1}, {-5, -1, -13, 14, -2, -1},
                {-3, -1, -13, 14, -4, -1}, {-3, -1, -13, 14, -2, -1}};
            fl1.AddReal(8, FKS::Type_t::QCD, real9, rlist9, rlist9);

            FKS::RegionList rlist10 = {FKS::Region(5, 4), FKS::Region(5, 0)};
            FKS::RegionList rlist10all = {FKS::Region(5, 4), FKS::Region(5, 0),
                                          FKS::Region(4, 0)};
            std::vector<FKS::PDGList> real10 = {
                {-1, 21, -13, 14, -2, 21}, {-1, 21, -13, 14, -4, 21},
                {-5, 21, -13, 14, -4, 21}, {-5, 21, -13, 14, -2, 21},
                {-3, 21, -13, 14, -4, 21}, {-3, 21, -13, 14, -2, 21}};
            fl1.AddReal(22, FKS::Type_t::QCD, real10, rlist10, rlist10all);

            FKS::RegionList rlist11 = {FKS::Region(5, -1)};
	    FKS::RegionList rlist11all = {FKS::Region(5, -1), FKS::Region(5, -2), 
	      FKS::Region(4,-1), FKS::Region(4,-2)};
            std::vector<FKS::PDGList> real11 = {
                {21, 21, -13, 14, -2, 1}, {21, 21, -13, 14, -4, 1},
                {21, 21, -13, 14, -4, 5}, {21, 21, -13, 14, -2, 5},
                {21, 21, -13, 14, -4, 3}, {21, 21, -13, 14, -2, 3}};
            fl1.AddReal(35, FKS::Type_t::QCD, real11, rlist11, rlist11all);

            FKS::RegionList rlist12 = {FKS::Region(5, -2)};
	    FKS::RegionList rlist12all = {FKS::Region(5, -1), FKS::Region(5, -2)};
            std::vector<FKS::PDGList> real12 = {
                {-1, -1, -13, 14, -2, -1}, {-1, -1, -13, 14, -4, -1},
                {-5, -5, -13, 14, -4, -5}, {-5, -5, -13, 14, -2, -5},
                {-3, -3, -13, 14, -4, -3}, {-3, -3, -13, 14, -2, -3}};
            fl1.AddReal(2, FKS::Type_t::QCD, real12, rlist12, rlist12all);
        }
        if (EW) {
            FKS::RegionList rlist = {FKS::Region(5, 4), FKS::Region(5, 2),
                                     FKS::Region(5, -1)};
            FKS::RegionList rlistall;
            if(CutOnReal) {
                rlistall = rlist;
            } else {
                rlistall = {FKS::Region(5, 4),  FKS::Region(5, 2),
                            FKS::Region(5, -1), FKS::Region(4, -2)};
            }
            std::vector<FKS::PDGList> real = {{-1, 21, -13, 14, -2, 22},
                                              {-1, 21, -13, 14, -4, 22},
                                              {-5, 21, -13, 14, -4, 22},
                                              {-5, 21, -13, 14, -2, 22},
                                              {-3, 21, -13, 14, -4, 22},
                                              {-3, 21, -13, 14, -2, 22}};
            fl1.AddReal(45, FKS::Type_t::EW, real, rlist, rlistall);
        }

        list.push_back(fl1);
    }

    return list;
}
