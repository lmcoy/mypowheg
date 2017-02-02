#include "zdata.h"

#include "config/file.h"

#include "cuts.h"
#include "matrixelement.h"
#include "me/parameters_sm.h"
#include "myhistograms.h"
#include "scales.h"

using namespace UserProcess;

static FKS::ProcessList generateProcesses(bool QCD, bool EW);

ZData::~ZData() {
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

int ZData::ProcessInit(Config::File &cfile) {
    auto icuts = std::shared_ptr<DrellYanCuts>(new DrellYanCuts());
    if (cfile.GetDoubleInterval("cuts", "mll", &icuts->mllmin,
                                &icuts->mllmax) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }

    if (cfile.GetDoubleInterval("cuts", "pT", &icuts->pTmin, &icuts->pTmax) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    if (cfile.GetDoubleInterval("cuts", "eta", &icuts->EtaMin,
                                &icuts->EtaMax) !=
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

    auto me = std::make_shared<DrellYanME>();
    me->Init(Params, AlphaScheme::Gmu, .0);

    MatrixElement = me;

    hists = new MyHistograms(3);
    hists->Init1D(0, hist_mll[0], hist_mll[1], hist_mll[2]);
    hists->Init1D(1, hist_pt[0], hist_pt[1], hist_pt[2]);
    hists->Init1D(2, hist_y[0], hist_y[1], hist_y[2]);

    Process = generateProcesses(RadiationType.QCD, RadiationType.EW);

    // add Z resonance
    Resonance.pdg = 23;
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

void ZData::ProcessPrint() const {
    auto dycuts = std::static_pointer_cast<DrellYanCuts>(cuts);
    printf("/******************************************************************"
           "***\n"
           " *                               cuts                              "
           "  *\n"
           " ******************************************************************"
           "***/\n");
    printf(" %-6s = (%g, %g)\n", "mll", dycuts->mllmin, dycuts->mllmax);
    printf(" %-6s = (%g, %g)\n", "pT", dycuts->pTmin, dycuts->pTmax);
    printf(" %-6s = (%g, %g)\n", "eta", dycuts->EtaMin, dycuts->EtaMax);
    printf(" %-6s = (%g, %g)\n", "y", dycuts->Ymin, dycuts->Ymax);
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
    using R = DrellYanME::SubProcesses;
    FKS::ProcessList list;
    FKS::ColorFlow color1 = {0, 501, 0, 0};
    FKS::ColorFlow color2 = {501, 0, 0, 0};
    FKS::FlavourConfig fl1(0, {{-2, 2, -13, 13}}, {{-2, 2, -4, 4}}, color1,
                           color2);
    if (EW) {
        FKS::ResonanceList resonances;
        size_t isr_res = resonances.Add(
            FKS::Resonance({{2, 3}}, 91.1876, 2.4952, 2.0 / 3.0));
        size_t fsr_res =
            resonances.Add(FKS::Resonance({{2, 3, 4}}, 91.1876, 2.4952, 1.0));
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 2, fsr_res));
        rlist.push_back(FKS::Region(4, 3, fsr_res));
        rlist.push_back(FKS::Region(4, 0, isr_res));
        fl1.AddRealDY(R::UXU_MUXMU_A, FKS::Type_t::EW, {{-2, 2, -13, 13, 22}},
                      {{-2, 2, -4, 4}}, rlist, resonances);
    }
    if (QCD) {
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 0));
        fl1.AddRealDY(R::UXU_MUXMU_G, FKS::Type_t::QCD, {{-2, 2, -13, 13, 21}},
                      {{-2, 2, -4, 4}}, rlist);
        fl1.AddRealDY(R::UXG_MUXMU_UX, FKS::Type_t::QCD,
                      {{-2, 21, -13, 13, -2}}, {{-2, 0, -4, 0}}, rlist);
        fl1.AddRealDY(R::GU_MUXMU_U, FKS::Type_t::QCD, {{21, 2, -13, 13, 2}},
                      {{0, 2, 0, 4}}, rlist);
    }
    list.push_back(fl1);

    color1 = {501, 0, 0, 0};
    color2 = {0, 501, 0, 0};
    FKS::FlavourConfig fl2(1, {{2, -2, -13, 13}}, {{2, -2, 4, -4}}, color1,
                           color2);
    if (EW) {
        FKS::ResonanceList resonances;
        size_t isr_res = resonances.Add(
            FKS::Resonance({{2, 3}}, 91.1876, 2.4952, 2.0 / 3.0));
        size_t fsr_res =
            resonances.Add(FKS::Resonance({{2, 3, 4}}, 91.1876, 2.4952, 1.0));
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 2, fsr_res));
        rlist.push_back(FKS::Region(4, 3, fsr_res));
        rlist.push_back(FKS::Region(4, 0, isr_res));
        fl2.AddRealDY(R::UUX_MUXMU_A, FKS::Type_t::EW, {{2, -2, -13, 13, 22}},
                      {{2, -2, 4, -4}}, rlist, resonances);
    }
    if (QCD) {
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 0));
        fl2.AddRealDY(R::UUX_MUXMU_G, FKS::Type_t::QCD, {{2, -2, -13, 13, 21}},
                      {{2, -2, 4, -4}}, rlist);
        fl2.AddRealDY(R::UG_MUXMU_U, FKS::Type_t::QCD, {{2, 21, -13, 13, 2}},
                      {{2, 0, 4, 0}}, rlist);
        fl2.AddRealDY(R::GUX_MUXMU_UX, FKS::Type_t::QCD,
                      {{21, -2, -13, 13, -2}}, {{0, -2, 0, -4}}, rlist);
    }
    list.push_back(fl2);

    color1 = {0, 501, 0, 0};
    color2 = {501, 0, 0, 0};
    FKS::FlavourConfig fl3(2, {{-1, 1, -13, 13}}, {{-1, 1, -3, 3, -5, 5}},
                           color1, color2);
    if (EW) {
        FKS::ResonanceList resonances;
        size_t isr_res = resonances.Add(
            FKS::Resonance({{2, 3}}, 91.1876, 2.4952, 1.0 / 3.0));
        size_t fsr_res =
            resonances.Add(FKS::Resonance({{2, 3, 4}}, 91.1876, 2.4952, 1.0));
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 2, fsr_res));
        rlist.push_back(FKS::Region(4, 3, fsr_res));
        rlist.push_back(FKS::Region(4, 0, isr_res));
        fl3.AddRealDY(R::DXD_MUXMU_A, FKS::Type_t::EW, {{-1, 1, -13, 13, 22}},
                      {{-1, 1, -3, 3, -5, 5}}, rlist, resonances);
    }
    if (QCD) {
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 0));
        fl3.AddRealDY(R::DXD_MUXMU_G, FKS::Type_t::QCD, {{-1, 1, -13, 13, 21}},
                      {{-1, 1, -3, 3, -5, 5}}, rlist);
        fl3.AddRealDY(R::DXG_MUXMU_DX, FKS::Type_t::QCD,
                      {{-1, 21, -13, 13, -1}}, {{-1, 0, -3, 0, -5, 0}}, rlist);
        fl3.AddRealDY(R::GD_MUXMU_D, FKS::Type_t::QCD, {{21, 1, -13, 13, 1}},
                      {{0, 1, 0, 3, 0, 5}}, rlist);
    }
    list.push_back(fl3);

    color1 = {501, 0, 0, 0};
    color2 = {0, 501, 0, 0};
    FKS::FlavourConfig fl4(3, {{1, -1, -13, 13}}, {{1, -1, 3, -3, 5, -5}},
                           color1, color2);
    if (EW) {
        FKS::ResonanceList resonances;
        size_t isr_res = resonances.Add(
            FKS::Resonance({{2, 3}}, 91.1876, 2.4952, 1.0 / 3.0));
        size_t fsr_res =
            resonances.Add(FKS::Resonance({{2, 3, 4}}, 91.1876, 2.4952, 1.0));
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 2, fsr_res));
        rlist.push_back(FKS::Region(4, 3, fsr_res));
        rlist.push_back(FKS::Region(4, 0, isr_res));
        fl4.AddRealDY(R::DDX_MUXMU_A, FKS::Type_t::EW, {{1, -1, -13, 13, 22}},
                      {{1, -1, 3, -3, 5, -5}}, rlist, resonances);
    }
    if (QCD) {
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 0));
        fl4.AddRealDY(R::DDX_MUXMU_G, FKS::Type_t::QCD, {{1, -1, -13, 13, 21}},
                      {{1, -1, 3, -3, 5, -5}}, rlist);
        fl4.AddRealDY(R::DG_MUXMU_D, FKS::Type_t::QCD, {{1, 21, -13, 13, 1}},
                      {{1, 0, 3, 0, 5, 0}}, rlist);
        fl4.AddRealDY(R::GDX_MUXMU_DX, FKS::Type_t::QCD,
                      {{21, -1, -13, 13, -1}}, {{0, -1, 0, -3, 0, -5}}, rlist);
    }
    list.push_back(fl4);

    return list;
}
