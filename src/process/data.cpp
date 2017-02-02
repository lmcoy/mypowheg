#include "process/data.h"

#include <complex>
#include <iostream>
#include "config/file.h"

using namespace UserProcess;

Data::Data() {
    hists = 0;
    cuts = 0;
    Params = 0;
    Params_as = 0;
}

int Data::Init(const char *filename) {
    Config::File cfile;
    Config::File::ReadError cerror = cfile.ReadFromFile(filename);
    if (cerror != Config::File::ReadError::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    int ret = readConfig(cfile);
    if (ret == 1) {
        return 1;
    }
    ret =  ProcessInit(cfile);
    if (ret == 1) {
        return 1;
    }
    if (!MatrixElement) {
        return 1;
    }
    if (Process.size() == 0) {
        return 1;
    }
    if (!cuts) {
        return 1;
    }
    return 0;
}

int Data::readConfig(Config::File & cfile) {
    RecombinationParam rparam;
    if (cfile.GetDouble("recombination", "dR", &rparam.dR) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        exit(1);
    }

    double sqrtS = -1.0;
    if (cfile.GetDouble("config", "SqrtS", &sqrtS) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    {
        long verb = 0;
        auto c = cfile.GetInt("config", "Verbose", &verb);
        switch(c) {
            case Config::File::Error::NoError:
                break;
            case Config::File::Error::KeyNotFound:
                verb = 0;
                break;
            default:
                verb = 0;
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (verb > 0) {
            Verbose = verb;
        }
    }
    {
        long xsec = 0;
        auto c = cfile.GetInt("config", "CalculateXSec", &xsec);
        switch(c) {
            case Config::File::Error::NoError:
                break;
            case Config::File::Error::KeyNotFound:
                xsec = 0;
                break;
            default:
                xsec = 0;
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (xsec > 0) {
            CalculateXSec = true;
        }
    }

    auto &iparams = IntParams;
    if (cfile.GetInt("integration.setup", "iterations",
                     &iparams.iterations_setup) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    if (cfile.GetInt("integration.setup", "nevents", &iparams.nevents_setup) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    if (cfile.GetInt("integration", "iterations", &iparams.iterations) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    if (cfile.GetInt("integration", "nevents", &iparams.nevents) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    {
        long ignorevirt = 0;
        auto c = cfile.GetInt("integration.setup", "IgnoreVirtual", &ignorevirt);
        switch(c) {
            case Config::File::Error::NoError:
                break;
            case Config::File::Error::KeyNotFound:
                ignorevirt = 0;
                break;
            default:
                ignorevirt = 0;
            std::cerr << cfile.ErrorMsg() << "\n";
            return 1;
        }
        if (ignorevirt > 0) {
            IgnoreVirtualInSetup = true;
        }
    }
    if (cfile.GetInt("xsec.setup", "iterations",
                     &iparams.iterations_xsec_setup) !=
        Config::File::Error::NoError) {
        iparams.iterations_xsec_setup = iparams.iterations_setup;
    }
    if (cfile.GetInt("xsec.setup", "nevents", &iparams.nevents_xsec_setup) !=
        Config::File::Error::NoError) {
        iparams.nevents_xsec_setup = iparams.nevents_setup;
    }
    if (cfile.GetInt("xsec", "iterations", &iparams.iterations_xsec) !=
        Config::File::Error::NoError) {
        iparams.iterations_xsec = iparams.iterations;
    }
    if (cfile.GetInt("xsec", "nevents", &iparams.nevents_xsec) !=
        Config::File::Error::NoError) {
        iparams.nevents_xsec = iparams.nevents;
    }

    long seed1;
    if (cfile.GetInt("random", "seed1", &seed1) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }

    long seed2;
    if (cfile.GetInt("random", "seed2", &seed2) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }

    long remnXi;
    if (cfile.GetInt("RadiationSampling", "RemnantXi", &remnXi) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }

    long realXi;
    if (cfile.GetInt("RadiationSampling", "RealXi", &realXi) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }

    long realY;
    if (cfile.GetInt("RadiationSampling", "RealY", &realY) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }

    long realPhi;
    if (cfile.GetInt("RadiationSampling", "RealPhi", &realPhi) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }

    long RadQCD;
    if (cfile.GetInt("RadiationType", "QCD", &RadQCD) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    if (RadQCD < 0 || RadQCD >1) {
        std::cerr << "RadiationType.QCD must be 0 to disable QCD radiation or "
                     "1 to enable it. RadiationType.QCD = " << RadQCD
                  << "is invalid\n";
        return 1;
    }
    if (RadQCD == 1) {
        RadiationType.QCD = true;
    }
    long RadEW;
    if (cfile.GetInt("RadiationType", "EW", &RadEW) !=
        Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    if (RadEW < 0 || RadEW >1) {
        std::cerr << "RadiationType.EW must be 0 to disable EW radiation or "
                     "1 to enable it. RadiationType.EW = " << RadEW
                  << "is invalid\n";
        return 1;
    }
    if (RadEW == 1) {
        RadiationType.EW = true;
    }

    if (cfile.GetInt("PDF", "LHAID", &Lhaid) != Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }

    long genN = 0;
    if (cfile.GetInt("EventGeneration", "nevents", &genN) != Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }

    long genNperIt = 0;
    if (cfile.GetInt("EventGeneration", "nperiteration", &genNperIt) != Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    long seed3 = 0;
    if (cfile.GetInt("EventGeneration", "seed", &seed3) != Config::File::Error::NoError) {
        std::cerr << cfile.ErrorMsg() << "\n";
        return 1;
    }
    if (genN <= 0) {
        std::cerr << "invalid number of events: " << genN << "\n";
        return 1;
    }
    GenEvent.N = genN;
    if(genNperIt <= 0) {
        std::cerr << "invalid number of events per iteration: " << genNperIt << "\n";
        return 1;
    }
    GenEvent.NperIt = genNperIt;
    GenEvent.seed3 = seed3;
    {
        long negevents = 0;
        auto c = cfile.GetInt("EventGeneration", "NegativeEvents", &negevents);
        switch(c) {
            case Config::File::Error::NoError:
                break;
            case Config::File::Error::KeyNotFound:
                negevents = 0;
                break;
            default:
                negevents = 0;
            std::cerr << cfile.ErrorMsg() << "\n";
            return 1;
        }
        if (negevents > 0) {
            NegativeEvents = true;
        }
    }
    {
        long bornunweighting = 0;
        auto c = cfile.GetInt("EventGeneration", "BornUnweighting",
                              &bornunweighting);
        switch(c) {
            case Config::File::Error::NoError:
                break;
            case Config::File::Error::KeyNotFound:
                bornunweighting = 0;
                break;
            default:
                bornunweighting = 0;
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (bornunweighting > 0) {
            BornUnweighting = true;
        }
    }
    {
        long genevents = 1;
        auto c = cfile.GetInt("EventGeneration", "GenerateEvents", &genevents);
        switch(c) {
            case Config::File::Error::NoError:
                break;
            case Config::File::Error::KeyNotFound:
                genevents = 1;
                break;
            default:
                genevents = 1;
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (genevents == 0) {
            GenerateEvents = false;
        }
    }
    {
        long bornonly = 0;
        auto c = cfile.GetInt("EventGeneration", "BornOnly", &bornonly);
        switch(c) {
            case Config::File::Error::NoError:
                break;
            case Config::File::Error::KeyNotFound:
                bornonly = 0;
                break;
            default:
                bornonly = 0;
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (bornonly > 0) {
            BornOnly = true;
        } else {
            BornOnly = false;
        }
    }

    {
        long cutreal = 0;
        auto c = cfile.GetInt("cuts", "CutReal", &cutreal);
        switch(c) {
            case Config::File::Error::NoError:
                break;
            case Config::File::Error::KeyNotFound:
            case Config::File::Error::CategoryNotFound:
                cutreal = 0;
                break;
            default:
                cutreal = 0;
            std::cerr << cfile.ErrorMsg() << "\n";
            return 1;
        }
        if (cutreal > 0) {
            CutOnReal = true;
        }
    }

    {
        std::string approx = "";
        auto c = cfile.GetString("RadiationType", "Approximation", &approx);
        switch (c) {
        case Config::File::Error::NoError:
            break;
        case Config::File::Error::KeyNotFound:
            approx = "";
            break;
        default:
            approx = "";
            std::cerr << cfile.ErrorMsg() << "\n";
            return 1;
        }
        if (approx == "") {
            RadiatePhoton = true;
            OnlyVirtualEW = false;
        } else if (approx == "NoPhotonRadiation") {
            RadiatePhoton = false;
            OnlyVirtualEW = false;
        } else if (approx == "NeglectRealEW") {
            RadiatePhoton = false;
            OnlyVirtualEW = true;
        } else {
            std::cerr << "error: unknown approximation \"" << approx
                      << "\" (valid approximations: NoPhotonRadiation, "
                         "NeglectRealEW)\n";
            return 1;
        }

        if (!RadiationType.EW && (!RadiatePhoton || OnlyVirtualEW)) {
            std::cerr << "warning: approximation useless "
                         "because there are no EW corrections anyway\n";
        }
    }

    {
        long rad_qed = 1;
        auto c = cfile.GetInt("RadiationType", "RadiationQED", &rad_qed);
        switch (c) {
            case Config::File::Error::NoError:
                break;
            case Config::File::Error::KeyNotFound:
            case Config::File::Error::CategoryNotFound:
                rad_qed = 1;
                break;
            default:
                rad_qed = 1;
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (rad_qed > 0) {
            RadiatePhoton = true;
        } else {
            RadiatePhoton = false;
        }
    }

    {
        long rad_qcd = 1;
        auto c = cfile.GetInt("RadiationType", "RadiationQCD", &rad_qcd);
        switch (c) {
            case Config::File::Error::NoError:
                break;
            case Config::File::Error::KeyNotFound:
            case Config::File::Error::CategoryNotFound:
                rad_qcd = 1;
                break;
            default:
                rad_qcd = 1;
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (rad_qcd > 0) {
            RadiateQCD = true;
        } else {
            RadiateQCD = false;
        }
    }

    {
        std::string splitME = "";
        auto c = cfile.GetString("RadiationType", "SplitME", &splitME);
        switch (c) {
        case Config::File::Error::NoError:
            break;
        case Config::File::Error::KeyNotFound:
            splitME = "";
            break;
        default:
            splitME = "";
            std::cerr << cfile.ErrorMsg() << "\n";
            return 1;
        }
        if (splitME == "" || splitME == "no" || splitME == "0") {
            modBbar = false;
            modSudakov = false;
        } else if (splitME == "yes" || splitME == "1") {
            modBbar = true;
            modSudakov = true;
        } else {
            std::cerr << "error: unknown value for SplitME \"" << splitME
                      << "\" (valid values (yes/no). If SplitME is not "
                         "specified, 'no' is used as default.)\n";
            return 1;
        }
    }

    {
        long modS = 0;
        auto c = cfile.GetInt("RadiationType", "ModS", &modS);
        switch (c) {
            case Config::File::Error::NoError:
                break;
            case Config::File::Error::KeyNotFound:
                modS = 0;
                break;
            default:
                modS = 0;
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (modS < 0 || modS > 2) {
            std::cerr << "error: unknown value for ModS: " << modS
                      << "(expected: 0, 1, 2)\n";
        }
        useResonanesInS = modS;
    }

    {
        std::string inter = "";
        auto c = cfile.GetString("RadiationType", "NoInterference", &inter);
        switch (c) {
            case Config::File::Error::NoError:
                break;
            case Config::File::Error::KeyNotFound:
                inter = "";
                break;
            default:
                inter = "";
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (inter == "" || inter == "no" || inter == "0") {
            noInterference = false;
        } else if (inter == "yes" || inter == "1") {
            noInterference = true;
        } else {
            std::cerr << "error: unknown value for NoInterference \"" << inter
                      << "\" (valid values (yes/no). If NoInterference is not "
                         "specified, 'no' is used as default.)\n";
            return 1;
        }
    }

    {
        std::string pdf = "";
        auto c = cfile.GetString("RadiationType", "QCDPDFScheme", &pdf);
        switch (c) {
            case Config::File::Error::NoError:
                break;
            case Config::File::Error::KeyNotFound:
                pdf = "";
                break;
            default:
                pdf = "";
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (pdf == "" || pdf == "MSbar") {
            PDFRenorm.QCD = FKS::PDFRenorm::MSbar;
        } else if (pdf == "DIS" ) {
            PDFRenorm.QCD = FKS::PDFRenorm::DIS;
        } else {
            std::cerr << "error: unknown value for QCDPDFScheme \"" << pdf
                      << "\" (valid values (MSbar/DIS). \n";
            return 1;
        }
    }

    {
        std::string pdf = "";
        auto c = cfile.GetString("RadiationType", "EWPDFScheme", &pdf);
        switch (c) {
            case Config::File::Error::NoError:
                break;
            case Config::File::Error::KeyNotFound:
                pdf = "";
                break;
            default:
                pdf = "";
                std::cerr << cfile.ErrorMsg() << "\n";
                return 1;
        }
        if (pdf == "" || pdf == "DIS") {
            PDFRenorm.EW = FKS::PDFRenorm::DIS;
        } else if (pdf == "MSbar" ) {
            PDFRenorm.EW = FKS::PDFRenorm::MSbar;
        } else {
            std::cerr << "error: unknown value for EWPDFScheme \"" << pdf
                      << "\" (valid values (MSbar/DIS). \n";
            return 1;
        }
    }
    {
        long guessv = 0;
        auto c = cfile.GetInt("Unweighting", "GuessVirtual", &guessv);
        switch (c) {
        case Config::File::Error::NoError:
            break;
        case Config::File::Error::KeyNotFound:
            guessv = 1;
            break;
        default:
            guessv = 0;
            std::cerr << cfile.ErrorMsg() << "\n";
            return 1;
        }
        switch (guessv) {
        case 0:
            Unweighting.GuessVirtual = false;
            break;
        case 1:
            Unweighting.GuessVirtual = true;
            break;
        default:
            std::cerr << "error: unknown value for Unweighting.GuessVirtual: "
                         "expected 0 for disable or 1 for enable\n";
            return 1;
        }
    }
    {
        long l = 0;
        auto c = cfile.GetInt("Unweighting", "MaxMultiple", &l);
        switch (c) {
        case Config::File::Error::NoError:
            break;
        case Config::File::Error::KeyNotFound:
            l = 3;
            break;
        default:
            l = 3;
            std::cerr << cfile.ErrorMsg() << "\n";
            return 1;
        }
        if (l < 1) {
            std::cerr
                << "error: Unweighting.MaxMultiple has to be >= 1: using 1\n";
            l = 1;
        }
        Unweighting.MaxMultiple = l;
    }

    {
        double d = -1.0;
        if (cfile.GetDouble("Radiation", "kT2min", &d) !=
            Config::File::Error::NoError) {
            std::cerr << cfile.ErrorMsg() << "\n";
            return 1;
        }
        if (d <= 0.0) {
            std::cerr << "error: Radiation.kT2min has to be > 0.0\n";
            return 1;
        }
        RadiationParameter.kT2min = d;
    }
    Recomb = rparam;

    SqrtS = sqrtS;

    Seed1 = (int)seed1;
    Seed2 = (int)seed2;

    NRealXi = (int)realXi;
    NRealY = (int)realY;
    NRealPhi = (int)realPhi;
    NRemnXi = (int)remnXi;

    return 0;
}

void Data::Reset() {
    BornMEStatus = BornMEStatus_t::None;
    hists->Reset();
}

Data::~Data() {}

namespace {

void print_pdfren(const char label[], FKS::PDFRenorm pdfren) {

    char dis[] = "DIS";
    char msbar[] = "MS bar";
    switch (pdfren) {
    case FKS::PDFRenorm::DIS:
        printf(" %-15s = %s\n", label, dis);
        break;
    case FKS::PDFRenorm::MSbar:
        printf(" %-15s = %s\n", label, msbar);
        break;
    }
}

} // namespace

void Data::Print() const {
    printf("/******************************************************************"
           "***\n"
           " *                        Approximations                           "
           "  *\n"
           " ******************************************************************"
           "***/\n");
    printf(" internal parameters:\n");
    printf(" %-15s = %d\n", "OnlyVirtualEW", OnlyVirtualEW);
    printf(" %-15s = %d\n", "RadiatePhoton", RadiatePhoton);
    printf(" %-15s = %d\n", "RadiateQCD", RadiateQCD);
    printf("\n\n");

    printf("/******************************************************************"
           "***\n"
           " *                      Resonance Handling                         "
           "  *\n"
           " ******************************************************************"
           "***/\n");
    printf(" split matrix element:\n");
    printf(" %-15s = %d\n", "mod. Bbar", modBbar);
    printf(" %-15s = %d\n", "mod. Sudakov", modSudakov);
    printf("\n");

    printf(" modified S functions (Drell Yan):\n");
    printf(" %-15s = %d\n", "mod S", useResonanesInS);
    printf("\n");

    printf(" disable interference terms: %d\n", noInterference);
    printf("\n\n");

    printf("/******************************************************************"
           "***\n"
           " *                     PDF renormalization                         "
           "  *\n"
           " ******************************************************************"
           "***/\n");
    print_pdfren("QCD scheme", PDFRenorm.QCD);
    print_pdfren("EW scheme", PDFRenorm.EW);
    printf("\n\n");

    printf("/******************************************************************"
           "***\n"
           " *                     Unweighting                                 "
           "  *\n"
           " ******************************************************************"
           "***/\n");
    printf( " %-15s = %d\n", "guess virtual", Unweighting.GuessVirtual );
    printf( " %-15s = %d\n", "MaxMultiple", Unweighting.MaxMultiple );
    printf("\n\n");

    printf("/******************************************************************"
           "***\n"
           " *                     Radiation                                   "
           "  *\n"
           " ******************************************************************"
           "***/\n");
    printf( " %-15s = %g\n", "kT2min", RadiationParameter.kT2min);
    printf("\n\n");

    printf("/******************************************************************"
           "***\n"
           " *                     NLO                                         "
           "  *\n"
           " ******************************************************************"
           "***/\n");
    printf( " %-15s = %d\n", "apply cuts to real", CutOnReal);
    printf("\n\n");

    printf("/******************************************************************"
           "***\n"
           " *                     Events                                      "
           "  *\n"
           " ******************************************************************"
           "***/\n");
    printf( " %-15s = %d\n", "negative events", NegativeEvents);
    printf( " %-15s = %d\n", "born unweighting", BornUnweighting);
    printf( " %-15s = %d\n", "born only", BornOnly);
    printf( " %-15s = %d\n", "ignore V in grid setup", IgnoreVirtualInSetup);
    printf("\n\n");

    ProcessPrint(); 
}
