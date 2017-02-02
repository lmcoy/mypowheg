#ifndef PROCESS_DATA_H_
#define PROCESS_DATA_H_

#include <memory>

#include "fks/histlist.h"
#include "fks/process.h"
#include "fks/radiationregion.h"
#include "fks/remnants.h"
#include "fks/scales.h"
#include "lhe/eventbuffer.h"
#include "pdf/pdfinterface.h"
#include "physics/alphas.h"
#include "powheg/resonance.h"
#include "process/cuts.h"
#include "process/histograms.h"
#include "process/matrixelement.h"
#include "random/rnd.h"
#include "util/histogram.h"

namespace Config {
class File;
}

namespace FKS {
class Param;
class Param_as;
}

namespace UserProcess {

struct PDFRenorm_t {
    FKS::PDFRenorm QCD = FKS::PDFRenorm::MSbar;
    FKS::PDFRenorm EW = FKS::PDFRenorm::DIS;
};

struct RecombinationParam {
    double dR = 0.1;
};

struct GenEventStats {
    int N = 0;
    int N_ENORM = 0;    ///< Norm for upper bounding function too small
    int N_ERADVAR = 0;  ///< radiation variables are out of bounds
    int N_NEG = 0;      ///< btilde is negative
    int N_MAX = 0;      ///< maximum of btilde is wrong
    int N_RADERROR = 0; ///< Error in RadiationGenerator
    int N_REJECT = 0;   ///< Rejected Events
    int N_BORN = 0;     ///< born events
    int N_REAL = 0;     ///< events with radiation
    int N_REJECTWOVIRTUAL = 0;
    int N_REJECTBORN = 0;
    int N_2TIMES = 0; ///< event written 2 times
    int N_3TIMES = 0; ///< event written 3 times
    int N_4TIMES = 0; ///< event written 4 times
    int N_NTIMES = 0; ///< event written >4 times
    int N_WRONGV = 0; ///< number of events where the guessed V was too small.
};

struct EventGen_t {
    int seed3;
    int N;
    int NperIt;
};

struct BtildeState_t {
    double Max = 0.0;
    double MaxBorn = 0.0;
    double MaxVoverB = 0.0;
    bool UpdateMax = true;
    HistList hists;

    BtildeState_t() {}
    ~BtildeState_t() {}

    void Reset() {
        Max = 0.0;
        MaxBorn = 0.0;
        MaxVoverB = 0.0;
    }

    void FillMax(double x) {
        auto hist = hists[0];
        if (hist == 0) { return; }
        hist->AddValue(x, 1.0);
    }

    void InitHists() {
        double xmin = 1.0;
        double xmax = 1e8;
        if (Max > 0.0) {
            xmax = 1.1 * Max;
            xmin = Max / 300.0;
        }
        hists.Append(300, xmin, xmax);
    }

    void MergeHistograms(Util::DataBuffer *buffer) {
        HistList histlist;
        histlist.ReadFromBuffer(buffer);
        hists.Add(histlist);
    }

    void WriteHistograms(Util::DataBuffer *buffer) {
        Util::DataBuffer b(Util::DataBuffer::GETSIZE);
        hists.WriteToBuffer(&b);
        int size = b.GetDataSize() + 100;

        buffer->ClearAndAllocate(size);
        hists.WriteToBuffer(buffer);
    }

    bool WriteHistogramsToFile(const char * filename) {
        std::ofstream of(filename, std::ios::out);
        if (!of.good()) {
            return false;
        }
        WriteHistogramsToStream(of);
        of.close();
        return true;
    }

    void WriteHistogramsToStream(std::ostream &of) {
        of << "set title \"Bbar\"\n";
        of << "set xlabel \"Bbar\"\n";
        of << "set ylabel \"#\"\n";
        of << "a = " << hists[0]->GetUnderflow() << "\n";
        of << "b = " << hists[0]->GetNumEntries() << "\n";
        of << "csum(x)=(a=a+x,a)/b\n";
        of << "plot \"-\" using 1:(csum($2)) with lines ti \"cumulative distribution/total\"\n";
        hists[0]->WriteToStream(of, false);
        of << "e\n\n";
    }
};

class IScales {
  public:
    virtual ~IScales() {}
    virtual double Factorization(const Phasespace::Phasespace &ps) const = 0;
    virtual double Renorm(const Phasespace::Phasespace &ps) const = 0;
};

enum class BornMEStatus_t { None, CrossSection, GenerateEvents };

struct Data {
    struct {
        long iterations_setup;
        long nevents_setup;
        long iterations;
        long nevents;
        long iterations_xsec_setup;
        long nevents_xsec_setup;
        long iterations_xsec;
        long nevents_xsec;
    } IntParams;
    FKS::ProcessList Process;
    int ProcessID;
    FKS::RadiationRegionList RadiationRegions;

    IScales *Scales = 0;
    FKS::Scales PowhegScales;
    BornMEStatus_t BornMEStatus = BornMEStatus_t::None;
    double counterterm;
    double SqrtS;
    RecombinationParam Recomb;
    Histograms *hists = 0;
    ICutsPtr cuts;
    FKS::Param *Params = 0;
    FKS::Param_as *Params_as = 0;
    int Seed1;
    int Seed2;
    int NRealXi;
    int NRealY;
    int NRealPhi;
    int NRemnXi;

    double x_random[20] = {0.0};
    double n_random = 0;

    IMatrixElementPtr MatrixElement;

    struct {
        double kT2min = 0.8; ///< minimal kT2 for radiation
    } RadiationParameter;

    struct {
        /**
         * If GuessedVirtual is true, the virtual part of Btilde is guessed
         * first to be F*born where F is a constant and F*born has to be an
         * upper bound for the virtual part. The virtual part is only calculated
         * if the guessed btilde passed the unweighting.
         * F is given bei BtildeState.MaxVoverB.
         */
        bool GuessVirtual = true;
        /**
         * MaxMultiple specifies how often an event can be written to the event
         * file. If Btilde is larger then BtildeState.Max during event
         * generation, the event will be written up to MaxMultiple times to the
         * event file. It should not happen too often that an event is written
         * more than once! Note: The events are not the they only share a common
         * born configuration.
         */
        int MaxMultiple = 3;
    } Unweighting;

    struct {
        /**
         * @brief Enable/Disable QCD
         */
        bool QCD = false;
        /**
         * @brief Enable/Disable EW
         */
        bool EW = false;
    } RadiationType;

    /**
     * If RadiatePhoton is false no photon radiation is generated in
     * Powheg::GenerateEvents().
     */
    bool RadiatePhoton = true;

    /**
     * If RadiateQCD is false no QCD radiation is generated in
     * Powheg::GenerateEvents().
     */
    bool RadiateQCD = true;

    /**
     * If OnlyVirtualEW is true the real terms for EW corrections are neglected,
     * i.e. POWHEG Btilde = B + V. This is an approximation to the full NLO
     * POWHEG result because the real contribution is only expected to give
     * percent level effects while the main correction including Sudakov logs
     * are in the virtual part. RadiatePhoton should be set to false if
     * OnlyVirtualEW is set to true.
     */
    bool OnlyVirtualEW = false;

    /**
     * split the matrix element in ISR, FSR and interference term. This
     * splitting is introduced for DrellYan EW because of the Z resonance. If
     * this flag is set to true, modSudakov should also be set to true to be
     * consistent.
     */
    bool modBbar = false;
    /**
     * introduce the same matrix element splitting as for modBbar. If set to
     * true, modBbar should also be set to true to be consistent.
     */
    bool modSudakov = false;

    /**
     * Use the modified S functions proposed in 1509.09071. You have to add the
     * resonances to the FKS::Region and a RegionList to FKS::FlavourConfig.
     * useResonancesInS can have the following values:
     *    0: no modification
     *    1: standard Breit-Wigner
     *    2: multiply Breit-Wigner by charge for EW
     */
    int useResonanesInS = 0;

    /**
     * Don't use any interference diagrams. This is used to examine the
     * difference between a split matrix element and the resonance aware S
     * functions. The split matrix element is the best thing you can do for
     * matrix elements without interference terms. We want to see how good the
     * resonance S functions do the mapping.
     */
    bool noInterference = false;

    /**
     * If true, cuts are also applied on the real phase space. This is important
     * if you want to calculate the real cross section. For event generation set
     * it to false because you want to be as inclusive as possible during
     * unweighting. All cuts and recombination has to be applied to the events
     * (except for a generation cut to make the cross section finite).
     */
    bool CutOnReal = false;

    /**
     * If true, events with negative weights are allowed. If false, events with
     * negative weights are neglected.
     */
    bool NegativeEvents = false;

    /**
     * Every event weight is multiplied with this factor. If NegativeEvents is
     * true, the weight should be set to the cross section, where we use the
     * absolute value of integrand.
     */
    double EventWeight = 1.0;

    /**
     * Determines if the virtual matrix element is ignored when the grid for
     * event generation is generated.
     */
    bool IgnoreVirtualInSetup = false;

    int Verbose = 0;

    /**
     * If false, no events are generated.
     */
    bool GenerateEvents = true;

    bool CalculateXSec = false;

    struct {
        int Ninit = 50000;
        int N = 100000;
        /**
         * Set Norm of upper bounding function such that the percentage of
         * smaller values is Ratio.
         */
        double Ratio = 0.995;
    } UpperBoundingParams;

    PDFRenorm_t PDFRenorm;

    Random::RNG rng;
    LHE::EventBuffer EventBuffer;

    GenEventStats GenEventStatistics;
    BtildeState_t BtildeState;

    Powheg::Resonance Resonance;
    EventGen_t GenEvent;

    /*
     * If true, the virtual part in Powheg::Btilde is set to 0.0. This can be
     * useful if the virtual part is computationally expensive and you want to
     * ignore it in the grid setup.
     */
    bool IgnoreVirtualInBtilde = false;

    /**
     * If true, the event generation is unweighted on the born matrix element
     * first. This should exclude born configurations with very small
     * contributions to the cross section.
     */
    bool BornUnweighting = false;

    /**
     * If true, B is used instead of B~, i.e. NLO corrections are disabled.
     */
    bool BornOnly = false;

    Data();
    virtual ~Data();
    int Init(const char *filename);
    void Reset();
    void Print() const;

    void InitEventGeneration() {
        BtildeState.UpdateMax = false;
        GenEventStatistics.N = 0;
        GenEventStatistics.N_ENORM = 0;
        GenEventStatistics.N_ERADVAR = 0;
        GenEventStatistics.N_NEG = 0;
        GenEventStatistics.N_MAX = 0;
        GenEventStatistics.N_RADERROR = 0;
        GenEventStatistics.N_REJECT = 0;
        GenEventStatistics.N_BORN = 0;
        GenEventStatistics.N_REAL = 0;
        GenEventStatistics.N_REJECTWOVIRTUAL = 0;
        GenEventStatistics.N_REJECTBORN = 0;
        GenEventStatistics.N_2TIMES = 0;
        GenEventStatistics.N_3TIMES = 0;
        GenEventStatistics.N_4TIMES = 0;
        GenEventStatistics.N_NTIMES = 0;
        GenEventStatistics.N_WRONGV = 0;
    };

    struct PDFS {
        static const unsigned long First = 0x1;
        static const unsigned long Second = 0x2;
        static const unsigned long Third = 0x4;
        static const unsigned long Fourth = 0x8;
        static const unsigned long Fifth = 0x10;
        static const unsigned long All = 0xffff;
    };
    void UsePDF(unsigned long f) { UsedPDF = f; }
    unsigned long UsedPDF = PDFS::All;
    std::shared_ptr<PDF::Interface> pdf;
    long int Lhaid = 0;

    std::shared_ptr<Physics::IAlphaS> AlphaS;

    bool UseBornCache = false;
    IMatrixElement::Result BornCache;

  protected:
    virtual int ProcessInit(Config::File &) = 0;
    virtual void ProcessPrint() const = 0;

  private:
    int readConfig(Config::File &);
};

} // end Process

#endif
