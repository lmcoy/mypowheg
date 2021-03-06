#ifndef FKS_PROCESS_H_
#define FKS_PROCESS_H_

#include <array>
#include <string>
#include <vector>

#include "fks/regions.h"
#include "fks/splitting.h"

namespace Phasespace {
class Phasespace;
}

namespace FKS {

typedef std::vector<int> PDGList;
typedef std::vector<std::array<int, 2> > PDFList;
typedef std::vector<PDGList> AllPDGs;
typedef std::vector<int> ColorFlow;

/**
 * @brief Correction Type
 */
enum class Type_t {
    EW, ///< Electroweak corrections
    QCD ///< QCD corrections
};

/**
 * @brief Born Process
 *
 * Drell-Yan: u u~ > mu+ mu-
 * Flavours: 2, -2, -13, 13
 * PDF:      2, -2,  4, -4 (use u and c quarks)
 * ID:       depends on impl.
 */
struct Born_t {
    PDGList Flavours;                 ///< Born flavours
    std::vector<PDGList> AllFlavours; ///< all born flavours;
    PDFList PDF; ///< List of PDF sets which are used for the inital state.
    int ID;      ///< Number which identifies the matrix element in BornME
    ColorFlow Color[2];
};

/**
 * @brief Breit-Wigner resonance
 *
 * If a process has a resonance, we can model it by a Breit-Wigner peak which is
 * defined by the resonance mass and width. The Resonance class models such a
 * peak and also specifies the decay products of the resonance.
 */
class Resonance {
  private:
    double Width = 1.0;
    double Mass = 0.0;
    double Q2 = 1.0;
    std::vector<int> Daughters;

    /** If true, BreitWigner returns always 1.0. */
    bool Disabled = false;

  public:
    /**
     * @brief Create a fake resonance
     *
     * Fake resonance means a that there is no Breit-Wigner Peak. It will always
     * return 1.
     */
    Resonance() : Disabled(true) {}

    /**
     * @brief Constructor
     *
     * @param daughters The indices of the decay products. These indices are
     * used
     *                  to identify the momentum position in the phase-space.
     * @param mass mass of the resonance
     * @param width width of the resonance
     */
    Resonance(std::vector<int> &&daughters, double mass, double width, double Q)
        : Width(width), Mass(mass), Q2(Q * Q), Daughters(std::move(daughters)) {
    }

    /**
     * @brief Breit-Wigner propagator
     *
     * BreitWigner returns the value of a Breit-Wigner propagator for the
     * phase-space point ps.
     */
    double BreitWigner(const Phasespace::Phasespace &ps, bool useCharge) const;
};

/**
 * @brief ResonanceList
 *
 * ResonanceList is a wrapper of std::vector<Resonance>. It always adds a 0th
 * element with a fake resonance.
 */
class ResonanceList {
  public:
    ResonanceList() : list(1, Resonance()) {}

    /**
     * @brief Append Resonance
     *
     * Add appends the Resonance res to the list and returns its index in the
     * list.
     */
    size_t Add(const Resonance &res) {
        list.push_back(res);
        return list.size() - 1;
    }

    const Resonance &operator[](size_t i) const { return list[i]; }

    Resonance &operator[](size_t i) { return list[i]; }

    size_t size() const { return list.size(); }

    std::vector<Resonance>::iterator begin() { return list.begin(); }

    std::vector<Resonance>::const_iterator begin() const {
        return list.begin();
    }

    std::vector<Resonance>::iterator end() { return list.end(); }

    std::vector<Resonance>::const_iterator end() const { return list.end(); }

  private:
    std::vector<Resonance> list;
};

/**
 * @brief Real Process
 */
struct Real_t {
    PDGList Flavours; ///< Real flavours
    PDFList PDF;      ///< List of PDF sets which are used for the inital state.
    RegionList Regions; ///< FKS regions
    Type_t Type;        ///< correction type
    int ID; ///< Number which identifies the matrix element in RealME
    ResonanceList Resonances;
    std::vector<AllPDGs> AllFlavours;
    RegionList AllRegions;
};

/**
 * @brief Collinear Remnant
 */
struct Remnant_t {
    Remnant_t(const FKS::Splitting &sp) : Splitting(sp) {}
    FKS::Splitting Splitting;
    std::vector<int> PDF;
    Type_t Type;
};

/**
 * @brief process information.
 *
 * FlavourConfig represents the flavour structures of a NLO process.
 * It contains the born flavour structure Born, a list of real flavour
 * structures and the collinear remnants which repesent the renormalization of
 * the pdfs due to factorization.
 */
class FlavourConfig {
  public:
    Born_t Born;                     ///< Born process
    std::vector<Real_t> Real;        ///< List of real processes
    std::vector<Remnant_t> Remnant1; ///< List of remnants for the 1st initial
                                     /// state
    std::vector<Remnant_t> Remnant2; ///< List of remnants for the 2nd initial
                                     /// state
    std::vector<double> Scales;      ///< List of scales which are used to scale
                                     /// BornME
    bool QCD = false;                ///< True, if a QCD real process is in Real
    bool EW = false;                 ///< True, if a EW real process is in Real

    /**
     * @brief Construct
     *
     * Creates a new FlavourConfig for a born process.
     *
     * @deprecated {@see AddReal()}
     *
     * @param id number which identifies the matrix element in BornME
     * @param pdgs  pdgs of the born process
     * @param pdfs  pdfs for the initial state partons
     * @param color1 color flow tags for particles in pdgs (color)
     * @param color2 color flow tags for particles in pdgs (anti-color)
     * @param scales scales for the matrix element. There has to be a scale for
     *               every pdf configuration. Therefore, scales has to have a
     *               length of pdfs.size()/2 or 0. If the length is zero a scale
     *               of 1 is assumed.
     */
    FlavourConfig(int id, const PDGList &pdgs, const std::vector<int> &pdfs,
                  const ColorFlow &color1, const ColorFlow &color2,
                  const std::vector<double> &scales = {});

    /**
     * @brief Construct
     *
     * Creates a new FlavourConfig for a born process. A born process can
     * have more than one pdg configuration. This is useful for processes with
     * the same matrix element but different pdfs. The code uses the first
     * process in `pdgs` as main process and adds the pdf information of the
     * other processes in `pdgs`.
     *
     * @param id number which identifies the matrix element in
     *           UserProcess::IMatrixElement::Born()
     * @param pdgs  pdgs of the born process
     * @param color1 color flow tags for particles in pdgs (color)
     * @param color2 color flow tags for particles in pdgs (anti-color)
     * @param scales scales for the matrix element. There has to be a scale for
     *               every element in `pdgs`. If the length is zero a scale
     *               of 1 is assumed.
     */
    FlavourConfig(int id, const std::vector<PDGList> &pdgs, int dummy,
                  const ColorFlow &color1, const ColorFlow &color2,
                  const std::vector<double> &scales = {});
    /**
     * AddRealDY appends a new real process to FlavourConfig.
     *
     * @deprecated {This version of AddReal is deprecated because it is not
     * general enough. This functions assumes that we can simply substitute
     * new pdfs to intial state partons to obtain a new process. This is in
     * general not true. Consider a very simple example: g u > Z u. We have
     * to change the intial state quark and the final state quark pdg to
     * have a correct real process.}
     *
     * @param id number which identifies the matrix element in RealME
     * @param pdgs pdgs for the real process
     * @param pdfs pdfs for the initial state partons
     * @param regions singular regions of the matrix element
     */
    void AddRealDY(int id, Type_t type, const PDGList &pdgs,
                   const std::vector<int> &pdfs, const RegionList &regions,
                   const ResonanceList &resonances = ResonanceList());

    /**
     * AddReal appends a new real process to FlavourConfig. A real process can
     * have more than one pdg configuration. This is useful for processes with
     * the same matrix element but different pdfs. The code uses the first
     * process in `pdgs` as main process and adds the pdf information of the
     * other processes in `pdgs`.
     *
     * @param id number which identifies the matrix element in
     *           UserProcess::IMatrixElement::Real().
     * @param pdgs pdgs for the real process.
     * @param regions singular regions of the matrix element
     * @param resonances resonances of the process. This is used to modify the S
     *                   functions as described in 1509.09071 eq. (2.9).
     */
    void AddReal(int id, Type_t type, const std::vector<PDGList> &pdgs,
                 const RegionList &regions, const RegionList &allregions,
                 const ResonanceList &resonances = ResonanceList());
    /**
     * Print prints information about FlavourConfig to stdout.
     */
    void Print() const;

  private:
    void AddRemnantInitial1(Type_t type, FKS::Splitting split,
                            const PDGList &pdgs) {
        Remnant_t remn(split);
        remn.Type = type;
        remn.Splitting = split;
        remn.PDF = pdgs;
        Remnant1.push_back(remn);
    }

    void AddRemnantInitial2(Type_t type, FKS::Splitting split,
                            const PDGList &pdgs) {
        Remnant_t remn(split);
        remn.Type = type;
        remn.PDF = pdgs;
        Remnant2.push_back(remn);
    }
};

typedef std::vector<FlavourConfig> ProcessList;

} // end namespace FKS

#endif
