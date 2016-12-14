#ifndef PHYSICS_PDGCODE_H_
#define PHYSICS_PDGCODE_H_

#include <map>
#include <string>
#include <vector>

namespace Physics {

/**
pdg codes
*/
typedef int PDGCode;

class PDG {
  public:
    ///@name Quarks
    ///@{
    static const PDGCode DQuark = 1;
    static const PDGCode UQuark = 2;
    static const PDGCode SQuark = 3;
    static const PDGCode CQuark = 4;
    static const PDGCode BQuark = 5;
    static const PDGCode TQuark = 6;
    ///@}

    ///@name AntiQuarks
    ///@{
    static const PDGCode AntiDQuark = -1;
    static const PDGCode AntiUQuark = -2;
    static const PDGCode AntiSQuark = -3;
    static const PDGCode AntiCQuark = -4;
    static const PDGCode AntiBQuark = -5;
    static const PDGCode AntiTQuark = -6;
    ///@}

    ///@name Leptons
    ///@{
    static const PDGCode Electron = 11;
    static const PDGCode ElectronNeutrino = 12;
    static const PDGCode Muon = 13;
    static const PDGCode MuonNeutrino = 14;
    static const PDGCode Tau = 15;
    static const PDGCode TauNeutrino = 16;
    ///@}

    ///@name AntiLeptons
    ///@{
    static const PDGCode Positron = -11;
    static const PDGCode AntiElectronNeutrino = -12;
    static const PDGCode AntiMuon = -13;
    static const PDGCode AntiMuonNeutrino = -14;
    static const PDGCode AntiTau = -15;
    static const PDGCode AntiTauNeutrino = -16;
    ///@}

    ///@name Gauge & Higgs bosons
    ///@{
    static const PDGCode Gluon = 21;
    static const PDGCode Photon = 22;
    static const PDGCode Z = 23;
    static const PDGCode WPlus = 24;
    static const PDGCode WMinus = -24;
    static const PDGCode Higgs = 25;
    static const PDGCode ZPrime = 32;
    static const PDGCode ZPrimePrime = 33;
    static const PDGCode WPrime = 34;
    static const PDGCode HeavyHiggs = 35;
    static const PDGCode PseudoscalarHiggs = 36;
    static const PDGCode HPlus = 37;
    static const PDGCode HMinus = -37;
    ///@}

    ///@name string functions
    ///@{
    static std::string Name(PDGCode pdg);
    static std::string LatexName(PDGCode pdg);
    static std::string CodesToName(const std::vector<int> &codes);
    static std::string Type(PDGCode pdg);
    ///@}

    ///@name particle type functions
    ///@{
    static bool IsQuark(PDGCode pdg);
    static bool IsLepton(PDGCode pdg);
    static bool IsChargedLepton(PDGCode pdg);
    static bool IsChargedFermion(PDGCode pdg);
    static bool IsNeutrino(PDGCode pdg);
    static bool IsGluon(PDGCode pdg) {
        if (pdg == 21 || pdg == 0)
            return true;
        return false;
    }
    static bool IsPhoton(PDGCode pdg) {
        if (pdg == 22) {
            return true;
        }
        return false;
    }
    ///@}
    
    static double Charge(PDGCode pdg);
    static double ChargeAbs(PDGCode pdg);

  private:
    static const std::map<PDGCode, std::pair<std::string, std::string> > names;
};

} // namespace

#endif
