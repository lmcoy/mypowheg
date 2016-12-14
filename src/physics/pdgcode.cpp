#include "physics/pdgcode.h"

#include "util/stringutil.h"
#include <cassert>

namespace Physics {

// map of string representations of pdg codes
const std::map<PDGCode, std::pair<std::string, std::string> > PDG::names{
    { 0, { "g", "g" } },
    { 2, { "u", "u" } },
    { 1, { "d", "d" } },
    { 3, { "s", "s" } },
    { 4, { "c", "c" } },
    { 5, { "b", "b" } },
    { 6, { "t", "t" } },
    { -2, { "u~", "\\bar{u}" } },
    { -1, { "d~", "\\bar{d}" } },
    { -3, { "s~", "\\bar{s}" } },
    { -4, { "c~", "\\bar{c}" } },
    { -5, { "b~", "\\bar{b}" } },
    { -6, { "t~", "\\bar{t}" } },
    { 11, { "e-", "e^{-}" } },
    { -11, { "e+", "e^{+}" } },
    { 13, { "mu-", "\\mu^{-}" } },
    { -13, { "mu+", "\\mu^{+}" } },
    { 15, { "tau-", "\\tau^{-}" } },
    { -15, { "tau+", "\\tau^{+}" } },
    { 12, { "ve", "\\nu_{e}" } },
    { -12, { "ve~", "\\bar{\\nu}_{e}" } },
    { 14, { "vm", "\\nu_{\\mu}" } },
    { -14, { "vm~", "\\bar{\\nu}_{\\mu}" } },
    { 16, { "vt", "\\nu_{\\tau}" } },
    { -16, { "vt~", "\\bar{\\nu}_{\\tau}" } },
    { 21, { "g", "g" } },
    { 22, { "a", "\\gamma" } },
    { 23, { "Z", "Z" } },
    { 24, { "W+", "W^{+}" } },
    { -24, { "W-", "W^{-}" } },
    { 25, { "h", "h" } }
};

/**
Name returns a string representation of particle pdg.
*/
std::string PDG::Name(PDGCode pdg) {
    auto x = PDG::names.find(pdg);
    if (x != PDG::names.end()) {
        return x->second.first;
    }
    return Strings::Format("unknown particle (pdg=%d)", (int)pdg);
}

/**
LatexName returns a string representation of particle pdg in Latex notation.
*/
std::string PDG::LatexName(PDGCode pdg) {
    auto x = PDG::names.find(pdg);
    if (x != PDG::names.end()) {
        return x->second.second;
    }
    return Strings::Format("\\text{pdg}=%d", (int)pdg);
}

/**
CodesToName returns a string represenation of the (2 -> n) process.
*/
std::string PDG::CodesToName(const std::vector<int> &codes) {
    std::string r;
    for (size_t i = 0; i < codes.size(); i++) {
        r.append(PDG::Name(codes[i]));
        if (i < codes.size() - 1) {
            r.append(" ");
        }
        if (i == 1) {
            r += "-> ";
        }
    }
    return r;
}

/**
IsQuark returns true if pdg is the pdg number for a quark or antiquark
*/
bool PDG::IsQuark(PDGCode pdg) {
    switch (pdg) {
    case PDG::UQuark:
    case PDG::DQuark:
    case PDG::SQuark:
    case PDG::CQuark:
    case PDG::TQuark:
    case PDG::BQuark:
    case PDG::AntiUQuark:
    case PDG::AntiDQuark:
    case PDG::AntiSQuark:
    case PDG::AntiCQuark:
    case PDG::AntiTQuark:
    case PDG::AntiBQuark:
        return true;
    }
    return false;
}
/**
IsLepton returns true if the particle pdg is a lepton.
*/
bool PDG::IsLepton(PDGCode pdg) {
    return PDG::IsChargedLepton(pdg) || PDG::IsNeutrino(pdg);
}

/**
IsChargedLepton returns true if the particle pdg is a charged lepton (\f$ l^\pm
\f$).
*/
bool PDG::IsChargedLepton(PDGCode pdg) {
    switch (pdg) {
    case PDG::Electron:
    case PDG::Muon:
    case PDG::Tau:
    case PDG::Positron:
    case PDG::AntiMuon:
    case PDG::AntiTau:
        return true;
    }
    return false;
}

bool PDG::IsChargedFermion(PDGCode pdg) {
    if (PDG::IsChargedLepton(pdg)) {
        return true;
    }
    if (PDG::IsQuark(pdg)) {
        return true;
    }
    return false;
}

/**
IsNeutrino returns ture if the particle pdg is a neutrino (\f$ \nu,\bar{\nu}
\f$).
*/
bool PDG::IsNeutrino(PDGCode pdg) {
    switch (pdg) {
    case PDG::ElectronNeutrino:
    case PDG::MuonNeutrino:
    case PDG::TauNeutrino:
    case PDG::AntiElectronNeutrino:
    case PDG::AntiMuonNeutrino:
    case PDG::AntiTauNeutrino:
        return true;
    }
    return false;
}

double PDG::Charge(PDGCode pdg) {
    // charge of leptons and quarks (array index == pdg particle number)
    static double charge[26] = { 0.0,       -1.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0,
                                 2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0, -1.0 / 3.0,
                                 2.0 / 3.0, 0.0,        0.0,       -1.0,
                                 0.0,       -1.0,       0.0,       -1.0,
                                 0.0,       -1.0,       0.0,       0.0,
                                 0.0,       0.0,        0.0,       0.0,
                                 1.0,       0.0 };
    if (pdg > 0 && pdg <= 25) {
        return charge[pdg];
    }
    if( pdg < 0 && pdg >= -25) {
        return -charge[-pdg];
    }
    assert( 0 && "pdg not implemented" );
    return 0.0;
}

double PDG::ChargeAbs(PDGCode pdg) {
    // charge of leptons and quarks (array index == pdg particle number)
    static double charge[26] = { 0.0,       1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0,
                                 2.0 / 3.0, 1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0,
                                 2.0 / 3.0, 0.0,       0.0,       1.0,
                                 0.0,       1.0,       0.0,       1.0,
                                 0.0,       1.0,       0.0,       0.0,
                                 0.0,       0.0,       0.0,       0.0,
                                 1.0,       0.0 };
    int apdg = abs(pdg);
    if (apdg >= 0 && apdg <= 25) {
        return charge[apdg];
    }
    assert( 0 && "pdg not implemented" );
    return 0.0;
}

std::string PDG::Type(PDGCode pdg) {
    if (pdg >= 1 && pdg <= 6) {
        return "q";
    }
    if (pdg <= -1 && pdg >= -6) {
        return "q~";
    }
    if (pdg == 11 || pdg == 13 || pdg == 15) {
        return "l-";
    }
    if (pdg == -11 || pdg == -13 || pdg == -15) {
        return "l+";
    }
    if (pdg == 12 || pdg == 14 || pdg == 16) {
        return "nu";
    }
    if (pdg == -12 || pdg == -14 || pdg == -16) {
        return "nu~";
    }
    return PDG::Name(pdg);
}

} // namespace Physics
