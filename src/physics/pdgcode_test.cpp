#include "gtest/gtest.h"

#include "pdgcode.h"

using namespace Physics;

TEST(Name, Unknown) {
    EXPECT_EQ(PDG::Name(9000001), "unknown particle (pdg=9000001)");
    EXPECT_EQ(PDG::Name(-9000001), "unknown particle (pdg=-9000001)");
}

TEST(Name, Quarks) {
    EXPECT_EQ(PDG::Name(PDG::UQuark), "u");
    EXPECT_EQ(PDG::Name(PDG::DQuark), "d");
    EXPECT_EQ(PDG::Name(PDG::SQuark), "s");
    EXPECT_EQ(PDG::Name(PDG::CQuark), "c");
    EXPECT_EQ(PDG::Name(PDG::BQuark), "b");
    EXPECT_EQ(PDG::Name(PDG::TQuark), "t");
}

TEST(Name, AntiQuarks) {
    EXPECT_EQ(PDG::Name(PDG::AntiUQuark), "u~");
    EXPECT_EQ(PDG::Name(PDG::AntiDQuark), "d~");
    EXPECT_EQ(PDG::Name(PDG::AntiSQuark), "s~");
    EXPECT_EQ(PDG::Name(PDG::AntiCQuark), "c~");
    EXPECT_EQ(PDG::Name(PDG::AntiBQuark), "b~");
    EXPECT_EQ(PDG::Name(PDG::AntiTQuark), "t~");
}

TEST(Name, Leptons) {
    EXPECT_EQ(PDG::Name(PDG::Electron), "e-");
    EXPECT_EQ(PDG::Name(PDG::Muon), "mu-");
    EXPECT_EQ(PDG::Name(PDG::Tau), "tau-");
    EXPECT_EQ(PDG::Name(PDG::ElectronNeutrino), "ve");
    EXPECT_EQ(PDG::Name(PDG::MuonNeutrino), "vm");
    EXPECT_EQ(PDG::Name(PDG::TauNeutrino), "vt");
}

TEST(Name, AntiLeptons) {
    EXPECT_EQ(PDG::Name(PDG::Positron), "e+");
    EXPECT_EQ(PDG::Name(PDG::AntiMuon), "mu+");
    EXPECT_EQ(PDG::Name(PDG::AntiTau), "tau+");
    EXPECT_EQ(PDG::Name(PDG::AntiElectronNeutrino), "ve~");
    EXPECT_EQ(PDG::Name(PDG::AntiMuonNeutrino), "vm~");
    EXPECT_EQ(PDG::Name(PDG::AntiTauNeutrino), "vt~");
}

TEST(LatexName, Unknown) {
    EXPECT_EQ(PDG::LatexName(9000001), "\\text{pdg}=9000001");
    EXPECT_EQ(PDG::LatexName(-9000001), "\\text{pdg}=-9000001");
}

TEST(LatexName, Quarks) {
    EXPECT_EQ(PDG::LatexName(PDG::UQuark), "u");
    EXPECT_EQ(PDG::LatexName(PDG::DQuark), "d");
    EXPECT_EQ(PDG::LatexName(PDG::SQuark), "s");
    EXPECT_EQ(PDG::LatexName(PDG::CQuark), "c");
    EXPECT_EQ(PDG::LatexName(PDG::BQuark), "b");
    EXPECT_EQ(PDG::LatexName(PDG::TQuark), "t");
}

TEST(LatexName, AntiQuarks) {
    EXPECT_EQ(PDG::LatexName(PDG::AntiUQuark), "\\bar{u}");
    EXPECT_EQ(PDG::LatexName(PDG::AntiDQuark), "\\bar{d}");
    EXPECT_EQ(PDG::LatexName(PDG::AntiSQuark), "\\bar{s}");
    EXPECT_EQ(PDG::LatexName(PDG::AntiCQuark), "\\bar{c}");
    EXPECT_EQ(PDG::LatexName(PDG::AntiBQuark), "\\bar{b}");
    EXPECT_EQ(PDG::LatexName(PDG::AntiTQuark), "\\bar{t}");
}

TEST(LatexName, Leptons) {
    EXPECT_EQ(PDG::LatexName(PDG::Electron), "e^{-}");
    EXPECT_EQ(PDG::LatexName(PDG::Muon), "\\mu^{-}");
    EXPECT_EQ(PDG::LatexName(PDG::Tau), "\\tau^{-}");
    EXPECT_EQ(PDG::LatexName(PDG::ElectronNeutrino),
              "\\nu_{e}");
    EXPECT_EQ(PDG::LatexName(PDG::MuonNeutrino),
              "\\nu_{\\mu}");
    EXPECT_EQ(PDG::LatexName(PDG::TauNeutrino),
              "\\nu_{\\tau}");
}

TEST(LatexName, AntiLeptons) {
    EXPECT_EQ(PDG::LatexName(PDG::Positron), "e^{+}");
    EXPECT_EQ(PDG::LatexName(PDG::AntiMuon), "\\mu^{+}");
    EXPECT_EQ(PDG::LatexName(PDG::AntiTau), "\\tau^{+}");
    EXPECT_EQ(PDG::LatexName(PDG::AntiElectronNeutrino),
              "\\bar{\\nu}_{e}");
    EXPECT_EQ(PDG::LatexName(PDG::AntiMuonNeutrino),
              "\\bar{\\nu}_{\\mu}");
    EXPECT_EQ(PDG::LatexName(PDG::AntiTauNeutrino),
              "\\bar{\\nu}_{\\tau}");
}

TEST(IsQuark, Quark) {
    EXPECT_TRUE(PDG::IsQuark(PDG::UQuark));
    EXPECT_TRUE(PDG::IsQuark(PDG::SQuark));
    EXPECT_TRUE(PDG::IsQuark(PDG::CQuark));
    EXPECT_TRUE(PDG::IsQuark(PDG::DQuark));
    EXPECT_TRUE(PDG::IsQuark(PDG::BQuark));
    EXPECT_TRUE(PDG::IsQuark(PDG::TQuark));
    // anti quarks
    EXPECT_TRUE(PDG::IsQuark(PDG::AntiUQuark));
    EXPECT_TRUE(PDG::IsQuark(PDG::AntiSQuark));
    EXPECT_TRUE(PDG::IsQuark(PDG::AntiCQuark));
    EXPECT_TRUE(PDG::IsQuark(PDG::AntiDQuark));
    EXPECT_TRUE(PDG::IsQuark(PDG::AntiBQuark));
    EXPECT_TRUE(PDG::IsQuark(PDG::AntiTQuark));
}

TEST(IsQuark, NoQuarks) {
    EXPECT_FALSE(PDG::IsQuark(PDG::Gluon));
    EXPECT_FALSE(PDG::IsQuark(PDG::Electron));
    EXPECT_FALSE(PDG::IsQuark(PDG::Photon));
    EXPECT_FALSE(PDG::IsQuark(PDG::Higgs));
}

TEST(IsChargedLepton, ChargedLeptons) {
    EXPECT_TRUE(PDG::IsChargedLepton(PDG::Electron));
    EXPECT_TRUE(PDG::IsChargedLepton(PDG::Positron));
    EXPECT_TRUE(PDG::IsChargedLepton(PDG::Muon));
    EXPECT_TRUE(PDG::IsChargedLepton(PDG::AntiMuon));
    EXPECT_TRUE(PDG::IsChargedLepton(PDG::AntiTau));
    EXPECT_TRUE(PDG::IsChargedLepton(PDG::Tau));
}

TEST(IsChargedLepton, NoChargedLeptons) {
    EXPECT_FALSE(PDG::IsChargedLepton(PDG::Gluon));
    EXPECT_FALSE(PDG::IsChargedLepton(PDG::ElectronNeutrino));
    EXPECT_FALSE(PDG::IsChargedLepton(PDG::Photon));
    EXPECT_FALSE(PDG::IsChargedLepton(PDG::Higgs));
}

TEST(IsNeutrino, Neutrinos) {
    EXPECT_TRUE(PDG::IsNeutrino(PDG::ElectronNeutrino));
    EXPECT_TRUE(PDG::IsNeutrino(PDG::AntiElectronNeutrino));
    EXPECT_TRUE(PDG::IsNeutrino(PDG::MuonNeutrino));
    EXPECT_TRUE(PDG::IsNeutrino(PDG::AntiMuonNeutrino));
    EXPECT_TRUE(PDG::IsNeutrino(PDG::AntiTauNeutrino));
    EXPECT_TRUE(PDG::IsNeutrino(PDG::TauNeutrino));
}

TEST(IsNeutrino, NoNeutrinos) {
    EXPECT_FALSE(PDG::IsNeutrino(PDG::Gluon));
    EXPECT_FALSE(PDG::IsNeutrino(PDG::Electron));
    EXPECT_FALSE(PDG::IsNeutrino(PDG::Photon));
    EXPECT_FALSE(PDG::IsNeutrino(PDG::Higgs));
}

TEST(IsLepton, Leptons) {
    EXPECT_TRUE(PDG::IsLepton(PDG::Electron));
    EXPECT_TRUE(PDG::IsLepton(PDG::Positron));
    EXPECT_TRUE(PDG::IsLepton(PDG::Muon));
    EXPECT_TRUE(PDG::IsLepton(PDG::AntiMuon));
    EXPECT_TRUE(PDG::IsLepton(PDG::AntiTau));
    EXPECT_TRUE(PDG::IsLepton(PDG::Tau));
    EXPECT_TRUE(PDG::IsLepton(PDG::ElectronNeutrino));
    EXPECT_TRUE(PDG::IsLepton(PDG::AntiElectronNeutrino));
    EXPECT_TRUE(PDG::IsLepton(PDG::MuonNeutrino));
    EXPECT_TRUE(PDG::IsLepton(PDG::AntiMuonNeutrino));
    EXPECT_TRUE(PDG::IsLepton(PDG::AntiTauNeutrino));
    EXPECT_TRUE(PDG::IsLepton(PDG::TauNeutrino));
}

TEST(IsLepton, NoLeptons) {
    EXPECT_FALSE(PDG::IsLepton(PDG::Gluon));
    EXPECT_FALSE(PDG::IsLepton(PDG::Photon));
    EXPECT_FALSE(PDG::IsLepton(PDG::Higgs));
}

TEST(Charge, Quarks) {
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::UQuark), 2.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::CQuark), 2.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::TQuark), 2.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::DQuark), -1.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::SQuark), -1.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::BQuark), -1.0 / 3.0);
}

TEST(Charge, AntiQuarks) {
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::AntiUQuark), -2.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::AntiCQuark), -2.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::AntiTQuark), -2.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::AntiDQuark), 1.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::AntiSQuark), 1.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::AntiBQuark), 1.0 / 3.0);
}

TEST(Charge, ChargedLeptons) {
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::Electron), -1.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::Muon), -1.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::Tau), -1.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::Positron), 1.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::AntiMuon), 1.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::AntiTau), 1.0);
}

TEST(Charge, Neutrinos) {
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::ElectronNeutrino), 0.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::AntiElectronNeutrino), 0.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::MuonNeutrino), 0.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::AntiMuonNeutrino), 0.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::TauNeutrino), 0.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::AntiTauNeutrino), 0.0);
}

TEST(Charge, GaugeBosons) {
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::Photon), 0.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::Gluon), 0.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::Z), 0.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::WPlus), 1.0);
    EXPECT_DOUBLE_EQ(PDG::Charge(PDG::WMinus), -1.0);
}

TEST(Charge, Higgs) { EXPECT_DOUBLE_EQ(PDG::Charge(PDG::Higgs), 0.0); }

TEST(ChargeAbs, Quarks) {
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::UQuark), 2.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::CQuark), 2.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::TQuark), 2.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::DQuark), 1.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::SQuark), 1.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::BQuark), 1.0 / 3.0);
}

TEST(ChargeAbs, AntiQuarks) {
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::AntiUQuark), 2.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::AntiCQuark), 2.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::AntiTQuark), 2.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::AntiDQuark), 1.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::AntiSQuark), 1.0 / 3.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::AntiBQuark), 1.0 / 3.0);
}

TEST(ChargeAbs, ChargeAbsdLeptons) {
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::Electron), 1.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::Muon), 1.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::Tau), 1.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::Positron), 1.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::AntiMuon), 1.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::AntiTau), 1.0);
}

TEST(ChargeAbs, Neutrinos) {
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::ElectronNeutrino), 0.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::AntiElectronNeutrino), 0.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::MuonNeutrino), 0.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::AntiMuonNeutrino), 0.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::TauNeutrino), 0.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::AntiTauNeutrino), 0.0);
}

TEST(ChargeAbs, GaugeBosons) {
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::Photon), 0.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::Gluon), 0.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::Z), 0.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::WPlus), 1.0);
    EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::WMinus), 1.0);
}

TEST(ChargeAbs, Higgs) { EXPECT_DOUBLE_EQ(PDG::ChargeAbs(PDG::Higgs), 0.0); }
