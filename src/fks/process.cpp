#include "fks/process.h"
#include "phasespace/phasespace.h"
#include "physics/pdgcode.h"
#include <cassert>
#include <cassert>
#include <iostream>

using namespace FKS;

namespace {
bool IsValidISRSplitting(Type_t type, int born, int real, int rad) {
    auto isquark = [](int pdg) -> bool {
        return pdg != 0 && pdg > -6 && pdg < 6;
    };
    auto isgluon = [](int pdg) -> bool { return pdg == 0 || pdg == 21; };
    auto isphoton = [](int pdg) -> bool { return pdg == 22; };
    bool isquark_born = isquark(born);
    bool isquark_real = isquark(real);
    bool isquark_rad = isquark(rad);
    bool isg_born = (type == Type_t::QCD) ? isgluon(born) : isphoton(born);
    bool isg_real = (type == Type_t::QCD) ? isgluon(real) : isphoton(real);
    bool isg_rad = (type == Type_t::QCD) ? isgluon(rad) : isphoton(rad);
    if (isquark_real && isquark_rad && isg_born && real == rad) {
        return true;
    }
    if (isquark_real && isquark_born && isg_rad && real == born) {
        return true;
    }
    if (isg_real && isquark_born && isquark_rad && rad == -born) {
        return true;
    }
    if (type == Type_t::QCD && isg_real && isg_born && isg_rad) {
        return true;
    }
    return false;
}
} // namespace

static bool is_equvivalent(const Real_t &a, const Real_t &b) {
    if (a.ID != b.ID) {
        return false;
    }
    if (a.AllFlavours.size() != b.AllFlavours.size()) {
        return false;
    }
    if (a.Regions.size() != b.Regions.size()) {
        return false;
    }
    for (size_t i = 0; i < a.AllFlavours.size(); i++) {
        if (a.AllFlavours[i][0] != b.AllFlavours[i][0]) {
            return false;
        }
        if (a.AllFlavours[i][1] != b.AllFlavours[i][1]) {
            return false;
        }
    }
    for (size_t i = 0; i < a.Regions.size(); i++) {
        if (a.Regions[i].I != b.Regions[i].I) {
            return false;
        }
        if (a.Regions[i].J != b.Regions[i].J) {
            return false;
        }
    }
    return true;
}

void FlavourConfig::Print() const {
    printf("# born:\n");
    printf("%s\n", Physics::PDG::CodesToName(Born.Flavours).c_str());
    printf("pdfs:\n");
    for (const auto &pdf : Born.PDF) {
        printf("* %s", Physics::PDG::Name(pdf[0]).c_str());
        printf(" %s\n", Physics::PDG::Name(pdf[1]).c_str());
    }
    printf("## real:\n");
    for (const auto &real : Real) {
        printf("%s\n", Physics::PDG::CodesToName(real.Flavours).c_str());
        printf("  pdfs:\n");
        for (const auto &pdf : real.PDF) {
            printf("  * %s", Physics::PDG::Name(pdf[0]).c_str());
            printf(" %s\n", Physics::PDG::Name(pdf[1]).c_str());
        }
        printf("  regions:\n");
        for (const auto &region : real.Regions) {
            printf("  * I = %d, J = %d\n", region.I, region.J);
        }
    }
    printf("## collinear remnants\n");
    printf("splitting of parton 1\n");
    for (const auto &remn : Remnant1) {
        int btype = remn.Splitting.BornPDG;
        int rtype = remn.Splitting.RealPDG;
        printf("  * %s", Physics::PDG::Type(rtype).c_str());
        printf("-> %s", Physics::PDG::Type(btype).c_str());
        printf(" pdfs =");
        for (size_t i = 0; i < remn.PDF.size(); i += 2) {
            printf("%d %d, ", remn.PDF[i], remn.PDF[i + 1]);
        }
        printf("\n");
    }
    printf("splitting of parton 2\n");
    for (const auto &remn : Remnant2) {
        int btype = remn.Splitting.BornPDG;
        int rtype = remn.Splitting.RealPDG;
        printf("  * %s", Physics::PDG::Type(rtype).c_str());
        printf("-> %s", Physics::PDG::Type(btype).c_str());
        printf(" pdfs =");
        for (size_t i = 0; i < remn.PDF.size(); i += 2) {
            printf("%d %d, ", remn.PDF[i], remn.PDF[i + 1]);
        }
        printf("\n");
    }
}

FlavourConfig::FlavourConfig(int id, const PDGList &pdgs,
                             const std::vector<int> &pdfs,
                             const ColorFlow &color1, const ColorFlow &color2,
                             const std::vector<double> &scales) {
    Born.ID = id;
    Born.Flavours = pdgs;
    Born.PDF.reserve(pdfs.size() / 2);
    for (size_t i = 0; i < pdfs.size(); i += 2) {
        Born.PDF.push_back({{pdfs[i], pdfs[i + 1]}});
    }
    Born.AllFlavours = std::vector<PDGList>(Born.PDF.size(), pdgs);
    for (size_t i = 0; i < Born.AllFlavours.size(); i++) {
        Born.AllFlavours[i][0] = Born.PDF[i][0];
        Born.AllFlavours[i][1] = Born.PDF[i][1];
    }
    if (scales.size() == 0) {
        Scales = std::vector<double>(Born.PDF.size(), 1.0);
    } else {
        assert(Born.PDF.size() == scales.size() &&
               "need a scale for every pdf");
        Scales = scales;
    }
    assert(color1.size() == Born.Flavours.size() &&
           "need a color for every particle");
    assert(color2.size() == Born.Flavours.size() &&
           "need a color for every particle");
    Born.Color[0] = color1;
    Born.Color[1] = color2;
}

FlavourConfig::FlavourConfig(int id, const std::vector<PDGList> &pdgs,
                             int dummy, const ColorFlow &color1,
                             const ColorFlow &color2,
                             const std::vector<double> &scales) {
    Born.ID = id;
    Born.Flavours = pdgs[0];
    Born.AllFlavours = pdgs;
    for (const auto &pdg : pdgs) {
        Born.PDF.push_back({{pdg[0], pdg[1]}});
    }

    if (scales.size() == 0) {
        Scales = std::vector<double>(Born.PDF.size(), 1.0);
    } else {
        assert(Born.PDF.size() == scales.size() &&
               "need a scale for every pdf");
        Scales = scales;
    }
    assert(color1.size() == Born.Flavours.size() &&
           "need a color for every particle");
    assert(color2.size() == Born.Flavours.size() &&
           "need a color for every particle");
    Born.Color[0] = color1;
    Born.Color[1] = color2;
}

void FlavourConfig::AddRealDY(int id, Type_t type, const PDGList &pdgs,
                              const std::vector<int> &pdfs,
                              const RegionList &regions,
                              const ResonanceList &resonances) {
    Real_t real;
    real.ID = id;
    real.Type = type;
    if (type == FKS::Type_t::QCD) {
        QCD = true;
    }
    if (type == FKS::Type_t::EW) {
        EW = true;
    }
    real.Flavours = pdgs;

    for (size_t i = 0; i < pdfs.size(); i += 2) {
        real.PDF.push_back({{pdfs[i], pdfs[i + 1]}});
    }
    std::vector<PDGList> all(real.PDF.size(), pdgs);
    for (size_t i = 0; i < all.size(); i++) {
        all[i][0] = real.PDF[i][0];
        all[i][1] = real.PDF[i][1];
        auto rad_index = pdgs.size() - 1;
        int rad = all[i][rad_index];
        if (rad != 0 && rad > -6 && rad < 6) {
            if (all[i][0] == 0 || all[i][0] == 21) {
                all[i][rad_index] = -Born.AllFlavours[i][0];
            }
            if (all[i][1] == 0 || all[i][1] == 21) {
                all[i][rad_index] = -Born.AllFlavours[i][1];
            }
        }
    }
    real.Regions = regions;
    real.AllRegions = regions;
    real.Resonances = resonances;
    real.AllFlavours.push_back(all);
    Real.push_back(real);

    int pdg_rad = pdgs[pdgs.size() - 1];
    int born0 = Born.Flavours[0];
    int real0 = pdgs[0];

    if (IsValidISRSplitting(type, born0, real0, pdg_rad)) {
        AddRemnantInitial1(type, FKS::Splitting(real0, born0), pdfs);
    }
    int born1 = Born.Flavours[1];
    int real1 = pdgs[1];
    if (IsValidISRSplitting(type, born1, real1, pdg_rad)) {
        AddRemnantInitial2(type, FKS::Splitting(real1, born1), pdfs);
    }
}

void FlavourConfig::AddReal(int id, Type_t type,
                            const std::vector<PDGList> &pdgs,
                            const RegionList &regions,
                            const RegionList &allregions,
                            const ResonanceList &resonances) {
    Real_t real;
    real.ID = id;
    real.Type = type;
    if (type == FKS::Type_t::QCD) {
        QCD = true;
    }
    if (type == FKS::Type_t::EW) {
        EW = true;
    }
    real.AllFlavours.push_back(pdgs);
    real.Flavours = pdgs[0];
    std::vector<int> pdfs;
    pdfs.reserve(2 * pdgs.size());
    for (size_t i = 0; i < pdgs.size(); i++) {
        real.PDF.push_back({{pdgs[i][0], pdgs[i][1]}});
        pdfs.push_back(pdgs[i][0]);
        pdfs.push_back(pdgs[i][1]);
    }
    assert(regions.size() > 0);
    real.Regions = regions;
    real.AllRegions = allregions;
    real.Resonances = resonances;
    bool inserted = false;
    for (auto &r : Real) {
        if (is_equvivalent(real, r)) {
            r.AllFlavours.push_back(pdgs);
            inserted = true;
            break;
        }
    }
    if (!inserted) {
        Real.push_back(real);
    }

    for (const auto &reg : regions) {
        if (reg.J < 2 && reg.I != 4) {
            int pdg_rad = real.Flavours[reg.I];
            int born0 = Born.Flavours[0];
            int real0 = real.Flavours[0];
            if (reg.J == 0 || reg.J == -1) {
                if (IsValidISRSplitting(type, born0, real0, pdg_rad)) {
                    AddRemnantInitial1(type, FKS::Splitting(real0, born0),
                                       pdfs);
                }
            }
            int born1 = Born.Flavours[1];
            int real1 = real.Flavours[1];
            if (reg.J == 0 || reg.J == -2) {
                if (IsValidISRSplitting(type, born1, real1, pdg_rad)) {
                    AddRemnantInitial2(type, FKS::Splitting(real1, born1),
                                       pdfs);
                }
            }
        }
    }
}

double Resonance::BreitWigner(const Phasespace::Phasespace &ps,
                              bool useCharge) const {
    if (Disabled) {
        return 1.0;
    }
    size_t len = Daughters.size();
    assert(len >= 2 && "need at least two daughters");
    auto p = ps.Momenta[Daughters[0]];

    for (size_t i = 1; i < len; i++) {
        p = p.Plus(ps.Momenta[Daughters[i]]);
    }
    double s = p.Dot(p);

    double M2 = Mass * Mass;
    double G2 = Width * Width;
    assert(Mass > 0.0);
    assert(Width > 0.0);

    double tmp = (s / M2 - 1.0);
    double res = 1.0 / (tmp * tmp + G2 / M2);
    if (useCharge) {
        res *= Q2;
    }
    return res;
}
