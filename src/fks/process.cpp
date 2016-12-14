#include "fks/process.h"
#include <cassert>
#include <iostream>
#include <cassert>
#include "physics/pdgcode.h"
#include "phasespace/phasespace.h"

using namespace FKS;

namespace {
bool IsValidISRSplitting(Type_t type, int born, int real, int rad) {
    auto isquark = [](int pdg)->bool {
        return pdg != 0 && pdg > -6 && pdg < 6;
    };
    auto isgluon = [](int pdg)->bool {
        return pdg == 0 || pdg == 21;
    };
    auto isphoton = [](int pdg)->bool {
        return pdg == 22;
    };
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
    for (const auto &remn: Remnant1) {
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
    for (const auto &remn: Remnant2) {
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
                             const std::vector<double> &scales) {
    Born.ID = id;
    Born.Flavours = pdgs;
    Born.PDF.reserve(pdfs.size() / 2);
    for (size_t i = 0; i < pdfs.size(); i += 2) {
        Born.PDF.push_back({ { pdfs[i], pdfs[i + 1] } });
    }
    if( scales.size() == 0 ) {
        Scales = std::vector<double>(Born.PDF.size(), 1.0);
    } else {
        assert(Born.PDF.size() == scales.size() &&
               "need a scale for every pdf");
        Scales = scales;
    }
}

void FlavourConfig::AddReal(int id, Type_t type, const PDGList &pdgs,
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
        real.PDF.push_back({ { pdfs[i], pdfs[i + 1] } });
    }
    real.Regions = regions;
    real.Resonances = resonances;
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
