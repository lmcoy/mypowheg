#include "fks/phasespaces.h"

#include "phasespace/recombination.h"
#include "phasespace/realphasespace.h"
#include "process/data.h"

using namespace FKS;

namespace {

Phasespace::Phasespace recombine_lab(FKS::Type_t type, const int *pdgs,
                                     const Phasespace::Phasespace &ps_real,
                                     double dR) {
    Phasespace::Phasespace ps_recombined;
    if (type == FKS::Type_t::EW) {
        ps_recombined = Phasespace::Recombine(pdgs, ps_real, dR);
    } else {
        ps_recombined = ps_real;
    }

    Phasespace::Phasespace ps_recombined_lab;
    ps_recombined_lab.SetToLabFromCMS(&ps_recombined);
    return ps_recombined_lab;
}

} // end namespace

void Phasespaces::Generate(const Phasespace::Phasespace &ps, int ifks, int jfks,
                              double xi, double y, double phi) {

    double xi_soft = 1e-10;
    double y1_coll = 0.999999;
    double y2_coll = -y1_coll;
    Phasespace::GenRealPhasespace(&Real, &ps, ifks, jfks, xi, y, phi);
    Phasespace::GenRealPhasespace(&Soft, &ps, ifks, jfks, xi_soft, y, phi);
    Phasespace::GenRealPhasespace(&Collinear1, &ps, ifks, jfks, xi, y1_coll,
                                  phi);
    Phasespace::GenRealPhasespace(&Collinear2, &ps, ifks, jfks, xi, y2_coll,
                                  phi);
    Born = ps;
}

Phasespaces Phasespaces::Recombined(FKS::Type_t type, const int *pdgs, 
                         const UserProcess::RecombinationParam &recomb) const {
    Phasespaces result;
    // recombine lepton & photon if they are close to each other
    result.Real = recombine_lab(type, pdgs, Real, recomb.dR);
    result.Collinear1 = recombine_lab(type, pdgs, Collinear1, recomb.dR);
    result.Collinear2 = recombine_lab(type, pdgs, Collinear2, recomb.dR);

    return result;
}
