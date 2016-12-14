#ifndef FINDNORM_H_CYKTF2PD
#define FINDNORM_H_CYKTF2PD

#include "fks/process.h"
#include "phasespace/phasespace.h"
#include "process/data.h"
#include "process/cuts.h"
#include "powheg/upperbounding.h"

namespace Powheg {

static int findNorm(bool search_bmax,
                    const Phasespace::Phasespace &ps, double x1, double x2,
                    double x3, double wgt, double *result,
                    UserProcess::Data *userdata) {

    Phasespace::Phasespace ps_lab;
    ps_lab.SetToLabFromCMS(&ps);
    userdata->PowhegScales.mu = userdata->Scales->Renorm(ps);
    userdata->PowhegScales.muF = userdata->Scales->Factorization(ps);
    userdata->PowhegScales.Q2 =
        userdata->PowhegScales.mu * userdata->PowhegScales.mu;
    size_t failed_cuts = 0;
    if (userdata->RadiationRegions.size() == 0) {
        *result = 0.0;
        return 0;
    }
    for (auto & radreg : userdata->RadiationRegions) {
        const int * pdgs = radreg.FlavourConfig->Born.Flavours.data();
        bool passcuts =
            userdata->cuts->ApplyCuts(ps.N + 2, pdgs, ps_lab.Momenta.data());
        if (!passcuts) {
            failed_cuts += 1;
            continue;
        }
        FKS::Type_t type = radreg.Type; 
        Powheg::GetNormForBoundingFcts(type, search_bmax, radreg, ps, x1, x2,
                                       x3, userdata);

    }
    *result = 0.0;
    auto s = userdata->RadiationRegions.size();
    if (s > 1 && failed_cuts >= s / 2) {
        // if we sample too many events that don't pass cuts, we have only very
        // few events in our Norm histograms. Therefore, we draw an additional
        // phase space point when cuts were not passed for too many radiation
        // regions.
        return -1;
    }
    return 0;
}

int InitFindNorm(const Phasespace::Phasespace &ps, double x1, double x2, double x3,
             double wgt, double *result, UserProcess::Data *userdata) {
    return findNorm(true, ps, x1, x2, x3, wgt, result, userdata);
}

int FindNorm(const Phasespace::Phasespace &ps, double x1, double x2, double x3,
             double wgt, double *result, UserProcess::Data *userdata) {
    return findNorm(false, ps, x1, x2, x3, wgt, result, userdata);
}


} // end namespace Powheg

#endif /* end of include guard: FINDNORM_H_CYKTF2PD */ 
