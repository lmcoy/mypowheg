#include "generateevents.h"

#include <iostream>
#include <map>

#include "fks/param.h"
#include "fks/xsec.h"
#include "lhe/event.h"
#include "phasespace/phasespace.h"
#include "phasespace/realphasespace.h"
#include "powheg/generateradiation.h"
#include "powheg/pickelement.h"
#include "powheg/reshuffle_momenta.h"
#include "powheg/resonance.h"
#include "powheg/unweighting.h"
#include "process/data.h"
#include "process/matrixelement.h"
#include "random/rnd.h"
#include "color.h"

#define EVENT_REJECT 1
#define EVENT_SUCCESS 2
#define EVENT_ERROR 3
#define EVENT_NEG 4

namespace Powheg {

namespace {

enum class Type { QCD, QED };

int GenRadiation(const BornConfig &bornconfig, const Phasespace::Phasespace &ps,
                 double wgt, UserProcess::Data *userdata,
                 Phasespace::Phasespace *ps_out, int nn, int *pdg_out,
                 Radiation *rad_out) {
    assert(nn >= ps.N + 2);

    auto pdf = bornconfig.PDFindex;
    auto fl = &userdata->Process[bornconfig.ProcessIndex];
    Radiation radiation;
    double kT2minG = userdata->RadiationParameter.kT2min;
    radiation.kT2 = 0.0;
    // TODO: Set alpha_s for born. alpha_s should be arbitrary because we use
    // Params_as in R and B. Therefore, it cancels in R/B except the radiation
    // alpha. This is set later.
    auto bme = userdata->MatrixElement->Born(fl->Born.ID, ps,
                                             userdata->PowhegScales.mu,
                                             userdata->Params, false, false);
    double B = bme.M2;
    bool isborn = true;
    int n_used_radreg = 0;
    for (const auto &r : userdata->RadiationRegions) {
        if (fl->Born.ID != r.FlavourConfig->Born.ID) {
            continue;
        }
        Radiation rad;
        RadiationType type = RadiationType::BORN;
        if (!userdata->RadiatePhoton && r.Type == FKS::Type_t::EW) {
            continue;
        }
        if (!userdata->RadiateQCD && r.Type == FKS::Type_t::QCD) {
            continue;
        }
        n_used_radreg += 1;

        double kT2min = 1.0;
        switch (r.Type) {
        case FKS::Type_t::QCD:
            kT2min = std::max(kT2minG, radiation.kT2);
            type = Powheg::QCD::GenerateRadiation(
                pdf, B, r, ps, kT2min, userdata, &userdata->rng, &rad);
            break;
        case FKS::Type_t::EW:
            kT2min = std::max(kT2minG, radiation.kT2);
            type = Powheg::QED::GenerateRadiation(
                pdf, B, r, ps, kT2min, userdata, &userdata->rng, &rad);
            break;
        }

        if (type == RadiationType::ERROR) {
            userdata->GenEventStatistics.N_RADERROR += 1;
            return EVENT_ERROR;
        }
        if (type == RadiationType::ENORM) {
            userdata->GenEventStatistics.N_ENORM += 1;
            return EVENT_ERROR;
        }
        if (type == RadiationType::ERADVAR) {
            userdata->GenEventStatistics.N_ERADVAR += 1;
            return EVENT_ERROR;
        }

        if (type == RadiationType::BORN) {
            continue;
        }
        if (rad.kT2 > radiation.kT2) {
            radiation = rad;
            isborn = false;
        }
    }

    *rad_out = radiation;
    // born event
    if (isborn) {
        ps_out->SetToLabFromCMS(&ps);
        pdg_out[0] = fl->Born.PDF[pdf][0];
        pdg_out[1] = fl->Born.PDF[pdf][1];
        if (pdg_out[0] == 0) {
            pdg_out[0] = 21;
        }
        if (pdg_out[1] == 0) {
            pdg_out[1] = 21;
        }
        for (int i = 0; i < ps.N; i++) {
            pdg_out[i + 2] = fl->Born.AllFlavours[pdf][i + 2];
        }
        userdata->GenEventStatistics.N_BORN += 1;
        if (userdata->RadiationRegions.size() > 0 && n_used_radreg > 0) {
            // tried to generate radiation down to scale kT2minG but didn't
            // radiate. Therefore, the scale for the parton shower has to be
            // kT2minG.
            rad_out->kT2 = kT2minG;
        }
        return EVENT_SUCCESS;
    }

    Phasespace::Phasespace ps_real;
    Phasespace::GenRealPhasespace(&ps_real, &ps, radiation.i, radiation.j,
                                  radiation.xi, radiation.y, radiation.phi);
    ps_out->SetToLabFromCMS(&ps_real);
    pdg_out[0] = radiation.Real->PDF[pdf][0];
    pdg_out[1] = radiation.Real->PDF[pdf][1];
    if (pdg_out[0] == 0) {
        pdg_out[0] = 21;
    }
    if (pdg_out[1] == 0) {
        pdg_out[1] = 21;
    }
    size_t n_r = radiation.Real->AllFlavours.size();
    int rind = 0;
    if (n_r > 1) {
        double rn[n_r];
        for (size_t i = 0; i < n_r; i++) {
            rn[i] = 1.0;
        }
        rind = pick_element(n_r, rn, userdata->rng.Random());
    }
    for (int i = 0; i < ps_out->N; i++) {
        pdg_out[i + 2] = radiation.Real->AllFlavours[rind][pdf][i + 2];
    }
    userdata->GenEventStatistics.N_REAL += 1;
    return EVENT_SUCCESS;
}

void make_event(const Phasespace::Phasespace &ps_out, const Radiation &rad,
                int *pdgs, UserProcess::Data *params,
                const BornConfig &bornconfig) {
    LHE::Event event;
    double kT = sqrt(rad.kT2);
    event.Scale = kT;
    int N = ps_out.N + 2;
    if (kT == 0.0) {
        event.Scale = params->PowhegScales.muF;
    }
    event.Alpha = params->Params->alpha;
    double aS = params->AlphaS->AlphaS(event.Scale * event.Scale);
    event.AlphaS = aS;
    event.Weight = 1.0;
    if (bornconfig.negative) {
        event.Weight *= -1.0;
    }
    event.Weight *= params->EventWeight;
    event.ID = 0;

    event.N = N;
    Resonance resonance = params->Resonance;

    int color1[10] = {0};
    int color2[10] = {0};
    int pdgs_b[10];
    auto fl = &params->Process[bornconfig.ProcessIndex];
    size_t n_born = fl->Born.Flavours.size();
    assert(n_born <= 10);
    auto pdf = bornconfig.PDFindex;
    for(size_t i = 0; i < n_born; i++) {
        color1[i] = fl->Born.Color[0][i];
        color2[i] = fl->Born.Color[1][i];
        pdgs_b[i] = fl->Born.AllFlavours[pdf][i];
    }

    Powheg::Color(rad.i, rad.j, pdgs_b, pdgs, color1, color2, rad.y);

    // write initial state partons
    for (int i = 0; i < 2; i++) {
        auto &particle = event.Particles[i];
        particle.PDG = pdgs[i];
        particle.Status = -1;
        particle.Color1 = color1[i];
        particle.Color2 = color2[i];
        particle.Momentum = ps_out.Momenta[i];
    }
    if (rad.j >= 2) {
        // test if radiation from resonance daughter
        int mother = rad.j;
        int rpar = rad.i;
        if(rad.i < rad.j) {
            mother = rad.i;
            rpar = rad.j;
        }
        if(mother == resonance.ID[0] || mother == resonance.ID[1]) {
            resonance.ID[2] = rpar;
        }
    }
    // write resonance
    int res = 0;
    bool has_resonance = false;
    if (resonance.ID[0] != 0) {
        has_resonance = true;
        event.N += 1;
        auto &particle = event.Particles[2];
        int i = resonance.ID[0];
        int i2 = resonance.ID[1];
        int i3 = resonance.ID[2];
        particle.PDG = resonance.pdg;
        particle.Status = 2;
        particle.Mother1 = 1;
        particle.Mother2 = 1;
        const Math::FourMomentum &P1 = ps_out.Momenta[i];
        const Math::FourMomentum &P2 = ps_out.Momenta[i2];
        particle.Momentum = P1.Plus(P2);
        if (i3 != 0) {
            particle.Momentum.Add(ps_out.Momenta[i3]);
        }
        res = 1;
    }
    // write final state particles
    for (int i = 2; i < N; i++) {
        int lhe_index = i + res;
        auto &particle = event.Particles[lhe_index];
        particle.PDG = pdgs[i];
        particle.Status = 1;
        particle.Color1 = color1[i];
        particle.Color2 = color2[i];
        particle.Momentum = ps_out.Momenta[i];
        if (i == resonance.ID[0] || i== resonance.ID[1] || i == resonance.ID[2] ) {
                particle.Mother1 = 3;
                particle.Mother2 = 3;
        } else {
                particle.Mother1 = 1;
                particle.Mother2 = 2;
        }
    }

    if (has_resonance) {
        reshuffle_momenta(&event, resonance, params);
    }

    auto ret = params->EventBuffer.Append(event);
    switch (ret) {
    case decltype(ret)::SUCCESS:
        break;
    case decltype(ret)::FULL:
        std::cerr << "warning: event buffer full\n";
        break;
    case decltype(ret)::MALLOCERROR:
        std::cerr << "warning: problem in remalloc. event not written "
                     "to buffer\n";
        break;
    case decltype(ret)::INTERNALERROR:
        std::cerr << "error: internal buffer error\n";
        exit(1);
    }
}

} // end namespace

int GenerateEvents(const Phasespace::Phasespace &ps, double x1, double x2,
                   double x3, double wgt, double *ff,
                   UserProcess::Data *params) {
    params->GenEventStatistics.N += 1;
    auto bornconfig = unweighting(ps, x1, x2, x3, wgt, params);
    if (bornconfig.status == BornConfig::Status::NegativeBtilde) {
        params->GenEventStatistics.N_NEG += 1;
        *ff = 0.0;
        return -1;
    } else if (bornconfig.status == BornConfig::Status::Rejected) {
        params->GenEventStatistics.N_REJECT += 1;
        *ff = 0.0;
        return -1;
    } else if (bornconfig.status ==
               BornConfig::Status::RejectedWithGuessedVirtual) {
        params->GenEventStatistics.N_REJECTWOVIRTUAL += 1;
        *ff = 0.0;
        return -1;
    } else if (bornconfig.status ==
               BornConfig::Status::RejectedWithBorn) {
        params->GenEventStatistics.N_REJECTBORN += 1;
        *ff = 0.0;
        return -1;
    } else if (bornconfig.status == BornConfig::Status::UnderestimatedVirtual) {
        std::cerr << "warning: underestimated virtual matrix element. Should "
                     "not happen too often!\n";
        params->GenEventStatistics.N_WRONGV += 1;
        *ff = 0.0;
        return -1;
    }

    Phasespace::Phasespace ps_out;
    int pdgs[10];
    // set scales
    params->PowhegScales.mu = params->Scales->Renorm(ps);
    params->PowhegScales.muF = params->Scales->Factorization(ps);
    params->PowhegScales.Q2 = params->PowhegScales.mu * params->PowhegScales.mu;

    if (bornconfig.n > params->Unweighting.MaxMultiple) {
        params->GenEventStatistics.N_MAX += 1;
        std::cerr << "warning: maximum for unweighting is too small. btilde ~ "
                  << bornconfig.n << " * Max\n";
        *ff = 0.0;
        return -1;
    }
    int n = 0;
    for (int i = 0; i < bornconfig.n; i++) {
        Radiation rad;
        int status = Powheg::GenRadiation(bornconfig, ps, wgt, params, &ps_out,
                                          10, pdgs, &rad);
        if (status == EVENT_SUCCESS) {
            make_event(ps_out, rad, pdgs, params, bornconfig);
            n += 1;
        }
    }
    if (n >= 1) {
        // at least one event
        switch (n) {
        case 1:
            break;
        case 2:
            params->GenEventStatistics.N_2TIMES += 1;
            break;
        case 3:
            params->GenEventStatistics.N_3TIMES += 1;
            break;
        case 4:
            params->GenEventStatistics.N_4TIMES += 1;
            break;
        default:
            params->GenEventStatistics.N_NTIMES += 1;
            break;
        }
        // *ff = bornconfig.btilde / wgt;
        double neg = bornconfig.negative ? -1.0 : 1.0;
        *ff = params->BtildeState.Max / wgt * neg;
        return n - 1;
    }

    // no event written
    *ff = 0.0;
    return -1;
}

} // end namespace Powheg
