#include "generateevents.h" 

#include <iostream>
#include <map>

#include "phasespace/phasespace.h"
#include "phasespace/realphasespace.h"
#include "process/data.h"
#include "process/matrixelement.h"
#include "powheg/pickelement.h"
#include "powheg/generateradiation.h"
#include "fks/param.h"
#include "lhe/event.h"
#include "random/rnd.h"
#include "fks/xsec.h"
#include "powheg/resonance.h"
#include "powheg/reshuffle_momenta.h"
#include "powheg/unweighting.h"

#define EVENT_REJECT 1
#define EVENT_SUCCESS 2
#define EVENT_ERROR 3
#define EVENT_NEG 4

namespace Powheg {

namespace {

enum class Type {
    QCD,
    QED
};

int GenRadiation(const BornConfig & bornconfig, const Phasespace::Phasespace &ps, double wgt,
              UserProcess::Data *userdata, Phasespace::Phasespace *ps_out,
              int nn, int *pdg_out, Radiation *rad_out) {
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
            pdg_out[i + 2] = fl->Born.Flavours[i + 2];
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
    Phasespace::GenRealPhasespace(&ps_real, &ps, radiation.j, radiation.xi,
                                  radiation.y, radiation.phi);
    ps_out->SetToLabFromCMS(&ps_real);
    pdg_out[0] = radiation.Real->PDF[pdf][0];
    pdg_out[1] = radiation.Real->PDF[pdf][1];
    if (pdg_out[0] == 0) {
        pdg_out[0] = 21;
    }
    if (pdg_out[1] == 0) {
        pdg_out[1] = 21;
    }
    for (int i = 0; i < ps_out->N; i++) {
        pdg_out[i + 2] = radiation.Real->Flavours[i + 2];
    }
    userdata->GenEventStatistics.N_REAL += 1;
    return EVENT_SUCCESS;
}

void assign_color_drellyan(bool qcd, LHE::Particle *particle, int pdg, int i) {
    if (qcd) {                          // QCD
        if (i < 2) {                    // inital state
            if (pdg >= -5 && pdg < 0) { // initial state qbar
                particle->Color1 = 0;
                particle->Color2 = 501;
                return;
            }
            if (pdg > 0 && pdg <= 5) { // initial state q
                particle->Color1 = 502;
                particle->Color2 = 0;
                return;
            }
            if (pdg == 21) { // initial state gluon
                particle->Color1 = 501;
                particle->Color2 = 502;
                return;
            }
        } else {
            if (pdg >= -5 && pdg < 0) { // final state qbar
                particle->Color1 = 0;
                particle->Color2 = 502;
                return;
            }
            if (pdg > 0 && pdg <= 5) { // final state q
                particle->Color1 = 501;
                particle->Color2 = 0;
                return;
            }
            if (pdg == 21) { // final state gluon
                particle->Color1 = 502;
                particle->Color2 = 501;
                return;
            }
        }
    }
    if (pdg >0 && pdg <= 5) {
        particle->Color1 = 501;
        particle->Color2 = 0;
    }
    if (pdg >= -5 && pdg < 0) {
        particle->Color1 = 0;
        particle->Color2 = 501;
    }
}

void make_event(const Phasespace::Phasespace &ps_out, const Radiation &rad,
                int *pdgs, UserProcess::Data *params) {
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
        event.ID = 0;

        event.N = N;
        const Resonance &resonance =
            (rad.j >= 2) ? params->ResonanceFSR : params->ResonanceISR;

        bool qcd =false;
        for (int i = 0; i < N; i++) {
            if(pdgs[i] == 21) {
                qcd = true;
            }
        }

        bool has_resonance = false;
        for (int i = 0, lhe_i = 0; i < N; i++, lhe_i++) {
            if (i >= 2 && i == resonance.ID[0]) {
                lhe_i += 1;
                event.N += 1;
                has_resonance = true;
            }
            LHE::Particle & particle = event.Particles[lhe_i];
            particle.PDG= pdgs[i];
            if (i < 2) {
                particle.Status = -1;
            } else {
                particle.Status = 1;
                if (i == resonance.ID[0] || i == resonance.ID[1] ||
                    i == resonance.ID[2]) {
                    // the actual mother particle is the resonance ID[0].
                    // However, the LHE counting starts with 1, therefore, we
                    // have to add +1 to get the correct value.
                    particle.Mother1 = resonance.ID[0] + 1;
                    particle.Mother2 = resonance.ID[0] + 1;
                } else {
                    particle.Mother1 = 1;
                    particle.Mother2 = 2;
                }
            }
            particle.Color1 = 0;
            particle.Color2 = 0;
            particle.Momentum = ps_out.Momenta[i];
            // TODO: think of a proper way to introduce color! A method is
            // explained in the POWHEG paper but this also deals with complex
            // cases.
            assign_color_drellyan(qcd, &particle, pdgs[i], i);
        }
        if (event.N == N + 1) {
            int i = resonance.ID[0];
            int i2 = resonance.ID[1];
            int i3 = resonance.ID[2];
            LHE::Particle & particle = event.Particles[i];
            particle.PDG = resonance.pdg;
            particle.Mother1 = 1;
            particle.Mother2 = 2;
            particle.Status = 2;
            const Math::FourMomentum &P1 = ps_out.Momenta[i];
            const Math::FourMomentum &P2 = ps_out.Momenta[i2];
            particle.Momentum = P1.Plus(P2);
            if (i3 != 0) {
                particle.Momentum.Add(ps_out.Momenta[i3]);
            }
        }
        if (has_resonance) {
            reshuffle_momenta(&event, resonance, params);
        }

        auto ret = params->EventBuffer.Append(event);
        switch(ret) {
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
    } else if (bornconfig.status == BornConfig::Status::RejectedWithGuessedVirtual) {
        params->GenEventStatistics.N_REJECTWOVIRTUAL += 1;
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
        int status = Powheg::GenRadiation(bornconfig, ps, wgt, params, &ps_out, 10,
                                       pdgs, &rad);
        if (status == EVENT_SUCCESS) {
            make_event(ps_out, rad, pdgs, params);
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
        *ff = params->BtildeState.Max / wgt;
        return n - 1;
    }

    // no event written
    *ff = 0.0;
    return -1;
}

} // end namespace Powheg
