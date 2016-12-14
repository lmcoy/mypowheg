#include <cmath>

#include "powheg/unweighting.h"

#include "phasespace/phasespace.h"
#include "process/data.h"
#include "powheg/btilde.h"

using namespace Powheg;

BornConfig Powheg::unweighting(const Phasespace::Phasespace &ps, double x1,
                               double x2, double x3, double wgt,
                               UserProcess::Data *userdata) {
    BornConfig bconfig;
    Btilde_t btilde;
    btilde.CalcWOVirtual(ps, x1, x2, x3, wgt, userdata);
    double totg = -1.0;
    if (userdata->Unweighting.GuessVirtual) {
        totg = btilde.TotalWithGuessedVirtual();
    }

    double u = userdata->rng.Random();

    double Max = userdata->BtildeState.Max;
    if (totg < u * Max && totg > 0.0 && totg < Max) {
        // reject
        bconfig.status = BornConfig::Status::RejectedWithGuessedVirtual;
        return bconfig;
    }

    btilde.CalcVirtual(ps, x1, x2, x3, wgt, userdata);
    double tot = btilde.Total();
    if (totg > 0.0 && tot > totg) {
        bconfig.status = BornConfig::Status::UnderestimatedVirtual;
        return bconfig;
    }

    bconfig.btilde = tot;
    if (tot < 0.0) {
        bconfig.status = BornConfig::Status::NegativeBtilde;
        return bconfig;
    }
    if (tot == 0.0) {
        bconfig.status = BornConfig::Status::Rejected;
        return bconfig;
    }
    auto meta = btilde.PickElement(userdata->rng.Random());
    bconfig.ProcessIndex = meta.first;
    bconfig.PDFindex = meta.second;

    // if tot > Max unweighting does not work properly anymore. This function
    // returns how often we would have to write this event.
    if (tot > Max) {
        double r = tot / Max;
        int n = (int)trunc(r);
        double d = r - trunc(r);
        if (userdata->rng.Random() < d) {
            n += 1;
        }
        // accept event
        bconfig.status = BornConfig::Status::Accepted;
        bconfig.n = n;
        return bconfig;
    }
    if (tot < u * Max) {
        // reject
        bconfig.status = BornConfig::Status::Rejected;
        return bconfig;
    }
    // accept event
    bconfig.status = BornConfig::Status::Accepted;
    bconfig.n = 1;
    return bconfig;
}

