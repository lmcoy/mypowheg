#include "powheg/upperbounding.h"

#include <iostream>

#include "libconfig.h"

#include "powheg/roverb.h"
#include "fks/ximax.h"
#include "fks/process.h"
#include "fks/sfunctions.h"
#include "fks/param.h"
#include "phasespace/phasespace.h"
#include "math/math.h"
#include "process/data.h"
#include "process/matrixelement.h"
#include "fks/radiationregion.h"
#include "powheg/alphas.h"

#include "phasespace/realphasespace.h"

namespace {

double xi_from_x1(const Phasespace::Phasespace & ps, int jfks, double x1, double y) {
    if (jfks >= 2) {
        return x1 * FKS::XiMaxFSR(sqrt(ps.X1 * ps.X2 * ps.S), ps.Momenta[jfks].E());
    }
    return x1 * FKS::XiMaxISR(ps.X1, ps.X2, y);
}

double Born(int bflv, const Phasespace::Phasespace & ps, UserProcess::Data *userdata) {

    const FKS::Scales &scales = userdata->PowhegScales;
    if (userdata->BornMEStatus != UserProcess::BornMEStatus_t::CrossSection) {
        InitBornME(userdata->Params, userdata->Params_as, AlphaScheme::Gmu,
                   scales.mu * scales.mu, 0.0, &userdata->counterterm);
        userdata->BornMEStatus = UserProcess::BornMEStatus_t::CrossSection;
    }
    BornMEOut bme;
    BornME(bflv, ps, sqrt(scales.Q2), userdata->Params, userdata->Params_as,
           &bme, false, false, userdata->counterterm);

    return bme.M2;
}

} // end namespace

namespace Powheg {

void GetNormForBoundingFcts(FKS::Type_t type, bool init, FKS::RadiationRegion &radregion,
                  const Phasespace::Phasespace &ps, double x1, double x2,
                  double x3, UserProcess::Data *userdata) {
    double y = x2 * 2.0 - 1.0;
    double xi = xi_from_x1(ps, radregion.Region.J, x1, y);
    double phi = x3 * 2.0 * Math::Pi;

    double U = 0.0;
    double alpha = 1.0/132.;
    double scale = userdata->PowhegScales.muF;
    if(type == FKS::Type_t::QCD) {
        // set alpha_s
        double kT2 = 0.0;
        double s = ps.X1 * ps.X2 * ps.S;
        if (radregion.Region.J < 2) { // ISR
            kT2 = s / (4.0 * (1.0 - xi)) * xi * xi * (1.0 - y * y);
        } else { // FSR
            kT2 = s / 2.0 * xi * xi * (1.0 - y);
        }
        double kT = sqrt(kT2);
        if (kT2 < userdata->RadiationParameter.kT2min) {
            return;
        }
        RadiationAlphaS ralphas(userdata->AlphaS);
        alpha = ralphas.SimpleAlphaS(kT2);
        assert(alpha > 0.0);

        userdata->Params_as->Set(alpha);
        scale = kT;
        if(kT < 3.0) {
            // if the kT is too small, we get problems because the pdfs have
            // numerical problems with too small scales. We exclude those points
            // from the search for the norm because they are really rare in
            // event generation. If there are too many of those points, you see
            // it in the event generation process because the norm is not
            // correct.
            return;
        }
    } else {
        double as = userdata->pdf->AlphaS(userdata->PowhegScales.mu);
        userdata->Params_as->Set(as);
        alpha = userdata->Params->alpha;
    }
    if (radregion.Region.J < 2) { // ISR
        U = xi * (1.0 - y * y) / alpha;
    } else { // FSR
        U = xi * (1.0 - y) / alpha;
    }
    const auto & fl = radregion.FlavourConfig;
    double B = Born(fl->Born.ID, ps, userdata);
    double lumi = 0.0;
    const auto &pdgs = fl->Born.PDF;
    for (auto &pdg : pdgs) {
        double f1 = userdata->pdf->Xfx(ps.X1, scale, pdg[0]) / ps.X1;
        double f2 = userdata->pdf->Xfx(ps.X2, scale, pdg[1]) / ps.X2;
        if (f1 * f2 > 0.0) {
            lumi += f1 * f2;
        }
    }

    int jfks = radregion.Region.J;
    assert(jfks >= 0);

    Phasespace::Phasespace ps_real;
    Phasespace::GenRealPhasespace(&ps_real, &ps, jfks, xi, y, phi);

    Util::StaticMatrix32 lumi_ratio(0, 0, 0);
    LumiRatio(radregion, ps, ps_real, userdata, scale, -1, &lumi_ratio);

    auto roverb = RoverB(B, radregion, ps, ps_real, alpha, userdata);

    assert((size_t)lumi_ratio.Rows() <= radregion.NPDF &&
           "too many born pdfs! increase NPDF");
    for (int pdf = 0; pdf < lumi_ratio.Rows(); pdf++) {
        double result = 0.0;
        for (int rflavour = 0; rflavour < lumi_ratio.Cols(); rflavour++) {
            double L = lumi_ratio.Get(pdf, rflavour);
            if (L < 0.0) {
                L = 0.0;
            }
            result += L * roverb[rflavour];
        }
        double N = result * U;

        if (!init) {
            LIB_ASSERT(
                radregion.NormHist[pdf],
                "there is no histogram allocated for pdf=%d. Did you supply "
                "an Nmax in radregion.CreateHistograms?",
                pdf);
            radregion.NormHist[pdf]->AddValue(N, 1.0);
        } else {
            // TODO: add option in userdata to set these values
            N *= 2;
            radregion.InitHist[pdf]->AddValue(N, 1.0);
        }
    }
}

} // end namespace Powheg
