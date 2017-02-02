#include "powheg/upperbounding.h"

#include <iostream>

#include "libconfig.h"

#include "fks/param.h"
#include "fks/process.h"
#include "fks/radiationregion.h"
#include "fks/sfunctions.h"
#include "fks/ximax.h"
#include "math/math.h"
#include "phasespace/phasespace.h"
#include "powheg/alphas.h"
#include "powheg/roverb.h"
#include "process/data.h"
#include "process/matrixelement.h"

#include "phasespace/realphasespace.h"

namespace {

double xi_from_x1(const Phasespace::Phasespace &ps, int jfks, double x1,
                  double y) {
    if (jfks >= 2) {
        return x1 *
               FKS::XiMaxFSR(sqrt(ps.X1 * ps.X2 * ps.S), ps.Momenta[jfks].E());
    }
    return x1 * FKS::XiMaxISR(ps.X1, ps.X2, y);
}

} // end namespace

namespace Powheg {

void GetNormForBoundingFcts(FKS::Type_t type, bool init,
                            FKS::RadiationRegion &radregion,
                            const Phasespace::Phasespace &ps, double x1,
                            double x2, double x3, UserProcess::Data *userdata) {
    double y = x2 * 2.0 - 1.0;
    int jfks = radregion.Region.J;
    int ifks = radregion.Region.I;
    int mother = (ifks > jfks) ? jfks : ifks;
    double xi = xi_from_x1(ps, mother, x1, y);
    double phi = x3 * 2.0 * Math::Pi;

    double U = 0.0;
    double alpha = 1.0 / 132.;
    double scale = userdata->PowhegScales.muF;
    if (type == FKS::Type_t::QCD) {
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

        userdata->Params->SetAlphaS(alpha);
        scale = kT;
        if (kT < 3.0) {
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
        userdata->Params->SetAlphaS(as);
        alpha = userdata->Params->alpha;
    }
    if (radregion.Region.J < 2) { // ISR
        U = xi * (1.0 - y * y) / alpha;
    } else { // FSR
        U = xi * (1.0 - y) / alpha;
    }
    const auto &fl = radregion.FlavourConfig;
    auto bme = userdata->MatrixElement->Born(fl->Born.ID, ps,
                                             sqrt(userdata->PowhegScales.Q2),
                                             userdata->Params, false, false);
    double B = bme.M2;
    double lumi = 0.0;
    const auto &pdgs = fl->Born.PDF;
    for (auto &pdg : pdgs) {
        double f1 = userdata->pdf->Xfx(ps.X1, scale, pdg[0]) / ps.X1;
        double f2 = userdata->pdf->Xfx(ps.X2, scale, pdg[1]) / ps.X2;
        if (f1 * f2 > 0.0) {
            lumi += f1 * f2;
        }
    }

    Phasespace::Phasespace ps_real;
    Phasespace::GenRealPhasespace(&ps_real, &ps, ifks, jfks, xi, y, phi);

    Util::StaticMatrix128 lumi_ratio(0, 0, 0);
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
            radregion.InitHist[pdf]->AddValue(N, 1.0);
        }
    }
}

} // end namespace Powheg
