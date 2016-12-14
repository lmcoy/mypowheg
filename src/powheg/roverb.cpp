#include "powheg/roverb.h"

#include <cassert>
#include <array>
#include <algorithm>

#include "pdf/pdfinterface.h"
#include "fks/ximax.h"
#include "fks/process.h"
#include "fks/sfunctions.h"
#include "fks/param.h"
#include "phasespace/phasespace.h"
#include "phasespace/realphasespace.h"
#include "math/math.h"
#include "process/data.h"
#include "process/matrixelement.h"
#include "fks/radiationregion.h"

namespace {

#define PDFCACHE_ENABLE

class PDFCache {
  public:
    PDFCache() {}
    double Get(const PDF::InterfacePtr &  pdf, double x, double scale, int pdg) {
#ifdef PDFCACHE_ENABLE
        for (const auto &v : cache) {
            if (v.pdg == 0xff) {
                break;
            }
            if (x == v.x && scale == v.muF && pdg == v.pdg) {
                return v.fx;
            }
        }
        PDFValue nv;
        nv.pdg = pdg;
        nv.x = x;
        nv.muF = scale;
        double xfx = pdf->Xfx(x, scale, pdg);
        nv.fx = xfx / x;

        cache[pos] = nv;
        pos += 1;
        if (pos >= 12) {
            pos = 0;
        }

        return nv.fx;
#else
        return pdf->Xfx(x, scale, pdg) / x;
#endif
    }

  private:
#ifdef PDFCACHE_ENABLE
    struct PDFValue {
        int pdg = 0xff;
        double x = 0.0;
        double muF = 0.0;
        double fx = 0.0;
    };
    std::array<PDFValue, 12> cache;
    size_t pos = 0;
#endif
};

} // end namespace

namespace Powheg {

void LumiRatio(const FKS::RadiationRegion &radreg,
               const Phasespace::Phasespace &ps_born,
               const Phasespace::Phasespace &ps_real,
               UserProcess::Data *userdata, double scale, int usepdf,
               Util::StaticMatrix32 *output) {
    assert(ps_born.N + 1 == ps_real.N);
    const auto fl = radreg.FlavourConfig;
    size_t pdf_len = fl->Born.PDF.size();
    size_t n_real = radreg.RealFlavour.size();
    output->Reset(pdf_len, n_real, 0.0);
    if (ps_real.X1 >= 1.0 || ps_real.X2 >= 1.0) {
        // TODO: It is possible that we have an unphysical real phase space
        // because the xi_max which is used in the sudakov for ISR radiation is
        // not the real xi_max. In the derivation of xi_max it is only used that
        // x_1 * x_2 < 1 but this doesn't require x1,x2 < 1.
        // The actual xi_max depends on y and is given in the derivation of the
        // initial state phase space generation (c.f. FKS::XiMaxISR).
        // Using this leads to a more complicated integral in the sudakov
        // because of the integration boundaries. Therefore, we want to avoid
        // this.
        // We simply use the wrong xi_max and return 0 for the R over B ratio.
        // Therefore, these unphysical points are never used. This is equivalent
        // to introducing theta(xi_max(y)-xi) to the R/B sudakov and integration
        // xi up to the wrong xi_max. I assume that this is correct but I'm not
        // 100% sure.
        return;
    }

    PDFCache pdfcache;

    bool FSR = radreg.Region.J >= 2;
    assert(output);
    for (size_t pdf = 0; pdf < pdf_len; pdf++) {
        if (usepdf >= 0 && (size_t)usepdf != pdf) {
            continue;
        }
        int pdg_born_1 = fl->Born.PDF[pdf][0];
        int pdg_born_2 = fl->Born.PDF[pdf][1];

        for (size_t i = 0; i < n_real; i++) {
            const auto real = radreg.RealFlavour[i];
            assert(real->PDF.size() == pdf_len);
            int pdg_real_1 = real->PDF[pdf][0];
            int pdg_real_2 = real->PDF[pdf][1];
            bool ratio1 = !(FSR && pdg_born_1 == pdg_real_1);
            bool ratio2 = !(FSR && pdg_born_2 == pdg_real_2);

            double denom1 = 1.0;
            double denom2 = 1.0;
            if (!ratio1 && !ratio2) {
                output->Set(pdf, i, 1.0);
                continue;
            }
            if (ratio1) {
                denom1 =
                    pdfcache.Get(userdata->pdf, ps_born.X1, scale, pdg_born_1);
            }
            if (ratio2) {
                denom2 =
                    pdfcache.Get(userdata->pdf, ps_born.X2, scale, pdg_born_2);
            }
            if (denom1 == 0.0 || denom2 == 0.0) {
                output->Set(pdf, i, -1.0);
                continue;
            }

            double f1 = 1.0;
            double f2 = 1.0;
            if (ratio1) {
                f1 =
                    pdfcache.Get(userdata->pdf, ps_real.X1, scale, pdg_real_1) /
                    denom1;
            }
            if (ratio2 && f1 > 0.0) {
                f2 =
                    pdfcache.Get(userdata->pdf, ps_real.X2, scale, pdg_real_2) /
                    denom2;
            }
            output->Set(pdf, i, f1 * f2);
        }
    }
}

std::array<double, 8> RoverB(double B, const FKS::RadiationRegion &radreg,
                             const Phasespace::Phasespace &ps,
                             const Phasespace::Phasespace &ps_real,
                             double rad_alpha, UserProcess::Data *userdata) {
    std::array<double, 8> output;
    output.fill(0.0);
    assert(ps_real.N == ps.N + 1);

    double J_rad = ps_real.Jacobian / ps.Jacobian;

    double flux_real = 1.0 / (2.0 * ps_real.X1 * ps_real.X2 * ps_real.S);
    double flux_born = 1.0 / (2.0 * ps.X1 * ps.X2 * ps.S);

    size_t n_real = radreg.RealFlavour.size();
    assert(n_real < 8);
    if (ps_real.X1 >= 1.0 || ps_real.X2 >= 1.0) {
        // TODO: It is possible that we have an unphysical real phase space
        // because the xi_max which is used in the sudakov for ISR radiation is
        // not the real xi_max. In the derivation of xi_max it is only used that
        // x_1 * x_2 < 1 but this doesn't require x1,x2 < 1.
        // The actual xi_max depends on y and is given in the derivation of the
        // initial state phase space generation (c.f. FKS::XiMaxISR).
        // Using this leads to a more complicated integral in the sudakov
        // because of the integration boundaries. Therefore, we want to avoid
        // this.
        // We simply use the wrong xi_max and return 0 for the R over B ratio.
        // Therefore, these unphysical points are never used. This is equivalent
        // to introducing theta(xi_max(y)-xi) to the R/B sudakov and integration
        // xi up to the wrong xi_max. I assume that this is correct but I'm not
        // 100% sure.
        return output;
    }

    double pre = flux_real / flux_born * J_rad / B;

    bool mod = userdata->modSudakov;
    double alpha = 1.0;
    switch (radreg.Type) {
    case FKS::Type_t::EW:
        alpha = userdata->Params->alpha;
        break;
    case FKS::Type_t::QCD:
        alpha = userdata->Params->alphaS();
        break;
    }

    for(size_t j = 0; j < n_real; j++) {
        int rflavour = radreg.RealFlavour[j]->ID;
        const FKS::RegionList &rlist = radreg.RealFlavour[j]->Regions;
        const auto &real = *radreg.RealFlavour[j];
        double S = 1.0;
        if (rlist.size() > 1) {
            S = FKS::SFunction(ps_real, radreg.Region, real,
                               userdata->useResonanesInS);
        }
        if (!mod && S == 0.0) {
            continue;
        }
        bool nointer = userdata->noInterference;
        double R_isr = 0.0;
        double R_fsr = 0.0;
        if (mod || nointer) {
            using Diag = UserProcess::IMatrixElement::Diagrams;
            R_isr = userdata->MatrixElement->Real(
                rflavour, ps_real, userdata->Params, Diag::ONLYISR);
            R_fsr = userdata->MatrixElement->Real(
                rflavour, ps_real, userdata->Params, Diag::ONLYFSR);
        }
        double R_int = 0.0;
        if (!nointer) {
            double R = userdata->MatrixElement->Real(rflavour, ps_real,
                                                     userdata->Params);
            R_int = R - R_isr - R_fsr;
        }

        double ff = 0.0;
        if (radreg.Region.J < 2) {
            if (mod) {
                ff = pre * (R_isr + S * R_int);
            } else {
                ff = pre * S * (R_isr + R_fsr + R_int);
            }
        } else {
            FKS::RegionList rlist_wo_isr;
            if (mod) {
                auto real_wo_isr = real;
                real_wo_isr.Regions.erase(
                    std::remove_if(real_wo_isr.Regions.begin(),
                                   real_wo_isr.Regions.end(),
                                   [](FKS::Region r) { return r.J < 2; }),
                    real_wo_isr.Regions.end());
                double S_fsr =
                    FKS::SFunction(ps_real, radreg.Region, real_wo_isr,
                                   userdata->useResonanesInS);
                ff = pre * (S_fsr * R_fsr + S * R_int);
            } else {
                ff = pre * S * (R_fsr + R_isr + R_int);
            }
        }

        output[j] = ff * rad_alpha / alpha;
    }
    return output;
}

} // end namespace powheg

