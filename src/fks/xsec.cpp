#include "fks/xsec.h"

#include <cassert>
#include <iostream>
#include <memory>

#include "fks/fks_g.h"
#include "fks/limits.h"
#include "fks/luminosity.h"
#include "fks/param.h"
#include "fks/phasespaces.h"
#include "fks/process.h"
#include "fks/remnants.h"
#include "fks/scales.h"
#include "fks/subtraction.h"
#include "fks/virtual.h"
#include "fks/ximax.h"
#include "math/math.h"
#include "pdf/pdfcache.h"
#include "phasespace/phasespace.h"
#include "process/cuts.h"
#include "process/data.h"
#include "process/matrixelement.h"

namespace {
using namespace FKS;

const double CrossSectionToPb = 3.8937944929e8;

enum class Type { QED, QCD };

struct CutsPassed {
    bool Real;
    bool Collinear1;
    bool Collinear2;
};

CutsPassed realCuts(const FKS::Phasespaces &PS_recomb_lab,
                    const UserProcess::ICutsPtr &cuts, int jfks,
                    const int *pdgs_r) {
    CutsPassed results;
    int N = PS_recomb_lab.Real.N + 2;

    results.Real =
        cuts->ApplyCuts(N, pdgs_r, PS_recomb_lab.Real.Momenta.data());
    results.Collinear1 =
        cuts->ApplyCuts(N, pdgs_r, PS_recomb_lab.Collinear1.Momenta.data());

    results.Collinear2 = false;
    if (jfks < 2) {
        results.Collinear2 =
            cuts->ApplyCuts(N, pdgs_r, PS_recomb_lab.Collinear2.Momenta.data());
    }
    return results;
}

typedef std::array<FKS::PartonLuminosity, FKS::MAXPDF> LumiReal;

/**
 * @brief lumi for the real cross section part
 *
 * realLuminosity returns the luminosities needed for the real subtracted part
 * of the cross section. It returns an array which contains the luminosity for
 * the subprocesses (e.g. (uu~, cc~))
 *
 * @param pdf         PDF cache
 * @param usedpdf     bit mask, specifies if sub process is used (0x1 use first
 *                    pdf,..., 0xfffff use all pdfs)
 * @param ps          born phase space
 * @param ps_real     real phase space
 * @param pdgs        list of pdgs
 * @param muF         factorization scale
 * @param x           x integration variabel
 * @param Xi          max. values for x
 */
LumiReal realLuminosity(PDF::Cache *pdf, unsigned long usedpdf,
                        const Phasespace::Phasespace &ps,
                        const Phasespace::Phasespace &ps_real,
                        const std::vector<std::array<int, 2>> &pdgs, double muF,
                        double x, const FKS::Xi &Xi) {
    LumiReal result;

    size_t len = pdgs.size();
    for (size_t i = 0; i < len; i++) {
        if (!(usedpdf & (1u << i))) {
            continue;
        }
        const auto &pdfs = pdgs[i];
        int fl1 = pdfs[0];
        int fl2 = pdfs[1];
        result[i] = FKS::Luminosity(pdf, fl1, fl2, x, ps.X1, ps.X2, ps_real.X1,
                                    ps_real.X2, Xi, muF);
    }

    return result;
}

FKS::MatrixElement realME(const FKS::Real_t &real, const FKS::Phasespaces &PS,
                          const FKS::Region &region, double x,
                          const FKS::Xi &Xi, double y, double phi,
                          const UserProcess::IMatrixElement::Result &bme,
                          const CutsPassed &cutspassed, bool born_passcuts,
                          const UserProcess::Data *userdata) {
    FKS::MatrixElement ME;
    Util::Matrix2 b(1, 0.0);
    b.SetLen(1);
    b.Set(0, 0, bme.M2);
    const auto typ = real.Type;
    if (cutspassed.Real) {
        double xi = x * Xi.Max;
        switch (typ) {
        case Type_t::EW:
            ME.Real = FKS::QED::SxG(real, PS, bme.M2, b, bme.SpinCorr, region,
                                    xi, y, phi, userdata);
            break;
        case Type_t::QCD:
            ME.Real = FKS::QCD::SxG(real, PS, bme.M2, bme.ColorCorr,
                                    bme.SpinCorr, region, xi, y, phi, userdata);
            break;
        }
    }

    FKS::MatrixElement lim;
    if (born_passcuts || cutspassed.Collinear1 || cutspassed.Collinear2) {
        switch (typ) {
        case Type_t::EW:
            lim = FKS::QED::Limits(real, PS, bme.M2, b, bme.SpinCorr, region, x,
                                   Xi, y, phi, userdata);
            break;
        case Type_t::QCD:
            lim =
                FKS::QCD::Limits(real, PS, bme.M2, bme.ColorCorr, bme.SpinCorr,
                                 region, x, Xi, y, phi, userdata);
            break;
        }
    }
    if (born_passcuts) {
        ME.Soft = lim.Soft;
        ME.SoftCollinear1 = lim.SoftCollinear1;
        ME.SoftCollinear2 = lim.SoftCollinear2;
    }
    if (cutspassed.Collinear1) {
        ME.Collinear1 = lim.Collinear1;
    }
    if (cutspassed.Collinear2) {
        ME.Collinear2 = lim.Collinear2;
    }
    return ME;
}

double sum_result(const Result &result) {
    double r = 0.0;
    for (size_t i = 0; i < result.size(); i++) {
        r += result[i];
    }
    return r;
}

void add_result(Result *r1, const Result &r2) {
    for (size_t i = 0; i < r1->size(); i++) {
        (*r1)[i] += r2[i];
    }
}

struct Results {
    Results() {
        Born1.fill(0.0);
        Born2.fill(0.0);
        Real1.fill(0.0);
        Real2.fill(0.0);
        Coll1.fill(0.0);
        Coll2.fill(0.0);
    }
    Result Born1;
    Result Born2;
    Result Real1;
    Result Real2;
    Result Coll1;
    Result Coll2;
    Result Combine() {
        Result result;
        result.fill(0.0);
        size_t len = Born1.size();
        for (size_t i = 0; i < len; i++) {
            result[i] += Born1[i] + Real1[i] + Coll1[i];
            result[i] += Born2[i] + Real2[i] + Coll2[i];
        }
        return result;
    }
};

Results realSubtraction(const Phasespace::Phasespace &ps_s, double s_r,
                        int jfks, const LumiReal &lumi,
                        const FKS::MatrixElement &ME, double x,
                        const FKS::Xi &Xi, double y) {
    Results results;

    // multiply only by the born jacobian det. The jacobian for the n+1
    // particle is applied in RealISR.
    int i_max = lumi.size();
    double xi = x * Xi.Max;
    for (int i = 0; i < i_max; i++) {
        if (jfks >= 2) {
            FKS::SubtractionTerms terms = FKS::RealFSR(
                lumi[i].Real, ME, x, Xi.Max, y, ps_s.Momenta[jfks].E(), s_r);
            results.Real1[i] = terms.Real;
            results.Born1[i] = terms.Soft + terms.SoftCollinear;
            results.Coll1[i] = terms.Collinear;
        } else {
            if (xi > 1e-9 && fabs(y) < 1.0 - 1e-7) {
                // need this cut to avoid numerical failure of subtraction
                auto terms1 = FKS::RealISR2(lumi[i], ME, 1, x, Xi, y);
                results.Real1[i] = terms1.Real;
                results.Born1[i] = terms1.Soft + terms1.SoftCollinear;
                results.Coll1[i] = terms1.Collinear;

                auto terms2 = FKS::RealISR2(lumi[i], ME, -1, x, Xi, y);
                results.Real2[i] = terms2.Real;
                results.Born2[i] = terms2.Soft + terms2.SoftCollinear;
                results.Coll2[i] = terms2.Collinear;
            }
        }
    }
    return results;
}

void realFillHists(Histograms *hists, const Results &results, double wgt,
                   const Phasespace::Phasespace &ps_lab,
                   const FKS::Phasespaces &ps_rec_lab, bool born_passcuts,
                   const CutsPassed &cutspassed,
                   const UserProcess::Data *data) {
    Result result_born1 = results.Born1;
    Result result_born2 = results.Born2;
    Result result_real1 = results.Real1;
    Result result_real2 = results.Real2;
    Result result_coll1 = results.Coll1;
    Result result_coll2 = results.Coll2;
    const auto &scale = data->Process[data->ProcessID].Scales;
    for (size_t i = 0; i < scale.size(); i++) {
        result_born1[i] *= scale[i];
        result_born2[i] *= scale[i];
        result_real1[i] *= scale[i];
        result_real2[i] *= scale[i];
        result_coll1[i] *= scale[i];
        result_coll2[i] *= scale[i];
    }
    double r_born1 = sum_result(result_born1);
    double r_born2 = sum_result(result_born2);
    double r_real1 = sum_result(result_real1);
    double r_real2 = sum_result(result_real2);
    double r_coll1 = sum_result(result_coll1);
    double r_coll2 = sum_result(result_coll2);

    if (born_passcuts) {
        hists->Fill(&ps_lab, wgt * (r_born1 + r_born2));
    }
    if (cutspassed.Real) {
        hists->Fill(&ps_rec_lab.Real, wgt * (r_real1 + r_real2));
    }
    if (cutspassed.Collinear1) {
        hists->Fill(&ps_rec_lab.Collinear1, wgt * r_coll1);
    }
    if (cutspassed.Collinear2) {
        hists->Fill(&ps_rec_lab.Collinear2, wgt * r_coll2);
    }
}

typedef std::array<FKS::LumRemnants, FKS::MAXPDF> LumiRemnants;

/**
 * @brief Luminosity for coll. remnants
 *
 * remnantPDF returns the luminosity for the collinear remnant part of the cross
section.
 */
LumiRemnants remnantPDF(PDF::Cache *pdf, unsigned long usedpdf, double x1,
                        double x2, const std::vector<int> &remnantpdf,
                        double muF, double x, double xi_max_coll, int y) {
    assert(y == 1 || y == -1);
    LumiRemnants lumir;

    size_t len = remnantpdf.size();
    for (size_t j = 0; j < len; j += 2) {
        int i = j / 2;
        if (!(usedpdf & (1u << i))) {
            continue;
        }
        int flv1 = remnantpdf[j];
        int flv2 = remnantpdf[j + 1];
        auto l = FKS::LuminosityRemnants(pdf, flv1, flv2, y, x1, x2, x,
                                         xi_max_coll, muF);
        lumir[i].Born = l.Born;
        lumir[i].Remnant = l.Remnant;
    }

    return lumir;
}

typedef std::array<double, FKS::MAXPDF> LumiBorn;

class Integrand {
  public:
    Integrand(std::shared_ptr<PDF::Interface> p) : pdf(p) {}
    void Init(const Phasespace::Phasespace &ps, bool Virtual, double wgt,
              UserProcess::Data *userdata);

    bool PSCuts(const UserProcess::ICutsPtr &) const;
    Result Born(const Phasespace::Phasespace &ps, const LumiBorn &lumi,
                Histograms *hists) const;
    Result Virtual(const Phasespace::Phasespace &ps, const LumiBorn &lumi,
                   Histograms *hists) const;
    Result CollRemnant(const Phasespace::Phasespace &ps, double x1,
                       const UserProcess::Data *userdata,
                       Histograms *hists) const;
    Result Real(const Phasespace::Phasespace &ps, bool, double x1, double x2,
                double x3, const UserProcess::Data *userdata,
                Histograms *hists) const;

    LumiBorn LuminosityBorn(unsigned long, double, double) const;

  private:
    Result real(const Phasespace::Phasespace &ps,
                const UserProcess::Data *userdata, double x, double y,
                double phi, bool born_passcuts, const Real_t &rflavour,
                const FKS::Region &region, double) const;

    Result remnant(Type_t type, const FKS::Splitting &splitting,
                   const LumiRemnants &lumir, const Phasespace::Phasespace &ps,
                   double x, double xi_max_coll, double alpha_rem,
                   PDFRenorm pdfren) const;

  private:
    UserProcess::Data *data;
    double alpha_ew;
    double alpha_s;
    int *pdgs_b;
    Phasespace::Phasespace ps_lab;
    mutable PDF::Cache pdf;
    UserProcess::IMatrixElement::Result bme;
    double wgt;
    FKS::Scales scales;
};

LumiBorn Integrand::LuminosityBorn(unsigned long usedpdf, double x1,
                                   double x2) const {
    // evalute pdf for the born phase space integration
    LumiBorn lumi;
    lumi.fill(0.0);
    const auto &fl = data->Process[data->ProcessID];
    size_t len = fl.Born.PDF.size();
    for (size_t i = 0; i < len; i++) {
        if (!(usedpdf & (1u << i))) {
            continue;
        }
        const auto &pdfs = fl.Born.PDF[i];
        int fl1 = pdfs[0];
        int fl2 = pdfs[1];
        double f1 = pdf.Get(x1, scales.muF, fl1);
        double f2 = pdf.Get(x2, scales.muF, fl2);
        lumi[i] = f1 * f2;
    }
    return lumi;
}

void Integrand::Init(const Phasespace::Phasespace &ps, bool Virtual,
                     double nwgt, UserProcess::Data *userdata) {
    wgt = nwgt;
    alpha_ew = userdata->Params->alpha;
    data = userdata;

    int proc = userdata->ProcessID;
    pdgs_b = userdata->Process[proc].Born.Flavours.data();

    ps_lab.SetToLabFromCMS(&ps);

    // define scales
    scales.muF = userdata->Scales->Factorization(ps);
    scales.mu = userdata->Scales->Renorm(ps);
    // use Q = mu! The Vfin term in BornME normalized to Q = mu
    scales.Q2 = scales.mu * scales.mu;

    alpha_s = userdata->AlphaS->AlphaS(scales.mu * scales.mu);
    userdata->Params->SetAlphaS(alpha_s);

    // calculate born matrix element
    if (userdata->BornMEStatus != UserProcess::BornMEStatus_t::CrossSection) {
        userdata->BornMEStatus = UserProcess::BornMEStatus_t::CrossSection;
    }
    int id = userdata->Process[proc].Born.ID;
    bool QCD = userdata->Process[proc].QCD & Virtual;
    bool EW = userdata->Process[proc].EW & Virtual;
    if (Virtual) {
        userdata->UseBornCache = false;
    }
    if (userdata->UseBornCache) {
        bme = userdata->BornCache;
    } else {
        bme = userdata->MatrixElement->Born(id, ps, scales.mu, userdata->Params,
                                            QCD, EW);
        userdata->BornCache = bme;
    }
    userdata->UseBornCache = false;
}

bool Integrand::PSCuts(const UserProcess::ICutsPtr &cuts) const {
    int N = ps_lab.N + 2;
    return cuts->ApplyCuts(N, pdgs_b, ps_lab.Momenta.data());
}

Result Integrand::Born(const Phasespace::Phasespace &ps, const LumiBorn &lumi,
                       Histograms *hists) const {
    Result result;
    double s = ps.X1 * ps.X2 * ps.S;
    double dsigma = 1.0 / (2.0 * s) * bme.M2;

    double result_born = dsigma * CrossSectionToPb;
    for (size_t i = 0; i < lumi.size(); i++) {
        result[i] = result_born * lumi[i];
    }

    Result result_hist = result;
    const auto &scale = data->Process[data->ProcessID].Scales;
    for (size_t i = 0; i < scale.size(); i++) {
        result_hist[i] *= scale[i] * ps.Jacobian;
    }
    
    hists->Fill(&ps_lab, wgt * sum_result(result_hist));
    return result;
}

Result Integrand::Virtual(const Phasespace::Phasespace &ps,
                          const LumiBorn &lumi, Histograms *hists) const {
    Result result;
    double v = 0.0;
    int N = ps.N + 2;
    double s = ps.X1 * ps.X2 * ps.S;
    bool QCD = data->Process[data->ProcessID].QCD;
    bool EW = data->Process[data->ProcessID].EW;
    if (EW) {
        Util::Matrix2 born_qed(1, bme.M2); // color correlated in QCD
        v += FKS::QED::Virtual(N, pdgs_b, ps.Momenta.data(), bme.M2, born_qed,
                               bme.VfinEW, alpha_ew, sqrt(s), scales);
    }
    if (QCD) {
        v += FKS::QCD::Virtual(N, pdgs_b, ps.Momenta.data(), bme.M2,
                               bme.ColorCorr, bme.VfinQCD, alpha_s, sqrt(s),
                               scales);
    }

    double flux = 1.0 / (2.0 * s);
    double result_v = v * flux * CrossSectionToPb;
    for (size_t i = 0; i < lumi.size(); i++) {
        result[i] = result_v * lumi[i];
    }
    Result result_hist = result;
    const auto &scale = data->Process[data->ProcessID].Scales;
    for (size_t i = 0; i < scale.size(); i++) {
        result_hist[i] *= scale[i] * ps.Jacobian;
    }
    hists->Fill(&ps_lab, wgt * sum_result(result_hist));
    return result;
}

Result Integrand::CollRemnant(const Phasespace::Phasespace &ps, double x1,
                              const UserProcess::Data *userdata,
                              Histograms *hists) const {
    Result result;
    result.fill(0.0);
    double xi_max_coll1 = FKS::XiMaxCollinearISR(ps.X1, ps.X2, 1);
    double xi_max_coll2 = FKS::XiMaxCollinearISR(ps.X1, ps.X2, -1);
    double alpha_rem = 0.0;
    int NumXi = userdata->NRemnXi;
    double denom = 1.0 / static_cast<double>(NumXi);

    int proc_id = userdata->ProcessID;
    const auto &fl = userdata->Process[proc_id];
    for (int i = 0; i < NumXi; i++) {
        double x1_i = (x1 + static_cast<double>(i)) * denom;
        for (const auto &remn : fl.Remnant1) {
            FKS::PDFRenorm pdfren = FKS::PDFRenorm::MSbar;
            if (remn.Type == Type_t::QCD) {
                alpha_rem = alpha_s;
                pdfren = userdata->PDFRenorm.QCD;
                if (!fl.QCD) {
                    continue;
                }
                pdfren = userdata->PDFRenorm.QCD;
            } else {
                if (userdata->OnlyVirtualEW) {
                    continue;
                }
                alpha_rem = alpha_ew;
                pdfren = userdata->PDFRenorm.EW;
                if (!fl.EW) {
                    continue;
                }
            }

            auto lumi = remnantPDF(&pdf, userdata->UsedPDF, ps.X1, ps.X2,
                                   remn.PDF, scales.muF, x1_i, xi_max_coll1, 1);
            auto res_remn = remnant(remn.Type, remn.Splitting, lumi, ps, x1_i,
                                    xi_max_coll1, alpha_rem, pdfren);
            for (size_t i = 0; i < result.size(); i++) {
                result[i] += res_remn[i];
            }
        }

        for (const auto &remn : fl.Remnant2) {
            FKS::PDFRenorm pdfren = FKS::PDFRenorm::MSbar;
            if (remn.Type == Type_t::QCD) {
                alpha_rem = alpha_s;
                pdfren = userdata->PDFRenorm.QCD;
                if (!fl.QCD) {
                    continue;
                }
            } else {
                if (userdata->OnlyVirtualEW) {
                    continue;
                }
                alpha_rem = alpha_ew;
                pdfren = userdata->PDFRenorm.EW;
                if (!fl.EW) {
                    continue;
                }
            }
            auto lumi =
                remnantPDF(&pdf, userdata->UsedPDF, ps.X1, ps.X2, remn.PDF,
                           scales.muF, x1_i, xi_max_coll2, -1);
            auto res_remn = remnant(remn.Type, remn.Splitting, lumi, ps, x1_i,
                                    xi_max_coll2, alpha_rem, pdfren);

            for (size_t i = 0; i < result.size(); i++) {
                result[i] += res_remn[i];
            }
        }
    }
    for (size_t i = 0; i < result.size(); i++) {
        result[i] *= denom;
    }
    Result result_hist = result;
    const auto &scale = data->Process[data->ProcessID].Scales;
    for (size_t i = 0; i < scale.size(); i++) {
        result_hist[i] *= scale[i] * ps.Jacobian;
    }
    hists->Fill(&ps_lab, wgt * sum_result(result_hist));
    return result;
}

Result Integrand::Real(const Phasespace::Phasespace &ps, bool born_passcuts,
                       double x1, double x2, double x3,
                       const UserProcess::Data *userdata,
                       Histograms *hists) const {
    Result result;
    result.fill(0.0);
    int procid = userdata->ProcessID;
    auto &fl = userdata->Process[procid];
    int NumXi = userdata->NRealXi;
    int NumY = userdata->NRealY;
    int NumPhi = userdata->NRealPhi;
    double denom_i = 1.0 / static_cast<double>(NumXi);
    double denom_j = 1.0 / static_cast<double>(NumY);
    double denom_k = 1.0 / static_cast<double>(NumPhi);
    double denom = denom_i * denom_j * denom_k;

    for (int i = 0; i < NumXi; i++) {
        for (int j = 0; j < NumY; j++) {
            for (int k = 0; k < NumPhi; k++) {
                double x1_i = (x1 + static_cast<double>(i)) * denom_i;
                double x2_j = (x2 + static_cast<double>(j)) * denom_j;
                double x3_k = (x3 + static_cast<double>(k)) * denom_k;
                // set variables for the radiated particle
                double y = 2.0 * x2_j - 1.0;
                double phi = 2.0 * Math::Pi * x3_k;
                for (const auto &rflavour : fl.Real) {
                    if (userdata->OnlyVirtualEW &&
                        rflavour.Type == FKS::Type_t::EW) {
                        continue;
                    }
                    double multiplicity =
                        static_cast<double>(rflavour.AllFlavours.size());
                    assert(multiplicity > 0.0);
                    for (const auto &region : rflavour.Regions) {
                        double w = 4.0 * Math::Pi * denom * CrossSectionToPb *
                                   multiplicity * ps.Jacobian;
                        Result r =
                            real(ps, userdata, x1_i, y, phi, born_passcuts,
                                 rflavour, region, w * wgt);
                        for (size_t k = 0; k < result.size(); k++) {
                            result[k] += multiplicity * r[k];
                        }
                    }
                }
            }
        }
    }
    for (size_t k = 0; k < result.size(); k++) {
        // the factor 4pi is the jacobian det of the change from y in [-1,1]
        // and phi in [0,2pi) to [0,1].
        result[k] *= 4.0 * Math::Pi * denom * CrossSectionToPb;
    }
    return result;
}

Result Integrand::real(const Phasespace::Phasespace &ps,
                       const UserProcess::Data *userdata, double x, double y,
                       double phi, bool born_passcuts, const Real_t &rflavour,
                       const FKS::Region &region, double rwgt) const {
    const UserProcess::RecombinationParam &recomb = userdata->Recomb;
    const int *pdgs_r = rflavour.Flavours.data();
    Histograms *hists = userdata->hists;
    int jfks = region.J;
    int ifks = region.I;

    FKS::Xi Xi;
    if (jfks >= 2) {
        int mother = (ifks > jfks) ? jfks : ifks;
        Xi.Max =
            FKS::XiMaxFSR(sqrt(ps.X1 * ps.X2 * ps.S), ps.Momenta[mother].E());
    } else {
        Xi.Max = FKS::XiMaxISR(ps.X1, ps.X2, y);
        Xi.Max_Coll1 = FKS::XiMaxCollinearISR(ps.X1, ps.X2, 1);
        Xi.Max_Coll2 = FKS::XiMaxCollinearISR(ps.X1, ps.X2, -1);
    }

    double xi = Xi.Max * x;

    // generate real phase space
    FKS::Phasespaces PS;
    PS.Generate(ps, ifks, jfks, xi, y, phi);
    auto PS_recomb_lab = PS.Recombined(rflavour.Type, pdgs_r, recomb);
    CutsPassed cutspassed;
    cutspassed.Real = true;
    cutspassed.Collinear1 = true;
    cutspassed.Collinear2 = (jfks < 2);
    if (userdata->CutOnReal) {
        cutspassed = realCuts(PS_recomb_lab, userdata->cuts, jfks, pdgs_r);
    }

    LumiReal lumi;
    if (born_passcuts || cutspassed.Real || cutspassed.Collinear1 ||
        cutspassed.Collinear2) {
        auto pdgs = rflavour.PDF;
        lumi = realLuminosity(&pdf, userdata->UsedPDF, ps, PS.Real, pdgs,
                              scales.muF, x, Xi);
    }

    FKS::MatrixElement ME = realME(rflavour, PS, region, x, Xi, y, phi, bme,
                                   cutspassed, born_passcuts, userdata);

    double s_r = PS.Real.X1 * PS.Real.X2 * PS.Real.S;
    Results results = realSubtraction(PS.Soft, s_r, jfks, lumi, ME, x, Xi, y);

    realFillHists(hists, results, rwgt, ps_lab, PS_recomb_lab,
                  born_passcuts, cutspassed, data);
    return results.Combine();
}

/**
 * @brief remnant part of cross section
 *
 * remnant returns the collinear remnant part of the cross section. The return
 * type is an array which elements are the pdf contribution as defined by lumir.
 *
 * @param splitting      splitting of the inital state
 * @param lumir          array with pdf contributions (e.g. ( uu~, cc~ ) )
 * @param ps             born phase space
 * @param x              x integration variable
 * @param xi_max_coll    max value of xi in collinear limit
 * @param alpha_rem      coupling constant (alpha_s or alpha_ew)
 * @param scales         scales of the process
 * @param bornme         born matrix element
 */
Result Integrand::remnant(FKS::Type_t type, const FKS::Splitting &splitting,
                          const LumiRemnants &lumir,
                          const Phasespace::Phasespace &ps, double x,
                          double xi_max_coll, double alpha_rem,
                          PDFRenorm pdfren) const {
    Result results;
    results.fill(0.0);
    double s_b = ps.X1 * ps.X2 * ps.S;
    int i_max = lumir.size();
    double xi = x;
    for (int i = 0; i < i_max; i++) {
        double r = 0.0;
        if (lumir[i].Born == 0.0 && lumir[i].Remnant == 0.0) {
            continue;
        }
        switch (type) {
        case Type_t::EW:
            r = FKS::QED::Remnant(scales, splitting, xi, xi_max_coll, s_b,
                                  alpha_rem, lumir[i], bme.M2, pdfren);
            break;
        case Type_t::QCD:
            r = FKS::QCD::Remnant(scales, splitting, xi, xi_max_coll, s_b,
                                  alpha_rem, lumir[i], bme.M2, pdfren);
            break;
        }
        results[i] = r * CrossSectionToPb;
    }
    return results;
}

template <bool Born, bool Virtual, bool Real, bool Remnant>
Result montecarlo_integrand(const Phasespace::Phasespace &ps, double x1,
                            double x2, double x3, double wgt,
                            UserProcess::Data *userdata) {
    assert(x2 > 0.0 && x2 < 1.0 && "will get numerical problems");
    assert(x1 > 0.0 && "will get numerical problems");
    Result result;
    result.fill(0.0);

    if (Born) {
        userdata->UseBornCache = false;
    }
    Integrand integrand(userdata->pdf);
    integrand.Init(ps, Virtual, wgt, userdata);

    bool born_passcuts = integrand.PSCuts(userdata->cuts);

    if (!born_passcuts && (!Real || !userdata->CutOnReal)) {
        return result;
    }
    /*********************************************************************
     *                          born + virtual                           *
     *********************************************************************/
    if ((Born || Virtual || Remnant) && born_passcuts) {
        LumiBorn lumi =
            integrand.LuminosityBorn(userdata->UsedPDF, ps.X1, ps.X2);

        if (Born) {
            Result result_b = integrand.Born(ps, lumi, userdata->hists);
            add_result(&result, result_b);
        }

        if (Virtual) {
            // calculate part which is integrated over the born phase space
            // (without born contributions).
            Result result_v = integrand.Virtual(ps, lumi, userdata->hists);
            add_result(&result, result_v);
        }

        if (Remnant) {
            // collinear remnants
            Result result_c =
                integrand.CollRemnant(ps, x1, userdata, userdata->hists);
            add_result(&result, result_c);
        }
    }

    /*********************************************************************
     *                               real                                *
     *********************************************************************/
    if (Real) {
        Result result_r = integrand.Real(ps, born_passcuts, x1, x2, x3,
                                         userdata, userdata->hists);
        add_result(&result, result_r);
    }

    for (size_t i = 0; i < result.size(); i++) {
        result[i] *= ps.Jacobian;
    }
    const auto &scale = userdata->Process[userdata->ProcessID].Scales;
    for (size_t i = 0; i < scale.size(); i++) {
        result[i] *= scale[i];
    }

    return result;
}

} // namespace

//
namespace FKS {

Result XSecFullByPDF(const Phasespace::Phasespace &ps, double x1, double x2,
                     double x3, double wgt, UserProcess::Data *params) {
    return montecarlo_integrand<true, true, true, true>(ps, x1, x2, x3, wgt,
                                                        params);
}

Result XSecBornByPDF(const Phasespace::Phasespace &ps, double x1, double x2,
                     double x3, double wgt, UserProcess::Data *params) {
    return montecarlo_integrand<true, false, false, false>(ps, x1, x2, x3, wgt,
                                                           params);
}

Result XSecRealByPDF(const Phasespace::Phasespace &ps, double x1, double x2,
                     double x3, double wgt, UserProcess::Data *params) {
    return montecarlo_integrand<false, false, true, false>(ps, x1, x2, x3, wgt,
                                                           params);
}

Result XSecVirtualByPDF(const Phasespace::Phasespace &ps, double x1, double x2,
                        double x3, double wgt, UserProcess::Data *params) {
    return montecarlo_integrand<false, true, false, false>(ps, x1, x2, x3, wgt,
                                                           params);
}
Result XSecRemnantByPDF(const Phasespace::Phasespace &ps, double x1, double x2,
                        double x3, double wgt, UserProcess::Data *params) {
    return montecarlo_integrand<false, false, false, true>(ps, x1, x2, x3, wgt,
                                                           params);
}

int XSecFull(const Phasespace::Phasespace &ps, double x1, double x2, double x3,
             double wgt, double *out, UserProcess::Data *params) {
    auto result = montecarlo_integrand<true, true, true, true>(ps, x1, x2, x3,
                                                               wgt, params);
    *out = sum_result(result);
    return 0;
}

int XSecReal(const Phasespace::Phasespace &ps, double x1, double x2, double x3,
             double wgt, double *out, UserProcess::Data *params) {
    auto result = montecarlo_integrand<false, false, true, false>(
        ps, x1, x2, x3, wgt, params);
    *out = sum_result(result);
    return 0;
}

int XSecVirtual(const Phasespace::Phasespace &ps, double x1, double x2,
                double x3, double wgt, double *out, UserProcess::Data *params) {
    auto result = montecarlo_integrand<false, true, false, false>(
        ps, x1, x2, x3, wgt, params);
    *out = sum_result(result);
    return 0;
}

int XSecBorn(const Phasespace::Phasespace &ps, double x1, double x2, double x3,
             double wgt, double *out, UserProcess::Data *params) {
    auto result = montecarlo_integrand<true, false, false, false>(
        ps, x1, x2, x3, wgt, params);
    *out = sum_result(result);
    return 0;
}

int XSecRemnant(const Phasespace::Phasespace &ps, double x1, double x2,
                double x3, double wgt, double *out, UserProcess::Data *params) {
    auto result = montecarlo_integrand<false, false, false, true>(
        ps, x1, x2, x3, wgt, params);
    *out = sum_result(result);
    return 0;
}

int XSecVirtualPlusRemnant(const Phasespace::Phasespace &ps, double x1,
                           double x2, double x3, double wgt, double *out,
                           UserProcess::Data *params) {
    auto result = montecarlo_integrand<false, true, false, true>(ps, x1, x2, x3,
                                                                 wgt, params);
    *out = sum_result(result);
    return 0;
}

} // end namespace FKS
