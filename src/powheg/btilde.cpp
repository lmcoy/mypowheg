#include <map>

#include "powheg/btilde.h"

#include "phasespace/phasespace.h"
#include "process/data.h"
#include "fks/xsec.h"

namespace Powheg {

int Btilde(const Phasespace::Phasespace &ps, double x1, double x2, double x3,
              double wgt, double *out, UserProcess::Data *params) {

    double result = 0.0;
    double result_born = 0.0;
    double voverbmax = 0.0;
    bool onlyborn = params->BornOnly;
    for (size_t i = 0; i < params->Process.size(); i++) {
        params->ProcessID = i;
        params->UseBornCache = false;
        auto born = FKS::XSecBornByPDF(ps, x1, x2, x3, wgt, params);
        FKS::Result real;
        real.fill(0.0);
        FKS::Result remn;
        remn.fill(0.0);
        if (!onlyborn) {
            params->UseBornCache = true;
            real = FKS::XSecRealByPDF(ps, x1, x2, x3, wgt, params);
            params->UseBornCache = true;
            remn = FKS::XSecRemnantByPDF(ps, x1, x2, x3, wgt, params);
        }
        FKS::Result virt;
        virt.fill(0.0);
        if (!params->IgnoreVirtualInBtilde && !onlyborn) {
            params->UseBornCache = false;
            virt = FKS::XSecVirtualByPDF(ps, x1, x2, x3, wgt, params);
        }
        params->UseBornCache = false;
        assert(born.size() == virt.size() && born.size()== remn.size());
        assert(born.size() == real.size());
        for(size_t n = 0; n < real.size(); n++) {
            result += real[n];
        }
        for(size_t n = 0; n < born.size(); n++) {
            result += born[n];
            result_born += born[n];
            result += virt[n];
            result += remn[n];
            if (born[n] > 0.0 && virt[n] / born[n] > voverbmax) {
                voverbmax = virt[n] / born[n];
            }
        }
    }
    *out = result;
    if (params->BtildeState.UpdateMax) {
        params->BtildeState.FillMax(wgt * result);
        double *max = &params->BtildeState.Max;
        if (*max == 0.0 || (wgt * result > *max && wgt * result < *max * 2.0)) {
            *max = wgt * result;
        }
        double maxb = params->BtildeState.MaxBorn;
        if (maxb == 0.0 || (wgt * result_born > maxb && wgt * result_born < maxb * 2.0)) {
            params->BtildeState.MaxBorn = wgt * result_born;
        }
        if (voverbmax > 0.0 && voverbmax > params->BtildeState.MaxVoverB) {
            params->BtildeState.MaxVoverB = voverbmax;
        }
    }

    return 0;
}

double Btilde_t::Born(const Phasespace::Phasespace &ps, double x1, double x2,
                      double x3, double wgt, UserProcess::Data *params) {
    double xsec = 0.0;
    for (size_t i = 0; i < params->Process.size(); i++) {
        params->ProcessID = i;
        auto born = FKS::XSecBornByPDF(ps, x1, x2, x3, wgt, params);
        for (size_t pdf = 0; pdf < born.size(); pdf++) {
            xsec += born[pdf];
        }
    }
    return xsec;
}

void Btilde_t::CalcWOVirtual(const Phasespace::Phasespace &ps, double x1,
                             double x2, double x3, double wgt,
                             UserProcess::Data *params) {
    bool onlyborn = params->BornOnly;
    bool negative = false;
    for (size_t i = 0; i < params->Process.size(); i++) {
        params->ProcessID = i;
        auto born = FKS::XSecBornByPDF(ps, x1, x2, x3, wgt, params);
        FKS::Result real;
        real.fill(0.0);
        FKS::Result remn;
        remn.fill(0.0);
        if (!onlyborn) {
            params->UseBornCache = true;
            real = FKS::XSecRealByPDF(ps, x1, x2, x3, wgt, params);
            params->UseBornCache = true;
            remn = FKS::XSecRemnantByPDF(ps, x1, x2, x3, wgt, params);
        }
        params->UseBornCache = false;
        for (size_t pdf = 0; pdf < born.size(); pdf++) {
            double xsec = born[pdf] + real[pdf] + remn[pdf];
            double F = params->BtildeState.MaxVoverB;
            if (onlyborn) {
                F = 0.0;
            }
            double xsec_guessed_v = xsec + born[pdf] * F;
            if (xsec < 0.0) {
                // if xsec < 0, xsec_guessed_v could be >0 due to an
                // overestimated virtual part. 
                negative = true;
            }
            if (xsec != 0.0) {
                Append(i, pdf, xsec * wgt, xsec_guessed_v * wgt);
            }
        }
    }
    if (negative) {
        sum_guessed = -1.0;
    }
}

void Btilde_t::CalcVirtual(const Phasespace::Phasespace &ps, double x1,
                           double x2, double x3, double wgt,
                           UserProcess::Data *params) {
    if (params->BornOnly) {
        has_virtual = true;
        return;
    }
    std::map<Meta, double> cache;
    for (size_t i = 0; i < meta.size(); i++) {
        auto prev = cache.find(meta[i]);
        if (prev != cache.end()) {
            data[i] += prev->second;
            sum += prev->second;
            continue;
        }
        auto processID = meta[i].first;
        auto pdfIndex = meta[i].second;
        params->ProcessID = processID;
        params->UseBornCache = false;
        auto virt = FKS::XSecVirtualByPDF(ps, x1, x2, x3, wgt, params);
        for (size_t p = 0; p < virt.size(); p++) {
            if (p == pdfIndex) {
                data[i] += virt[p] * wgt;
                sum += virt[p] * wgt;
            } else {
                cache[Meta(processID, p)] = virt[p] * wgt;
            }
        }
    }
    has_virtual = true;
}

} // end namespace Powheg
