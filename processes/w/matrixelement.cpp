#include "matrixelement.h"

#include <iostream>

#include "phasespace/phasespace.h"

#include "me/parameters_sm.h"

#include "me/virtual.h"
#include "me/wg.h"
#include "me/wq.h"
#include "me/wa.h"
#include "me/wqx.h"

void MatrixElementW::Init(const FKS::Param *param, AlphaScheme scheme,
                          double inveps) {
    auto p = static_cast<const Parameters_sm *>(param);
    W_InitDimReg(*p, inveps);
}

UserProcess::IMatrixElement::Result
MatrixElementW::Born(int process, const Phasespace::Phasespace &ps, double Q,
                     const FKS::Param *param, bool QCD, bool EW) {
    using R = SubProcesses;
    int perm[4] = { 0,1,2,3 };
    switch (process) {
    case R::QXQ_MUPNU:
        perm[0] = 0;
        perm[1] = 1;
        break;
    case R::QQX_MUPNU:
        perm[0] = 1;
        perm[1] = 0;
        break;
    default:
        assert(0 && "unknown born matrix element");
    }
    double M2 = 0.0;
    double VfinQCD = 0.0;
    double VfinEW = 0.0;
    if (!QCD && !EW) {
        W_MatrixElement(&M2, &VfinQCD, ps.Momenta.data(), perm, param->alpha,
                        param->alphaS(), Q, Type::None, Scheme::Alpha0);
    }
    if (QCD) {
        W_MatrixElement(&M2, &VfinQCD, ps.Momenta.data(), perm, param->alpha,
                        param->alphaS(), Q, Type::QCD, Scheme::Alpha0);
    }
    if (EW) {
        W_MatrixElement(&M2, &VfinEW, ps.Momenta.data(), perm, param->alpha,
                        param->alphaS(), Q, Type::EW, Scheme::Gmu);
    }
    //double born = M2;
    //double s = ps.X1 * ps.X2 *ps.S;
    //double tmp = log(Q*Q/s);
    //double vfin = 4.0 / 3.0 * (M_PI * M_PI - 8.0 - tmp * (tmp + 3.0)) * born;
    IMatrixElement::Result me;
    me.M2 = M2;
    me.VfinQCD = VfinQCD;
    me.VfinEW = VfinEW;
    me.ColorCorr.SetLen(4);
    me.ColorCorr.Set(0, 1, 4.0 / 3.0 * M2);
    me.ColorCorr.Set(1, 0, 4.0 / 3.0 * M2);
    return me;
}

double MatrixElementW::Real(int process, const Phasespace::Phasespace &ps,
                            const FKS::Param *param, Diagrams diagrams) {
    using R = SubProcesses;
    const Parameters_sm *p = static_cast<const Parameters_sm *>(param);
    if (diagrams == Diagrams::ONLYFSR) {
        switch (process) {
        case R::QXG_MUPNU_QX:
        case R::QXQ_MUPNU_G:
        case R::GQX_MUPNU_QX:
        case R::GQ_MUPNU_Q:
        case R::QQX_MUPNU_G:
        case R::QG_MUPNU_Q:
            // QCD has only ISR diagrams
            return 0.0;
        }
    }
    switch (process) {
    // --- QCD ----------------------------------------------
    case R::QXG_MUPNU_QX: {
        Wqx me;
        int perm[] = { 1, 0, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    case R::QXQ_MUPNU_G: {
        Wg me;
        int perm[] = { 1, 0, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    case R::GQX_MUPNU_QX: {
        Wqx me;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    case R::GQ_MUPNU_Q: {
        Wq me;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    case R::QQX_MUPNU_G: {
        Wg me;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    case R::QG_MUPNU_Q: {
        Wq me;
        int perm[] = { 1, 0, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    // --- EW -----------------------------------------------
    case R::QXQ_MUPNU_A: {
        Wa me;
        int perm[] = { 1, 0, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    case R::QQX_MUPNU_A: {
        Wa me;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    }
    assert(-1);
    return -1.0;
}
