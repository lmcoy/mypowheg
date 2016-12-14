#include "matrixelement.h"

#include "phasespace/phasespace.h"

#include "me/parameters_sm.h"

#include "me/virtual_qed.h"
#include "me/uux_mupmuma.h"
#include "me/ddx_mupmuma.h"

#include "me/uux_mupmumg.h"
#include "me/gu_mupmumu.h"
#include "me/gux_mupmumux.h"

#include "me/ddx_mupmumg.h"
#include "me/gd_mupmumd.h"
#include "me/gdx_mupmumdx.h"

void DrellYanME::Init(const FKS::Param *param, AlphaScheme scheme,
                      double inveps) {
    auto p = static_cast<const Parameters_sm *>(param);
    double c = 0.0;
    InitVfin_qxq_lxl(*p, scheme, inveps, &c);
    counterterm = c;
    initialized = true;
}

UserProcess::IMatrixElement::Result
DrellYanME::Born(int process, const Phasespace::Phasespace &ps, double Q,
                 const FKS::Param *param, bool QCD, bool EW) {

    assert(initialized && "call Init() first!");
    int perm[4] = {0};
    perm[2] = 2;
    perm[3] = 3;
    int proc_nb = 0;
    assert(process >= 0 && process < 4);
    switch (process) {
    case 0:
        perm[0] = 0;
        perm[1] = 1;
        proc_nb = 2;
        break;
    case 1:
        perm[0] = 1;
        perm[1] = 0;
        proc_nb = 2;
        break;
    case 2:
        perm[0] = 0;
        perm[1] = 1;
        proc_nb = 1;
        break;
    case 3:
        perm[0] = 1;
        perm[1] = 0;
        proc_nb = 1;
        break;
    default:
        break;
    }

    double M2 = 0.0;
    double VfinQCD = 0.0;
    double VfinEW = 0.0;
    Vfin_qxq_lxl(&M2, &VfinQCD, ps.Momenta.data(), perm, CorrectionType::None,
                 proc_nb, param->alpha, param->alphaS(), Q * Q, 0.0);
    if (QCD) {
        Vfin_qxq_lxl(&M2, &VfinQCD, ps.Momenta.data(), perm,
                     CorrectionType::QCD, proc_nb, param->alpha,
                     param->alphaS(), Q * Q, 0.0);
    }
    if (EW) {
        double me_tmp = 0.0;
        Vfin_qxq_lxl(&me_tmp, &VfinEW, ps.Momenta.data(), perm,
                     CorrectionType::EW, proc_nb, param->alpha, param->alphaS(),
                     Q * Q, counterterm);
        assert(me_tmp - M2 < 1e-8);
    }

    IMatrixElement::Result me;
    me.M2 = M2;
    me.VfinQCD = VfinQCD;
    me.VfinEW = VfinEW;
    me.ColorCorr.SetLen(4);
    me.ColorCorr.Set(0, 1, 4.0 / 3.0 * M2);
    me.ColorCorr.Set(1, 0, 4.0 / 3.0 * M2);
    return me;
}

double DrellYanME::Real(int process, const Phasespace::Phasespace &ps,
                          const FKS::Param *param,
                          Diagrams diagrams) {
    using R = SubProcesses;
    const Parameters_sm *p = static_cast<const Parameters_sm *>(param);
    if (diagrams == Diagrams::ONLYFSR) {
        switch (process) {
        case R::UXU_MUXMU_G:
        case R::UXG_MUXMU_UX:
        case R::GU_MUXMU_U:
        case R::UUX_MUXMU_G:
        case R::UG_MUXMU_U:
        case R::GUX_MUXMU_UX:
        case R::DXD_MUXMU_G:
        case R::DXG_MUXMU_DX:
        case R::GD_MUXMU_D:
        case R::DDX_MUXMU_G:
        case R::DG_MUXMU_D:
        case R::GDX_MUXMU_DX:
            // QCD has only ISR diagrams
            return 0.0;
        }
    }
    switch (process) {
    case R::UXU_MUXMU_G: {
        ME_uux_mupmumg me;
        int perm[] = { 1, 0, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    case R::UXG_MUXMU_UX: {
        ME_gux_mupmumux me;
        int perm[] = { 1, 0, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    case R::GU_MUXMU_U: {
        ME_gu_mupmumu me2;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me2.Calculate(ps, perm, *p);
    }
    case R::UUX_MUXMU_G: {
        ME_uux_mupmumg me;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    case R::UG_MUXMU_U: {
        ME_gu_mupmumu me;
        int perm[] = { 1, 0, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    case R::GUX_MUXMU_UX: {
        ME_gux_mupmumux me;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    // down type quarks
    case R::DXD_MUXMU_G: {
        ME_ddx_mupmumg me;
        int perm[] = { 1, 0, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    case R::DXG_MUXMU_DX: {
        ME_gdx_mupmumdx me;
        int perm[] = { 1, 0, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    case R::GD_MUXMU_D: {
        ME_gd_mupmumd me2;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me2.Calculate(ps, perm, *p);
    }
    case R::DDX_MUXMU_G: {
        ME_ddx_mupmumg me;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    case R::DG_MUXMU_D: {
        ME_gd_mupmumd me;
        int perm[] = { 1, 0, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    case R::GDX_MUXMU_DX: {
        ME_gdx_mupmumdx me;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me.Calculate(ps, perm, *p);
    }
    // ---- EW ----
    case R::UXU_MUXMU_A: {
        ME_uux_mupmuma me;
        int perm[] = { 1, 0, 2, 3, 4 };
        int flag = ME_uux_mupmuma::DEFAULT;
        switch (diagrams) {
        case Diagrams::ONLYISR:
            flag = ME_uux_mupmuma::ONLYISR;
            break;
        case Diagrams::ONLYFSR:
            flag = ME_uux_mupmuma::ONLYFSR;
            break;
        case Diagrams::ALL:
            flag = ME_uux_mupmuma::DEFAULT;
            break;
        }

        return me.Calculate(ps, perm, *p, flag);
    }
    case R::UUX_MUXMU_A: {
        ME_uux_mupmuma me;
        int perm[] = { 0, 1, 2, 3, 4 };
        int flag = ME_uux_mupmuma::DEFAULT;
        switch (diagrams) {
        case Diagrams::ONLYISR:
            flag = ME_uux_mupmuma::ONLYISR;
            break;
        case Diagrams::ONLYFSR:
            flag = ME_uux_mupmuma::ONLYFSR;
            break;
        case Diagrams::ALL:
            flag = ME_uux_mupmuma::DEFAULT;
            break;
        }
        return me.Calculate(ps, perm, *p, flag);
    }
    // down type quarks
    case R::DXD_MUXMU_A: {
        ME_ddx_mupmuma me;
        int perm[] = { 1, 0, 2, 3, 4 };
        int flag = ME_ddx_mupmuma::DEFAULT;
        switch (diagrams) {
        case Diagrams::ONLYISR:
            flag = ME_ddx_mupmuma::ONLYISR;
            break;
        case Diagrams::ONLYFSR:
            flag = ME_ddx_mupmuma::ONLYFSR;
            break;
        case Diagrams::ALL:
            flag = ME_ddx_mupmuma::DEFAULT;
            break;
        }
        return me.Calculate(ps, perm, *p, flag);
    }
    case R::DDX_MUXMU_A: {
        ME_ddx_mupmuma me;
        int perm[] = { 0, 1, 2, 3, 4 };
        int flag = ME_ddx_mupmuma::DEFAULT;
        switch (diagrams) {
        case Diagrams::ONLYISR:
            flag = ME_ddx_mupmuma::ONLYISR;
            break;
        case Diagrams::ONLYFSR:
            flag = ME_ddx_mupmuma::ONLYFSR;
            break;
        case Diagrams::ALL:
            flag = ME_ddx_mupmuma::DEFAULT;
            break;
        }
        return me.Calculate(ps, perm, *p, flag);
    }
    }
    assert(-1);
    return -1.0;
}
