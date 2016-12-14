#include "processes.h"
#include "process/matrixelement.h"
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

struct R {
    // QCD
    static constexpr int UXU_MUXMU_G = 0;
    static constexpr int UXG_MUXMU_UX = 1;
    static constexpr int GU_MUXMU_U =2;
    static constexpr int UUX_MUXMU_G = 3;
    static constexpr int UG_MUXMU_U = 4;
    static constexpr int GUX_MUXMU_UX = 5;
    static constexpr int DXD_MUXMU_G = 6;
    static constexpr int DXG_MUXMU_DX = 7;
    static constexpr int GD_MUXMU_D = 8;
    static constexpr int DDX_MUXMU_G = 9;
    static constexpr int DG_MUXMU_D = 10;
    static constexpr int GDX_MUXMU_DX = 11;
    // EW
    static constexpr int UXU_MUXMU_A = 12;
    static constexpr int UUX_MUXMU_A = 13;
    static constexpr int DXD_MUXMU_A = 14;
    static constexpr int DDX_MUXMU_A = 15;
};

FKS::ProcessList GenerateProcesses(bool QCD, bool EW) {
    FKS::ProcessList list;
    FKS::FlavourConfig fl1(0, { { -2, 2, -13, 13 } }, { { -2, 2, -4, 4 } });
    if (EW) {
        FKS::ResonanceList resonances;
        size_t isr_res = resonances.Add(
            FKS::Resonance({ { 2, 3 } }, 91.1876, 2.4952, 2.0 / 3.0));
        size_t fsr_res = resonances.Add(
            FKS::Resonance({ { 2, 3, 4 } }, 91.1876, 2.4952, 1.0));
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 2, fsr_res));
        rlist.push_back(FKS::Region(4, 3, fsr_res));
        rlist.push_back(FKS::Region(4, 0, isr_res));
        fl1.AddReal(R::UXU_MUXMU_A, FKS::Type_t::EW, { { -2, 2, -13, 13, 22 } },
                    { { -2, 2, -4, 4 } }, rlist, resonances);

    }
    if (QCD) {
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 0));
        fl1.AddReal(R::UXU_MUXMU_G, FKS::Type_t::QCD,
                    { { -2, 2, -13, 13, 21 } }, { { -2, 2, -4, 4 } }, rlist);
        fl1.AddReal(R::UXG_MUXMU_UX, FKS::Type_t::QCD,
                    { { -2, 21, -13, 13, -2 } }, { { -2, 0, -4, 0 } }, rlist);
        fl1.AddReal(R::GU_MUXMU_U, FKS::Type_t::QCD, { { 21, 2, -13, 13, 2 } },
                    { { 0, 2, 0, 4 } }, rlist);
    }
    list.push_back(fl1);

    FKS::FlavourConfig fl2(1, { { 2, -2, -13, 13 } }, { { 2, -2, 4, -4 } });
    if (EW) {
        FKS::ResonanceList resonances;
        size_t isr_res = resonances.Add(
            FKS::Resonance({ { 2, 3 } }, 91.1876, 2.4952, 2.0 / 3.0));
        size_t fsr_res = resonances.Add(
            FKS::Resonance({ { 2, 3, 4 } }, 91.1876, 2.4952, 1.0));
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 2, fsr_res));
        rlist.push_back(FKS::Region(4, 3, fsr_res));
        rlist.push_back(FKS::Region(4, 0, isr_res));
        fl2.AddReal(R::UUX_MUXMU_A, FKS::Type_t::EW, { { 2, -2, -13, 13, 22 } },
                    { { 2, -2, 4, -4 } }, rlist, resonances);
    }
    if (QCD) {
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 0));
        fl2.AddReal(R::UUX_MUXMU_G, FKS::Type_t::QCD,
                    { { 2, -2, -13, 13, 21 } }, { { 2, -2, 4, -4 } }, rlist);
        fl2.AddReal(R::UG_MUXMU_U, FKS::Type_t::QCD, { { 2, 21, -13, 13, 2 } },
                    { { 2, 0, 4, 0 } }, rlist);
        fl2.AddReal(R::GUX_MUXMU_UX, FKS::Type_t::QCD,
                    { { 21, -2, -13, 13, -2 } }, { { 0, -2, 0, -4 } }, rlist);
    }
    list.push_back(fl2);

    FKS::FlavourConfig fl3(2, { { -1, 1, -13, 13 } },
                           { { -1, 1, -3, 3, -5, 5 } });
    if (EW) {
        FKS::ResonanceList resonances;
        size_t isr_res = resonances.Add(
            FKS::Resonance({ { 2, 3 } }, 91.1876, 2.4952, 1.0 / 3.0));
        size_t fsr_res = resonances.Add(
            FKS::Resonance({ { 2, 3, 4 } }, 91.1876, 2.4952, 1.0));
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 2, fsr_res));
        rlist.push_back(FKS::Region(4, 3, fsr_res));
        rlist.push_back(FKS::Region(4, 0, isr_res));
        fl3.AddReal(R::DXD_MUXMU_A, FKS::Type_t::EW, { { -1, 1, -13, 13, 22 } },
                    { { -1, 1, -3, 3, -5, 5 } }, rlist, resonances);
    }
    if (QCD) {
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 0));
        fl3.AddReal(R::DXD_MUXMU_G, FKS::Type_t::QCD,
                    { { -1, 1, -13, 13, 21 } }, { { -1, 1, -3, 3, -5, 5 } },
                    rlist);
        fl3.AddReal(R::DXG_MUXMU_DX, FKS::Type_t::QCD,
                    { { -1, 21, -13, 13, -1 } }, { { -1, 0, -3, 0, -5, 0 } },
                    rlist);
        fl3.AddReal(R::GD_MUXMU_D, FKS::Type_t::QCD, { { 21, 1, -13, 13, 1 } },
                    { { 0, 1, 0, 3, 0, 5 } }, rlist);
    }
    list.push_back(fl3);

    FKS::FlavourConfig fl4(3, { { 1, -1, -13, 13 } },
                           { { 1, -1, 3, -3, 5, -5 } });
    if (EW) {
        FKS::ResonanceList resonances;
        size_t isr_res = resonances.Add(
            FKS::Resonance({ { 2, 3 } }, 91.1876, 2.4952, 1.0 / 3.0));
        size_t fsr_res = resonances.Add(
            FKS::Resonance({ { 2, 3, 4 } }, 91.1876, 2.4952, 1.0));
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 2, fsr_res));
        rlist.push_back(FKS::Region(4, 3, fsr_res));
        rlist.push_back(FKS::Region(4, 0, isr_res));
        fl4.AddReal(R::DDX_MUXMU_A, FKS::Type_t::EW, { { 1, -1, -13, 13, 22 } },
                    { { 1, -1, 3, -3, 5, -5 } }, rlist, resonances);
    }
    if (QCD) {
        FKS::RegionList rlist;
        // add all real flavour structures
        rlist.push_back(FKS::Region(4, 0));
        fl4.AddReal(R::DDX_MUXMU_G, FKS::Type_t::QCD,
                    { { 1, -1, -13, 13, 21 } }, { { 1, -1, 3, -3, 5, -5 } },
                    rlist);
        fl4.AddReal(R::DG_MUXMU_D, FKS::Type_t::QCD, { { 1, 21, -13, 13, 1 } },
                    { { 1, 0, 3, 0, 5, 0 } }, rlist);
        fl4.AddReal(R::GDX_MUXMU_DX, FKS::Type_t::QCD,
                    { { 21, -1, -13, 13, -1 } }, { { 0, -1, 0, -3, 0, -5 } },
                    rlist);
    }
    list.push_back(fl4);

    return list;
}

void InitBornME(const FKS::Param *param, const FKS::Param_as *param_as,
                AlphaScheme scheme, double mu2, double inveps, double *counterterm) {
    const Parameters_sm *p = dynamic_cast<const Parameters_sm *>(param);
    InitVfin_qxq_lxl(*p, scheme, inveps, counterterm);
}
// Calculate the born matrix element, the finite part of the virtual part and
// the color correlated born amplitude.
int BornME(int process, const Phasespace::Phasespace &ps, double Q,
           const FKS::Param *param, const FKS::Param_as *param_aS,
           BornMEOut *me, bool QCD, bool EW, double counterterm) {
    int perm[4] = {0};
    perm[2] = 2;
    perm[3] = 3;
    int proc_nb;
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
        return -1;
    }

    Vfin_qxq_lxl(&me->M2, &me->VfinQCD, ps.Momenta.data(), perm,
                 CorrectionType::None, proc_nb, param->alpha, param_aS->aS,
                 Q * Q, 0.0);
    if (QCD) {
        Vfin_qxq_lxl(&me->M2, &me->VfinQCD, ps.Momenta.data(), perm,
                     CorrectionType::QCD, proc_nb, param->alpha, param_aS->aS, Q*Q, 0.0);
    }
    if (EW) {
        double me_tmp = 0.0;
        Vfin_qxq_lxl(&me_tmp, &me->VfinEW, ps.Momenta.data(), perm,
                     CorrectionType::EW, proc_nb, param->alpha, param_aS->aS, Q*Q, counterterm);
        assert(me_tmp - me->M2 < 1e-8);
    }
    // double Q2 = Q * Q;
    // double s = ps.S * ps.X1 * ps.X2;
    // double tmp = log(Q2 / s);
    me->ColorCorr.SetLen(4);
    me->ColorCorr.Set(0, 1, 4.0 / 3.0 * me->M2);
    me->ColorCorr.Set(1, 0, 4.0 / 3.0 * me->M2);
        // me->Vfin = 4.0 / 3.0 * (M_PI * M_PI - 8.0 - tmp * (tmp + 3.0)) * me->M2;
    return 0;
}

// Calculate the real matrix element
double RealME(int process, const Phasespace::Phasespace &ps,
              const FKS::Param *param, const FKS::Param_as *param_aS,
              Diagrams diagrams) {
    const Parameters_sm *p = dynamic_cast<const Parameters_sm *>(param);
    const Parameters_alphaS * p_as = (Parameters_alphaS*)param_aS;
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
        return me.Calculate(ps, perm, *p, *p_as);
    }
    case R::UXG_MUXMU_UX: {
        ME_gux_mupmumux me;
        int perm[] = { 1, 0, 2, 3, 4 };
        return me.Calculate(ps, perm, *p, *p_as);
    }
    case R::GU_MUXMU_U: {
        ME_gu_mupmumu me2;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me2.Calculate(ps, perm, *p, *p_as);
    }
    case R::UUX_MUXMU_G: {
        ME_uux_mupmumg me;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me.Calculate(ps, perm, *p, *p_as);
    }
    case R::UG_MUXMU_U: {
        ME_gu_mupmumu me;
        int perm[] = { 1, 0, 2, 3, 4 };
        return me.Calculate(ps, perm, *p, *p_as);
    }
    case R::GUX_MUXMU_UX: {
        ME_gux_mupmumux me;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me.Calculate(ps, perm, *p, *p_as);
    }
    // down type quarks
    case R::DXD_MUXMU_G: {
        ME_ddx_mupmumg me;
        int perm[] = { 1, 0, 2, 3, 4 };
        return me.Calculate(ps, perm, *p, *p_as);
    }
    case R::DXG_MUXMU_DX: {
        ME_gdx_mupmumdx me;
        int perm[] = { 1, 0, 2, 3, 4 };
        return me.Calculate(ps, perm, *p, *p_as);
    }
    case R::GD_MUXMU_D: {
        ME_gd_mupmumd me2;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me2.Calculate(ps, perm, *p, *p_as);
    }
    case R::DDX_MUXMU_G: {
        ME_ddx_mupmumg me;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me.Calculate(ps, perm, *p, *p_as);
    }
    case R::DG_MUXMU_D: {
        ME_gd_mupmumd me;
        int perm[] = { 1, 0, 2, 3, 4 };
        return me.Calculate(ps, perm, *p, *p_as);
    }
    case R::GDX_MUXMU_DX: {
        ME_gdx_mupmumdx me;
        int perm[] = { 0, 1, 2, 3, 4 };
        return me.Calculate(ps, perm, *p, *p_as);
    }
    // ---- EW ----
    case R::UXU_MUXMU_A: {
        ME_uux_mupmuma me;
        int perm[] = { 1, 0, 2, 3, 4 };
        int flag = ME_uux_mupmuma::DEFAULT;
        switch(diagrams) { 
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

        return me.Calculate(ps, perm, *p, *p_as, flag);
    }
    case R::UUX_MUXMU_A: {
        ME_uux_mupmuma me;
        int perm[] = { 0, 1, 2, 3, 4 };
        int flag = ME_uux_mupmuma::DEFAULT;
        switch(diagrams) { 
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
        return me.Calculate(ps, perm, *p, *p_as, flag);
    }
    // down type quarks
    case R::DXD_MUXMU_A: {
        ME_ddx_mupmuma me;
        int perm[] = { 1, 0, 2, 3, 4 };
        int flag = ME_ddx_mupmuma::DEFAULT;
        switch(diagrams) { 
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
        return me.Calculate(ps, perm, *p, *p_as, flag);
    }
    case R::DDX_MUXMU_A: {
        ME_ddx_mupmuma me;
        int perm[] = { 0, 1, 2, 3, 4 };
        int flag = ME_ddx_mupmuma::DEFAULT;
        switch(diagrams) { 
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
        return me.Calculate(ps, perm, *p, *p_as, flag);
    }
    }
    assert(-1.0);
    return -1.0;
}

            

