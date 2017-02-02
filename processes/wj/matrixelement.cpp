#include "matrixelement.h"

#include <iostream>

#include "phasespace/phasespace.h"

#include "me/ddx_dux.h"
#include "me/dxdx_uxdx.h"
#include "me/gdx_gux.h"
#include "me/gg_dux.h"
#include "me/gu_gd.h"
#include "me/parameters_sm.h"
#include "me/uc_us.h"
#include "me/ucx_dcx.h"
#include "me/ud_dd.h"
#include "me/udx_ccx.h"
#include "me/udx_ddx.h"
#include "me/udx_gg.h"
#include "me/udx_uux.h"
#include "me/usx_ucx.h"
#include "me/uu_ud.h"
#include "me/uux_dux.h"
#include "me/uux_scx.h"
#include "me/uxdx_uxux.h"
#include "me/uxsx_uxcx.h"
#include "me/udx_ga.h"
#include "me/gu_da.h"
#include "me/gdx_uxa.h"
#include "me/virtual.h"
#include "me/wg.h"
#include "me/wq.h"
#include "me/wqx.h"

void MatrixElementWj::Init(const FKS::Param *param, AlphaScheme scheme,
                           double inveps) {
    auto p = static_cast<const Parameters_sm *>(param);
    Wj_InitDimReg(*p, inveps);
}

UserProcess::IMatrixElement::Result
MatrixElementWj::Born(int process, const Phasespace::Phasespace &ps, double Q,
                      const FKS::Param *param, bool QCD, bool EW) {
    const Parameters_sm *p = static_cast<const Parameters_sm *>(param);
    IMatrixElement::Result result;

    constexpr double NC = 3;

    switch (process) {
    case 0: {
        Wq me;
        int perm[] = {0, 1, 2, 3, 4};
        auto spin = me.Calculate(ps, perm, *p);
        double born = spin.Born();
        result.SpinCorr = spin;
        result.ColorCorr.SetLen(5);
        result.ColorCorr.Set(0, 1, NC / 2.0 * born);
        result.ColorCorr.Set(1, 0, NC / 2.0 * born);
        result.ColorCorr.Set(4, 0, NC / 2.0 * born);
        result.ColorCorr.Set(0, 4, NC / 2.0 * born);
        result.ColorCorr.Set(4, 1, -1 / (2 * NC) * born);
        result.ColorCorr.Set(1, 4, -1 / (2 * NC) * born);
        result.M2 = born;
    } break;
    case 1: {
        Wqx me;
        int perm[] = {0, 1, 2, 3, 4};
        auto spin = me.Calculate(ps, perm, *p);
        double born = spin.Born();
        result.SpinCorr = spin;
        result.ColorCorr.SetLen(5);
        result.ColorCorr.Set(0, 1, NC / 2.0 * born);
        result.ColorCorr.Set(1, 0, NC / 2.0 * born);
        result.ColorCorr.Set(4, 0, NC / 2.0 * born);
        result.ColorCorr.Set(0, 4, NC / 2.0 * born);
        result.ColorCorr.Set(4, 1, -1 / (2 * NC) * born);
        result.ColorCorr.Set(1, 4, -1 / (2 * NC) * born);
        result.M2 = born;
    } break;
    case 2: {
        Wq me;
        int perm[] = {1, 0, 2, 3, 4};
        auto spin = me.Calculate(ps, perm, *p);
        double born = spin.Born();
        result.SpinCorr = spin;
        result.ColorCorr.SetLen(5);
        result.ColorCorr.Set(0, 1, NC / 2.0 * born);
        result.ColorCorr.Set(1, 0, NC / 2.0 * born);
        result.ColorCorr.Set(4, 0, -1 / (2 * NC) * born);
        result.ColorCorr.Set(0, 4, -1 / (2 * NC) * born);
        result.ColorCorr.Set(4, 1, NC / 2.0 * born);
        result.ColorCorr.Set(1, 4, NC / 2.0 * born);
        result.M2 = born;
    } break;
    case 3: {
        Wg me;
        int perm[] = {0, 1, 2, 3, 4};
        auto spin = me.Calculate(ps, perm, *p);
        double born = spin.Born();
        result.SpinCorr = spin;
        result.ColorCorr.SetLen(5);
        result.ColorCorr.Set(0, 1, -1 / (2 * NC) * born);
        result.ColorCorr.Set(1, 0, -1 / (2 * NC) * born);
        result.ColorCorr.Set(4, 1, NC / 2.0 * born);
        result.ColorCorr.Set(1, 4, NC / 2.0 * born);
        result.ColorCorr.Set(4, 0, NC / 2.0 * born);
        result.ColorCorr.Set(0, 4, NC / 2.0 * born);
        result.M2 = born;
    } break;
    case 4: {
        Wqx me;
        int perm[] = {1, 0, 2, 3, 4};
        auto spin = me.Calculate(ps, perm, *p);
        double born = spin.Born();
        result.SpinCorr = spin;
        result.ColorCorr.SetLen(5);
        result.ColorCorr.Set(0, 1, NC / 2.0 * born);
        result.ColorCorr.Set(1, 0, NC / 2.0 * born);
        result.ColorCorr.Set(4, 1, NC / 2.0 * born);
        result.ColorCorr.Set(1, 4, NC / 2.0 * born);
        result.ColorCorr.Set(4, 0, -1 / (2 * NC) * born);
        result.ColorCorr.Set(0, 4, -1 / (2 * NC) * born);
        result.M2 = born;
    } break;
    case 5: {
        Wg me;
        int perm[] = {1, 0, 2, 3, 4};
        auto spin = me.Calculate(ps, perm, *p);
        double born = spin.Born();
        result.SpinCorr = spin;
        result.ColorCorr.SetLen(5);
        result.ColorCorr.Set(0, 1, -1 / (2 * NC) * born);
        result.ColorCorr.Set(1, 0, -1 / (2 * NC) * born);
        result.ColorCorr.Set(4, 1, NC / 2.0 * born);
        result.ColorCorr.Set(1, 4, NC / 2.0 * born);
        result.ColorCorr.Set(4, 0, NC / 2.0 * born);
        result.ColorCorr.Set(0, 4, NC / 2.0 * born);
        result.M2 = born;
    } break;
    }

    if (QCD) {
        int perm[5] = {0, 1, 3, 2, 4};
        double born = 0.0;
        double virt = 0.0;
        int proc = 0;
        switch (process) {
        case 0:
            perm[0] = 1;
            perm[1] = 0;
            proc = 2;
            break;
        case 1:
            perm[0] = 1;
            perm[1] = 0;
            proc = 4;
            break;
        case 2:
            perm[0] = 0;
            perm[1] = 1;
            proc = 2;
            break;
        case 3:
            perm[0] = 0;
            perm[1] = 1;
            proc = 0;
            break;
        case 4:
            perm[0] = 0;
            perm[1] = 1;
            proc = 4;
            break;
        case 5:
            perm[0] = 1;
            perm[1] = 0;
            proc = 0;
            break;
        default:
            assert(0 && "unknown born process");
        }
        Wj_MatrixElement(proc, &born, &virt, ps.Momenta.data(), perm,
                         param->alpha, param->alphaS(), Q, Type::QCD,
                         Scheme::Alpha0);
        result.VfinQCD = virt;

        if (result.M2 > 0.0 && fabs(1.0 - born / result.M2) > 1e-6) {
            printf("process: %d mismatch in born matrix element %g <-> %g "
                   "(diff = %g)\n",
                   process, born, result.M2, 1.0 - born / result.M2);
        }
    }

    if (EW) {
        int perm[5] = {0, 1, 3, 2, 4};
        double born = 0.0;
        double virt = 0.0;
        int proc = 0;
        switch (process) {
        case 0:
            perm[0] = 1;
            perm[1] = 0;
            proc = 2;
            break;
        case 1:
            perm[0] = 1;
            perm[1] = 0;
            proc = 4;
            break;
        case 2:
            perm[0] = 0;
            perm[1] = 1;
            proc = 2;
            break;
        case 3:
            perm[0] = 0;
            perm[1] = 1;
            proc = 0;
            break;
        case 4:
            perm[0] = 0;
            perm[1] = 1;
            proc = 4;
            break;
        case 5:
            perm[0] = 1;
            perm[1] = 0;
            proc = 0;
            break;
        default:
            assert(0 && "unknown born process");
        }
        Wj_MatrixElement(proc, &born, &virt, ps.Momenta.data(), perm,
                         param->alpha, param->alphaS(), Q, Type::EW,
                         Scheme::Gmu);
        result.VfinEW = virt;
    }

    return result;
}

double MatrixElementWj::Real(int process, const Phasespace::Phasespace &ps,
                             const FKS::Param *param, Diagrams diagrams) {

    const Parameters_sm *p = static_cast<const Parameters_sm *>(param);
    switch (process) {
    case 0: {
        UUX_SCX me; // 12
        int perm[] = {1, 0, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 1: {
        UUX_DUX me; // 4
        int perm[] = {1, 0, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 2: {
        DXDX_UXDX me; // 9
        int perm[] = {0, 1, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 3: {
        UC_US me; // 10
        int perm[] = {0, 1, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 4: {
        UDX_DDX me; // 6
        int perm[] = {0, 1, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 5: {
        UDX_DDX me; // 6
        int perm[] = {1, 0, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 6: {
        UDX_UUX me; // 5
        int perm[] = {1, 0, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 7: {
        DDX_DUX me; // 7
        int perm[] = {0, 1, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 8: {
        UXSX_UXCX me; // 16
        int perm[] = {1, 0, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 9: {
        UXDX_UXUX me; // 8
        int perm[] = {1, 0, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 10: {
        UUX_DUX me; // 4
        int perm[] = {0, 1, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 11: {
        GDX_GUX me; // 19
        int perm[] = {0, 1, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 12: {
        UXSX_UXCX me; // 16
        int perm[] = {0, 1, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 13: {
        UDX_CCX me; // 14
        int perm[] = {0, 1, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 14: {
        GDX_GUX me; // 19
        int perm[] = {1, 0, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 15: {
        UDX_DDX me; // 6
        int perm[] = {1, 0, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 16: {
        USX_UCX me; // 15
        int perm[] = {1, 0, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 17: {
        USX_UCX me; // 15
        int perm[] = {0, 1, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 18: {
        UU_UD me; // 2
        int perm[] = {0, 1, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 19: {
        UCX_DCX me; // 13
        int perm[] = {1, 0, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 20: {
        UDX_DDX me; // 6
        int perm[] = {0, 1, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 21: {
        UDX_GG me; // 1
        int perm[] = {1, 0, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 22: {
        GDX_GUX me; // 19
        int perm[] = {1, 0, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 23: {
        GU_GD me; // 18
        int perm[] = {0, 1, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 24: {
        DDX_DUX me; // 7
        int perm[] = {1, 0, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 25: {
        UDX_UUX me; // 5
        int perm[] = {0, 1, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 26: {
        GG_DUX me; // 20
        int perm[] = {0, 1, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 27: {
        UC_US me; // 10
        int perm[] = {1, 0, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 28: {
        UD_DD me; // 3
        int perm[] = {1, 0, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 29: {
        UUX_SCX me; // 12
        int perm[] = {0, 1, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 30: {
        GDX_GUX me; // 19
        int perm[] = {0, 1, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 31: {
        GU_GD me; // 18
        int perm[] = {1, 0, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 32: {
        UDX_GG me; // 1
        int perm[] = {0, 1, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 33: {
        UD_DD me; // 3
        int perm[] = {0, 1, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 34: {
        UXDX_UXUX me; // 8
        int perm[] = {0, 1, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 35: {
        GG_DUX me; // 20
        int perm[] = {0, 1, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 36: {
        GU_GD me; // 18
        int perm[] = {1, 0, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 37: {
        UDX_CCX me; // 14
        int perm[] = {1, 0, 2, 3, 5, 4};
        return me.Calculate(ps, perm, *p);
    }
    case 38: {
        UCX_DCX me; // 13
        int perm[] = {0, 1, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 39: {
        GU_GD me; // 18
        int perm[] = {0, 1, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    // --- EW ---------------------------------------
    case 40: {
        UDX_GA me;
        int perm[] = {0, 1, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 41: {
        UDX_GA me;
        int perm[] = {1, 0, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 42: {
        GU_DA me;
        int perm[] = {0, 1, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 43: {
        GU_DA me;
        int perm[] = {1, 0, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 44: {
        GDX_UXA me;
        int perm[] = {0, 1, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    case 45: {
        GDX_UXA me;
        int perm[] = {1, 0, 2, 3, 4, 5};
        return me.Calculate(ps, perm, *p);
    }
    default:
        assert(0 && "matrix element id not implemented");
    }

    return -1.0;
}
