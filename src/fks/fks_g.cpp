#include <algorithm>

#include "fks_g.h"

#include "fks/process.h"
#include "phasespace/phasespace.h"
#include "fks/sfunctions.h"
#include "process/matrixelement.h"
#include "fks/param.h"
#include "fks/ximax.h"
#include "process/data.h"
#include "util/matrix.h"
#include "fks/limits.h"
#include "fks/phasespaces.h"

namespace {

enum class Type { QED, QCD };

template <Type type, typename... Args>
double softlimit(Args... args) {
    switch (type) {
    case Type::QED:
        return FKS::QED::SoftLimit(args...);
    case Type::QCD:
        return FKS::QCD::SoftLimit(args...);
    }
    assert(0.0);
    return 0.0;
}

template <Type type, typename... Args>
double collinearLimitFSR(Args... args) {
    switch (type) {
    case Type::QED:
        return FKS::QED::CollinearLimitFSR(args...);
    case Type::QCD:
        return FKS::QCD::CollinearLimitFSR(args...);
    }
    assert(0.0);
    return 0.0;
}

template <Type type, typename... Args>
double collinearLimitISR(Args... args) {
    switch (type) {
    case Type::QED:
        return FKS::QED::CollinearLimitISR(args...);
    case Type::QCD:
        return FKS::QCD::CollinearLimitISR(args...);
    }
    assert(0.0);
    return 0.0;
}

template <Type type, typename... Args>
double softCollinearLimitFSR(Args... args) {
    switch (type) {
    case Type::QED:
        return FKS::QED::SoftCollinearLimitFSR(args...);
    case Type::QCD:
        return FKS::QCD::SoftCollinearLimitFSR(args...);
    }
    assert(0.0);
    return 0.0;
}

template <Type type, typename... Args>
double softCollinearLimitISR(Args... args) {
    switch (type) {
    case Type::QED:
        return FKS::QED::SoftCollinearLimitISR(args...);
    case Type::QCD:
        return FKS::QCD::SoftCollinearLimitISR(args...);
    }
    assert(0.0);
    return 0.0;
}

template <Type type>
double SxG_ISR(const FKS::Real_t &real, const FKS::Phasespaces &PS,
               double bornme, const Util::Matrix2 &collcorr, double s_b,
               const FKS::Region &region, double xi, double y, double phi,
               const UserProcess::Data *userdata) {
    static const double delta = 1e-6;
    int proc = userdata->ProcessID;

    const int *pdgs_b = userdata->Process[proc].Born.Flavours.data();
    const int *pdgs_r = real.Flavours.data();

    const FKS::Param *param = userdata->Params;
    double alpha = param->alpha;
    if(type == Type::QCD) {
        alpha = param->alphaS();
    }
    if (1.0 - y < delta) {
        if (xi < 1e-6) {
            return softCollinearLimitISR<type>(pdgs_r, pdgs_b, 1, s_b, alpha,
                                               bornme);
        }
        return collinearLimitISR<type>(pdgs_r, pdgs_b, xi, 1, s_b, alpha,
                                       bornme);
    }
    if (1.0 + y < delta) {
        if (xi < 1e-6) {
            return softCollinearLimitISR<type>(pdgs_r, pdgs_b, -1, s_b, alpha,
                                               bornme);
        }
        return collinearLimitISR<type>(pdgs_r, pdgs_b, xi, -1, s_b, alpha,
                                       bornme);
    }
    bool mod = userdata->modBbar;
    bool nointer = userdata->noInterference;

    int pdg_rad = pdgs_r[region.I];
    if (xi < 1e-8) {
        // ps = ps_real in soft limit
        double S =
            FKS::SFunction(PS.Soft, region, real, userdata->useResonanesInS);
        double soft_isr = 0.0;
        double soft_fsr = 0.0;
        const Phasespace::Phasespace &ps = PS.Born;
        int N = PS.Real.N; // N is the number of final states because of a
                           // stupid decision
        if (mod || nointer) {
            int pdgs_isr[N + 1];
            int pdgs_fsr[N + 1];
            for (int k = 0; k < 2; k++) {
                pdgs_isr[k] = pdgs_r[k];
                pdgs_fsr[k] = 22;
            }
            for (int k = 2; k < N + 1; k++) {
                pdgs_isr[k] = 22;
                pdgs_fsr[k] = pdgs_r[k];
            }
            soft_isr =
                softlimit<type>(N + 1, ps.Momenta.data(), pdgs_isr, pdg_rad,
                                s_b, 0, alpha, collcorr, y, phi);
            soft_fsr =
                softlimit<type>(N + 1, ps.Momenta.data(), pdgs_fsr, pdg_rad,
                                s_b, 0, alpha, collcorr, y, phi);
        }
        double soft_int = 0.0;
        if (!nointer) {
            double soft =
                softlimit<type>(N + 1, ps.Momenta.data(), pdgs_r, pdg_rad, s_b,
                                0, alpha, collcorr, y, phi);
            soft_int = soft - soft_isr - soft_fsr;
        }

        if (mod) {
            return soft_isr + S * soft_int;
        }
        return S * (soft_isr + soft_fsr + soft_int);
    }

    double real_isr = 0.0;
    double real_fsr = 0.0;
    if (mod || nointer) {
        using Diag = UserProcess::IMatrixElement::Diagrams;
        real_isr = userdata->MatrixElement->Real(real.ID, PS.Real, param,
                                                 Diag::ONLYISR);
        real_fsr = userdata->MatrixElement->Real(real.ID, PS.Real, param,
                                                 Diag::ONLYFSR);
    }
    double real_int = 0.0;
    if (!nointer) {
        double realme = userdata->MatrixElement->Real(real.ID, PS.Real, param);
        real_int = realme - real_isr - real_fsr;
    }

    double pre = (1.0 - y * y) * xi * xi;
    double S = FKS::SFunction(PS.Real, region, real, userdata->useResonanesInS);

    if (mod) {
        return pre * (real_isr + S * real_int);
    }
    return pre * S * (real_isr + real_fsr + real_int);
}

template <Type type>
double SxG_FSR(const FKS::Real_t &real, const FKS::Phasespaces &PS,
               double bornme, const Util::Matrix2 &collcorr, double s_r,
               const FKS::Region &region, double xi, double y, double phi,
               const UserProcess::Data *userdata) {
    int proc = userdata->ProcessID;

    const int *pdgs_b = userdata->Process[proc].Born.Flavours.data();
    const int *pdgs_r = real.Flavours.data();
    
    const FKS::Param *param = userdata->Params;
    double alpha = param->alpha;
    if(type == Type::QCD) {
        alpha = param->alphaS();
    }

    int pdg_b = pdgs_b[region.J];
    int pdg_r = pdgs_r[region.J];
    if (y > 1.0 - 1e-5) {
        if (xi < 1e-6) {
            return softCollinearLimitFSR<type>(pdg_r, pdg_b, s_r, alpha,
                                               bornme);
        }
        return collinearLimitFSR<type>(pdg_r, pdg_b, xi, s_r, alpha, bornme);
    }

    bool mod = userdata->modBbar;
    bool nointer = userdata->noInterference;
    int pdg_rad = pdgs_r[region.I];
    if (xi < 1e-6) {
        double soft_isr = 0.0;
        double soft_fsr = 0.0;
        int N = PS.Real.N; // N is the number of final states because of a
                           // stupid decision
        const Phasespace::Phasespace &ps = PS.Born;
        if (mod || nointer) {
            int pdgs_isr[N + 1];
            int pdgs_fsr[N + 1];
            for (int k = 0; k < 2; k++) {
                pdgs_isr[k] = pdgs_r[k];
                pdgs_fsr[k] = 22;
            }
            for (int k = 2; k < N + 1; k++) {
                pdgs_isr[k] = 22;
                pdgs_fsr[k] = pdgs_r[k];
            }

            soft_isr =
                softlimit<type>(N + 1, ps.Momenta.data(), pdgs_isr, pdg_rad,
                                s_r, region.J, alpha, collcorr, y, phi);
            soft_fsr =
                softlimit<type>(N + 1, ps.Momenta.data(), pdgs_fsr, pdg_rad,
                                s_r, region.J, alpha, collcorr, y, phi);
        }
        double soft_int = 0.0;
        if (!nointer) {
            double soft =
                softlimit<type>(N + 1, ps.Momenta.data(), pdgs_r, pdg_rad, s_r,
                                region.J, alpha, collcorr, y, phi);
            soft_int = soft - soft_isr - soft_fsr;
        }
        // ps = ps_real in soft limit
        double S =
            FKS::SFunction(PS.Soft, region, real, userdata->useResonanesInS);

        if (mod) {
            auto real_wo_isr = real;
            real_wo_isr.Regions.erase(
                std::remove_if(real_wo_isr.Regions.begin(),
                               real_wo_isr.Regions.end(),
                               [](FKS::Region r) { return r.J < 2; }),
                real_wo_isr.Regions.end());
            double S_fsr = FKS::SFunction(PS.Soft, region, real_wo_isr,
                                          userdata->useResonanesInS);
            return S_fsr * soft_fsr + S * soft_int;
        }
        return S * (soft_fsr + soft_isr + soft_int);
    }

    double real_isr = 0.0;
    double real_fsr = 0.0;
    if (mod || nointer) {
        using Diag = UserProcess::IMatrixElement::Diagrams;
        real_fsr = userdata->MatrixElement->Real(real.ID, PS.Real, param,
                                                 Diag::ONLYFSR);
        real_isr = userdata->MatrixElement->Real(real.ID, PS.Real, param,
                                                 Diag::ONLYISR);
    }
    double real_int = 0.0;
    if (!nointer) {
        double realme =
            userdata->MatrixElement->Real(real.ID, PS.Real, param);
        real_int = realme - real_fsr - real_isr;
    }
    double pre = (1.0 - y) * xi * xi;
    double S = FKS::SFunction(PS.Real, region, real, userdata->useResonanesInS);
    if (mod) {
        auto real_wo_isr = real;
        real_wo_isr.Regions.erase(
            std::remove_if(real_wo_isr.Regions.begin(),
                           real_wo_isr.Regions.end(),
                           [](FKS::Region r) { return r.J < 2; }),
            real_wo_isr.Regions.end());
        double S_fsr = FKS::SFunction(PS.Real, region, real_wo_isr,
                                      userdata->useResonanesInS);
        return pre * (S_fsr * real_fsr + S * real_int);
    }
    return S * pre * (real_isr + real_fsr + real_int);
}

} // end namespace

namespace FKS {
namespace QED {

FKS::MatrixElement Limits(const FKS::Real_t &real, const FKS::Phasespaces &PS,
                          double bornme, const Util::Matrix2 &collcorr,
                          const FKS::Region &region, double x,
                          const FKS::Xi &Xi, double y, double phi,
                          const UserProcess::Data *userdata) {
    FKS::MatrixElement lim;
    if (region.J >= 2) {
        const auto &ps_real = PS.Real;
        double s_r = ps_real.X1 * ps_real.X2 * ps_real.S;
        double xi = x * Xi.Max;
        lim.Soft = SxG_FSR<Type::QED>(real, PS, bornme, collcorr, s_r, region,
                                      0.0, y, phi, userdata);
        lim.Collinear1 = SxG_FSR<Type::QED>(real, PS, bornme, collcorr, s_r,
                                            region, xi, 1.0, phi, userdata);
        lim.SoftCollinear1 = SxG_FSR<Type::QED>(
            real, PS, bornme, collcorr, s_r, region, 0.0, 1.0, phi, userdata);
    } else {
        const auto &ps= PS.Born;
        double s_b = ps.X1 * ps.X2 * ps.S;
        double xi1 = x * Xi.Max_Coll1;
        double xi2 = x * Xi.Max_Coll2;
        lim.Soft = SxG_ISR<Type::QED>(real, PS, bornme, collcorr, s_b, region,
                                      0.0, y, phi, userdata);
        lim.Collinear1 =
            SxG_ISR<Type::QED>(real, PS, bornme, collcorr, s_b, region,
                               xi1, 1.0, phi, userdata);
        lim.Collinear2 = SxG_ISR<Type::QED>(real, PS, bornme, collcorr, s_b,
                                            region, xi2, -1.0, phi, userdata);
        lim.SoftCollinear1 = SxG_ISR<Type::QED>(
            real, PS, bornme, collcorr, s_b, region, 0.0, 1.0, phi, userdata);
        lim.SoftCollinear2 = SxG_ISR<Type::QED>(
            real, PS, bornme, collcorr, s_b, region, 0.0, -1.0, phi, userdata);
    }
    return lim;
}

double SxG(const FKS::Real_t &real,
           const FKS::Phasespaces &PS, double bornme,
           const Util::Matrix2 &collcorr, const FKS::Region &region, double xi,
           double y, double phi, 
           const UserProcess::Data *userdata) {
    const auto &ps_real = PS.Real;
    double s_r = ps_real.X1 * ps_real.X2 * ps_real.S;
    if (region.J >= 2) {
        return SxG_FSR<Type::QED>(real, PS, bornme, collcorr, s_r, region, xi,
                                  y, phi, userdata);
    }
    const auto &ps = PS.Born;
    double s_b = ps.X1 * ps.X2 * ps.S;
    return SxG_ISR<Type::QED>(real, PS, bornme, collcorr, s_b, region, xi, y,
                              phi, userdata);
}

} // end namespace QED

namespace QCD {

FKS::MatrixElement Limits(const FKS::Real_t &real,
                          const FKS::Phasespaces &PS,
                          double bornme, const Util::Matrix2 &collcorr,
                          const FKS::Region &region, double x,
                          const FKS::Xi &Xi, double y, double phi,
                          const UserProcess::Data *userdata) {
    FKS::MatrixElement lim;
    if (region.J >= 2) {
        double xi = x * Xi.Max;
        const auto &ps_real = PS.Real;
        double s_r = ps_real.X1 * ps_real.X2 * ps_real.S;
        lim.Soft = SxG_FSR<Type::QCD>(real, PS, bornme, collcorr, s_r, region,
                                      0.0, y, phi, userdata);
        lim.Collinear1 = SxG_FSR<Type::QCD>(real, PS, bornme, collcorr, s_r,
                                            region, xi, 1.0, phi, userdata);
        lim.SoftCollinear1 = SxG_FSR<Type::QCD>(
            real, PS, bornme, collcorr, s_r, region, 0.0, 1.0, phi, userdata);
    } else {
        const auto &ps = PS.Born;
        double s_b = ps.X1 * ps.X2 * ps.S;
        double xi1 = x * Xi.Max_Coll1;
        double xi2 = x * Xi.Max_Coll2;
        lim.Soft = SxG_ISR<Type::QCD>(real, PS, bornme, collcorr, s_b, region,
                                      0.0, y, phi, userdata);
        lim.Collinear1 = SxG_ISR<Type::QCD>(real, PS, bornme, collcorr, s_b,
                                            region, xi1, 1.0, phi, userdata);
        lim.Collinear2 = SxG_ISR<Type::QCD>(real, PS, bornme, collcorr, s_b,
                                            region, xi2, -1.0, phi, userdata);
        lim.SoftCollinear1 = SxG_ISR<Type::QCD>(
            real, PS, bornme, collcorr, s_b, region, 0.0, 1.0, phi, userdata);
        lim.SoftCollinear2 = SxG_ISR<Type::QCD>(
            real, PS, bornme, collcorr, s_b, region, 0.0, -1.0, phi, userdata);
    }
    return lim;
}

double SxG(const FKS::Real_t &real, const FKS::Phasespaces &PS, double bornme,
           const Util::Matrix2 &collcorr, const FKS::Region &region, double xi,
           double y, double phi, const UserProcess::Data *userdata) {
    const auto &ps_real = PS.Real;
    double s_r = ps_real.X1 * ps_real.X2 * ps_real.S;
    if (region.J >= 2) {
        return SxG_FSR<Type::QCD>(real, PS, bornme, collcorr, s_r, region,
                                  xi, y, phi, userdata);
    }
    const auto &ps = PS.Born;
    double s_b = ps.X1 * ps.X2 * ps.S;
    return SxG_ISR<Type::QCD>(real, PS, bornme, collcorr, s_b, region, xi,
                              y, phi, userdata);
}

} // end namespace QCD

} // end namespace FKS

