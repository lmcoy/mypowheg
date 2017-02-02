#include <algorithm>

#include "fks_g.h"

#include "fks/limits.h"
#include "fks/param.h"
#include "fks/phasespaces.h"
#include "fks/process.h"
#include "fks/sfunctions.h"
#include "fks/ximax.h"
#include "phasespace/phasespace.h"
#include "process/data.h"
#include "process/matrixelement.h"
#include "process/matrixelement.h"
#include "util/matrix.h"

namespace {

enum class Type { QED, QCD };

template <Type type, typename... Args> double softlimit(Args... args) {
    switch (type) {
    case Type::QED:
        return FKS::QED::SoftLimit(args...);
    case Type::QCD:
        return FKS::QCD::SoftLimit(args...);
    }
    assert(0.0);
    return 0.0;
}

template <Type type, typename... Args> double collinearLimitFSR(Args... args) {
    switch (type) {
    case Type::QED:
        return FKS::QED::CollinearLimitFSR(args...);
    case Type::QCD:
        return FKS::QCD::CollinearLimitFSR(args...);
    }
    assert(0.0);
    return 0.0;
}

template <Type type, typename... Args> double collinearLimitISR(Args... args) {
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
               double bornme, const Util::Matrix2 &collcorr,
               const UserProcess::SpinCorrelated &spincorr, double s_b,
               const FKS::Region &region, double xi, double y, double phi,
               const UserProcess::Data *userdata) {
    static constexpr double delta_y = 1e-5;
    static constexpr double delta_xi = 1e-5;
    static constexpr double delta_xi_2 = 1e-4;
    int proc = userdata->ProcessID;

    const int *pdgs_b = userdata->Process[proc].Born.Flavours.data();
    int Nborn = userdata->Process[proc].Born.Flavours.size();
    const int *pdgs_r = real.Flavours.data();

    const FKS::Param *param = userdata->Params;
    double alpha = param->alpha;
    if (type == Type::QCD) {
        alpha = param->alphaS();
    }
    if (1.0 - y < delta_y) {
        if (region.J == -2) {
            // the S function is 0.0 in this limit
            return 0.0;
        }
        if (xi < delta_xi_2) {
            return softCollinearLimitISR<type>(Nborn, pdgs_r, pdgs_b, 1, phi,
                                               s_b, alpha, bornme, spincorr);
        }
        return collinearLimitISR<type>(Nborn, pdgs_r, pdgs_b, xi, 1, phi, s_b,
                                       alpha, bornme, spincorr);
    }
    if (1.0 + y < delta_y) {
        if (region.J == -1) {
            // the S function is 0.0 in this limit
            return 0.0;
        }
        if (xi < delta_xi_2) {
            return softCollinearLimitISR<type>(Nborn, pdgs_r, pdgs_b, -1, phi,
                                               s_b, alpha, bornme, spincorr);
        }
        return collinearLimitISR<type>(Nborn, pdgs_r, pdgs_b, xi, -1, phi, s_b,
                                       alpha, bornme, spincorr);
    }
    bool mod = userdata->modBbar;
    bool nointer = userdata->noInterference;

    int pdg_rad = pdgs_r[region.I];
    if (xi < delta_xi) {
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
    double S = FKS::SFunction(PS.Real, region, real, userdata->useResonanesInS);
    double real_int = 0.0;
    if (!nointer && S > 1e-6) {
        double realme = userdata->MatrixElement->Real(real.ID, PS.Real, param);
        real_int = realme - real_isr - real_fsr;
    }

    double pre = (1.0 - y * y) * xi * xi;

    if (mod) {
        return pre * (real_isr + S * real_int);
    }
    return pre * S * (real_isr + real_fsr + real_int);
}

template <Type type>
double SxG_FSR(const FKS::Real_t &real, const FKS::Phasespaces &PS,
               double bornme, const Util::Matrix2 &collcorr,
               const UserProcess::SpinCorrelated &spincorr, double s_r,
               const FKS::Region &region, double xi, double y, double phi,
               const UserProcess::Data *userdata) {
    constexpr double delta_y = 1e-6;
    constexpr double delta_xi = 1e-8;
    constexpr double delta_xi_2 = 1e-6;
    int proc = userdata->ProcessID;

    const int *pdgs_b = userdata->Process[proc].Born.Flavours.data();
    const int *pdgs_r = real.Flavours.data();

    const FKS::Param *param = userdata->Params;
    double alpha = param->alpha;
    if (type == Type::QCD) {
        alpha = param->alphaS();
    }

    int mother_index = region.J;
    if (region.J > region.I) {
        mother_index = region.I;
    }

    int pdg_b = pdgs_b[mother_index];
    int pdg_r = pdgs_r[region.J];
    int pdg_rad = pdgs_r[region.I];
    if (y > 1.0 - delta_y) {
        if (xi < delta_xi_2) {
            return softCollinearLimitFSR<type>(pdg_r, pdg_b, PS.Born,
                                               mother_index, phi, alpha, bornme,
                                               spincorr);
        }
        double S = 1.0;
        if (type == Type::QCD) {
            if (pdg_r == 21 || pdg_r == 0) {
                if (pdg_rad == 21 || pdg_rad == 0) {
                    double Ej = PS.Real.Momenta[region.J].E();
                    double Ei = PS.Real.Momenta[region.I].E();
                    // S function is not 1 for collinear limit of g -> gg
                    // splitting. We have to add the h factor.
                    S = 1.0 / (1.0 + Ei / Ej);
                }
            }
        }
        return S * collinearLimitFSR<type>(pdg_r, pdg_b, PS.Born, mother_index,
                                           xi, phi, alpha, bornme, spincorr);
    }

    bool mod = userdata->modBbar;
    bool nointer = userdata->noInterference;
    if (xi < delta_xi) {
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
                                s_r, mother_index, alpha, collcorr, y, phi);
            soft_fsr =
                softlimit<type>(N + 1, ps.Momenta.data(), pdgs_fsr, pdg_rad,
                                s_r, mother_index, alpha, collcorr, y, phi);
        }
        double soft_int = 0.0;
        if (!nointer) {
            double soft =
                softlimit<type>(N + 1, ps.Momenta.data(), pdgs_r, pdg_rad, s_r,
                                mother_index, alpha, collcorr, y, phi);
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
    double S = FKS::SFunction(PS.Real, region, real, userdata->useResonanesInS);
    if (!nointer && S > 1e-6) {
        double realme = userdata->MatrixElement->Real(real.ID, PS.Real, param);
        real_int = realme - real_fsr - real_isr;
    }
    double pre = (1.0 - y) * xi * xi;

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
                          const UserProcess::SpinCorrelated &spincorr,
                          const FKS::Region &region, double x,
                          const FKS::Xi &Xi, double y, double phi,
                          const UserProcess::Data *userdata) {
    FKS::MatrixElement lim;
    if (region.J >= 2) {
        const auto &ps_real = PS.Real;
        double s_r = ps_real.X1 * ps_real.X2 * ps_real.S;
        double xi = x * Xi.Max;
        lim.Soft = SxG_FSR<Type::QED>(real, PS, bornme, collcorr, spincorr, s_r,
                                      region, 0.0, y, phi, userdata);
        lim.Collinear1 =
            SxG_FSR<Type::QED>(real, PS, bornme, collcorr, spincorr, s_r,
                               region, xi, 1.0, phi, userdata);
        lim.SoftCollinear1 =
            SxG_FSR<Type::QED>(real, PS, bornme, collcorr, spincorr, s_r,
                               region, 0.0, 1.0, phi, userdata);
    } else {
        const auto &ps = PS.Born;
        double s_b = ps.X1 * ps.X2 * ps.S;
        double xi1 = x * Xi.Max_Coll1;
        double xi2 = x * Xi.Max_Coll2;
        lim.Soft = SxG_ISR<Type::QED>(real, PS, bornme, collcorr, spincorr, s_b,
                                      region, 0.0, y, phi, userdata);
        lim.Collinear1 =
            SxG_ISR<Type::QED>(real, PS, bornme, collcorr, spincorr, s_b,
                               region, xi1, 1.0, phi, userdata);
        lim.Collinear2 =
            SxG_ISR<Type::QED>(real, PS, bornme, collcorr, spincorr, s_b,
                               region, xi2, -1.0, phi, userdata);
        lim.SoftCollinear1 =
            SxG_ISR<Type::QED>(real, PS, bornme, collcorr, spincorr, s_b,
                               region, 0.0, 1.0, phi, userdata);
        lim.SoftCollinear2 =
            SxG_ISR<Type::QED>(real, PS, bornme, collcorr, spincorr, s_b,
                               region, 0.0, -1.0, phi, userdata);
    }
    return lim;
}

double SxG(const FKS::Real_t &real, const FKS::Phasespaces &PS, double bornme,
           const Util::Matrix2 &collcorr,
           const UserProcess::SpinCorrelated &spincorr,
           const FKS::Region &region, double xi, double y, double phi,
           const UserProcess::Data *userdata) {
    const auto &ps_real = PS.Real;
    double s_r = ps_real.X1 * ps_real.X2 * ps_real.S;
    if (region.J >= 2) {
        return SxG_FSR<Type::QED>(real, PS, bornme, collcorr, spincorr, s_r,
                                  region, xi, y, phi, userdata);
    }
    const auto &ps = PS.Born;
    double s_b = ps.X1 * ps.X2 * ps.S;
    return SxG_ISR<Type::QED>(real, PS, bornme, collcorr, spincorr, s_b, region,
                              xi, y, phi, userdata);
}

} // end namespace QED

namespace QCD {

FKS::MatrixElement Limits(const FKS::Real_t &real, const FKS::Phasespaces &PS,
                          double bornme, const Util::Matrix2 &collcorr,
                          const UserProcess::SpinCorrelated &spincorr,
                          const FKS::Region &region, double x,
                          const FKS::Xi &Xi, double y, double phi,
                          const UserProcess::Data *userdata) {
    FKS::MatrixElement lim;
    if (region.J >= 2) {
        double xi = x * Xi.Max;
        const auto &ps_real = PS.Real;
        double s_r = ps_real.X1 * ps_real.X2 * ps_real.S;
        lim.Soft = SxG_FSR<Type::QCD>(real, PS, bornme, collcorr, spincorr, s_r,
                                      region, 0.0, y, phi, userdata);
        lim.Collinear1 =
            SxG_FSR<Type::QCD>(real, PS, bornme, collcorr, spincorr, s_r,
                               region, xi, 1.0, phi, userdata);
        lim.SoftCollinear1 =
            SxG_FSR<Type::QCD>(real, PS, bornme, collcorr, spincorr, s_r,
                               region, 0.0, 1.0, phi, userdata);
    } else {
        const auto &ps = PS.Born;
        double s_b = ps.X1 * ps.X2 * ps.S;
        double xi1 = x * Xi.Max_Coll1;
        double xi2 = x * Xi.Max_Coll2;
        lim.Soft = SxG_ISR<Type::QCD>(real, PS, bornme, collcorr, spincorr, s_b,
                                      region, 0.0, y, phi, userdata);
        lim.Collinear1 =
            SxG_ISR<Type::QCD>(real, PS, bornme, collcorr, spincorr, s_b,
                               region, xi1, 1.0, phi, userdata);
        lim.Collinear2 =
            SxG_ISR<Type::QCD>(real, PS, bornme, collcorr, spincorr, s_b,
                               region, xi2, -1.0, phi, userdata);
        lim.SoftCollinear1 =
            SxG_ISR<Type::QCD>(real, PS, bornme, collcorr, spincorr, s_b,
                               region, 0.0, 1.0, phi, userdata);
        lim.SoftCollinear2 =
            SxG_ISR<Type::QCD>(real, PS, bornme, collcorr, spincorr, s_b,
                               region, 0.0, -1.0, phi, userdata);
    }
    return lim;
}

double SxG(const FKS::Real_t &real, const FKS::Phasespaces &PS, double bornme,
           const Util::Matrix2 &collcorr,
           const UserProcess::SpinCorrelated &spincorr,
           const FKS::Region &region, double xi, double y, double phi,
           const UserProcess::Data *userdata) {
    const auto &ps_real = PS.Real;
    double s_r = ps_real.X1 * ps_real.X2 * ps_real.S;
    if (region.J >= 2) {
        return SxG_FSR<Type::QCD>(real, PS, bornme, collcorr, spincorr, s_r,
                                  region, xi, y, phi, userdata);
    }
    const auto &ps = PS.Born;
    double s_b = ps.X1 * ps.X2 * ps.S;
    return SxG_ISR<Type::QCD>(real, PS, bornme, collcorr, spincorr, s_b, region,
                              xi, y, phi, userdata);
}

} // end namespace QCD

} // end namespace FKS
