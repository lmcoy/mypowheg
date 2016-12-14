#include "fks/remnants.h"
#include "fks/splitting.h"
#include "fks/scales.h"
#include "fks/luminosity.h"
#include "math/math.h"

#include <cmath>

#include "physics/pdgcode.h"

using namespace FKS;

namespace {

double K_DIS_QED(int realp, int bornp, double xi) {
    char rf = '0';
    char bf = '0';
    if (Physics::PDG::IsChargedFermion(realp)) {
        rf = 'f';
    }
    if (Physics::PDG::IsChargedFermion(bornp)) {
        bf = 'f';
    }
    if (rf == 'f' && bf == 'f') {
        double Q = Physics::PDG::Charge(realp);
        double z = 1.0 -xi;
        double ret = (1.0 + z * z) / (1.0 - z) * (log((1.0 - z) / z) - 0.75) + 2.25 +
               1.25 * z;
        return Q * Q * ret;
    }
    assert(0 && "not implemented");
    return 0.0;
}

double K_DIS_QCD(int realp, int bornp, double xi) {
    assert(0 && "DIS for QCD not implemented");
    return 0.0;
}

enum class T {
    QCD, QED
};

template <T type>
double K_DIS(int realp, int bornp, double xi) {
    switch(type) {
        case T::QCD:
            return K_DIS_QCD(realp, bornp, xi);
        case T::QED:
            return K_DIS_QED(realp, bornp, xi);
    }
    return 0.0;
}

template <T type>
double PxXi(int r, int b, double xi) {
    switch(type) {
        case T::QCD:
            return FKS::QCD::splittingTimesXi(r, b, xi);
        case T::QED:
            return FKS::QED::splittingTimesXi(r, b, xi);
    }
    return 0.0;
}

template <T type>
double PxXiSoft(int r, int b) {
    switch(type) {
        case T::QCD:
            return FKS::QCD::splittingTimesXiSoft(r, b);
        case T::QED:
            return FKS::QED::splittingTimesXiSoft(r, b);
    }
    return 0.0;
}

template <T type>
double Peps(int r, int b, double xi) {
    switch(type) {
        case T::QCD:
            return FKS::QCD::splittingEps(r, b, xi);
        case T::QED:
            return FKS::QED::splittingEps(r, b, xi);
    }
    return 0.0;
}

template <T type>
double remnant(const Scales &scales, const Splitting &sp, double xi,
               double xi_max, double s_b, double alpha,
               const LumRemnants &lumir, double born, PDFRenorm pdfren) {

    int realp = sp.RealPDG;
    int bornp = sp.BornPDG;
    double theta = (xi < xi_max) ? 1.0 : 0.0;
    double splitting = PxXi<type>(realp, bornp, xi);
    double splitting_soft = PxXiSoft<type>(realp, bornp);
    double plus = 1.0 / xi * (log(s_b / scales.muF / scales.muF / (1.0 - xi)) +
                              2.0 * log(xi));
    double plus_soft =
        1.0 / xi * (log(s_b / scales.muF / scales.muF) + 2.0 * log(xi));

    double ME_rem = lumir.Remnant * born;
    double ME_born = lumir.Born * born;
    double ret = (plus * splitting - Peps<type>(realp, bornp, xi)) /
                     (1.0 - xi) * ME_rem * theta -
                 plus_soft * splitting_soft * ME_born;
    if (pdfren == PDFRenorm::DIS) { // DIS pdf renormalization
        double K = K_DIS<type>(realp, bornp, xi);
        ret -= K * (theta * ME_rem / (1.0 - xi) - ME_born);
    }
    ret *= 1.0 / (2.0 * s_b);
    return ret * alpha / (2.0 * Math::Pi);
}

}

double FKS::QED::Remnant(const Scales &scales, const Splitting &sp, double xi,
                         double xi_max, double s_b, double alpha,
                         const LumRemnants &lumir, double born,
                         PDFRenorm pdfrenorm) {
    return remnant<T::QED>(scales, sp, xi, xi_max, s_b, alpha, lumir, born,
                           pdfrenorm);
}

double FKS::QCD::Remnant(const Scales &scales, const Splitting &sp, double xi,
                         double xi_max, double s_b, double alpha_s,
                         const LumRemnants &lumir, double born,
                         PDFRenorm pdfrenorm) {
    return remnant<T::QCD>(scales, sp, xi, xi_max, s_b, alpha_s, lumir, born,
                           pdfrenorm);
}
