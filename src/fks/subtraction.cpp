#include "fks/subtraction.h"

#include "math/math.h"
#include "fks/ximax.h"
#include "fks/luminosity.h"
#include "fks/limits.h"

namespace FKS {

SubtractionTerms RealISR2(const PartonLuminosity &lum, const MatrixElement &me,
                          int colldir, double x, const Xi &xi, double y) {
    double xi_max_coll = 0.0;
    double xi_real = xi.Max * x;
    double xi_coll = 0.0;
    // constant factor from jacobian det of n+1 particle
    double prefact = 1.0 / (64.0 * Math::Pi_Cube);

    double lum_real = lum.Real;
    double lum_coll = 0.0;
    double colllimit = 0.0;
    double softcolllimit = 0.0;
    double y_pre = 0.0;
    assert(colldir == 1 || colldir == -1);
    if (colldir == 1) {
        lum_coll = lum.Collinear1;
        xi_coll = xi.Max_Coll1 * x;
        xi_max_coll = xi.Max_Coll1;
        colllimit = me.Collinear1;
        softcolllimit = me.SoftCollinear1;
        y_pre = 1.0 / (1.0 - y);
    } else {
        lum_coll = lum.Collinear2;
        xi_coll = xi.Max_Coll2 * x;
        xi_max_coll = xi.Max_Coll2;
        colllimit = me.Collinear2;
        softcolllimit = me.SoftCollinear2;
        y_pre = 1.0 / (1.0 + y);
    }

    double lum_born = lum.Born;

    double softlimit = me.Soft;

    double fks_real = lum_real * 1.0 / (1.0 - xi_real) * me.Real;
    double fks_coll = lum_coll * 1.0 / (1.0 - xi_coll) * colllimit;
    double fks_soft = lum_born * softlimit;
    double fks_softcoll = lum_born * softcolllimit;

    double pre = 0.25 * prefact * y_pre / x;
    double pre_end = 0.25 * y_pre * prefact;
    SubtractionTerms terms;
    terms.Real = pre * fks_real;
    terms.Soft = (-pre + pre_end * log(xi.Max)) * fks_soft;
    terms.Collinear = -pre * fks_coll;
    terms.SoftCollinear = (pre - pre_end * log(xi_max_coll)) * fks_softcoll;

    return terms;
}

SubtractionTerms RealFSR(double lumi, const MatrixElement &me, double x,
                         double ximax, double y, double len_kj_born, double s) {
    double sqrts = sqrt(s);
    double xi = x * ximax;
    double M2 = s - 2.0 * len_kj_born * sqrts;
    // flux factor cancels
    double common_J = 1.0 / (64.0 * Math::Pi_Cube) / len_kj_born;

    double k0 = sqrts * xi * 0.5;
    double len_kn_real =
        (s - M2 - 2 * sqrts * k0) / (2.0 * (sqrts - k0 * (1.0 - y)));
    double len_kn_soft = (s - M2) / (2.0 * sqrts);
    double len_kn_collinear = (s - M2 - 2.0 * sqrts * k0) / (2.0 * sqrts);

    double J_real = common_J * len_kn_real / (2.0 - xi * (1.0 - y));
    double J_soft = common_J * len_kn_soft / 2.0;
    double J_collinear = common_J * len_kn_collinear / 2.0;

    double G_real = J_real * me.Real;
    double G_soft = J_soft * me.Soft;
    double G_collinear = J_collinear * me.Collinear1;
    double G_softcoll = J_soft * me.SoftCollinear1;

    double pre = lumi * 1.0 / (1.0 - y) / xi;
    double pre_end = lumi * log(ximax) / (1.0 - y);
    SubtractionTerms terms;
    terms.Real = pre * G_real;
    terms.Soft = (-pre + pre_end) * G_soft;
    terms.Collinear = -pre * G_collinear;
    terms.SoftCollinear = (pre - pre_end) * G_softcoll;
    return terms;
}

} // namespace FKS
