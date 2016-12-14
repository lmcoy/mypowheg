#include "fks/luminosity.h"
#include "fks/ximax.h"

namespace FKS {

/**
Luminosity returns the parton luminosities for the real matrix elements.

It returns all end points for the parton luminosity f(x1)*f(x2) needed for the
plus distributions in:

(1/xi)_+ (1/(1-+y)_+ f(x1)f(x2) R

@param flavour1 flavour of the first initial state particle (in lhapdf
numbering)
@param flavour2 flavour of the second initial state particle (in lhapdf
numbering)
@param x integration variable in [0,1) for xi integration
@param x1_b momentum fraction of the underlying born phase space
@param x2_b momentum fraction of the underlying born phase space
@param x1_real momentum fraction of the n+1 particle phase space
@param x2_real momentum fraction of the n+1 particle phase space
@param colldir defines if the collinear limit is in y = 1 or y = -1 direction
@param xi_max_real maximal value for xi
@param xi_max_coll maximal value for xi in the collinear limit
*/
PartonLuminosity Luminosity(PDF::Cache *pdf, int flavour1, int flavour2,
                            double x, double x1_b, double x2_b, double x1_real,
                            double x2_real, const Xi &xi, double muF) {
    double xi_coll1 = xi.Max_Coll1 * x;
    double xi_coll2 = xi.Max_Coll2 * x;

    double x_coll1 = x1_b / (1.0 - xi_coll1);
    double x_coll2 = x2_b / (1.0 - xi_coll2);

    double lum_real =
        pdf->Get(x1_real, muF, flavour1) * pdf->Get(x2_real, muF, flavour2);
    double f1 = pdf->Get(x1_b, muF, flavour1);
    double f2 = pdf->Get(x2_b, muF, flavour2);
    
    double lum_coll1 = pdf->Get(x_coll1, muF, flavour1) * f2;
    double lum_coll2 = f1 * pdf->Get(x_coll2, muF, flavour2);

    double lum_born = f1 * f2;

    return PartonLuminosity(lum_born, lum_coll1, lum_coll2, lum_real);
}

LumRemnants LuminosityRemnants(PDF::Cache *pdf, int flavour1, int flavour2,
                               int colldir, double x1b, double x2b, double x,
                               double xi_max_col, double muF) {
    double f1b = pdf->Get(x1b, muF, flavour1);
    double f2b = pdf->Get(x2b, muF, flavour2);

    double xi = x;
    if (xi > xi_max_col) {
        return LumRemnants(0.0, f1b * f2b);
    }

    double x_rem = 1.0 / (1 - xi);
    if (colldir == 1) {
        x_rem *= x1b;
        double fc = pdf->Get(x_rem, muF, flavour1);
        return LumRemnants(fc * f2b, f1b * f2b);
    }
    x_rem *= x2b;
    double fc = pdf->Get(x_rem, muF, flavour2);
    return LumRemnants(fc * f1b, f1b * f2b);
}

} // end namespace FKS
