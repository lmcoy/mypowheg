#ifndef FKS_LUMINOSITY_H_
#define FKS_LUMINOSITY_H_

#include "libconfig.h"

#include "pdf/pdfcache.h"

namespace FKS {

struct Xi;

/**
PartonLuminosity combines the parton luminosites for the real phase space and
its limits.
*/
struct LIB_LOCAL PartonLuminosity {
    PartonLuminosity() : Born(0.0), Collinear1(0.0), Collinear2(0.0), Real(0.0) {}
    PartonLuminosity(double born, double collinear1, double collinear2, double real)
        : Born(born), Collinear1(collinear1), Collinear2(collinear2), Real(real) {
    }

    double Born;      ///< luminosity in the soft and soft collinear limit
    double Collinear1; ///< luminosity in the collinear limit y == 1
    double Collinear2; ///< luminosity in the collinear limit y == -1 (only for ISR)
    double Real;      ///< luminosity at real phase space point

    PartonLuminosity & Add(const PartonLuminosity & l) {
        Born += l.Born;
        Collinear1 += l.Collinear1;
        Collinear2 += l.Collinear2;
        Real += l.Real;
        return *this;
    }
};

LIB_LOCAL PartonLuminosity Luminosity(PDF::Cache *pdf, int flavour1,
                                      int flavour2, double x, double x1_b,
                                      double x2_b, double x1_real,
                                      double x2_real, const Xi &xi, double muF);

struct LIB_LOCAL LumRemnants {
    LumRemnants() : Remnant(0.0), Born(0.0) {}
    LumRemnants(double remn, double born) : Remnant(remn), Born(born) {}
    double Remnant;
    double Born;
};

LIB_LOCAL LumRemnants LuminosityRemnants(PDF::Cache *pdf, int flavour1,
                                         int flavour2, int colldir, double x1b,
                                         double x2b, double x,
                                         double xi_max_col, double muF);
} // end namespace FKS
#endif
