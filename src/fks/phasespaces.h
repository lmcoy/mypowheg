#ifndef FKS_FKSPHASESPACES_H_
#define FKS_FKSPHASESPACES_H_

#include "phasespace/phasespace.h"

#include "fks/process.h"

namespace UserProcess {
struct RecombinationParam;
}

namespace FKS {

struct Phasespaces {
    Phasespace::Phasespace Real;
    Phasespace::Phasespace Soft;
    Phasespace::Phasespace Collinear1;
    Phasespace::Phasespace Collinear2;
    Phasespace::Phasespace Born;

    void Generate(const Phasespace::Phasespace &ps, int ifks, int jfks,
                  double xi, double y, double phi);

    Phasespaces Recombined(FKS::Type_t type, const int *pdgs,
                           const UserProcess::RecombinationParam &recomb) const;
};

} // end namespace FKS

#endif
