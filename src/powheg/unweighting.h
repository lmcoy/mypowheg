#ifndef UNWEIGHTING_H_ZWLWID6X
#define UNWEIGHTING_H_ZWLWID6X

#include <cstdlib>

#include "libconfig.h"

namespace UserProcess {
class Data;
}

namespace Phasespace {
class Phasespace;
}

namespace Powheg {

class BornConfig {
  public:
    enum class Status {
        Accepted,
        Rejected,
        NegativeBtilde,
        RejectedWithGuessedVirtual,
        UnderestimatedVirtual
    };
    size_t PDFindex = 0;
    size_t ProcessIndex = 0;
    double btilde = 0.0;
    int n = 0;
    Status status;
};

BornConfig unweighting(const Phasespace::Phasespace &ps, double x1, double x2,
                       double x3, double wgt, UserProcess::Data *userdata);

} /* Powheg */

#endif /* end of include guard: UNWEIGHTING_H_ZWLWID6X */
