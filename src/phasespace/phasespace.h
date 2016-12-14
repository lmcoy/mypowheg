#ifndef PHASESPACE_PHASESPACE_H_
#define PHASESPACE_PHASESPACE_H_

#include <array>

#include "math/fourmomentum.h"
#include "phasespace/internal.h"

namespace Phasespace {

const int MAXMOM = 6;

class Phasespace {
  public:
    void SetToCMSFromLab(Phasespace const *const ps_lab);
    void SetToLabFromCMS(Phasespace const *const ps_cms);

    double X1; ///< momentum fraction of the 1st initial state particles
    double X2; ///< momentum fraction of the 2nd initial state particles
    double S;  ///< CMS energy of the protons
    std::array<double, MAXMOM-2> Masses; ///< masses of the final state particles
    std::array<Math::FourMomentum, MAXMOM> Momenta;
    double Jacobian; ///< jacobian determinante
    int N;
};

// ----------------------------------------------------------------------------
// Implementation of phasespace function.
//
// Note: They have to be implemented in the header because the compiler has to
// now the functions. If they would be in a cpp file, the compiler wouldn't know
// them if the program is compiled to a lib.
// ----------------------------------------------------------------------------

/**
SetToCMSFromLab sets this phasespace to ps_lab in the partonic CMS frame.

@param ps_lab phase space in lab frame
*/
inline void Phasespace::SetToCMSFromLab(Phasespace const *const ps_lab) {
    const double x1 = ps_lab->X1;
    const double x2 = ps_lab->X2;
    const double g = gamma(x1, x2);
    const double bg = -betagamma(x1, x2);

    N = ps_lab->N;
    X1 = x1;
    X2 = x2;
    S = ps_lab->S;
    Masses = ps_lab->Masses;
    Momenta = ps_lab->Momenta;
    for (int i = 0; i < N+2; i++) {
        boost_z(g, bg, &Momenta[i]);
    }
}

/**
SetToLabFromCMS sets this phasespace to ps_cms in the lab frame.

@param ps_cms phase space in partonic CMS.
*/
inline void Phasespace::SetToLabFromCMS(Phasespace const *const ps_cms) {
    const double x1 = ps_cms->X1;
    const double x2 = ps_cms->X2;
    const double g = gamma(x1, x2);
    const double bg = betagamma(x1, x2);

    N = ps_cms->N;
    X1 = x1;
    X2 = x2;
    S = ps_cms->S;
    Masses = ps_cms->Masses;
    Momenta = ps_cms->Momenta;
    for (int i = 0; i < N+2; i++) {
        boost_z(g, bg, &Momenta[i]);
    }
}

} // end namespace Phasespace

#endif
