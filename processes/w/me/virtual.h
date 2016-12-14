#ifndef VIRTUAL_H_JA38FH2P
#define VIRTUAL_H_JA38FH2P

#include "math/fourmomentum.h"
#include "parameters_sm.h"

enum class Type {
    None,
    EW,
    QCD
};

enum class Scheme {
    Alpha0,
    Gmu
};

void W_InitMassReg(const Parameters_sm &param, double lambda);

void W_InitDimReg(const Parameters_sm &param, double inveps);

void W_MatrixElement(double *born, double *virt,
                     const Math::FourMomentum *momenta, const int *perm,
                     double alpha, double alphas, double mu, Type type, Scheme scheme);

#endif

