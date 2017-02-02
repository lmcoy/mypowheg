#ifndef VIRTUAL_H_HAKSHLOA
#define VIRTUAL_H_HAKSHLOA

#include "math/fourmomentum.h"
#include "parameters_sm.h"

enum class Type { None, EW, QCD };

enum class Scheme { Alpha0, Gmu };

void Wj_InitDimReg(const Parameters_sm &param, double inveps);

void Wj_InitMassReg(const Parameters_sm &param, double lambda);

void Wj_MatrixElement(int proc, double *born, double *virt,
                      const Math::FourMomentum *momenta, const int *perm,
                      double alpha, double alphas, double mu, Type type,
                      Scheme scheme);

#endif
