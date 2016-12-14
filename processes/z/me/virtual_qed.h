#ifndef VIRTUAL_QED_H
#define VIRTUAL_QED_H

#include "math/fourmomentum.h"
#include "parameters_sm.h"
#include "process/matrixelement.h"
#include "../matrixelement.h"

enum class CorrectionType {
    None,
    QCD,
    EW
};

void Vfin_qxq_lxl(double *born, double *virt, const Math::FourMomentum *momenta,
                  const int *perm, CorrectionType type, int quark_flavour, double alpha, double alphas, double mu2,
                  double counterterm);

void InitVfin_qxq_lxl(const Parameters_sm &param, AlphaScheme scheme,
                      double inveps, double *counterterm);

#endif 
 
