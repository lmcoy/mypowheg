#ifndef BORNPSGENERATOR_H_AE458UJL
#define BORNPSGENERATOR_H_AE458UJL

#include "powheggen/phasespacegenerator.h"

#include "process/data.h"
#include "phasespace/twoparticlegenerator.h"

class BornPSGenerator : public PhasespaceGenerator {
  public:
    virtual Phasespace::Phasespace Gen(int n, double *x,
                                       UserProcess::Data *userdata);
    virtual int Dim() const { return 4; }
};

Phasespace::Phasespace BornPSGenerator::Gen(int n, double *x,
                                            UserProcess::Data *userdata) {
    double SqrtS = userdata->SqrtS;
    double y_b = x[2];
    Phasespace::TwoParticleGenerator gen;
    double xt[] = { x[0], x[1], y_b, x[3] };
    double masses[2] = { 0.0 };
    Phasespace::Phasespace ps_out;
    gen(&ps_out, 4, xt, SqrtS * SqrtS, 2, masses);
    return ps_out;
}

#endif
