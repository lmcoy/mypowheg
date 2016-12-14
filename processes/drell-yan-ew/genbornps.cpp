#include "genbornps.h"

#include "process/data.h"
#include "phasespace/phasespace.h"
#include "phasespace/twoparticlegenerator.h"

void GenBornPhasespace(Phasespace::Phasespace *ps_out, const int ndim,
                       const double x[], const UserProcess::Data *userdata) {
    double SqrtS = userdata->SqrtS;
    double y_b = x[2];
    Phasespace::TwoParticleGenerator gen;
    double xt[] = { x[0], x[1], y_b, x[6] };
    double masses[2] = { 0.0 };
    gen(ps_out, 4, xt, SqrtS * SqrtS, 2, masses);
}

