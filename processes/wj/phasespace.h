#ifndef PHASESPACE_SJD72MS7_H_
#define PHASESPACE_SJD72MS7_H_

#include <array>

#include "phasespace/phasespace.h"
#include "powheggen/phasespacegenerator.h"

#include "process/data.h"

Phasespace::Phasespace GenPhasespace(double S, const std::array<double, 7> &v);

class BornPSGenerator : public PhasespaceGenerator {
  public:
    virtual Phasespace::Phasespace Gen(int n, double *x,
                                       UserProcess::Data *userdata);
    virtual int Dim() const { return 7; }
};

#endif
