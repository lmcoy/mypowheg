#ifndef PHASESPACEGENERATOR_H_AUI98SU1
#define PHASESPACEGENERATOR_H_AUI98SU1

#include <memory>

#include "phasespace/phasespace.h"

namespace UserProcess {
struct Data;
}

class PhasespaceGenerator {
  public:
    /**
     * @brief Generate Born phase space
     *
     * Gen has to implement an n particle phase space generator. It is used to
     * generate the underlying born phase space.
     *
     * Gen has to construct the phase space from a set of kinematic variables x.
     * All variables have to bin the range [0,1).
     *
     * @param n number of kinematic variables in x
     * @param x kinematic varibales in [0,1)
     * @param userdata information about the process
     */
    virtual Phasespace::Phasespace Gen(int n, double *x,
                                       UserProcess::Data *userdata) = 0;

    /**
     * @brief number of kinematic variables
     *
     * Dim returns the number of kinematic variables that are needed to
     * construct the phase space.
     */
    virtual int Dim() const = 0;
};
typedef std::shared_ptr<PhasespaceGenerator> PhasespaceGeneratorPtr;

#endif
