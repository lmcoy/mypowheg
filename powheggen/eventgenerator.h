#ifndef EVENTGENERATOR_H_AD6DSK71
#define EVENTGENERATOR_H_AD6DSK71

#include <memory>

#include "powheg/generateevents.h"

#include "phasespacegenerator.h"
#include "integraltransformation.h"

#include "run.h"

struct VegasState;

using Func = decltype(&Powheg::GenerateEvents);

struct IntegrandParams {
    IntegralTransformationPtr trafo;
    PhasespaceGeneratorPtr psgen;
    Func func;
    UserProcess::Data *userdata;
};
typedef std::shared_ptr<IntegrandParams> IntegrandParamsPtr;

class EventGenerator {
  public:
    EventGenerator(IntegralTransformationPtr trafo,
                   PhasespaceGeneratorPtr psgen,
                   std::shared_ptr<UserProcess::Data> data);

    virtual ~EventGenerator();

    virtual int Init();
    virtual void Setup();
    virtual int GenerateEvents();

  protected:
    std::shared_ptr<UserProcess::Data> userdata;
    int rank;
    int num_procs;
    bool verbose = true;
    int NDIM;
    VegasState *vegasstate_ = 0;

  private:
    Run run;

    IntegrandParams GenIntegrandParams(Func func) {
        IntegrandParams p;
        p.trafo = trafo_;
        p.psgen = psgen_;
        p.userdata = userdata.get();
        p.func = func;
        return p;
    }

    void SearchNormForUpperBounding();
    void Integrate();

    IntegralTransformationPtr trafo_;
    PhasespaceGeneratorPtr psgen_;
};

#endif
