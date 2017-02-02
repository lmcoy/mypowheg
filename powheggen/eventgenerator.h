#ifndef EVENTGENERATOR_H_AD6DSK71
#define EVENTGENERATOR_H_AD6DSK71

#include <memory>

#include "powheg/generateevents.h"

#include "integraltransformation.h"
#include "phasespacegenerator.h"

#include "run.h"

struct VegasState;

using Func = decltype(&Powheg::GenerateEvents);

struct IntegrandParams {
    IntegralTransformationPtr trafo;
    PhasespaceGeneratorPtr psgen;
    Func func;
    UserProcess::Data *userdata;
    bool print_func_params = false;
    bool absint = false;
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
    virtual void Load();
    virtual int GenerateEvents();

    void DebugIntegrate(const std::vector<double> &);

  protected:
    std::shared_ptr<UserProcess::Data> userdata;
    int rank;
    int num_procs;
    bool verbose = false;
    int NDIM;
    VegasState *vegasstate_ = 0;

  private:
    const char *run_file = "run.dat";
    const char *grid_file = "vegas_grid.dat";
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
    void XSec();

    struct Integral {
        double integral = 0.0;
        double error = 0.0;
    };
    Integral xsec(IntegrandParams &xsec, double abserror, double relerror);

    IntegralTransformationPtr trafo_;
    PhasespaceGeneratorPtr psgen_;
};

#endif
