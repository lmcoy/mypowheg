#include <cassert>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <memory>

#include <mpi.h>

#include "powheggen/eventgenerator.h"

#include "phasespace.h"
#include "powheggen/integraltransformation.h"

#include "breitwignermapping.h"
#include "cuts.h"
#include "wjdata.h"

void fail(int ret) {
    MPI_Finalize();
    exit(ret);
}

int main(int argc, char *argv[]) {
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    MPI_Init(NULL, NULL);

    auto trafo = std::make_shared<BreitWignerMapping>();
    auto psgen = std::make_shared<BornPSGenerator>();
    EventGenerator generator(trafo, psgen, std::make_shared<WJData>());

    if (argc < 2) {
        std::cerr << "syntax: " << argv[0] << " x1 x2 ...\n";
        return -1;
    }
    std::vector<double> param;
    param.reserve(10);
    for (int i = 1; i < argc; i++) {
        double d = -1.0;
        if (Strings::ParseDouble(argv[i], &d) != 0) {
            std::cerr << "error in " << i << ": could not convert \"" << argv[i]
                      << "\" into number\n";
            return -1;
        }
        if (d < 0.0 || d > 1.0) {
            std::cerr << "error in " << i + 1 << ": " << d << " not in [0,1]\n";
            return -1;
        }
        param.push_back(d);
    }

    if (generator.Init() == -1) {
        fail(1);
    }
    generator.DebugIntegrate(param);

    MPI_Finalize();
    return 0;
}
