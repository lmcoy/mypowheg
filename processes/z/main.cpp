#include <cmath>
#include <iostream>
#include <ctime>
#include <cassert>
#include <cstdio>
#include <cfenv>
#include <memory>

#include <mpi.h>

#include "powheggen/eventgenerator.h"

#include "breitwignermapping.h"
#include "bornpsgenerator.h"

#include "zdata.h"
#include "cuts.h"

void fail(int ret) {
    MPI_Finalize();
    exit(ret);
}

int main(int argc, char *argv[]) {
    MPI_Init(NULL, NULL);

    auto trafo = std::make_shared<BreitWignerMapping>();
    auto psgen = std::make_shared<BornPSGenerator>();
    EventGenerator generator(trafo, psgen, std::make_shared<ZData>());

    if (generator.Init() == -1) {
        fail(1);
    }
    generator.Setup();

    generator.GenerateEvents();

    MPI_Finalize();
    return 0;
}

