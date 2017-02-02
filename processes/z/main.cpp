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

void print_usage(char * argv0) {
    printf( "usage: %s [-load]\n", argv0);
    printf("  -load is optinal. If -load is used the cross section and grid is "
           "loaded from a file.\n");
}

int main(int argc, char *argv[]) {
    MPI_Init(NULL, NULL);

    auto trafo = std::make_shared<BreitWignerMapping>();
    auto psgen = std::make_shared<BornPSGenerator>();
    EventGenerator generator(trafo, psgen, std::make_shared<ZData>());

    if (generator.Init() == -1) {
        fail(1);
    }
    if (argc == 1) {
        generator.Setup();
    }
    if (argc == 2) {
        if (strcmp(argv[1], "-load") == 0) {
            generator.Load();
        } else {
            print_usage(argv[0]);
            fail(1);
        }
    }
    if (argc > 2) {
        print_usage(argv[0]);
        fail(1);
    }

    generator.GenerateEvents();

    MPI_Finalize();
    return 0;
}

