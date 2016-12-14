#include <iostream>
#include <random>
#include <cstdint>
#include <chrono>
#include <thread>

#include <sys/types.h>
#include <unistd.h>

#include "random/rnd.h"

namespace Random {

int GetSeed() {
    uint pid = (uint)getpid();

    uint time =
        std::chrono::high_resolution_clock::now().time_since_epoch().count();

    std::random_device rd;
    uint rdv = rd();

    void *h = malloc(sizeof(int));
    uint addr = (int64_t)h;
    free(h);

    uint addrs = 0;
    addrs = (int64_t) & addrs;

    std::seed_seq seq{ time, pid, rdv, addr, addrs };

    std::vector<std::uint32_t> seeds(1);

    seq.generate(seeds.begin(), seeds.end());

    return seeds[0];
}

RNG::RNG(int seed) {
    r = new std::mt19937(seed);
    dis = new std::uniform_real_distribution<>(0.0, 1.0);
}

RNG::~RNG() {
    delete r;
    delete dis;
}

void RNG::Seed(int seed) { r->seed(seed); }

double RNG::Random() { return (*dis)(*r); }

void *init_rnd(void *args, int seed) { return new RNG(seed); }

void free_rnd(void *r) {
    RNG *rng = (RNG *)r;
    delete rng;
}

double get_rnd(void *r) {
    RNG *rng = (RNG *)r;
    return (*rng->dis)(*rng->r);
}

} // end namespace random
