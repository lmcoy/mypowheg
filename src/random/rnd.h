#ifndef RANDOM_RND_H_
#define RANDOM_RND_H_

#include <random>
/*********************************************************************
 *                     random number generation                      *
 *********************************************************************/

namespace Random {

/// GetSeed returns a not predictable seed (hopefully).
int GetSeed();

struct RNG {
    RNG(int seed = 0);
    ~RNG();

    void Seed(int seed);

    double Random();

    std::mt19937 *r;
    std::uniform_real_distribution<> *dis;
};

void *init_rnd(void *args, int seed);

void free_rnd(void *r);

double get_rnd(void *r);

} // namespace Random

#endif

