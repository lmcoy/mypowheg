#ifndef REGIONS_H_
#define REGIONS_H_

#include <vector>

namespace FKS {

/**
 * @brief singular region
 *
 * Region represents a singular FKS region.
 */
struct Region {
    Region(int i = -1, int j = -1, unsigned resonance_id = 0)
        : I(i), J(j), ResonanceID(resonance_id) {}
    /** @brief radiated particle */
    int I;
    /**
    @brief mother particle
    J is the index of the particle that radiated particle I. If J >= 2, J is the
    actual index. If J == 0, both intial state particles can be the mother
    particle. J == -1(-2) means that particle 0(1) is the mother particle.
    */
    int J;
    unsigned int ResonanceID;
};

typedef std::vector<Region> RegionList;

} // namespace FKS

#endif /* end of include guard: REGIONS_H_ */
