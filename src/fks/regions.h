#ifndef REGIONS_H_
#define REGIONS_H_

#include <vector>

namespace FKS { 

struct Region {
    Region(int i = -1, int j = -1, unsigned resonance_id = 0)
        : I(i), J(j), ResonanceID(resonance_id) {}
    int I, J;
    unsigned int ResonanceID;
};

typedef std::vector<Region> RegionList;

} // namespace FKS

#endif /* end of include guard: REGIONS_H_ */

