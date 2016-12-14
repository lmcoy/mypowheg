#ifndef RADIATIONREGION_H_NHHWPYQN
#define RADIATIONREGION_H_NHHWPYQN

#include <fstream>
#include <iostream>

#include "fks/regions.h"
#include "fks/process.h"
#include "fks/histlist.h"

namespace Util {
class Histogram;
class DataBuffer;
}

namespace FKS {

struct RadiationRegion {
    static const size_t NPDF = 3;
    RadiationRegion(const FKS::FlavourConfig *fl);
    ~RadiationRegion();

    bool CreateHistograms(int Nbins);

    void ComputeNorm();

    bool WriteNormHistBinaryToBuffer(Util::DataBuffer *buffer) const;
    bool MergeNormHistBinaryFromBuffer(Util::DataBuffer *buffer);
    int GetBinarySize() const;
    void WriteNormHistToStream(std::ostream &of) const;
    bool AppendNormHistToFile(const char *filename) const;

    FKS::Region Region;
    const FKS::FlavourConfig * FlavourConfig;

    std::vector<const Real_t *> RealFlavour;
    Type_t Type;
    double Norm[NPDF] = { 0.0 };
    Util::Histogram *NormHist[NPDF] = { 0 };
    HistList InitHist;

    double SearchNormHistUpperBound(size_t pdf);

    bool HasNormHist() const { return created_; }

  private:
    bool created_;
    int getind(int i, int j) { return j + i * 11; }
};

typedef std::vector<RadiationRegion> RadiationRegionList;

RadiationRegionList FindRadiationRegions(const FKS::ProcessList &l);

} /* FKS */

#endif /* end of include guard: RADIATIONREGION_H_NHHWPYQN */
