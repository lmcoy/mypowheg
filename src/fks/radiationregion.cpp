#include "fks/radiationregion.h"
#include <cmath>

#include "util/histogram.h"
#include "util/databuffer.h"
#include "libconfig.h"
#include "physics/pdgcode.h"

namespace FKS {

RadiationRegion::RadiationRegion(const FKS::FlavourConfig *fl)
    : FlavourConfig(fl) {
    double x[] = { 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6 };
    for (size_t i = 0; i < NPDF; i++) {
        InitHist.Append(11, x);
    }
}

RadiationRegion::~RadiationRegion() {
    for (size_t i = 0; i < NPDF; i++) {
        if (NormHist[i]) {
            delete NormHist[i];
        }
    }
}

bool RadiationRegion::CreateHistograms(int Nbins) {
    for (size_t i = 0; i < NPDF; i++) {
        double u = SearchNormHistUpperBound(i);
        if (u <= 0) {
            continue;
        }
        NormHist[i] = new Util::Histogram(Nbins, 0.0, u);
    }
    return true;
}

void RadiationRegion::ComputeNorm() {
    const char *born_pdgs =
        Physics::PDG::CodesToName(FlavourConfig->Born.Flavours).c_str();
    for (size_t i = 0; i < NPDF; i++) {
        if (NormHist[i] == 0) {
            break;
        }
        double started = false;
        int bins_with_zero = 0;
        double norm = -1.0;
        int bmax =
            static_cast<int>(static_cast<double>(NormHist[i]->N()) * 0.08);
        double m = 1e-5 * NormHist[i]->GetNumEntries();
        for (int bin = 0; bin < NormHist[i]->N(); bin++) {
            double y = NormHist[i]->GetBinContent(bin);
            if (y < m && started) {
                bins_with_zero += 1;
                if (bins_with_zero > bmax) {
                    norm = NormHist[i]->GetXUpper(bin - bmax);
                    break;
                }

            }
            if (y > 0.0) {
                started = true;
            }
        }
        LIB_ASSERT(norm > 0.0,
                   "born = %s, emitter = %d, pdf = %lu: It looks "
                   "like the distribution of norms is not close to "
                   "zero for large N in the histograms. You should "
                   "try to use more samples in InitFindNorm. Expected to find "
                   "at least %d bins without content but found only %d",
                   born_pdgs, Region.J, i, bmax, bins_with_zero);
        Norm[i] = norm;
        // double overflow = (double)NormHist[i]->GetOverflow();
        // double total = (double)NormHist[i]->GetNumEntries();
        // double p_over = overflow / total;
        // LIB_ASSERT(p_over < 0.002, "born = %d, emitter = %d, pdf = %lu: %.2f %% "
        //                            "of values are in the overflow bin. "
        //                            "Increase the number of points in the init "
        //                            "phase of FindNorm.",
        //            FlavourConfig->Born.ID, Region.J, i, 100.0 * p_over);
        // double mean = NormHist[i]->Mean();
        // double svariance = sqrt(NormHist[i]->Variance());
        // LIB_ASSERT(mean < NormHist[i]->GetXMax() / 10.0,
        //            "born = %s, emitter = %d, pdf = %lu: the mean of the norm "
        //            "histogram is larger than 10 %% of the whole range. This "
        //            "seems weird. Check what is going on! mean = %g, Nmax = %g",
        //            born_pdgs, Region.J, i, mean, NormHist[i]->GetXMax());
        // LIB_ASSERT(
        //     svariance < NormHist[i]->GetXMax() / 10.0,
        //     "born = %s, emitter = %d, pdf = %lu: the sqrt(variance) of the "
        //     "norm histogram is larger than 10 %% of the whole range, i.e. the "
        //     "width is pretty large. This will yield a bad performance. Check "
        //     "what is "
        //     "going on! "
        //     "sqrt(variance) = %g, Nmax = %g",
        //     born_pdgs, Region.J, i, svariance, NormHist[i]->GetXMax());
    }
}

bool
RadiationRegion::WriteNormHistBinaryToBuffer(Util::DataBuffer *buffer) const {
    uint64_t num_hists = 0;
    for (size_t i = 0; i < NPDF; i++) {
        if(NormHist[i] == 0) break;
        num_hists += 1;
    }
    buffer->AddUInt64((uint64_t)num_hists);
    int num = 0;
    for (size_t i = 0; i < num_hists; i++) {
        int success = NormHist[i]->WriteBinaryToBuffer(buffer);
        LIB_ASSERT(success == 0, "error while writing histogram %lu to buffer: %d",
                   i, success);
        num += 1;
    }

    return num != 0;
}

bool RadiationRegion::MergeNormHistBinaryFromBuffer(Util::DataBuffer *buffer) {
    uint64_t num = buffer->GetUInt64();
    LIB_ASSERT((uint64_t)NPDF >= num,
               "NPDF = %lu, but buffer has = %lu histograms", NPDF, num);
    uint64_t num_hists = 0;
    for (uint64_t i = 0; i < NPDF; i++) {
        if(NormHist[i] == 0) break;
        num_hists += 1;
    }
    LIB_ASSERT(num == num_hists,
               "RadiationRegion has %lu hists, but buffer has %lu", num_hists,
               num);
    for (size_t i = 0; i < num; i++) {
        Util::Histogram h(1, 0, 1);
        h.ReadBinaryFromBuffer(buffer);
        LIB_ASSERT(NormHist[i]->AddHistogram(&h),
                   "i=%lu: histograms do not have the same binning", i);
    }
    return true;
}

int RadiationRegion::GetBinarySize() const {
    Util::DataBuffer buffer(Util::DataBuffer::GETSIZE);
    WriteNormHistBinaryToBuffer(&buffer);
    return buffer.GetDataSize();
}

void RadiationRegion::WriteNormHistToStream(std::ostream &of) const {
    std::string born_pdgs =
        Physics::PDG::CodesToName(FlavourConfig->Born.Flavours);
    for (size_t i = 0; i < NPDF; i++) {
        if (!NormHist[i]) {
            break;
        }
        of << "set title \"born = " << born_pdgs << " pdf = " << i
           << " emitter = " << Region.J << "\"\n";
        of << "plot \"-\" using 1:2 with lines\n";
        NormHist[i]->WriteToStream(of, false);
        of << "e\n\n";
    }
}

bool RadiationRegion::AppendNormHistToFile(const char *filename) const {
    std::ofstream of(filename, std::ios::out | std::ios::app);
    if (!of.good()) {
        return false;
    }
    WriteNormHistToStream(of);
    of.close();
    return true;
}

double RadiationRegion::SearchNormHistUpperBound(size_t pdf) {
    Util::Histogram *hist = InitHist[pdf];

    if (hist->GetNumEntries() == 0) {
        return -1.0;
    }
    double threshold = 0.001 * hist->GetNumEntries();
    for (int bin = hist->N() - 1; bin >= 0; bin--) {
        double y = hist->GetBinContent(bin);
        if (y > threshold) {
            return hist->GetXUpper(bin);
        }
    }
    assert((double)hist->GetOverflow() < 4.0 * threshold &&
           "too many entries in overflow bin");
    return hist->GetXLower(0);
}

namespace {
int find_radregion(const RadiationRegionList &list,
                   const FKS::FlavourConfig &fl, FKS::Type_t type,
                   FKS::Region region) {
    size_t len = list.size();
    for (size_t r = 0; r < len; r++) {
        if (list[r].FlavourConfig != &fl) {
            continue;
        }
        if (list[r].Type != type) {
            continue;
        }
        if (list[r].Region.I == region.I && list[r].Region.J == region.J) {
            return r;
        }
    }
    return -1;
}
} // namespace

RadiationRegionList FindRadiationRegions(const FKS::ProcessList &l) {
    RadiationRegionList radlist;
    for (const auto &fl : l) {
        for (const auto &real : fl.Real) {
            for (const auto &region : real.Regions) {
                int pos = find_radregion(radlist, fl, real.Type, region);
                if (pos != -1) {
                    // region exists
                    radlist[pos].RealFlavour.push_back(&real);
                } else {
                    RadiationRegion radreg(&fl);
                    radreg.Region = region;
                    radreg.RealFlavour.push_back(&real);
                    radreg.Type = real.Type;
                    radlist.push_back(radreg);
                }
            }
        }
    }

    return radlist;
}

} /* FKS */
