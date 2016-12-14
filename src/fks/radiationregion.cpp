#include "fks/radiationregion.h"
#include <cmath>

#include "util/histogram.h"
#include "util/databuffer.h"
#include "libconfig.h"
#include "physics/pdgcode.h"

namespace FKS {

RadiationRegion::RadiationRegion(const FKS::FlavourConfig *fl)
    : FlavourConfig(fl) {
    double x[] = { 1e-07,   1e-06,  5e-06,  1e-05,  5e-05,  0.0001,
                   0.00025, 0.0005, 0.001,  0.0025, 0.005,  0.01,
                   0.025,   0.05,   0.1,    0.25,   0.5,    1,
                   2.5,     5.0,    10,     25.0,   50.0,   100,
                   500.0,   1000,   5000.0, 10000,  100000, 1000000 };
    for (size_t i = 0; i < NPDF; i++) {
        InitHist.Append(30, x);
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

void RadiationRegion::ComputeNorm(double Ratio) {
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
        double sum = NormHist[i]->GetUnderflow();
        double N = NormHist[i]->Sum();
        double n2 = 1e13;
        for (int bin = 0; bin < NormHist[i]->N(); bin++) {
            double y = NormHist[i]->GetBinContent(bin);
            sum += y;
            if (sum > N * Ratio && n2 > 1e12) {
                n2 = NormHist[i]->GetXUpper(bin);
            }
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
        if (norm < 0.0) {
            norm = 1e13;
        }
        norm = std::min(norm, n2);
        Norm[i] = norm;
        LIB_ASSERT(norm > 0.0,
                   "error: did not find norm for upper bounding function.");
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
               "NPDF = %lu, but buffer has = %llu histograms", NPDF, num);
    uint64_t num_hists = 0;
    for (uint64_t i = 0; i < NPDF; i++) {
        if(NormHist[i] == 0) break;
        num_hists += 1;
    }
    LIB_ASSERT(num == num_hists,
               "RadiationRegion has %llu hists, but buffer has %llu", num_hists,
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
        of << "set xlabel \"N\"\n";
        of << "set ylabel \"ratio\"\n";
        of << "a = " << NormHist[i]->GetUnderflow() << "\n";
        of << "b = " << NormHist[i]->GetNumEntries() << "\n";
        of << "csum(x)=(a=a+x,a)/b\n";
        of << "plot \"-\" using 1:(csum($2)) with lines ti \"cumulative distribution/total\"\n";
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
