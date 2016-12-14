#ifndef HISTLIST_H_WBFRFJPH
#define HISTLIST_H_WBFRFJPH

#include <vector>

#include "util/histogram.h"
#include "util/databuffer.h"
#include "libconfig.h"

class HistList {
  public:
    HistList() {}
    HistList(const HistList &o) {
        for (size_t i = 0; i < o.hists_.size(); i++) {
            auto hist = new Util::Histogram(*o.hists_[i]);
            hists_.push_back(hist);
        }
    }
    HistList &operator=(const HistList &) = delete;

    ~HistList() { dealloc(); }

    bool WriteToBuffer(Util::DataBuffer *buffer) const {
        buffer->AddUInt64((uint64_t)hists_.size());
        size_t n_written = 0;
        size_t num = hists_.size();
        for (size_t i = 0; i < num; i++) {
            int success = hists_[i]->WriteBinaryToBuffer(buffer);
            LIB_ASSERT(success == 0,
                       "error while writing histogram %lu to buffer: %d", i,
                       success);
            n_written += 1;
        }
        return n_written == hists_.size();
    }

    bool ReadFromBuffer(Util::DataBuffer *buffer) {
        uint64_t num = buffer->GetUInt64();

        dealloc();

        for (uint64_t i = 0; i < num; i++) {
            auto hist = new Util::Histogram(1, 0.0, 1.0);
            hist->ReadBinaryFromBuffer(buffer);
            hists_.push_back(hist);
        }
        return true;
    }

    void Append(int n, double xmin, double xmax) {
        hists_.push_back(new Util::Histogram(n, xmin, xmax));
    }

    void Append(int n, double *x) {
        hists_.push_back(new Util::Histogram(n, x));
    }

    bool Add(const HistList &o) {
        if (o.hists_.size() != hists_.size()) {
            return false;
        }
        bool status = true;
        for (size_t i = 0; i < hists_.size(); i++) {
            bool success = hists_[i]->AddHistogram(o.hists_[i]);
            LIB_ASSERT(success, "error: histogram %lu: wrong binning", i);
            if (!success) {
                status = false;
            }
        }
        return status;
    }

    Util::Histogram *operator[](int i) {
        if (i < 0) {
            return 0;
        }
        if ((size_t)i >= hists_.size()) {
            return 0;
        }
        return hists_[i];
    }

  private:
    void dealloc() {
        for (size_t i = 0; i < hists_.size(); i++) {
            if (hists_[i]) {
                delete hists_[i];
                hists_[i] = 0;
            }
        }
        hists_.clear();
    }

    std::vector<Util::Histogram *> hists_;
};

#endif /* end of include guard: HISTLIST_H_WBFRFJPH */
