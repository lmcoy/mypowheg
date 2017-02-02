#ifndef HISTOGRAMS_H
#define HISTOGRAMS_H

#include <vector>
#include <string>

#include "phasespace/phasespace.h"
#include "util/histogram.h"
#include "util/databuffer.h"

class Histograms {
  public:
    explicit Histograms(int n) : initialized_(0) {
        n_ = n;
        hists_ = new Util::Histogram *[n_];
        names.reserve(n);
        for (int i = 0; i < n_; ++i) {
            hists_[i] = 0;
            names.push_back(std::string(""));
        }
    }

    virtual ~Histograms() {
        for (int i = 0; i < n_; ++i) {
            if (hists_[i] != 0) {
                delete hists_[i];
            }
        }
        delete[] hists_;
    }

    virtual bool Fill(Phasespace::Phasespace const *const ps_lab, double wgt) {
        return false;
    }

    virtual bool Init1D(int n, int bins, double xmin, double xmax) {
        if (n < 0 || n >= n_) {
            return false;
        }
        if (hists_[n] != 0) {
            return false;
        }
        hists_[n] = new Util::Histogram(bins, xmin, xmax);
        initialized_ |= ((uint32_t)1u) << n;
        return true;
    }

    virtual int N() const { return n_; }

    virtual bool IsInit() const {
        uint32_t pattern = (((uint32_t)1) << n_) - 1;
        return pattern == initialized_;
    }

    virtual bool Write(int n, const char *filename, bool norm) {
        if (n < 0 || n >= n_) {
            return false;
        }
        if (hists_[n] == 0) {
            return false;
        }

        hists_[n]->WriteToFile(filename, norm);
        return true;
    }

    virtual bool WriteToStream(int n, std::ostream & ostr, bool norm) const {
        if (n < 0 || n >= n_) {
            assert(0);
            return false;
        }
        if (hists_[n] == 0) {
            assert(0);
            return false;
        }
        hists_[n]->WriteToStreamScaled(ostr, norm, 1.0);
        return true;
    }

    virtual void WriteBinaryToBuffer(Util::DataBuffer *buffer) const {
        if (!IsInit()) {
            return;
        }
        for (int i = 0; i < n_; ++i) {
            hists_[i]->WriteBinaryToBuffer(buffer);
        }
    }
    virtual void ReadBinaryFromBuffer(Util::DataBuffer *buffer) const {
        if (!IsInit()) {
            return;
        }
        for (int i = 0; i < n_; ++i) {
            hists_[i]->ReadBinaryFromBuffer(buffer);
        }
    }

    virtual void Add(const Histograms &h, double scale) {
        if (!IsInit()) {
            assert(0);
            return;
        }
        if (!h.IsInit()) {
            assert(0);
            return;
        }
        if (n_ != h.n_) {
            assert(0);
            return;
        }
        for (int i = 0; i < n_; ++i) {
            hists_[i]->AddScaledHistogram(h.hists_[i], scale);
        }
    }

    virtual void InitFrom(const Histograms &h) {
        if (n_ != h.n_) {
            return;
        }
        assert(h.IsInit() == true);
        initialized_ = h.initialized_;
        for (int i = 0; i < n_; ++i) {
            Init1D(i, h.hists_[i]->N(), h.hists_[i]->GetXMin(),
                   h.hists_[i]->GetXMax());
        }
    }

    void Reset() {
        for (int i = 0; i < n_; ++i) {
            hists_[i]->Reset();
        }
    }

    virtual Util::Histogram *GetHist(int n) {
        if (n >= n_) {
            return 0;
        }
        if (n < 0) {
            return 0;
        }
        return hists_[n];
    }

    virtual void SetName(int i, const std::string & name) {
        assert(i >= 0 && i < n_);
        names[i] = name;
    }

    virtual std::string GetName(int i) {
        assert(i >= 0 && i < n_);
        return names[i];
    }

  protected:
    Util::Histogram **hists_;
    std::vector<std::string> names;

  private:
    int n_;
    uint32_t initialized_;
};

#endif

