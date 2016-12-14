#ifndef UTIL_HISTOGRAM_H_
#define UTIL_HISTOGRAM_H_

#include <ostream>

#include "util/databuffer.h"

namespace Util {

class Histogram {
  public:
    Histogram(int n, double min, double max);
    Histogram(int n, double *boundaries);
    Histogram(const Histogram& o);
    Histogram & operator=(const Histogram & ) = delete;
    virtual ~Histogram();

    bool AddValue(double x, double value);
    bool SetValue(double x, double value);

    bool Get(double x, double *dest) const;

    int GetBin(double x) const;

    int N() const { return n_; }
    double GetXMin() const { return bound_[0]; }
    double GetXMax() const { return bound_[n_]; }

    double GetX(int bin) const;
    double GetXLower(int bin) const;
    double GetXUpper(int bin) const;

    double GetBinContent(int bin) const;

    bool WriteToStream(std::ostream &of, bool norm) const;
    bool WriteToStreamScaled(std::ostream &of, bool norm, double scale) const;
    bool WriteToFile(const char *filename, bool norm) const;

    int WriteBinaryToBuffer(Util::DataBuffer * buffer);
    bool ReadBinaryFromBuffer(Util::DataBuffer * buffer);

    bool AddHistogram(const Histogram * h);
    bool AddScaledHistogram(const Histogram * h, double scale);

    void Reset();
    void Scale(double factor);

    double Sum() const;
    double Sum2() const;

    int GetNumEntries() const;
    int GetOverflow() const { return overflow_; }
    int GetUnderflow() const { return underflow_; }

    double Mean() const;
    double Variance() const;

  private:
    int n_;
    double *bins_;
    double *bins_sqr_;
    double *bound_;

    int overflow_;
    int underflow_;
    double underflow_v_;
    double underflow_v2_;
    double overflow_v_;
    double overflow_v2_;
    int *entries_;
};

} // Util

#endif

