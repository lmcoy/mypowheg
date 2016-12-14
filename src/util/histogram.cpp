#include <cassert>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <cstring>

#include "util/histogram.h"
#include "util/stringutil.h"

using namespace Util;

Histogram::Histogram(int n, double min, double max) {
    assert(max > min);
    assert(n > 0);

    underflow_ = 0;
    overflow_ = 0;
    underflow_v_ = 0.0;
    overflow_v_ = 0.0;
    underflow_v2_ = 0.0;
    overflow_v2_ = 0.0;
    entries_ = 0;
    n_ = n;
    bins_ = new double[n];
    bins_sqr_ = new double[n];
    bound_ = new double[n + 1];
    entries_ = new int[n];

    double width = (max - min) / (double)n;
    for (int i = 0; i < n; ++i) {
        bins_[i] = 0.0;
        bins_sqr_[i] = 0.0;
        bound_[i] = min + ((double)i) * width;
        entries_[i] = 0;
    }
    bound_[n] = max;
}

Histogram::Histogram(int n, double * boundaries) {
    assert(n > 0);
    underflow_ = 0;
    overflow_ = 0;
    underflow_v_ = 0.0;
    overflow_v_ = 0.0;
    underflow_v2_ = 0.0;
    overflow_v2_ = 0.0;
    entries_ = 0;
    n_ = n - 1;
    bins_ = new double[n - 1];
    bins_sqr_ = new double[n - 1];
    bound_ = new double[n];
    entries_ = new int[n - 1];
    for (int i = 0; i < n - 1; i++) {
        bins_[i] = 0.0;
        bins_sqr_[i] = 0.0;
        bound_[i] = boundaries[i];
        entries_[i] = 0;
    }
    bound_[n - 1] = boundaries[n - 1];
}

Histogram::Histogram(const Histogram &o) {
    underflow_ = o.underflow_;
    overflow_ = o.overflow_;
    underflow_v_ = o.underflow_v_;
    overflow_v_ = o.overflow_v_;
    underflow_v2_ = o.underflow_v2_;
    overflow_v2_ = o.overflow_v2_;
    n_ = o.n_;
    bins_ = new double[n_];
    bins_sqr_ = new double[n_];
    bound_ = new double[n_ + 1];
    entries_ = new int[n_];
    memcpy(bins_, o.bins_, sizeof(double) * n_);
    memcpy(bins_sqr_, o.bins_sqr_, sizeof(double) * n_);
    memcpy(bound_, o.bound_, sizeof(double) * (n_ + 1));
    memcpy(entries_, o.entries_, sizeof(int) * n_);
}

Histogram::~Histogram() {
    delete[] bins_;
    delete[] bins_sqr_;
    delete[] bound_;
    delete[] entries_;
}

int Histogram::GetBin(double x) const {
    if (x < bound_[0]) {
        return -1;
    }
    if (x > bound_[n_]) {
        return -2;
    }
    int left = 0;
    int right = n_;
    while (left <= right) {
        int len = right - left;
        int mid = left + (right - left) / 2;
        if (bound_[mid] == x) {
            return mid - 1;
        }
        if (len == 1) {
            return left;
        }
        if (x > bound_[mid]) {
            left = mid;
        } else {
            right = mid;
        }
    }
    assert(0);
    return -1;
}

bool Histogram::Get(double x, double *dest) const {
    int index = GetBin(x);
    if (index < 0) {
        return false;
    }
    *dest = bins_[index];
    return true;
}

bool Histogram::AddValue(double x, double value) {
    int index = GetBin(x);
    if (index < -2) {
        return false;
    }
    switch (index) {
    case -1:
        underflow_ += 1;
        underflow_v_ += value;
        underflow_v2_ += value * value;
        return false;
    case -2:
        overflow_ += 1;
        overflow_v_ += value;
        overflow_v2_ += value * value;
        return false;
    default:
        bins_[index] += value;
        bins_sqr_[index] += value * value;
        entries_[index] += 1;
    }
    return true;
}

bool Histogram::SetValue(double x, double value) {
    int index = GetBin(x);
    if (index < 0) {
        return false;
    }
    bins_[index] = value;
    bins_sqr_[index] = value * value;
    entries_[index] = 1;
    return true;
}

double Histogram::GetX(int bin) const {
    assert(bin >= 0);
    assert(bin < n_);
    return (GetXUpper(bin) + GetXLower(bin)) / 2.0;
}

double Histogram::GetXLower(int bin) const {
    assert(bin >= 0);
    assert(bin < n_);
    return bound_[bin];
}

double Histogram::GetXUpper(int bin) const {
    assert(bin >= 0);
    assert(bin < n_);
    return bound_[bin + 1];
}

bool Histogram::WriteToStream(std::ostream &of, bool norm) const {
    return WriteToStreamScaled(of, norm, 1.0);
}

bool Histogram::WriteToStreamScaled(std::ostream &of, bool norm, double scale) const {
    double sum = 0.0;
    if (norm) {
        for (int i = 0; i < n_; ++i) {
            // double len = bound_[i + 1] - bound_[i];
            sum += bins_[i];
        }
        sum += overflow_v_ + underflow_v_;
    } else {
        sum = 1.0;
    }
    sum /= scale;
    of << "# histogram\n";
    double sum_2 = 0.0;
    for (int i = 0; i < n_; ++i) {
        double y = 0.0;
        if (sum != 0.0) {
            y = bins_[i] / sum;
        }
        of << Strings::Format("%10g    %20g\n", bound_[i], y);
        sum_2 += y;
    }
    of << Strings::Format("# sum = %g\n", sum_2);
    of << Strings::Format("# number of entries (with over/underflow) = %d\n",
                          GetNumEntries());
    of << Strings::Format("# underflow = %d, overflow = %d\n", underflow_,
                          overflow_);
    of << Strings::Format("# sum with over and underflow bin = %g\n",
                          scale * Sum());

    return true;
}

bool Histogram::WriteToFile(const char *filename, bool norm) const {
    std::ofstream of(filename, std::ios::out | std::ios::trunc);
    if (!of.good()) {
        return false;
    }
    WriteToStream(of, norm);
    of.close();
    return true;
}

bool Histogram::AddHistogram(const Histogram *h) {
    if (n_ != h->n_) {
        return false;
    }
    for (int i = 0; i < n_; ++i) {
        if (fabs(bound_[i] - h->bound_[i]) > 1e-15) {
            return false;
        }
    }
    overflow_ += h->overflow_;
    overflow_v_ += h->overflow_v_;
    underflow_ += h->underflow_;
    underflow_v_ += h->underflow_v_;
    overflow_v2_ += h->overflow_v2_;
    underflow_v2_ += h->underflow_v2_;
    for (int i = 0; i < n_; ++i) {
        bins_[i] += h->bins_[i];
        bins_sqr_[i] += h->bins_sqr_[i];
        entries_[i] += h->entries_[i];
    }
    return true;
}

bool Histogram::AddScaledHistogram(const Histogram *h, double scale) {
    if (n_ != h->n_) {
        return false;
    }
    for (int i = 0; i < n_; ++i) {
        if (fabs(bound_[i] - h->bound_[i]) > 1e-15) {
            return false;
        }
    }
    overflow_ += h->overflow_;
    overflow_v_ += scale * h->overflow_v_;
    underflow_ += h->underflow_;
    underflow_v_ += scale * h->underflow_v_;
    overflow_v2_ += scale * scale * h->overflow_v2_;
    underflow_v2_ += scale * scale * h->underflow_v2_;
    for (int i = 0; i < n_; ++i) {
        bins_[i] += scale * h->bins_[i];
        bins_sqr_[i] += scale * scale * h->bins_sqr_[i];
        entries_[i] += h->entries_[i];
    }
    return true;
}

void Histogram::Reset() {
    overflow_ = 0;
    underflow_ = 0;
    overflow_v_ = 0.0;
    underflow_v_ = 0.0;
    underflow_v2_ = 0.0;
    overflow_v2_ = 0.0;
    for (int i = 0; i < n_; i++) {
        bins_[i] = 0.0;
        bins_sqr_[i] = 0.0;
        entries_[i] = 0;
    }
}

void Histogram::Scale(double factor) {
    for (int i = 0; i < n_; i++) {
        bins_[i] *= factor;
        bins_sqr_[i] *= factor * factor;
    }
    underflow_v_ *= factor;
    overflow_v_ *= factor;
    underflow_v2_ *= factor * factor;
    overflow_v2_ *= factor * factor;
}

double Histogram::Sum() const {
    double sum = underflow_v_ + overflow_v_;
    for (int i = 0; i < n_; ++i) {
        sum += bins_[i];
    }
    return sum;
}

double Histogram::Sum2() const {
    double sum2 = underflow_v2_ + overflow_v2_;
    for (int i = 0; i < n_; ++i) {
        sum2 += bins_sqr_[i];
    }
    return sum2;
}

int Histogram::WriteBinaryToBuffer(Util::DataBuffer *dbuffer) {

    if (dbuffer->AddInt64((int64_t)n_) == -1) {
        return -1;
    }
    if (dbuffer->AddInt64((int64_t)overflow_) == -1) {
        return -2;
    }
    if (dbuffer->AddInt64((int64_t)underflow_) == -1) {
        return -3;
    }
    if (dbuffer->AddDouble(underflow_v_) == -1) {
        return -4;
    }
    if (dbuffer->AddDouble(overflow_v_) == -1) {
        return -5;
    }
    if (dbuffer->AddDouble(underflow_v2_) == -1) {
        return -6;
    }
    if (dbuffer->AddDouble(overflow_v2_) == -1) {
        return -7;
    }
    if (dbuffer->AddDoubleArray(n_, bins_) == -1) {
        return -8;
    }
    if (dbuffer->AddDoubleArray(n_, bins_sqr_) == -1) {
        return -9;
    }
    if (dbuffer->AddDoubleArray(n_ + 1, bound_) == -1) {
        return -10;
    }
    if (dbuffer->AddIntArray(n_, entries_) == -1) {
        return -11;
    }
    return 0;
}

bool Histogram::ReadBinaryFromBuffer(DataBuffer * buffer) {
    int s = (int)buffer->GetUInt64();
    if (s != n_) {
        delete[] bins_;
        delete[] bins_sqr_;
        delete[] entries_;
        delete[] bound_;
        bins_ = new double[s];
        bins_sqr_ = new double[s];
        entries_ = new int[s];
        bound_ = new double[s + 1];
    }
    n_ = s;
    overflow_ = (int)buffer->GetInt64();
    underflow_ = (int)buffer->GetInt64();
    underflow_v_ = buffer->GetDouble();
    overflow_v_ = buffer->GetDouble();
    underflow_v2_ = buffer->GetDouble();
    overflow_v2_ = buffer->GetDouble();
    buffer->GetDoubleArray(n_, bins_);
    buffer->GetDoubleArray(n_, bins_sqr_);
    buffer->GetDoubleArray(n_ + 1, bound_);
    buffer->GetIntArray(n_, entries_);

    return true;
}

int Histogram::GetNumEntries() const {
    int ret = overflow_ + underflow_;
    for (int i = 0; i < n_; i++) {
        ret += entries_[i];
    }
    return ret;
}

double Histogram::GetBinContent(int bin) const {
    assert(bin >= 0);
    assert(bin < n_);
    return bins_[bin];
}

double Histogram::Mean() const {
    double sum = 0.0;
    double sum_y = 0.0;
    for (int i = 0; i < n_; i++) {
        double x = GetX(i);
        double y = GetBinContent(i);
        sum += x * y;
        sum_y += fabs(y);
    }
    if (sum_y == 0.0) {
        return 0.0;
    }
    return sum / sum_y;
}

double Histogram::Variance() const {
    double sum = 0.0;
    double sum_y = 0.0;
    double sum_sqr = 0.0;
    for (int i = 0; i < n_; i++) {
        double x = GetX(i);
        double y = GetBinContent(i);

        sum_sqr += x * x * y;
        sum += x * y;
        sum_y += fabs(y);
    }
    if (sum_y == 0.0) {
        return 0.0;
    }
    return (sum_sqr - sum * sum / sum_y ) / sum_y;
}

