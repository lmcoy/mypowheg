#ifndef UTIL_STATICMATRIX_H_
#define UTIL_STATICMATRIX_H_

#include <cassert>

namespace Util {

template <int Capacity> class StaticMatrix {
  public:
    explicit StaticMatrix(int nrows, int ncols, double def)
        : nrows_(nrows), ncols_(ncols) {
        assert(nrows * ncols <= Capacity);
        for (int i = 0; i < Capacity; i++) {
            data_[i] = def;
        }
    }

    void Reset(int nrows, int ncols, double v) {
        assert(nrows * ncols <= Capacity);
        nrows_ = nrows;
        ncols_ = ncols;
        for (int i = 0; i < Capacity; i++) {
            data_[i] = v;
        }
    }

    double Get(int row, int col) const {
        assert(row >= 0 && row < nrows_);
        assert(col >= 0 && col < ncols_);
        return data_[row * ncols_ + col];
    }

    void Set(int row, int col, double value) {
        assert(row >= 0 && row < nrows_);
        assert(col >= 0 && col < ncols_);
        data_[row * ncols_ + col] = value;
    }

    void SetAll(double value) {
        for (int i = 0; i < ncols_ * nrows_; i++) {
            data_[i] = value;
        }
    }

    void SetFromMatrix(const StaticMatrix<Capacity> &m) {
        nrows_ = m.nrows_;
        ncols_ = m.ncols_;
        memcpy(data_, m.data_, sizeof(double) * Capacity);
    }

    void MulRow(int row, double factor) {
        assert(row >= 0 && row < nrows_);
        for (int i = 0; i < ncols_; i++) {
            data_[row * ncols_ + i] *= factor;
        }
    }

    void MulAll(double factor) {
        for (int i = 0; i < ncols_ * nrows_; i++) {
            data_[i] *= factor;
        }
    }

    int Rows() const { return nrows_; }

    int Cols() const { return ncols_; }

  private:
    int nrows_;
    int ncols_;
    double data_[Capacity];
};

typedef StaticMatrix<16> StaticMatrix16;
typedef StaticMatrix<32> StaticMatrix32;
typedef StaticMatrix<64> StaticMatrix64;
typedef StaticMatrix<128> StaticMatrix128;

} // end namespace Util

#endif
