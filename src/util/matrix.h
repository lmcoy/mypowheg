#ifndef UTIL_MATRIX_H_
#define UTIL_MATRIX_H_

#include <cassert>

namespace Util {

template <int M, int N> class Matrix {
  public:
    Matrix() {
        for (int i = 0; i < M * N; i++) {
            data_[i] = 0.0;
        }
    }
    double *operator[](int row) {
        assert(row >= 0 && row < M);
        return &(data_[row * N]);
    }
    const double *operator[](int row) const {
        assert(row >= 0 && row < M);
        return &(data_[row * N]);
    }

  private:
    double data_[M * N];
};

class Matrix2 {
  public:
    static const int CAPACITY = 8;
    explicit Matrix2(int n = 1, double value = 0.0) {
        assert(n <= CAPACITY);
        len_ = n;
        for (int i = 0; i < CAPACITY * CAPACITY; i++) {
            data_[i] = value;
        }
    }

    double Get(int i, int j) const {
        assert(i >= 0 && i < len_);
        assert(j >= 0 && j < len_);
        return data_[i * CAPACITY + j];
    }

    void Set(int i, int j, double value) {
        assert(i >= 0 && i < len_);
        assert(j >= 0 && j < len_);
        data_[i * CAPACITY + j] = value;
    }

    void SetAll(double value) {
        for (int i = 0; i < CAPACITY * CAPACITY; i++) {
            data_[i] = value;
        }
    }

    void Add(int i, int j, double value) {
        assert(i >= 0 && i < len_);
        assert(j >= 0 && j < len_);
        data_[i * CAPACITY + j] += value;
    }

    void MulAll(double factor) {
        for (int i = 0; i < CAPACITY * CAPACITY; i++) {
            data_[i] *= factor;
        }
    }

    int GetLen() const { return len_; }

    void SetLen(int n) {
        assert(n > 0 && n <= CAPACITY);
        len_ = n;
    }

    double * Data() {
        return data_;
    }

    int DataLen() const {
        return CAPACITY * CAPACITY;
    }

    static int IndexToRow(int index) {
        int column = IndexToColumn(index);
        return (index - column) / CAPACITY;
    }
    static int IndexToColumn(int index) {
        return index % CAPACITY;
    }

  private:
    double data_[CAPACITY * CAPACITY];
    int len_;
};

} // end namespace Math

#endif

