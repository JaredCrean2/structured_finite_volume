#ifndef STRUCTURED_FINITE_VOLUME_UTILS_MATRIX_H
#define STRUCTURED_FINITE_VOLUME_UTILS_MATRIX_H

#include "project_defs.h"
#include "vec.h"
#include <iostream>

namespace structured_fv {

struct Column
{
  UInt j;
};

struct Row
{
  UInt i;
};

template <typename T, UInt M, UInt N=M>
class Matrix
{
  public:
    constexpr Matrix(const FixedVec<T, M*N>& vals) :
      m_data(vals)
    {}

    constexpr Matrix(const T& val)
    {
      for (UInt i=0; i < M*N; ++i)
        m_data[i] = val;
    }

    constexpr Matrix() :
      m_data{}
    {}

    constexpr Matrix<T, M, N>& operator=(const FixedVec<T, M*N>& other)
    {
      for (UInt i=0; i < M*N; ++i)
        m_data[i] = other[i];

      return *this;
    }

    constexpr Matrix<T, M, N>& operator=(const T& val)
    {
      for (UInt i=0; i < M*N; ++i)
        m_data[i] = val;

      return *this;
    }

    constexpr T& operator()(UInt i, UInt j) {return m_data[getIdx(i, j)]; }

    constexpr const T& operator()(UInt i, UInt j) const {return m_data[getIdx(i, j)]; }

    constexpr FixedVec<T, N> operator()(const Row& row) const
    {
      FixedVec<T, N> vals;
      for (UInt j=0; j < N; ++j)
        vals[j] = operator()(row.i, j);

      return vals;
    }

    constexpr FixedVec<T, M> operator()(const Column& col) const
    {
      FixedVec<T, M> vals;
      for (UInt i=0; i < M; ++i)
        vals[i] = operator()(i, col.j);

      return vals;
    }    

    constexpr T* getData() { return m_data.data(); }

    constexpr const T* getData() const { return m_data.data(); }

    constexpr UInt extent0() const { return M; }

    constexpr UInt extent1() const { return N; }

    constexpr std::pair<UInt, UInt> size() const { return std::pair<UInt, UInt>{M, N}; }

  private:
    constexpr UInt getIdx(UInt i, UInt j) const
    {
      assert(i >= 0 && i < M);
      assert(j >= 0 && j < N);

      return i*N + j;
    }

    FixedVec<T, M*N> m_data;
};

template <typename T, UInt M, UInt N, UInt N2>
constexpr FixedVec<T, M> operator*(const Matrix<T, M, N>& mat, const FixedVec<T, N2>& x)
{
  static_assert(N == N2, "matrix and vector dimensions must agree");

  FixedVec<T, M> b{};
  for (UInt i=0; i < M; ++i)
    for (UInt j=0; j < N; ++j)
      b[i] += mat(i, j) * x[j];

  return b;
}

template <typename T, UInt M, UInt N, UInt N2, UInt K>
constexpr Matrix<T, M, K> operator*(const Matrix<T, M, N>& mat, const Matrix<T, N2, K>& x)
{
  static_assert(N == N2, "matrix dimensions must agree");

  Matrix<T, M, K> b;
  for (UInt i=0; i < M; ++i)
    for (UInt j=0; j < N; ++j)
      for (UInt k=0; k < K; ++k)
        b(i, k) += mat(i, j) * x(j, k);

  return b;
}

template <typename T, UInt M, UInt N>
std::ostream& operator<<(std::ostream& os, const Matrix<T, M, N>& mat)
{
  for (UInt i=0; i < M; ++i)
  {
    for (UInt j=0; j < N; ++j)
      os << mat(i, j) << (j < N-1 ? "," : "");

    os << (i < M-1 ? "\n" : "");
  }

  return os;
}

}

#endif