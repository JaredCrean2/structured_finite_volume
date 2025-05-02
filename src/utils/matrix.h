#ifndef STRUCTURED_FINITE_VOLUME_UTILS_MATRIX_H
#define STRUCTURED_FINITE_VOLUME_UTILS_MATRIX_H

#include "project_defs.h"
#include <array>
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
    constexpr Matrix(const std::array<T, M*N>& vals) :
      m_data(vals)
    {}

    constexpr Matrix() :
      m_data{}
    {}    

    constexpr T& operator()(UInt i, UInt j) {return m_data[getIdx(i, j)]; }

    constexpr const T& operator()(UInt i, UInt j) const {return m_data[getIdx(i, j)]; }

    constexpr std::array<T, N> operator()(const Row& row) const
    {
      std::array<T, N> vals;
      for (UInt j=0; j < N; ++j)
        vals[j] = operator()(row.i, j);

      return vals;
    }

    constexpr std::array<T, M> operator()(const Column& col) const
    {
      std::array<T, M> vals;
      for (UInt i=0; i < M; ++i)
        vals[i] = operator()(i, col.j);

      return vals;
    }    

    constexpr T* getData() { m_data.data(); }

    constexpr const T* getData() const { m_data.data(); }

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

    std::array<T, M*N> m_data;
};

template <typename T, UInt M, UInt N, size_t N2>
constexpr std::array<T, M> operator*(const Matrix<T, M, N>& mat, const std::array<T, N2>& x)
{
  static_assert(N == N2, "matrix and vector dimensions must agree");

  std::array<T, M> b{};
  for (UInt i=0; i < M; ++i)
    for (UInt j=0; j < N; ++j)
      b[i] += mat(i, j) * x[j];

  return b;
}

template <typename T, UInt M, UInt N>
std::ostream& operator<<(std::ostream& os, const Matrix<T, M, N>& mat)
{
  for (size_t i=0; i < M; ++i)
  {
    for (size_t j=0; j < N; ++j)
      os << mat(i, j) << (j < N-1 ? "," : "");

    os << (i < M-1 ? "\n" : "");
  }

  return os;
}

}

#endif