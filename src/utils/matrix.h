#ifndef STRUCTURED_FINITE_VOLUME_UTILS_MATRIX_H
#define STRUCTURED_FINITE_VOLUME_UTILS_MATRIX_H

#include "project_defs.h"
#include <array>
#include <iostream>

namespace structured_fv {

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