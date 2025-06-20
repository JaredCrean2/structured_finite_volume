#ifndef STRUCTURED_FINITE_VOLUME_UTILS_ARRAY3_H
#define STRUCTURED_FINITE_VOLUME_UTILS_ARRAY3_H

#include "project_defs.h"
#include "vec.h"

namespace structured_fv {


template <typename T, UInt M, UInt N=M, UInt P=M>
class Array3
{
  public:
    constexpr Array3(const FixedVec<T, M*N*P>& vals) :
      m_data(vals)
    {}

    constexpr Array3() :
      m_data{}
    {}    

    constexpr T& operator()(UInt i, UInt j, UInt k) {return m_data[getIdx(i, j, k)]; }

    constexpr const T& operator()(UInt i, UInt j, UInt k) const {return m_data[getIdx(i, j, k)]; }   

    constexpr T* getData() { return m_data.data(); }

    constexpr const T* getData() const { return m_data.data(); }

    constexpr UInt extent0() const { return M; }

    constexpr UInt extent1() const { return N; }

    constexpr UInt extent2() const { return P; }

    constexpr Vec3<UInt> size() const { return Vec3<UInt>{M, N, P}; }

  private:
    constexpr UInt getIdx(UInt i, UInt j, UInt k) const
    {
      assert(i >= 0 && i < M);
      assert(j >= 0 && j < N);

      return i*N*P + j*P + k;
    }

    FixedVec<T, M*N*P> m_data;
};

}

#endif