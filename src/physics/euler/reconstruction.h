#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_RECONSTRUCTION_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_RECONSTRUCTION_H

#include "reconstruction_enum.h"
#include "physics/common/slope_limiters.h"
#include "physics/euler/euler_flux.h"
#include "physics/euler/typedefs.h"
#include "euler_flux.h"

namespace structured_fv {
namespace euler {

class ReconstructionBase
{
  public:
    virtual ~ReconstructionBase() = default;
};

template <typename ReconstructionType>
constexpr bool IsReconstruction = std::is_base_of_v<ReconstructionBase, ReconstructionType>;


// applies a slope limiter to each component of the input vectors
template <typename SlopeLimiter>
class ReconstructionSimple final : public ReconstructionBase
{
  static_assert(common::IsSlopeLimiter<SlopeLimiter>);

  public:
    ReconstructionSimple(SlopeLimiter limiter) :
      m_limiter(limiter)
    {}

    //~ReconstructionSimple() override = default;

    // sign = +1 gives the solution at the right end of the cell,
    // sign = -1 gives the solution at the left end of the cell
    template <typename T>
    /*constexpr*/ Vec4<T> operator()(const Vec4<T>& q_im1, const Vec4<T>& q_i,
                                 const Vec4<T>& q_ip1, Int sign) const
    {
      constexpr T epsilon = 1e-15;
      Vec4<T> u_endpoint{};
      for (UInt j=0; j < 4; ++j)
      {
        T denom = q_ip1[j] - q_i[j];
        T epsilon2 = denom >= 0 ? epsilon : -epsilon;
        T r = (q_i[j] - q_im1[j])/(denom + epsilon2);
        T slope = 0.5*(q_ip1[j] - q_im1[j]);
        u_endpoint[j] = q_i[j] + 0.5*sign*m_limiter(r)*slope;
      }

      return u_endpoint;
    }

  private:
    SlopeLimiter m_limiter;
};


template <typename SlopeLimiter>
using ReconstructionConservative = ReconstructionSimple<SlopeLimiter>;


template <typename SlopeLimiter>
class ReconstructionPrimitive final : public ReconstructionBase
{
  public:
    ReconstructionPrimitive(SlopeLimiter limiter) :
      m_reconstruction(limiter)
    {}

    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& q_im1, const Vec4<T>& q_i,
                                 const Vec4<T>& q_ip1, Int sign) const
    {
      Vec4<T> prim_im1 = compute_primitive_variables(q_im1);
      Vec4<T> prim_i   = compute_primitive_variables(q_i);
      Vec4<T> prim_ip1 = compute_primitive_variables(q_ip1);
      Vec4<T> prim_endpoint =  m_reconstruction(prim_im1, prim_i, prim_ip1, sign);
      return compute_conservative_variables(prim_endpoint, PrimitiveVarTag());
    }

  private:
    ReconstructionSimple<SlopeLimiter> m_reconstruction;
};


template <typename SlopeLimiter>
constexpr Reconstruction get_enum(const ReconstructionConservative<SlopeLimiter>&)
{
  return Reconstruction::Conservative;
}

template <typename SlopeLimiter>
constexpr Reconstruction get_enum(const ReconstructionPrimitive<SlopeLimiter>&)
{
  return Reconstruction::Primitive;
}

template <typename ReconstructionType, 
          std::enable_if_t<IsReconstruction<ReconstructionType>, bool> = true>
std::ostream& operator<<(std::ostream& os, const ReconstructionType& recon)
{
  os << get_enum(recon);
  return os;
}

}
}

#endif