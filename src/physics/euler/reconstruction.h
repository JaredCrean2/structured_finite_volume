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

// sign = +1 gives the solution at the right end of the cell,
// sign = -1 gives the solution at the left end of the cell
template <typename T, typename SlopeLimiter>
constexpr Vec4<T> reconstruct(const SlopeLimiter& limiter, const Vec4<T>& q_im1,
                              const Vec4<T>& q_i,
                              const Vec4<T>& q_ip1, Int sign)
{
  static_assert(common::IsSlopeLimiter<SlopeLimiter>);

  constexpr T epsilon = 1e-15;
  Vec4<T> u_endpoint{};
  for (UInt j=0; j < 4; ++j)
  {
    T denom = q_ip1[j] - q_i[j];
    T epsilon2 = denom >= 0 ? epsilon : -epsilon;
    T r = (q_i[j] - q_im1[j])/(denom + epsilon2);
    T slope = 0.5*(q_ip1[j] - q_im1[j]);
    u_endpoint[j] = q_i[j] + 0.5*sign*limiter(r)*slope;
  }

  return u_endpoint;
}

template <typename T, typename SlopeLimiter>
/*constexpr*/ Vec4<T> reconstruct_jac(const SlopeLimiter& limiter,
                                  const Vec4<T>& q_im1, Vec4<T>& q_im1_jac,
                                  const Vec4<T>& q_i, Vec4<T>& q_i_jac,
                                  const Vec4<T>& q_ip1, Vec4<T>& q_ip1_jac,
                                  Int sign)
{
  static_assert(common::IsSlopeLimiter<SlopeLimiter>);

  constexpr T epsilon = 1e-15;
  Vec4<T> u_endpoint{};
  for (UInt j=0; j < 4; ++j)
  {
    T denom = q_ip1[j] - q_i[j];
    T denom_dot_ip1 = 1;
    T denom_dot_i  = -1;

    T epsilon2 = denom >= 0 ? epsilon : -epsilon;
    T denom_pert = denom + epsilon2;

    T r = (q_i[j] - q_im1[j])/denom_pert;
    T r_dot_im1 = -1/denom_pert;
    T r_dot_i = 1/denom_pert - (q_i[j] - q_im1[j])*denom_dot_i/(denom_pert*denom_pert);
    T r_dot_ip1 = -(q_i[j] - q_im1[j])*denom_dot_ip1/(denom_pert*denom_pert);


    T slope = 0.5*(q_ip1[j] - q_im1[j]);
    T slope_dot_im1 = -0.5;
    T slope_dot_ip1 =  0.5;
    auto [phi, phi_dot] = limiter(r, T(1));
    T phi_dot_im1 = phi_dot * r_dot_im1;
    T phi_dot_i   = phi_dot * r_dot_i;
    T phi_dot_ip1 = phi_dot * r_dot_ip1;


    u_endpoint[j] = q_i[j] + 0.5*sign*phi*slope;
    q_im1_jac[j] = 0.5*sign*(phi_dot_im1*slope + phi*slope_dot_im1);
    q_i_jac[j]   = 1 + 0.5*sign*(phi_dot_i*slope);
    q_ip1_jac[j] = 0.5*sign*(phi_dot_ip1*slope + phi*slope_dot_ip1);   
  }

  return u_endpoint;
}

/*
// applies a slope limiter to each component of the input vectors
template <typename SlopeLimiter>
class ReconstructionSimple final : public ReconstructionBase
{
  static_assert(common::IsSlopeLimiter<SlopeLimiter>);

  public:
    using Limiter = SlopeLimiter;

    ReconstructionSimple(SlopeLimiter limiter) :
      m_limiter(limiter)
    {}

    //~ReconstructionSimple() override = default;


    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& q_im1, const Vec4<T>& q_i,
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
    // 
    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& q_im1, Matrix<T, 4>& q_im1_jac,
                                 const Vec4<T>& q_i,   Matrix<T, 4>& q_i_jac,
                                 const Vec4<T>& q_ip1, Matrix<T, 4>& q_ip1_jac,
                                 Int sign) const
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
*/

//template <typename SlopeLimiter>
//using ReconstructionConservative = ReconstructionSimple<SlopeLimiter>;

template <typename SlopeLimiter>
class ReconstructionConservative final : public ReconstructionBase
{
  public:
    using Limiter = SlopeLimiter;

    ReconstructionConservative(SlopeLimiter limiter) :
      m_limiter(limiter)
    {}

    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& q_im1, const Vec4<T>& q_i,
                                 const Vec4<T>& q_ip1, Int sign) const
    {
      return reconstruct(m_limiter, q_im1, q_i, q_ip1, sign);
    }

    //TODO: the jacobian is diagonal for conservative variables, but the interface
    //      for other reconstructions is dense.  Maybe introduce a custom type
    //      with a multiplication operator?
    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& q_im1, Matrix<T, 4>& q_im1_jac,
                                 const Vec4<T>& q_i,   Matrix<T, 4>& q_i_jac,
                                 const Vec4<T>& q_ip1, Matrix<T, 4>& q_ip1_jac, 
                                 Int sign) const
    {
      Vec4<T> q_end, q_im1_dot, q_i_dot, q_ip1_dot;
      q_end = reconstruct_jac(m_limiter, q_im1, q_im1_dot,
                                         q_i,   q_i_dot,
                                         q_ip1, q_ip1_dot, sign);

      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
        {
          q_im1_jac(i, j) = 0;
          q_i_jac(i, j)   = 0;
          q_ip1_jac(i, j) = 0;
        }
      
      for (UInt i=0; i < 4; ++i)
      {
        q_im1_jac(i, i) = q_im1_dot[i];
        q_i_jac(i, i)   = q_i_dot[i];
        q_ip1_jac(i, i) = q_ip1_dot[i];
      }

      return q_end;
    }    

  private:
    SlopeLimiter m_limiter;
};

// takes in conservative variables, converts to primitive, applies slope
// limiter, then converts result back to conservative
template <typename SlopeLimiter>
class ReconstructionPrimitive final : public ReconstructionBase
{
  public:
    using Limiter = SlopeLimiter;

    ReconstructionPrimitive(SlopeLimiter limiter) :
      m_limiter(limiter)
    {}

    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& q_im1, const Vec4<T>& q_i,
                                 const Vec4<T>& q_ip1, Int sign) const
    {
      Vec4<T> prim_im1 = compute_primitive_variables(q_im1);
      Vec4<T> prim_i   = compute_primitive_variables(q_i);
      Vec4<T> prim_ip1 = compute_primitive_variables(q_ip1);
      Vec4<T> prim_endpoint =  reconstruct(m_limiter, prim_im1, prim_i, prim_ip1, sign);
      return compute_conservative_variables(prim_endpoint, PrimitiveVarTag());
    }

    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& q_im1, Matrix<T, 4>& q_im1_jac,
                                 const Vec4<T>& q_i,   Matrix<T, 4>& q_i_jac,
                                 const Vec4<T>& q_ip1, Matrix<T, 4>& q_ip1_jac, 
                                 Int sign) const
    {
      auto [prim_im1, prim_jac_im1] = compute_primitive_variables_jac(q_im1);
      auto [prim_i, prim_jac_i]     = compute_primitive_variables_jac(q_i);
      auto [prim_ip1, prim_jac_ip1] = compute_primitive_variables_jac(q_ip1);

      Vec4<T> prim_endpoint, prim_im1_dot, prim_i_dot, prim_ip1_dot;
      prim_endpoint = reconstruct_jac(m_limiter, prim_im1, prim_im1_dot,
                                                 prim_i,   prim_i_dot,
                                                 prim_ip1, prim_ip1_dot, sign);

      // overwrite prim_jacs with d(prim_endpoint)/dq
      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
        {
          prim_jac_im1(i, j) *= prim_im1_dot[i];
          prim_jac_i(i, j)   *= prim_i_dot[i];
          prim_jac_ip1(i, j) *= prim_ip1_dot[i];
        }

      const auto& prim_endpoint_dot_im1 = prim_jac_im1;
      const auto& prim_endpoint_dot_i   = prim_jac_i;
      const auto& prim_endpoint_dot_ip1 = prim_jac_ip1;
      const auto [q_endpoint, q_endpoint_jac] = compute_conservative_variables_jac(prim_endpoint);

      //TODO: the primitive variable jacobians are sparse
      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
        {
          q_im1_jac(i, j) = 0;
          q_i_jac(i, j)   = 0;
          q_ip1_jac(i, j) = 0;
          for (UInt k=0; k < 4; ++k)
          {
            q_im1_jac(i, j) += q_endpoint_jac(i, k) * prim_endpoint_dot_im1(k, j);
            q_i_jac(i, j)   += q_endpoint_jac(i, k) * prim_endpoint_dot_i(k, j);
            q_ip1_jac(i, j) += q_endpoint_jac(i, k) * prim_endpoint_dot_ip1(k, j);
          }
        }

      return q_endpoint;
    }

  private:
    SlopeLimiter m_limiter;
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