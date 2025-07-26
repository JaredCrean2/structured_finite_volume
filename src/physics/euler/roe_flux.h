#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_ROE_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_ROE_H

#include "euler_flux.h"
#include "numerical_flux_base.h"
#include "roe_state.h"


namespace structured_fv {
namespace euler {

// Roe flux with Hartens entropy fix.  Setting delta = 0 gives
// Roe's original scheme
class RoeFlux final : public NumericalFlux
{
  public:
    RoeFlux(Real delta=0.3) :
      m_delta(delta)
    {}
    
    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& qL, const Vec4<T>& qR, 
                                 const Vec2<Real>& normal) const
    {
      RoeAvgState<T> avg_state = compute_roe_avg(qL, qR);
      Vec4<T> q_avg = compute_conservative_variables(avg_state, RoeStateTag());

      Matrix<T, 4, 4> R, Rinv;
      Vec4<T> lambdas{};
      //TODO: it would be better of we could compute the eigendecomp
      //      directly from the Roe state
      compute_eigen_decomp(q_avg, normal, R, lambdas, Rinv);

      Vec4<T> tmp{};
      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
          tmp[i] += Rinv(i, j)*(qR[j] - qL[j]);

      for (UInt i=0; i < 4; ++i)
      {
        T lambda = smoothAbs(lambdas[i]);
        if (lambda < m_delta)
          lambda = (lambda*lambda + m_delta*m_delta)/(2*m_delta);
        tmp[i] *= lambda;
      }

      Vec4<T> flux{};
      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
          flux[i] += R(i, j) * tmp[j];

      auto fL = compute_euler_flux(qL, normal);
      auto fR = compute_euler_flux(qR, normal);

      for (UInt i=0; i < 4; ++i)
        flux[i] = 0.5*(fL[i] + fR[i] - flux[i]);
        
      return flux;
    }

    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& qL, const Vec4<T>& qR, 
                                 const Vec2<T>& normal,
                                 Matrix<T, 4>& flux_dotL, Matrix<T, 4>& flux_dotR) const
    {
      throw std::runtime_error("not implemented");
    }     

  private:
    Real m_delta;
};

}
}

#endif