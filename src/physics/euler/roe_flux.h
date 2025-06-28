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
    
    constexpr Vec4<Real> operator()(const Vec4<Real>& qL, const Vec4<Real>& qR, 
                                    const Vec2<Real>& normal) const
    {
      RoeAvgState<Real> avg_state = compute_roe_avg(qL, qR);
      Vec4<Real> q_avg = compute_conservative_variables(avg_state, RoeStateTag());

      Matrix<Real, 4, 4> R, Rinv;
      Vec4<Real> lambdas{};
      //TODO: it would be better of we could compute the eigendecomp
      //      directly from the Roe state
      compute_eigen_decomp(q_avg, normal, R, lambdas, Rinv);

      Vec4<Real> tmp{};
      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
          tmp[i] += Rinv(i, j)*(qR[j] - qL[j]);

      for (UInt i=0; i < 4; ++i)
      {
        Real lambda = std::abs(lambdas[i]);
        if (lambda < m_delta)
          lambda = (lambda*lambda + m_delta*m_delta)/(2*m_delta);
        tmp[i] *= lambda;
      }

      Vec4<Real> flux{};
      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
          flux[i] += R(i, j) * tmp[j];

      auto fL = compute_euler_flux(qL, normal);
      auto fR = compute_euler_flux(qR, normal);

      for (UInt i=0; i < 4; ++i)
        flux[i] = 0.5*(fL[i] + fR[i] - flux[i]);
        
      return flux;
    }

  private:
    Real m_delta;
};

}
}

#endif