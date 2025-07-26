#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_ROEHHT_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_ROEHHT_H

#include "euler_flux.h"
#include "numerical_flux_base.h"
#include "roe_state.h"

namespace structured_fv {
namespace euler {

// Roe flux with the Harten-Hyman entropy fix
class RoeHHFlux final : public NumericalFlux
{
  public:
    RoeHHFlux()
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

      Vec4<T> alphas{};
      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
          alphas[i] += Rinv(i, j)*(qR[j] - qL[j]);

      Vec4<T> q_kL = qL;
      Vec4<T> q_kR;
      Vec4<T> flux{};

      //TODO: getting R^T rather than R would be more efficient
      for (UInt i=0; i < 4; ++i)
      {
        q_kR = q_kL + alphas[i]*R(Column{i});
        T lambda_kL = compute_un(q_kL, normal);
        T lambda_kR = compute_un(q_kR, normal);
        if (i == 0 || i == 3)
        {
          Int sgn = i == 0 ? -1 : 1;
          T n_mag = std::sqrt(dot(normal, normal));
          lambda_kL += sgn*n_mag*compute_sos(q_kL);
          lambda_kR += sgn*n_mag*compute_sos(q_kR);
        }

        if (lambda_kL < 0 && lambda_kR > 0)
        {
          T beta = (lambda_kR - lambdas[i])/(lambda_kR - lambda_kL);
          T beta_c = 1 - beta;
          for (UInt j=0; j < 4; ++j)
          {
            flux[j] += (std::abs(lambda_kL)*beta  *alphas[i])*R(j, i);
            flux[j] += (std::abs(lambda_kR)*beta_c*alphas[i])*R(j, i);
          }
        } else
        {
          for (UInt j=0; j < 4; ++j)
            flux[j] += (std::abs(lambdas[i])*alphas[i])*R(j, i);
        }

        q_kL = q_kR;
      }


      auto fL = compute_euler_flux(qL, normal);
      auto fR = compute_euler_flux(qR, normal);

      for (UInt i=0; i < 4; ++i)
        flux[i] = 0.5*(fL[i] + fR[i] - flux[i]);
        
      return flux;
    }

    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& qL, const Vec4<T>& qR, 
                                 const Vec2<Real>& normal,
                                 Matrix<T, 4>& flux_dotL, Matrix<T, 4>& flux_dotR) const
    {
      throw std::runtime_error("not implemented");
    }     
};

}
}

#endif