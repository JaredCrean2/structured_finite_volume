#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_LFFLUX_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_LFFLUX_H

#include "euler_flux.h"
#include "numerical_flux_base.h"

namespace structured_fv {
namespace euler {

class LaxFriedrichFlux final : public NumericalFlux
{
  public:
    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& qL, const Vec4<T>& qR, 
                                 const Vec2<Real>& normal) const
    {
      Real n_mag = sqrt(dot(normal, normal));
      T aL = compute_sos(qL);
      T aR = compute_sos(qR);
      T unL = compute_un(qL, normal);
      T unR = compute_un(qR, normal);

      T lambda_max = std::max(smoothAbs(unL) + aL*n_mag, smoothAbs(unR) + aR*n_mag);
      Vec4<T> fL = compute_euler_flux(qL, normal);
      Vec4<T> fR = compute_euler_flux(qR, normal);

      Vec4<T> flux{0, 0, 0, 0};
      for (UInt i=0; i < 4; ++i)
        flux[i] = 0.5*(fL[i] + fR[i]) - 0.5*lambda_max*(qR[i] - qL[i]);   

      return flux;
    }

    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& qL, const Vec4<T>& qR, 
                                 const Vec2<Real>& normal,
                                 Matrix<T, 4>& flux_dotL, Matrix<T, 4>& flux_dotR) const
    {
      Real n_mag = sqrt(dot(normal, normal));
      auto [aL, aL_dot] = compute_sos_jac(qL);
      auto [aR, aR_dot] = compute_sos_jac(qR);
      auto [unL, unL_dot] = compute_un_jac(qL, normal);
      auto [unR, unR_dot] = compute_un_jac(qR, normal);

      auto [abs_unL, abs_unL_dot] = smoothAbs_dot(unL, 1);
      Vec4<T> abs_unL_dotL = abs_unL_dot * unL_dot;
      auto [abs_unR, abs_unR_dot] = smoothAbs_dot(unR, 1);
      Vec4<T> abs_unR_dotR = abs_unR_dot * unR_dot;

      T lambda_maxL = abs_unL + aL*n_mag;
      Vec4<T> lambda_maxL_dot = abs_unL_dotL + aL_dot * n_mag;

      T lambda_maxR = abs_unR + aR*n_mag;
      Vec4<T> lambda_maxR_dot = abs_unR_dotR + aR_dot * n_mag;

      T lambda_max = std::max(lambda_maxL, lambda_maxR);
      Int is_left = lambda_maxL >= lambda_maxR;
      Int is_right = 1 - is_left;
      Vec4<T> lambda_max_dotL = is_left * lambda_maxL_dot;
      Vec4<T> lambda_max_dotR = is_right * lambda_maxR_dot;

      auto [fL, fL_dotL] = compute_euler_flux_jac(qL, normal);
      auto [fR, fR_dotR] = compute_euler_flux_jac(qR, normal);

      Vec4<T> flux{0, 0, 0, 0};
      for (UInt i=0; i < 4; ++i)
      {
        flux[i] = 0.5*(fL[i] + fR[i]) - 0.5*lambda_max*(qR[i] - qL[i]);
        for (UInt j=0; j < 4; ++j)
        {
          flux_dotL(i, j) = 0.5*(fL_dotL(i, j) - lambda_max_dotL[j]*(qR[i] - qL[i]));
          flux_dotR(i, j) = 0.5*(fR_dotR(i, j) - lambda_max_dotR[j]*(qR[i] - qL[i]));
        }

        flux_dotL(i, i) += 0.5*lambda_max;
        flux_dotR(i, i) -= 0.5*lambda_max;
      }

      return flux;
    }    
};

}
}

#endif