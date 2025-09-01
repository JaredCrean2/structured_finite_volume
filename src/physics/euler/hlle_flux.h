#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_HLLE_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_HLLE_H

#include "euler_flux.h"
#include "numerical_flux_base.h"
#include "roe_state.h"

namespace structured_fv {
namespace euler {

class HLLEFlux final : public NumericalFlux
{
  public:
    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& qL, const Vec4<T>& qR, 
                                 const Vec2<Real>& normal) const
    {
      T n_mag = std::sqrt(dot(normal, normal));
      T aL = compute_sos(qL);
      T aR = compute_sos(qR);

      RoeAvgState<T> avg_state = compute_roe_avg(qL, qR);
      T aAvg = compute_sos(avg_state);

      T unL = compute_un(qL, normal);
      T unR = compute_un(qR, normal);
      T unAvg = compute_un(avg_state, normal);

      T s1 = std::min(unL - aL*n_mag, unAvg - aAvg*n_mag);
      T s2 = std::max(unR + aR*n_mag, unAvg + aAvg*n_mag);

      Vec4<T> fL = compute_euler_flux(qL, normal);
      Vec4<T> fR = compute_euler_flux(qR, normal);
      
      if (s1 < 0 && s2 > 0)
      {
        Vec4<T> flux{};
        for (UInt i=0; i < 4; ++i)
          flux[i] = (s2*fL[i] - s1*fR[i] + s1*s2*(qR[i] - qL[i]))/(s2 - s1);
                 
        return flux;
      } else if (s1 > 0)
      {
        return fL;
      } else  // s1 < s2 < 0
      {
        return fR;
      }
    }

    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& qL, const Vec4<T>& qR, 
                                 const Vec2<Real>& normal,
                                 Matrix<T, 4>& flux_dotL, Matrix<T, 4>& flux_dotR) const
    {
      Real n_mag2 = dot(normal, normal);
      Real n_mag = std::sqrt(n_mag2);
      auto [aL, aL_dotL] = compute_sos_jac(qL);
      auto [aR, aR_dotR] = compute_sos_jac(qR);

      RoeAvgState<T> avg_state;
      Matrix<T, 4> avg_state_dotL, avg_state_dotR;
      avg_state = compute_roe_avg_jac(qL, qR, avg_state_dotL, avg_state_dotR);
    
      auto [aAvg, aAvg_dot] = compute_sos_jac(avg_state);
      Vec4<T> aAvg_dotL{}, aAvg_dotR{};
      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
        {
          aAvg_dotL[i] += aAvg_dot[j] * avg_state_dotL(j, i);
          aAvg_dotR[i] += aAvg_dot[j] * avg_state_dotR(j, i);
        }

      auto [unL, unL_dotL] = compute_un_jac(qL, normal);
      auto [unR, unR_dotR] = compute_un_jac(qR, normal);
      auto [unAvg, unAvg_dot] = compute_un_jac(avg_state, normal);

      Vec4<T> unAvg_dotL{}, unAvg_dotR{};
      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
        {
          unAvg_dotL[i] += unAvg_dot[j] * avg_state_dotL(j, i);
          unAvg_dotR[i] += unAvg_dot[j] * avg_state_dotR(j, i);
        }
      
      T s1_L   = unL - aL*n_mag;
      T s1_avg = unAvg - aAvg*n_mag;
      T s1 = 0;
      Vec4<T> s1_dotL{}, s1_dotR{};
      if (s1_L <= s1_avg)
      {
        s1 = s1_L;
        s1_dotL = unL_dotL - aL_dotL * n_mag;
      } else
      {
        s1 = s1_avg;
        s1_dotL = unAvg_dotL - aAvg_dotL * n_mag;
        s1_dotR = unAvg_dotR - aAvg_dotR * n_mag;
      }

      T s2_R   = unR + aR*n_mag;
      T s2_avg = unAvg + aAvg*n_mag;
      T s2 = 0;
      Vec4<T> s2_dotL{}, s2_dotR{};
      if (s2_R >= s2_avg)
      {
        s2 = s2_R;
        s2_dotR = unR_dotR + aR_dotR*n_mag;
      } else
      {
        s2 = s2_avg;
        s2_dotL = unAvg_dotL + aAvg_dotL*n_mag;
        s2_dotR = unAvg_dotR + aAvg_dotR*n_mag;
      }

      Vec4<T> flux{};
      if (s1 < 0 && s2 > 0)
      {
        Vec4<T> fL, fR;
        std::tie(fL, flux_dotL) = compute_euler_flux_jac(qL, normal);
        std::tie(fR, flux_dotR) = compute_euler_flux_jac(qR, normal);
        for (UInt i=0; i < 4; ++i)
        {
          T num = (s2*fL[i] - s1*fR[i] + s1*s2*(qR[i] - qL[i]));
          T den = (s2 - s1);
          flux[i] = num/den;

          for (UInt j=0; j < 4; ++j)
          {
            T num_dotL = s2_dotL[j]*fL[i] + s2*flux_dotL(i, j) - s1_dotL[j]*fR[i] + 
                         (s1_dotL[j]*s2 + s1*s2_dotL[j])*(qR[i] - qL[i]) - s1*s2*del(i, j);
            T den_dotL = s2_dotL[j] - s1_dotL[j];


            T num_dotR = s2_dotR[j]*fL[i] - s1_dotR[j]*fR[i] - s1*flux_dotR(i, j) + 
                         (s1*s2_dotR[j] + s1_dotR[j]*s2)*(qR[i] - qL[i]) + s1*s2*del(i, j);
            T den_dotR = s2_dotR[j] - s1_dotR[j];
            
            flux_dotL(i, j) = num_dotL/den - num*den_dotL/(den*den);
            flux_dotR(i, j) = num_dotR/den - num*den_dotR/(den*den);
          }
        }
        return flux;
      } else if (s1 > 0)
      {
        std::tie(flux, flux_dotL) = compute_euler_flux_jac(qL, normal);
        flux_dotR = 0;
        return flux;
      } else  // s1 < s2 < 0
      {
        std::tie(flux, flux_dotR) = compute_euler_flux_jac(qR, normal);
        flux_dotL = 0;
        return flux;
      }        
      
    }  
};

}
}

#endif