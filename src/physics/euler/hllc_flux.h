#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_HLLC_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_HLLC_H

#include "euler_flux.h"
#include "numerical_flux_base.h"
#include "roe_state.h"


namespace structured_fv {
namespace euler {

class HLLCFlux final : public NumericalFlux
{
  public:
    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& qL, const Vec4<T>& qR, 
                                 const Vec2<Real>& normal) const
    {
      T n_mag2 = dot(normal, normal);
      T n_mag = std::sqrt(n_mag2);
      T aL = compute_sos(qL);
      T aR = compute_sos(qR);

      RoeAvgState avg_state = compute_roe_avg(qL, qR);
      T aAvg = compute_sos(avg_state);

      T unL = compute_un(qL, normal);
      T unR = compute_un(qR, normal);
      T unAvg = compute_un(avg_state, normal);

      T sL = std::min(unL - aL*n_mag, unAvg - aAvg*n_mag);
      T sR = std::max(unR + aR*n_mag, unAvg + aAvg*n_mag);


      //std::cout << "sL = " << sL << std::endl;
      //std::cout << "sR = " << sR << std::endl;

      if (sL > 0)
      {
        return compute_euler_flux(qL, normal);
      } else if (sR < 0)
      {
        return compute_euler_flux(qR, normal);
      } else
      {
        T pL = compute_pressure(qL);
        T pR = compute_pressure(qR);
        T deltaSL = sL - unL;
        T deltaSR = sR - unR;

        T s_star = ((pR - pL)*n_mag2 + qL[0]*unL*deltaSL - qR[0]*unR*deltaSR)/(qL[0]*deltaSL - qR[0]*deltaSR);
        //std::cout << "s_star = " << s_star << std::endl;

        Vec4<T> fL = compute_euler_flux(qL, normal);
        Vec4<T> fR = compute_euler_flux(qR, normal);

        const auto& qk = s_star > 0 ? qL : qR;
        const auto& fk = s_star > 0 ? fL : fR;
        const auto sk  = s_star > 0 ? sL : sR;
        const auto pk  = s_star > 0 ? pL : pR;
        const auto unk = s_star > 0 ? unL : unR;

        Vec4<T> d{0, normal[0]/n_mag, normal[1]/n_mag, s_star/n_mag};
        T pressure_term = sk*(pk*n_mag2 + qk[0]*(sk - unk)*(s_star - unk))/n_mag;

        Vec4<T> flux{0, 0, 0, 0};
        for (UInt i=0; i < 4; ++i)
          flux[i] = (s_star*(sk*qk[i] - fk[i]) + pressure_term*d[i])/(sk - s_star);

        return flux;
      }
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