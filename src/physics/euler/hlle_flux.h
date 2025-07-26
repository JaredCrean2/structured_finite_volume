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
      throw std::runtime_error("not implemented");
    }  
};

}
}

#endif