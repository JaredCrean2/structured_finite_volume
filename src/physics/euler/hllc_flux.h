#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_HLLC_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_HLLC_H

#include "euler_flux.h"
#include "numerical_flux_base.h"
#include "roe_state.h"
#include "utils/math.h"


namespace structured_fv {
namespace euler {

class HLLCFlux final : public NumericalFlux
{
  public:
    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& qL, const Vec4<T>& qR, 
                                 const Vec2<Real>& normal) const
    {
      bool do_output = imag(qR[0]) > 0;

      std::cout << "\ncomputing value" << std::endl;
      T n_mag2 = dot(normal, normal);
      T n_mag = std::sqrt(n_mag2);
      T aL = compute_sos(qL);
      T aR = compute_sos(qR);

      RoeAvgState<T> avg_state = compute_roe_avg(qL, qR);
      T aAvg = compute_sos(avg_state);

      T unL = compute_un(qL, normal);
      T unR = compute_un(qR, normal);
      T unAvg = compute_un(avg_state, normal);

      T sL = std::min(unL - aL*n_mag, unAvg - aAvg*n_mag);
      T sR = std::max(unR + aR*n_mag, unAvg + aAvg*n_mag);

      if (do_output)
      {
        PRINTVAR(unL - aL*n_mag);
        PRINTVAR(unAvg - aAvg*n_mag);
        PRINTVAR(sL);
      }

      if (do_output)
      {
        PRINTVAR(avg_state);
        PRINTVAR(unR);
        PRINTVAR(aR);
        PRINTVAR(unAvg);
        PRINTVAR(aAvg);

        T sR_R = unR + aR*n_mag;
        T sR_avg = unAvg + aAvg*n_mag;

        PRINTVAR(sR_R);
        PRINTVAR(imag(sR_R)/1e-80);
        PRINTVAR(sR_avg);
        PRINTVAR(imag(sR_avg)/1e-80);


        PRINTVAR(sL);
        PRINTVAR(sR);
        PRINTVAR(imag(unL)/1e-80);
        PRINTVAR(imag(aL)/1e-80);
        PRINTVAR(imag(unAvg)/1e-80);
        PRINTVAR(imag(aAvg)/1e-80);

        PRINTVAR(imag(sL)/1e-80);
        PRINTVAR(imag(sR)/1e-80);
      }

      if (sL > 0)
      {
        PRINTCONST("case 1");
        return compute_euler_flux(qL, normal);
      } else if (sR < 0)
      {
        PRINTCONST("case 2");
        return compute_euler_flux(qR, normal);
      } else
      {
        PRINTCONST("case 3");
        T pL = compute_pressure(qL);
        T pR = compute_pressure(qR);
        T deltaSL = sL - unL;
        T deltaSR = sR - unR;

        T s_star = ((pR - pL)*n_mag2 + qL[0]*unL*deltaSL - qR[0]*unR*deltaSR)/(qL[0]*deltaSL - qR[0]*deltaSR);
        if (do_output)
        {
          std::cout << "num_dot = " << imag(((pR - pL)*n_mag2 + qL[0]*unL*deltaSL - qR[0]*unR*deltaSR))/1e-80 << std::endl;
          std::cout << "den_dot = " << imag((qL[0]*deltaSL - qR[0]*deltaSR))/1e-80 << std::endl;
        }
        //std::cout << "s_star = " << s_star << std::endl;

        Vec4<T> fL = compute_euler_flux(qL, normal);
        Vec4<T> fR = compute_euler_flux(qR, normal);

        if (do_output)
          PRINTVAR(s_star);
        const auto& qk = s_star > 0 ? qL : qR;
        const auto& fk = s_star > 0 ? fL : fR;
        const auto sk  = s_star > 0 ? sL : sR;
        const auto pk  = s_star > 0 ? pL : pR;
        const auto unk = s_star > 0 ? unL : unR;

        Vec4<T> d{0, normal[0]/n_mag, normal[1]/n_mag, s_star/n_mag};
        T pressure_term = sk*(pk*n_mag2 + qk[0]*(sk - unk)*(s_star - unk))/n_mag;

        if (do_output)
        {
          PRINTVAR(imag(s_star)/1e-80);
          PRINTVAR(pressure_term);
          PRINTVAR(imag(pressure_term)/1e-80);
        }

        Vec4<T> flux{0, 0, 0, 0};
        for (UInt i=0; i < 4; ++i)
        {
          T num = (s_star*(sk*qk[i] - fk[i]) + pressure_term*d[i]);
          T den = (sk - s_star);

          if (i == 0 && do_output)
          {
            PRINTVAR(num);
            PRINTVAR(imag(num)/1e-80);
            PRINTVAR(den);
            PRINTVAR(imag(den)/1e-80);
          }
          flux[i] = (s_star*(sk*qk[i] - fk[i]) + pressure_term*d[i])/(sk - s_star);
        }

        return flux;
      }
    }

    template <typename T>
    constexpr Vec4<T> operator()(const Vec4<T>& qL, const Vec4<T>& qR, 
                                 const Vec2<Real>& normal,
                                 Matrix<T, 4>& flux_dotL, Matrix<T, 4>& flux_dotR) const
    {
      std::cout << "\ncomputing jacobian" << std::endl;
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

      //T unL = compute_un(qL, normal);
      auto [unL, unL_dotL] = compute_un_jac(qL, normal);
      //T unR = compute_un(qR, normal);
      auto [unR, unR_dotR] = compute_un_jac(qR, normal);
      //T unAvg = compute_un(avg_state, normal);
      auto [unAvg, unAvg_dot] = compute_un_jac(avg_state, normal);

      PRINTVAR(normal);
      PRINTVAR(unAvg_dot);
      PRINTVAR(avg_state_dotR);
      Vec4<T> unAvg_dotL{}, unAvg_dotR{};
      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
        {
          unAvg_dotL[i] += unAvg_dot[j] * avg_state_dotL(j, i);
          unAvg_dotR[i] += unAvg_dot[j] * avg_state_dotR(j, i);
        }

      //T sL = std::min(unL - aL*n_mag, unAvg - aAvg*n_mag);
      //T sR = std::max(unR + aR*n_mag, unAvg + aAvg*n_mag);

      PRINTVAR(avg_state)
      PRINTVAR(unL);
      PRINTVAR(aL);
      PRINTVAR(unAvg);
      PRINTVAR(aAvg);
      T sL_L   = unL - aL*n_mag;
      T sL_avg = unAvg - aAvg*n_mag;
      T sL = 0;
      Vec4<T> sL_dotL{}, sL_dotR{};
      if (sL_L <= sL_avg)
      {
        PRINTCONST("sL_L");
        sL = sL_L;
        sL_dotL = unL_dotL - aL_dotL * n_mag;
      } else
      {
        PRINTCONST("sL_avg");
        sL = sL_avg;
        sL_dotL = unAvg_dotL - aAvg_dotL * n_mag;
        sL_dotR = unAvg_dotR - aAvg_dotR * n_mag;
      }

      //PRINTVAR(unL_dotL[1]);
      //PRINTVAR(aL_dotL[1]);
      PRINTVAR(unAvg_dotR[1]);
      PRINTVAR(aAvg_dotR[1]);
      PRINTVAR(sL);
      PRINTVAR(sL_dotR[1]);

      T sR_R   = unR + aR*n_mag;
      T sR_avg = unAvg + aAvg*n_mag;
      T sR = 0;
      Vec4<T> sR_dotL{}, sR_dotR{};
      if (sR_R >= sR_avg)
      {
        PRINTCONST("sR_R");
        sR = sR_R;
        sR_dotR = unR_dotR + aR_dotR*n_mag;
      } else
      {
        PRINTCONST("sR_avg");
        sR = sR_avg;
        sR_dotL = unAvg_dotL + aAvg_dotL*n_mag;
        sR_dotR = unAvg_dotR + aAvg_dotR*n_mag;
      }

      PRINTVAR(sR_R);
      PRINTVAR(sL_dotR[0]);
      PRINTVAR(sR_dotR[0]);


      std::cout << "sL = " << sL << std::endl;
      std::cout << "sR = " << sR << std::endl;
      Vec4<T> flux{};
      if (sL > 0)
      {
        std::cout << "case 1" << std::endl;
        std::tie(flux, flux_dotL) = compute_euler_flux_jac(qL, normal);
        flux_dotR = 0;
        return flux;
      } else if (sR < 0)
      {
        std::cout << "case 2" << std::endl;
        std::tie(flux, flux_dotR) = compute_euler_flux_jac(qR, normal);
        flux_dotL = 0;
        return flux;
      } else
      {
        std::cout << "case 3" << std::endl;
        //TODO: maybe delay calculation of some derivatives until inside this branch

        //T pL = compute_pressure(qL);
        auto [pL, pL_dotL] = compute_pressure_jac(qL);
        //T pR = compute_pressure(qR);
        auto [pR, pR_dotR] = compute_pressure_jac(qR);
        T deltaSL = sL - unL;
        Vec4<T> deltaSL_dotL = sL_dotL - unL_dotL;
        Vec4<T> deltaSL_dotR = sL_dotR;

        T deltaSR = sR - unR;
        Vec4<T> deltaSR_dotL = sR_dotL;
        Vec4<T> deltaSR_dotR = sR_dotR - unR_dotR;

        //T s_star = ((pR - pL)*n_mag2 + qL[0]*unL*deltaSL - qR[0]*unR*deltaSR)/(qL[0]*deltaSL - qR[0]*deltaSR);

        T den = qL[0]*deltaSL - qR[0]*deltaSR;
        Vec4<T> den_dotL = qL[0]*deltaSL_dotL - qR[0]*deltaSR_dotL;
        den_dotL[0] += deltaSL;

        Vec4<T> den_dotR = qL[0]*deltaSL_dotR - qR[0]*deltaSR_dotR;
        den_dotR[0] += -deltaSR;

        //T s_star = ((pR - pL)*n_mag2 + qL[0]*unL*deltaSL - qR[0]*unR*deltaSR)/den;
        T num = (pR - pL)*n_mag2 + qL[0]*unL*deltaSL - qR[0]*unR*deltaSR;
        Vec4<T> num_dotL = -pL_dotL*n_mag2 + qL[0]*unL_dotL*deltaSL + qL[0]*unL*deltaSL_dotL - qR[0]*unR*deltaSR_dotL;
        num_dotL[0] += unL*deltaSL;

        Vec4<T> num_dotR = pR_dotR*n_mag2 + qL[0]*unL*deltaSL_dotR - qR[0]*unR_dotR*deltaSR - qR[0]*unR*deltaSR_dotR;
        num_dotR[0] += -unR*deltaSR;

        T s_star = num/den;
        Vec4<T> s_star_dotL = num_dotL/den - num*den_dotL/(den*den);
        Vec4<T> s_star_dotR = num_dotR/den - num*den_dotR/(den*den);

        //std::cout << "s_star = " << s_star << std::endl;

        auto [fL, fL_dotL] = compute_euler_flux_jac(qL, normal);
        auto [fR, fR_dotR] = compute_euler_flux_jac(qR, normal);

        PRINTVAR(s_star);
        PRINTVAR(num_dotR[0]);
        PRINTVAR(den_dotR[0]);
        if (s_star > 0)
        {
          Vec4<T> d{0, normal[0]/n_mag, normal[1]/n_mag, s_star/n_mag};
          T pressure_term = sL*(pL*n_mag2 + qL[0]*(sL - unL)*(s_star - unL))/n_mag;

          Vec4<T> pressure_term_dotL = sL*(pL_dotL*n_mag2 + qL[0]*(sL_dotL - unL_dotL)*(s_star - unL) + qL[0]*(sL - unL)*(s_star_dotL - unL_dotL));
          pressure_term_dotL[0] += sL*(sL - unL)*(s_star - unL);
          pressure_term_dotL += sL_dotL*(pL*n_mag2 + qL[0]*(sL - unL)*(s_star - unL));
          pressure_term_dotL /= n_mag;

          Vec4<T> pressure_term_dotR = sL*(qL[0]*sL_dotR*(s_star - unL) + qL[0]*(sL - unL)*(s_star_dotR)) +
                                       sL_dotR*(pL*n_mag2 + qL[0]*(sL - unL)*(s_star - unL));
          pressure_term_dotR /= n_mag;

          PRINTVAR(s_star_dotR[0]);
          PRINTVAR(pressure_term_dotR[0]);
          for (UInt i=0; i < 4; ++i)
          {
            T num = s_star*(sL*qL[i] - fL[i]) + pressure_term*d[i];
            T den = (sL - s_star);

            if (i == 0)
            {
              PRINTVAR(num);
              PRINTVAR(den);
            }
            flux[i] = num/den;
            for (UInt j=0; j < 4; ++j)
            {
              T num_dotL = s_star*(sL_dotL[j]*qL[i] + del(i, j)*sL - fL_dotL(i, j)) + s_star_dotL[j]*(sL*qL[i] - fL[i])
                           + pressure_term_dotL[j]*d[i] + pressure_term*del(i, 3)*s_star_dotL[j];
              T den_dotL = sL_dotL[j] - s_star_dotL[j];
              flux_dotL(i, j) = num_dotL/den - num*den_dotL/(den*den);

              T num_dotR = s_star*(sL_dotR[j]*qL[i]) + s_star_dotR[j]*(sL*qL[i] - fL[i]) +
                           pressure_term_dotR[j]*d[i]  + pressure_term*del(i, 3)*s_star_dotR[j]/n_mag;
              T den_dotR = sL_dotR[j] - s_star_dotR[j];

              if (i == 0 && j == 0)
              {
                PRINTVAR(num_dotR);
                PRINTVAR(den_dotR);
              }
              flux_dotR(i, j) = num_dotR/den - num*den_dotR/(den*den);
            }
          }
        } else
        {
          Vec4<T> d{0, normal[0]/n_mag, normal[1]/n_mag, s_star/n_mag};
          T pressure_term = sR*(pR*n_mag2 + qR[0]*(sR - unR)*(s_star - unR))/n_mag;

          Vec4<T> pressure_term_dotL = sR*(qR[0]*sR_dotL*s_star + qR[0]*(sR - unR)*s_star_dotL) +
                                       sR_dotL*(pR*n_mag2 + qR[0]*(sR - unR)*(s_star - unR));
          pressure_term_dotL /= n_mag;

          Vec4<T> pressure_term_dotR = sR*(pR_dotR*n_mag2 + qR[0]*(sR_dotR - unR_dotR)*(s_star - unR) + qR[0]*(sR - unR)*(s_star_dotR - unR_dotR));
          pressure_term_dotR[0] += sR*(sR - unR)*(s_star -  unR);
          pressure_term_dotR    += sR_dotR*(pR*n_mag2 + qR[0]*(sR - unR)*(s_star - unR));
          pressure_term_dotR    /= n_mag;

          for (UInt i=0; i < 4; ++i)
          {
            T num = s_star*(sR*qR[i] - fR[i]) + pressure_term*d[i];
            T den = sR - s_star;
            flux[i] = num/den;

            if (i == 0)
            {
              PRINTVAR(num);
              PRINTVAR(den);
            }            
            for (UInt j=0; j < 4; ++j)
            {
              T num_dotL = s_star*(sR_dotL[j]*qR[i]) + s_star_dotL[j]*(sR*qR[i] - fR[i]) + 
                           pressure_term_dotL[j]*d[i] + pressure_term*del(i, 3)*s_star_dotL[j]/n_mag;
              T den_dotL = sR_dotL[j] - s_star_dotL[j];
              flux_dotL(i, j) = num_dotL/den - num*den_dotL/(den*den);

              T num_dotR = s_star*(sR_dotR[j]*qR[i] + del(i, j)*sR - fR_dotR(i, j)) + s_star_dotR[j]*(sR*qR[i] - fR[i]) +
                           pressure_term_dotR[j]*d[i] + del(i, 3)*pressure_term*s_star_dotR[j]/n_mag;
              T den_dotR = sR_dotR[j] - s_star_dotR[j];

              if (i == 0 && j == 0)
              {
                PRINTVAR(num_dotR);
                PRINTVAR(den_dotR);
              }

              flux_dotR(i, j) = num_dotR/den - num*den_dotR/(den*den);
            }
          }        
        }
      }

      return flux;
    }   
};

}
}

#endif