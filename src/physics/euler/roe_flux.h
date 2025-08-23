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
                                 const Vec2<Real>& normal,
                                 Matrix<T, 4>& flux_dotL, Matrix<T, 4>& flux_dotR) const
    {
      Matrix<T, 4> avg_state_dotL, avg_state_dotR;
      RoeAvgState<T> avg_state = compute_roe_avg_jac(qL, qR, avg_state_dotL, avg_state_dotR);
      //Vec4<T> q_avg = compute_conservative_variables(avg_state, RoeStateTag());
      auto [q_avg, q_avg_dot] = compute_conservative_variables_jac(avg_state, RoeStateTag());
      //TODO could re-use the storage from avg_state_dotL,R?
      Matrix<T, 4> q_avg_dotL = q_avg_dot * avg_state_dotL;
      Matrix<T, 4> q_avg_dotR = q_avg_dot * avg_state_dotR;


      //TODO: maybe compute derivative of everything wrt q_avg, then only compute
      //      df/dq_avg * dq_avg/dqL at the last step

      Matrix<T, 4, 4> R, Rinv, lambdas_dot;
      Vec4<T> lambdas{};
      Array3<T, 4> R_dot, Rinv_dot;
      //TODO: it would be better of we could compute the eigendecomp
      //      directly from the Roe state
      compute_eigen_decomp_jac(q_avg, normal, R, R_dot, lambdas, lambdas_dot, Rinv, Rinv_dot);

      Array3<T, 4> R_dotL, R_dotR, Rinv_dotL, Rinv_dotR;
      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
          for (UInt k=0; k < 4; ++k)
            for (UInt p=0; p < 4; ++p)
            {
              R_dotL(i, j, p) += R_dot(i, j, k) * q_avg_dotL(k, p);
              R_dotR(i, j, p) += R_dot(i, j, k) * q_avg_dotR(k, p);
              Rinv_dotL(i, j, p) += Rinv_dot(i, j, k) * q_avg_dotL(k, p);
              Rinv_dotR(i, j, p) += Rinv_dot(i, j, k) * q_avg_dotR(k, p);
            }

      Matrix<T, 4> lambdas_dotL, lambdas_dotR;
      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
          for (UInt k=0; k < 4; ++k)
          {
            lambdas_dotL(i, j) += lambdas_dot(i, k) * q_avg_dotL(k, j);
            lambdas_dotR(i, j) += lambdas_dot(i, k) * q_avg_dotR(k, j);
          }
            
        

      Vec4<T> tmp{};
      Matrix<T, 4> tmp_dotL, tmp_dotR;
      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
        {
          tmp[i] += Rinv(i, j)*(qR[j] - qL[j]);

          for (UInt k=0; k < 4; ++k)
          {
            //TODO: if we only need sum_k Rinv_dotL(i, j, k), we could save some stack space
            tmp_dotL(i, k) += Rinv_dotL(i, j, k)*(qR[j] - qL[j]) - Rinv(i, j)*del(j, k);
            tmp_dotR(i, k) += Rinv_dotR(i, j, k)*(qR[j] - qL[j]) + Rinv(i, j)*del(j, k);
          }
        }

      for (UInt i=0; i < 4; ++i)
      {
        //T lambda = smoothAbs(lambdas[i]);
        auto [lambda, lambda_dot] = smoothAbs_dot(lambdas[i], 1);
        if (lambda < m_delta)
        {
          lambda_dot = 2*lambda*lambda_dot/(2*m_delta);
          lambda = (lambda*lambda + m_delta*m_delta)/(2*m_delta);
        }

        for (UInt j=0; j < 4; ++j)
        {
          Real lambda_dotL = lambda_dot * lambdas_dotL(i, j);
          Real lambda_dotR = lambda_dot * lambdas_dotR(i, j);
          tmp_dotL(i, j) = tmp_dotL(i, j)*lambda + tmp[i]*lambda_dotL;
          tmp_dotR(i, j) = tmp_dotR(i, j)*lambda + tmp[i]*lambda_dotR;
        }

        tmp[i] *= lambda;
      }
        

      Vec4<T> flux{};
      flux_dotL = 0;
      flux_dotR = 0;
      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
        {
          flux[i] += R(i, j) * tmp[j];
          for (UInt k=0; k < 4; ++k)
          {
            // again, only need sum_k R_dotL(i, j, k)
            flux_dotL(i, k) += R(i, j) * tmp_dotL(j, k) + R_dotL(i, j, k) * tmp[j];
            flux_dotR(i, k) += R(i, j) * tmp_dotR(j, k) + R_dotR(i, j, k) * tmp[j];
          }
        }

      auto [fL, fL_dotL] = compute_euler_flux_jac(qL, normal);
      auto [fR, fR_dotR] = compute_euler_flux_jac(qR, normal);

      for (UInt i=0; i < 4; ++i)
      {
        flux[i] = 0.5*(fL[i] + fR[i] - flux[i]);
        for (UInt j=0; j < 4; ++j)
        {
          flux_dotL(i, j) = 0.5*(fL_dotL(i, j) - flux_dotL(i, j));
          flux_dotR(i, j) = 0.5*(fR_dotR(i, j) - flux_dotR(i, j));
        }
      }
        
      return flux;
    }    

  private:
    Real m_delta;
};

}
}

#endif