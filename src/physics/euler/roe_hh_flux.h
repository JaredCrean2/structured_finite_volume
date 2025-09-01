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
          T n_mag = sqrt(dot(normal, normal));
          lambda_kL += sgn*n_mag*compute_sos(q_kL);
          lambda_kR += sgn*n_mag*compute_sos(q_kR);
        }

        if (lambda_kL < 0 && lambda_kR > 0)
        {
          T beta = (lambda_kR - lambdas[i])/(lambda_kR - lambda_kL);
          T beta_c = 1 - beta;
          for (UInt j=0; j < 4; ++j)
          {
            flux[j] += (smoothAbs(lambda_kL)*beta  *alphas[i])*R(j, i);
            flux[j] += (smoothAbs(lambda_kR)*beta_c*alphas[i])*R(j, i);
          }
        } else
        {        
          for (UInt j=0; j < 4; ++j)
            flux[j] += (smoothAbs(lambdas[i])*alphas[i])*R(j, i);
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
                      
      Vec4<T> alphas{};
      Matrix<T, 4> alphas_dotL, alphas_dotR;
      for (UInt i=0; i < 4; ++i)
        for (UInt j=0; j < 4; ++j)
        {
          alphas[i] += Rinv(i, j)*(qR[j] - qL[j]);

          for (UInt k=0; k < 4; ++k)
          {
            alphas_dotL(i, k) += Rinv_dotL(i, j, k)*(qR[j] - qL[j]);
            alphas_dotR(i, k) += Rinv_dotR(i, j, k)*(qR[j] - qL[j]);
          }

          alphas_dotL(i, j) -= Rinv(i, j);
          alphas_dotR(i, j) += Rinv(i, j);
        }
        
      Vec4<T> q_kL = qL;
      Vec4<T> q_kR;
      Matrix<T, 4> q_kL_dotL;
      Matrix<T, 4> q_kL_dotR;
      Matrix<T, 4> q_kR_dotL;
      Matrix<T, 4> q_kR_dotR;
      for (UInt i=0; i < 4; ++i)
        q_kL_dotL(i, i) = 1;

      Vec4<T> flux{};
      flux_dotL = 0;
      flux_dotR = 0;

      
      //TODO: getting R^T rather than R would be more efficient
      for (UInt i=0; i < 4; ++i)
      {
        for (UInt j=0; j < 4; ++j)
        {
          q_kR[j] = q_kL[j] + alphas[i]*R(j, i);
          for (UInt k=0; k < 4; ++k)
          {
            q_kR_dotL(j, k) = q_kL_dotL(j, k) + alphas[i]*R_dotL(j, i, k) + alphas_dotL(i, k)*R(j, i);
            q_kR_dotR(j, k) = q_kL_dotR(j, k) + alphas[i]*R_dotR(j, i, k) + alphas_dotR(i, k)*R(j, i);
          }
        }

        auto [lambda_kL, lambda_kL_dot] = compute_un_jac(q_kL, normal);
        auto [lambda_kR, lambda_kR_dot] = compute_un_jac(q_kR, normal);
        Vec4<T> lambda_kL_dotL{}, lambda_kL_dotR{};
        Vec4<T> lambda_kR_dotL{}, lambda_kR_dotR{};
        for (UInt j=0; j < 4; ++j)
          for (UInt k=0; k < 4; ++k)
          {
            lambda_kL_dotL[j] += lambda_kL_dot[k]*q_kL_dotL(k, j);
            lambda_kL_dotR[j] += lambda_kL_dot[k]*q_kL_dotR(k, j);
            lambda_kR_dotL[j] += lambda_kR_dot[k]*q_kR_dotL(k, j);
            lambda_kR_dotR[j] += lambda_kR_dot[k]*q_kR_dotR(k, j);
          }


        if (i == 0 || i == 3)
        {
          Int sgn = i == 0 ? -1 : 1;
          T n_mag = sqrt(dot(normal, normal));
          auto [aL, aL_dot] = compute_sos_jac(q_kL);
          auto [aR, aR_dot] = compute_sos_jac(q_kR);
          lambda_kL += sgn*n_mag*aL;
          lambda_kR += sgn*n_mag*aR;
          for (UInt j=0; j < 4; ++j)
          {
            T aL_dotLj = 0, aL_dotRj = 0;
            T aR_dotLj = 0, aR_dotRj = 0;
            for (UInt k=0; k < 4; ++k)
            {
              aL_dotLj += aL_dot[k]*q_kL_dotL(k, j);
              aL_dotRj += aL_dot[k]*q_kL_dotR(k, j);
              aR_dotLj += aR_dot[k]*q_kR_dotL(k, j);
              aR_dotRj += aR_dot[k]*q_kR_dotR(k, j);
            }

            lambda_kL_dotL[j] += sgn*n_mag*aL_dotLj;
            lambda_kL_dotR[j] += sgn*n_mag*aL_dotRj;
            lambda_kR_dotL[j] += sgn*n_mag*aR_dotLj;
            lambda_kR_dotR[j] += sgn*n_mag*aR_dotRj;
          }
        }

        if (lambda_kL < 0 && lambda_kR > 0)
        {
          T beta_num = (lambda_kR - lambdas[i]);
          T beta_den = (lambda_kR - lambda_kL);

          Vec4<T> beta_num_dotL{}, beta_num_dotR{};
          Vec4<T> beta_den_dotL{}, beta_den_dotR{};
          for (UInt j=0; j < 4; ++j)
          {
            beta_num_dotL[j] = lambda_kR_dotL[j] - lambdas_dotL(i, j);
            beta_num_dotR[j] = lambda_kR_dotR[j] - lambdas_dotR(i, j);
            beta_den_dotL[j] = lambda_kR_dotL[j] - lambda_kL_dotL[j];
            beta_den_dotR[j] = lambda_kR_dotR[j] - lambda_kL_dotR[j];
          }

          T beta = beta_num/beta_den;
          T beta_c = 1 - beta;

          Vec4<T> beta_dotL = beta_num_dotL/beta_den - beta_num*beta_den_dotL/(beta_den*beta_den);
          Vec4<T> beta_dotR = beta_num_dotR/beta_den - beta_num*beta_den_dotR/(beta_den*beta_den);
          
          for (UInt j=0; j < 4; ++j)
          {
            auto [abs_lambda_kL, abs_lambda_kL_dot] = smoothAbs_dot(lambda_kL, 1);
            auto [abs_lambda_kR, abs_lambda_kR_dot] = smoothAbs_dot(lambda_kR, 1);
            Vec4<T> abs_lambda_kL_dotL = abs_lambda_kL_dot * lambda_kL_dotL;
            Vec4<T> abs_lambda_kL_dotR = abs_lambda_kL_dot * lambda_kL_dotR;
            Vec4<T> abs_lambda_kR_dotL = abs_lambda_kR_dot * lambda_kR_dotL;
            Vec4<T> abs_lambda_kR_dotR = abs_lambda_kR_dot * lambda_kR_dotR;

            flux[j] += abs_lambda_kL*beta  *alphas[i]*R(j, i);
            flux[j] += abs_lambda_kR*beta_c*alphas[i]*R(j, i);

            
            for (UInt k=0; k < 4; ++k)
            {              
              flux_dotL(j, k) += abs_lambda_kL_dotL[k]*beta*alphas[i]*R(j, i) +
                                 abs_lambda_kL*beta_dotL[k]*alphas[i]*R(j, i) +
                                 abs_lambda_kL*beta*alphas_dotL(i, k)*R(j, i) +
                                 abs_lambda_kL*beta*alphas[i]*R_dotL(j, i, k);

              flux_dotR(j, k) += abs_lambda_kL_dotR[k]*beta*alphas[i]*R(j, i) +
                                 abs_lambda_kL*beta_dotR[k]*alphas[i]*R(j, i) +
                                 abs_lambda_kL*beta*alphas_dotR(i, k)*R(j, i) +
                                 abs_lambda_kL*beta*alphas[i]*R_dotR(j, i, k);

              // beta_c_dot = -beta_dot


              flux_dotL(j, k) += abs_lambda_kR_dotL[k]*beta_c*alphas[i]*R(j, i) +
                                 abs_lambda_kR*(-beta_dotL[k])*alphas[i]*R(j, i) +
                                 abs_lambda_kR*beta_c*alphas_dotL(i, k)*R(j, i) +              
                                 abs_lambda_kR*beta_c*alphas[i]*R_dotL(j, i, k);

              flux_dotR(j, k) += abs_lambda_kR_dotR[k]*beta_c*alphas[i]*R(j, i) +
                                 abs_lambda_kR*(-beta_dotR[k])*alphas[i]*R(j, i) +
                                 abs_lambda_kR*beta_c*alphas_dotR(i, k)*R(j, i) +
                                 abs_lambda_kR*beta_c*alphas[i]*R_dotR(j, i, k); 
                                 
            } 
          }
        } else
        {
          auto [abs_lambda, abs_lambda_dot] = smoothAbs_dot(lambdas[i], 1);
          for (UInt j=0; j < 4; ++j)
          {
            flux[j] += abs_lambda*alphas[i]*R(j, i);

            for (UInt k=0; k < 4; ++k)
            {
              T abs_lambda_dotLk = abs_lambda_dot * lambdas_dotL(i, k);
              T abs_lambda_dotRk = abs_lambda_dot * lambdas_dotR(i, k);

              flux_dotL(j, k) += abs_lambda_dotLk*alphas[i]*R(j, i) +
                                 abs_lambda*alphas_dotL(i, k)*R(j, i) +
                                 abs_lambda*alphas[i]*R_dotL(j, i, k);

              flux_dotR(j, k) += abs_lambda_dotRk*alphas[i]*R(j, i) +
                                 abs_lambda*alphas_dotR(i, k)*R(j, i) +
                                 abs_lambda*alphas[i]*R_dotR(j, i, k);
            }
          }
          
        }
        q_kL = q_kR;
        q_kL_dotL = q_kR_dotL;
        q_kL_dotR = q_kR_dotR;  // TODO: do some pointer swaps here instead of copy?
      }
      
      auto [fL, fL_dot] = compute_euler_flux_jac(qL, normal);
      auto [fR, fR_dot] = compute_euler_flux_jac(qR, normal);

      for (UInt i=0; i < 4; ++i)
      {
        flux[i] = 0.5*(fL[i] + fR[i] - flux[i]);

        for (UInt j=0; j < 4; ++j)
        {
          flux_dotL(i, j) = 0.5*(fL_dot(i, j) - flux_dotL(i, j));
          flux_dotR(i, j) = 0.5*(fR_dot(i, j) - flux_dotR(i, j));
        }
      }
    
      return flux; 
    }
      
};

}
}

#endif