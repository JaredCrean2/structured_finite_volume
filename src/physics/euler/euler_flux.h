#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_FLUX_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_FLUX_H

#include "typedefs.h"
#include "utils/array3.h"
#include "utils/matrix.h"
#include "utils/project_defs.h"

#include <iostream>

namespace structured_fv {
namespace euler {

// the solution vector is [rho, rho u, rho v, E]
// where E is the total energy per unit volume
// Some constituitive relations:
//   - e is the specific internal energy (energy per unit mass)
//   - rho e is the internal energy per unit volume
//   - 0.5(u^2 + v^2) is specific kinetic energy (kinetic energy per unit mass)
//   - 0.5 rho(u^2 + v^2) is kinetic energy per unit volume
//   - E = rho e + 0.5 rho (u^2 + v^2)
//   - e = cv T
//   - h = cp T =  e + p/rho (specific enthaly aka enthalpy per unit mass)
//   - gamma = cp/cv
//   - e = cv T = cv (p/(rho R)) = (cv/R)(p/rho) = p/((gamma-1) rho)
//   - this E = rho e + 0.5 rho (u^2 + v^2) = p/(gamma-1) + 0.5 rho(u^2 + v^2)
// Note that some authors define E to be total energy per unit mass and
// have the 4th component of the solution vector = rho E

template <typename T>
constexpr T compute_pressure(const Vec4<T>& sol)
{
  return Gamma_m1 * (sol[3] - 0.5*(sol[1]*sol[1] + sol[2]*sol[2])/sol[0]);
}

template <typename T>
constexpr ScalarVectorPair<T> compute_pressure_jac(const Vec4<T>& sol)
{
  Real velocity_squared = sol[1]*sol[1] + sol[2]*sol[2];
  T p =  Gamma_m1 * (sol[3] - 0.5*(velocity_squared)/sol[0]);
  Vec4<T> p_jac{};
  p_jac[0] =  Gamma_m1 * 0.5 * velocity_squared/(sol[0]*sol[0]);
  p_jac[1] = -Gamma_m1 * sol[1]/sol[0];
  p_jac[2] = -Gamma_m1 * sol[2]/sol[0];
  p_jac[3] =  Gamma_m1;

  return {p, p_jac};
}
/*
template <typename T>
constexpr std::pair<T, Vec4<T>> compute_pressure_jac(const Vec4<T>& sol)
{

  T p = Gamma_m1 * (sol[3] - 0.5*(sol[1]*sol[1] + sol[2]*sol[2])/sol[0]);
  Vec4<T> p_jac{(0.5*Gamma_m1)*(sol[1]*sol[1] + sol[2]*sol[2])/(sol[0]*sol[0]),
                 -Gamma_m1*sol[1]/sol[0],
                 -Gamma_m1*sol[2]/sol[0],
                  Gamma_m1};

  return {p, p_jac};
}
*/

template <typename T>
constexpr T compute_sos2(const Vec4<T>& sol)
{
  return Gamma * compute_pressure(sol)/sol[0];
}

template <typename T>
constexpr ScalarVectorPair<T> compute_sos2_jac(const Vec4<T>& sol)
{
  auto [p, p_jac] = compute_pressure_jac(sol);
  T sos2 =  Gamma * p/sol[0];

  p_jac[0] = Gamma * p_jac[0]/sol[0] - Gamma* p/(sol[0]*sol[0]);  
  for (UInt i=1; i < sol.size(); ++i)
    p_jac[i] = Gamma * p_jac[i]/sol[0];

  return {sos2, p_jac};
}

template <typename T>
constexpr T compute_sos(const Vec4<T>& sol)
{
  return sqrt(compute_sos2(sol));
}

template <typename T>
constexpr ScalarVectorPair<T> compute_sos_jac(const Vec4<T>& sol)
{
  auto [sos2, sos2_jac] = compute_sos2_jac(sol);
  T sos = sqrt(sos2);
  T fac = 1.0/(2*sos);
  for (UInt i=0; i < sol.size(); ++i)
     sos2_jac[i] *= fac;
    
  return {sos, sos2_jac};
}

template <typename T>
constexpr T compute_momentum_n(const Vec4<T>& sol, const Vec2<Real>& normal)
{
  return (sol[1]*normal[0] + sol[2]*normal[1]);
}

template <typename T>
constexpr ScalarVectorPair<T> compute_momentum_n_jac(const Vec4<T>& sol, const Vec2<Real>& normal)
{
  T rhou_n =  (sol[1]*normal[0] + sol[2]*normal[1]);
  Vec4<T> rhou_n_jac{0, normal[0], normal[1], 0};
  return {rhou_n, rhou_n_jac};
}

template <typename T>
constexpr T compute_un(const Vec4<T>& sol, const Vec2<Real>& normal)
{
  return compute_momentum_n(sol, normal)/sol[0];
}

template <typename T>
constexpr ScalarVectorPair<T> compute_un_jac(const Vec4<T>& sol, const Vec2<Real>& normal)
{
  auto [rhou_n, rhou_n_jac] = compute_momentum_n_jac(sol, normal);
  T u_n = rhou_n/sol[0];

  for (UInt i=0; i < sol.size(); ++i)
    rhou_n_jac[i] /= sol[0];

  rhou_n_jac[0] -= rhou_n/(sol[0]*sol[0]);
  return {u_n, rhou_n_jac};
}

template <typename T>
constexpr T compute_temperature(const Vec4<T>& sol)
{
  T rho_e = sol[3] - 0.5*(sol[1]*sol[1] + sol[2]*sol[2])/sol[0];
  return rho_e/(sol[0]*Cv);
}

template <typename T>
constexpr ScalarVectorPair<T> compute_temperature_jac(const Vec4<T>& sol)
{
  T momentum_squared = sol[1]*sol[1] + sol[2]*sol[2];
  T rho_e = sol[3] - 0.5*(momentum_squared)/sol[0];
  Vec4<T> rho_e_jac{0.5*momentum_squared/(sol[0]*sol[0]),
                   -sol[1]/sol[0],
                   -sol[2]/sol[0],
                   1};

  T temp = rho_e/(sol[0]*Cv);
  for (UInt i=0; i < sol.size(); ++i)
    rho_e_jac[i] /= sol[0]*Cv;
  rho_e_jac[0] -= rho_e/(sol[0]*sol[0]*Cv);

  return {temp, rho_e_jac};
}

template <typename T>
constexpr Vec4<T> compute_euler_flux(const Vec4<T>& sol, const Vec2<Real>& normal)
{
  T pressure = compute_pressure(sol);
  T momentum_n = (sol[1]*normal[0] + sol[2]*normal[1]);
  T velocity_n = momentum_n/sol[0];

  Vec4<T> flux{0, 0, 0, 0};
  flux[0] = momentum_n;
  flux[1] = sol[1]*velocity_n + pressure*normal[0];
  flux[2] = sol[2]*velocity_n + pressure*normal[1];
  flux[3] = (sol[3] + pressure)*velocity_n;
  
  return flux;
}

template <typename T>
constexpr VectorMatrixPair<T>
compute_euler_flux_jac(const Vec4<T>& sol, const Vec2<Real>& normal)
{
  auto [pressure, p_jac] = compute_pressure_jac(sol);
  T momentum_n = (sol[1]*normal[0] + sol[2]*normal[1]);
  T momentum_n_dot1 = normal[0];
  T momentum_n_dot2 = normal[1];

  T velocity_n = momentum_n/sol[0];
  T velocity_n_dot0 = -momentum_n/(sol[0]*sol[0]);
  T velocity_n_dot1 = momentum_n_dot1/sol[0];
  T velocity_n_dot2 = momentum_n_dot2/sol[0];

  Vec4<T> flux{0, 0, 0, 0};
  flux[0] = momentum_n;
  flux[1] = sol[1]*velocity_n + pressure*normal[0];
  flux[2] = sol[2]*velocity_n + pressure*normal[1];
  flux[3] = (sol[3] + pressure)*velocity_n;

  Matrix<T, 4, 4> flux_jac;
  flux_jac(0, 0) = 0;
  flux_jac(0, 1) = momentum_n_dot1;
  flux_jac(0, 2) = momentum_n_dot2;
  flux_jac(0, 3) = 0;

  flux_jac(1, 0) = sol[1]*velocity_n_dot0 + p_jac[0]*normal[0];
  flux_jac(1, 1) = velocity_n + sol[1]*velocity_n_dot1 + p_jac[1]*normal[0];
  flux_jac(1, 2) = sol[1]*velocity_n_dot2 + p_jac[2]*normal[0];
  flux_jac(1, 3) = p_jac[3]*normal[0];

  flux_jac(2, 0) = sol[2]*velocity_n_dot0 + p_jac[0]*normal[1];
  flux_jac(2, 1) = sol[2]*velocity_n_dot1 + p_jac[1]*normal[1];
  flux_jac(2, 2) = velocity_n + sol[2]*velocity_n_dot2 + p_jac[2]*normal[1];
  flux_jac(2, 3) = p_jac[3]*normal[1];

  flux_jac(3, 0) = (sol[3] + pressure)*velocity_n_dot0 + p_jac[0]*velocity_n;
  flux_jac(3, 1) = (sol[3] + pressure)*velocity_n_dot1 + p_jac[1]*velocity_n;
  flux_jac(3, 2) = (sol[3] + pressure)*velocity_n_dot2 + p_jac[2]*velocity_n;
  flux_jac(3, 3) = (1 + p_jac[3])*velocity_n;
  
  return {flux, flux_jac};
}

// Eigenvalue decomposition of df/dq = R Lambda R^-1 in a given normal direction

template <typename T>
constexpr void compute_eigen_decomp(const Vec4<T>& sol, const Vec2<Real>& normal,
                                    Matrix<T, 4>& R, Vec4<T>& lambda, Matrix<T, 4>& Rinv)
{
  Real n_mag = sqrt(dot(normal, normal));
  Vec3<Real> n_unit{normal[0]/n_mag, normal[1]/n_mag, 0};
  T u = sol[1]/sol[0];
  T v = sol[2]/sol[0];
  T U_squared_half = (u*u + v*v)/2.0;  //TODO: define op(Complex, int)
  T sos2 = compute_sos2(sol);
  T sos = sqrt(sos2);
  T H = U_squared_half + sos2/Gamma_m1;
  T Un = compute_un(sol, normal);


  R(0, 0) = 1;
  R(0, 1) = 1;
  R(0, 2) = 0;
  R(0, 3) = 1;

  R(1, 0) = u - n_unit[0]*sos;
  R(1, 1) = u;
  R(1, 2) = sol[0]*n_unit[1];
  R(1, 3) = u + n_unit[0]*sos;

  R(2, 0) = v - n_unit[1]*sos;
  R(2, 1) = v;
  R(2, 2) = -sol[0]*n_unit[0];
  R(2, 3) = v + n_unit[1]*sos;

  R(3, 0) = H - sos*Un/n_mag;  
  R(3, 1) = U_squared_half;
  R(3, 2) = sol[0]*(u*n_unit[1] - v*n_unit[0]);
  R(3, 3) = H + sos*Un/n_mag;

  lambda[0] = Un - sos*n_mag;
  lambda[1] = Un;
  lambda[2] = Un;
  lambda[3] = Un + sos*n_mag;


  T M2 = 2.0*U_squared_half/sos2;

  Rinv(0, 0) = (Gamma_m1/4)*M2 + Un/(2.0*sos*n_mag);
  Rinv(0, 1) = -(n_unit[0] + Gamma_m1*u/sos)/(2.0*sos);
  Rinv(0, 2) = -(n_unit[1] + Gamma_m1*v/sos)/(2.0*sos);
  Rinv(0, 3) = Gamma_m1/(2.0*sos2);  

  Rinv(1, 0) =  1.0 - 0.5*Gamma_m1*M2;
  Rinv(1, 1) =  Gamma_m1*u/sos2; 
  Rinv(1, 2) =  Gamma_m1*v/sos2;
  Rinv(1, 3) = -Gamma_m1/sos2;

  Rinv(2, 0) =  (v*n_unit[0] - u*n_unit[1])/sol[0];
  Rinv(2, 1) =  n_unit[1]/sol[0];
  Rinv(2, 2) = -n_unit[0]/sol[0];
  Rinv(2, 3) =  0;

  Rinv(3, 0) = (Gamma_m1/4)*M2 - Un/(2.0*sos*n_mag);
  Rinv(3, 1) = (n_unit[0] - Gamma_m1*u/sos)/(2.0*sos);
  Rinv(3, 2) = (n_unit[1] - Gamma_m1*v/sos)/(2.0*sos);
  Rinv(3, 3) = Gamma_m1/(2.0*sos2);
}

// For R_jac and Rinv, the 3rd dimension is the derivative dimension, ie.
// R[:, :, k] is the derivative wrt sol[k]
template <typename T>
constexpr void compute_eigen_decomp_jac(const Vec4<T>& sol, const Vec2<Real>& normal,
                                        Matrix<T, 4>& R, Array3<T, 4>& R_jac,
                                        Vec4<T>& lambda, Matrix<T, 4>& lambda_jac,
                                        Matrix<T, 4>& Rinv, Array3<T, 4>& Rinv_jac)
{  
  Real n_mag = sqrt(dot(normal, normal));
  Vec3<Real> n_unit{normal[0]/n_mag, normal[1]/n_mag, 0};
  T u = sol[1]/sol[0];
  Vec4<T> u_dot{-sol[1]/(sol[0]*sol[0]), 1.0/sol[0], 0, 0};

  T v = sol[2]/sol[0];
  Vec4<T> v_dot{-sol[2]/(sol[0]*sol[0]), 0, 1.0/sol[0], 0};

  T U_squared_half = (u*u + v*v)/2;
  Vec4<T> U_squared_half_dot{u*u_dot[0] + v*v_dot[0], u*u_dot[1], v*v_dot[2], 0};
  
  //T sos2 = compute_sos2(sol);
  auto [sos2, sos2_dot] = compute_sos2_jac(sol);

  T sos = sqrt(sos2);
  Vec4<T> sos_dot = sos2_dot / (2*sos);

  T H = U_squared_half + sos2/Gamma_m1;
  Vec4<T> H_dot = sos2_dot / Gamma_m1 + U_squared_half_dot;

  //T Un = compute_un(sol, normal);
  auto [Un, Un_dot] = compute_un_jac(sol, normal);

  R(0, 0) = 1;
  R(0, 1) = 1;
  R(0, 2) = 0;
  R(0, 3) = 1;
  for (UInt j=0; j < 4; ++j)
    for (UInt k=0; k < 4; ++k)
      R_jac(0, j, k) = 0;
    

  R(1, 0) = u - n_unit[0]*sos;
  R(1, 1) = u;
  R(1, 2) = sol[0]*n_unit[1];
  R(1, 3) = u + n_unit[0]*sos;
  for (UInt k=0; k < 4; ++k)
  {
    R_jac(1, 0, k) = u_dot[k] -n_unit[0]*sos_dot[k];
    R_jac(1, 1, k) = u_dot[k];
    R_jac(1, 2, k) = 0;
    R_jac(1, 3, k) = u_dot[k] + n_unit[0]*sos_dot[k];
  }
  R_jac(1, 2, 0) = n_unit[1];

  R(2, 0) = v - n_unit[1]*sos;
  R(2, 1) = v;
  R(2, 2) = -sol[0]*n_unit[0];
  R(2, 3) = v + n_unit[1]*sos;
  for (UInt k=0; k < 4; ++k)
  {
    R_jac(2, 0, k) = v_dot[k] - n_unit[1]*sos_dot[k];
    R_jac(2, 1, k) = v_dot[k];
    R_jac(2, 2, k) = 0;
    R_jac(2, 3, k) = v_dot[k] + n_unit[1]*sos_dot[k];
  }
  R_jac(2, 2, 0) = -n_unit[0];

  R(3, 0) = H - sos*Un/n_mag;  
  R(3, 1) = U_squared_half;
  R(3, 2) = sol[0]*(u*n_unit[1] - v*n_unit[0]);
  R(3, 3) = H + sos*Un/n_mag;

  for (UInt k=0; k < 4; ++k)
  {
    T term2 = (sos_dot[k]*Un + sos*Un_dot[k])/n_mag;;
    R_jac(3, 0, k) = H_dot[k] - term2;
    R_jac(3, 1, k) = U_squared_half_dot[k];
    R_jac(3, 2, k) = sol[0]*(u_dot[k]*n_unit[1] - v_dot[k]*n_unit[0]);
    R_jac(3, 3, k) = H_dot[k] + term2;
  }
  R_jac(3, 2, 0) += (u*n_unit[1] - v*n_unit[0]);

  lambda[0] = Un - sos*n_mag;
  lambda[1] = Un;
  lambda[2] = Un;
  lambda[3] = Un + sos*n_mag;

  for (UInt k=0; k < 4; ++k)
  {
    lambda_jac(0, k) = Un_dot[k] - sos_dot[k]*n_mag;
    lambda_jac(1, k) = Un_dot[k];
    lambda_jac(2, k) = Un_dot[k];
    lambda_jac(3, k) = Un_dot[k] + sos_dot[k]*n_mag;
  }


  T M2 = 2*U_squared_half/sos2;
  Vec4<T> M2_dot;
  for (UInt i=0; i < 4; ++i)
    M2_dot[i] = 2*(U_squared_half_dot[i]/sos2 - U_squared_half*sos2_dot[i]/(sos2*sos2));

  Rinv(0, 0) = (Gamma_m1/4)*M2 + Un/(2*sos*n_mag);
  Rinv(0, 1) = -(n_unit[0] + Gamma_m1*u/sos)/(2*sos);
  Rinv(0, 2) = -(n_unit[1] + Gamma_m1*v/sos)/(2*sos);
  Rinv(0, 3) = Gamma_m1/(2*sos2);

  for (UInt k=0; k < 4; ++k)
  {
    Rinv_jac(0, 0, k) = (Gamma_m1/4)*M2_dot[k] + Un_dot[k]/(2*sos*n_mag) - Un*sos_dot[k]/(2*sos2*n_mag);

    T t1     = -(n_unit[0] + Gamma_m1*u/sos);
    T t1_dot = -(Gamma_m1*(u_dot[k]/sos - u*sos_dot[k]/sos2));

    T t2     = -(n_unit[1] + Gamma_m1*v/sos);
    T t2_dot = -(Gamma_m1*(v_dot[k]/sos - v*sos_dot[k]/sos2));

    Rinv_jac(0, 1, k) = ( t1_dot/sos - t1*sos_dot[k]/sos2)/2;
    Rinv_jac(0, 2, k) = ( t2_dot/sos - t2*sos_dot[k]/sos2)/2;
    Rinv_jac(0, 3, k) = -Gamma_m1*sos2_dot[k]/(2*sos2*sos2);
  }

  Rinv(1, 0) =  1 - 0.5*Gamma_m1*M2;
  Rinv(1, 1) =  Gamma_m1*u/sos2;
  Rinv(1, 2) =  Gamma_m1*v/sos2;
  Rinv(1, 3) = -Gamma_m1/sos2;

  for (UInt k=0; k < 4; ++k)
  {
    Rinv_jac(1, 0, k) = -0.5*Gamma_m1*M2_dot[k];
    Rinv_jac(1, 1, k) = (Gamma_m1*u_dot[k])/sos2 - (Gamma_m1*u/(sos2*sos2))*sos2_dot[k];
    Rinv_jac(1, 2, k) = (Gamma_m1*v_dot[k])/sos2 - (Gamma_m1*v/(sos2*sos2))*sos2_dot[k];
    Rinv_jac(1, 3, k) = (Gamma_m1/(sos2*sos2))*sos2_dot[k];
  }

  Rinv(2, 0) =  (v*n_unit[0] - u*n_unit[1])/sol[0];
  Rinv(2, 1) =  n_unit[1]/sol[0];
  Rinv(2, 2) = -n_unit[0]/sol[0];
  Rinv(2, 3) =  0;

  for (UInt k=0; k < 4; ++k)
  {
    Rinv_jac(2, 0, k) = (v_dot[k]*n_unit[0] - u_dot[k]*n_unit[1])/sol[0];
    Rinv_jac(2, 1, k) = 0;
    Rinv_jac(2, 2, k) = 0;
    Rinv_jac(2, 3, k) = 0;
  }
  Rinv_jac(2, 0, 0) += -(v*n_unit[0] - u*n_unit[1])/(sol[0]*sol[0]);
  Rinv_jac(2, 1, 0) += -n_unit[1]/(sol[0]*sol[0]);
  Rinv_jac(2, 2, 0) +=  n_unit[0]/(sol[0]*sol[0]);

  Rinv(3, 0) = (Gamma_m1/4)*M2 - Un/(2*sos*n_mag);
  Rinv(3, 1) = (n_unit[0] - Gamma_m1*u/sos)/(2*sos);
  Rinv(3, 2) = (n_unit[1] - Gamma_m1*v/sos)/(2*sos);
  Rinv(3, 3) = Gamma_m1/(2*sos2);

  for (UInt k=0; k < 4; ++k)
  {
    Rinv_jac(3, 0, k) = (Gamma_m1/4)*M2_dot[k] - (Un_dot[k]/(2*sos*n_mag) - (Un/(2*sos2*n_mag))*sos_dot[k]);

    T t1 = n_unit[0] - Gamma_m1*u/sos;
    T t1_dot = -Gamma_m1*u_dot[k]/sos + (Gamma_m1*u*sos_dot[k]/sos2);

    T t2 = n_unit[1] - Gamma_m1*v/sos;
    T t2_dot = -Gamma_m1*v_dot[k]/sos + Gamma_m1*v*sos_dot[k]/sos2;

    Rinv_jac(3, 1, k) = t1_dot/(2*sos) - t1*sos_dot[k]/(2*sos2);
    Rinv_jac(3, 2, k) = t2_dot/(2*sos) - t2*sos_dot[k]/(2*sos2);
    Rinv_jac(3, 3, k) = -Gamma_m1*sos2_dot[k]/(2*sos2*sos2);
  }
}

struct PrimitiveVarTag {};
struct ConservativeVarTag {};

template <typename T>
constexpr Vec4<T> compute_conservative_variables(const Vec4<T>& primitive_vars, PrimitiveVarTag = PrimitiveVarTag{})
{
  T e = Cv*primitive_vars[3];
  T E = primitive_vars[0]*(e + 0.5*(primitive_vars[1]*primitive_vars[1] + primitive_vars[2]*primitive_vars[2]));

  return {primitive_vars[0], primitive_vars[0]*primitive_vars[1], primitive_vars[0]*primitive_vars[2], E};
}

template <typename T>
constexpr VectorMatrixPair<T> compute_conservative_variables_jac(const Vec4<T>& primitive_vars, PrimitiveVarTag = PrimitiveVarTag{})
{
  T e = Cv*primitive_vars[3];
  T total_energy = e + 0.5*(primitive_vars[1]*primitive_vars[1] + primitive_vars[2]*primitive_vars[2]);
  T E = primitive_vars[0]*total_energy;

  Vec4<T> cons{primitive_vars[0], primitive_vars[0]*primitive_vars[1], primitive_vars[0]*primitive_vars[2], E};
  Matrix<T, 4, 4> dqdp({1, 0, 0, 0,
                       primitive_vars[1], primitive_vars[0], 0, 0,
                       primitive_vars[2], 0, primitive_vars[0], 0,
                       total_energy, primitive_vars[0]*primitive_vars[1], primitive_vars[0]*primitive_vars[2], primitive_vars[0]*Cv});
  return {cons, dqdp};
}

template <typename T>
constexpr Vec4<T> compute_primitive_variables(const Vec4<T>& conservative_vars, ConservativeVarTag = ConservativeVarTag{})
{
  T u = conservative_vars[1]/conservative_vars[0];
  T v = conservative_vars[2]/conservative_vars[0];
  T e = conservative_vars[3]/conservative_vars[0] - 0.5*(u*u + v*v);
  T temp = e/Cv;

  return {conservative_vars[0], u, v, temp};
}

template <typename T>
constexpr VectorMatrixPair<T> compute_primitive_variables_jac(const Vec4<T>& conservative_vars, ConservativeVarTag = ConservativeVarTag{})
{
  T rho = conservative_vars[0];

  T u = conservative_vars[1]/rho;
  T du_drho = -conservative_vars[1]/(rho*rho);
  T du_drhou = 1.0/rho;

  T v = conservative_vars[2]/rho;
  T dv_drho = -conservative_vars[2]/(rho*rho);
  T dv_drhov = 1.0/rho;

  T e = conservative_vars[3]/rho - 0.5*(u*u + v*v);
  T de_drho = -conservative_vars[3]/(rho*rho) - u*du_drho - v*dv_drho ;
  T de_drhou = -u*du_drhou;
  T de_drhov = -v*dv_drhov;
  T de_dE = 1.0/rho;

  T temp = e/Cv;

  Vec4<T> prim{rho, u, v, temp};

  Matrix<T, 4, 4> dpdq({1, 0, 0, 0,
                        du_drho, du_drhou, 0, 0,
                        dv_drho, 0, dv_drhov, 0,
                        de_drho/Cv, de_drhou/Cv, de_drhov/Cv, de_dE/Cv});

  return {prim, dpdq};
}

// For now, rotate solution into the x direction and rotate the flux
// back to the wall normal direction.
// Eventually work out the eigen decomposition of the flux jacobian in the
// wall normal direction 
/*
class Rotator
{
  public:
    constexpr Rotator(const Vec2<Real>& normal) :
      m_normal(normal),
      m_x_normal{sqrt(normal[0]*normal[0] + normal[1]*normal[1]), 0}
    {}

    constexpr Vec4<Real> rotateForward(const Vec4<Real>& sol) const
    {
      Real mag = m_x_normal[0];
      return {sol[0], (m_normal[0]*sol[1] + m_normal[1]*sol[1])/mag,
              (-m_normal[1]*sol[1] + m_normal[0]*sol[2])/mag, sol[3]};
    }

    constexpr Vec4<Real> rotateBack(const Vec4<Real>& flux) const
    {
      Real mag = m_x_normal[0];
      return {flux[0], (m_normal[0]*flux[1] - m_normal[1]*flux[2])/mag,
              (m_normal[1]*flux[1] + m_normal[0]*flux[0])/mag, flux[3]};
    }

    constexpr const Vec2<Real>& getRotatedNormal() const { return m_x_normal; }

  private:
    const Vec2<Real> m_normal;
    const Vec2<Real> m_x_normal;
};
*/




}
}

#endif