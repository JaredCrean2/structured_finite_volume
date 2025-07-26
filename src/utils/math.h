#ifndef UTILS_MATH_H
#define UTILS_MATH_H


#include "project_defs.h"
#include "vec.h"

namespace structured_fv {
/*
  template <typename T>
  using Vec3 = FixedVec<T, 3>;

  template <typename T>
  using Vec2 = FixedVec<T, 2>;
*/
  static constexpr double PI = 3.141592653589793238462643383279502884;
  //using structured_fv::Real;
}

namespace std {

template <typename T>
using EnableIfIntegral = std::enable_if_t<std::is_integral_v<T>, bool>;

template <typename IntType, EnableIfIntegral<IntType> = true>
constexpr structured_fv::Complex operator+(const structured_fv::Complex& a, IntType b)
{
  return structured_fv::Complex(a.real() + b, a.imag());
}

template <typename IntType, EnableIfIntegral<IntType> = true>
constexpr structured_fv::Complex operator+(IntType b, const structured_fv::Complex& a)
{
  return a + b;
}

template <typename IntType, EnableIfIntegral<IntType> = true>
constexpr structured_fv::Complex operator-(const structured_fv::Complex& a, IntType b)
{
  return structured_fv::Complex(a.real() - b, a.imag());
}

template <typename IntType, EnableIfIntegral<IntType> = true>
constexpr structured_fv::Complex operator-(IntType b, const structured_fv::Complex& a)
{
  return structured_fv::Complex(b - a.real(), -a.imag());
}

template <typename IntType, EnableIfIntegral<IntType> = true>
constexpr structured_fv::Complex operator*(const structured_fv::Complex& a, IntType b)
{
  return structured_fv::Complex(b*a.real(), b*a.imag());
}

template <typename IntType, EnableIfIntegral<IntType> = true>
constexpr structured_fv::Complex operator*(IntType b, const structured_fv::Complex& a)
{
  return a*b;
}

template <typename IntType, EnableIfIntegral<IntType> = true>
constexpr structured_fv::Complex operator/(const structured_fv::Complex& a, IntType b)
{
  return a/static_cast<structured_fv::Real>(b);
}

template <typename IntType, EnableIfIntegral<IntType> = true>
constexpr structured_fv::Complex operator/(IntType b, const structured_fv::Complex& a)
{
  return static_cast<structured_fv::Real>(b)/a;
}


// these definitions make Complex numbers work as dual numbers
constexpr bool operator<(const structured_fv::Complex& lhs, const structured_fv::Complex& rhs)
{
  return lhs.real() < rhs.real();
}

constexpr bool operator<=(const structured_fv::Complex& lhs, const structured_fv::Complex& rhs)
{
  return lhs.real() <= rhs.real();
}

constexpr bool operator>(const structured_fv::Complex& lhs, const structured_fv::Complex& rhs)
{
  return lhs.real() > rhs.real();
}

constexpr bool operator>=(const structured_fv::Complex& lhs, const structured_fv::Complex& rhs)
{
  return lhs.real() >= rhs.real();
}

constexpr structured_fv::Complex max(const structured_fv::Complex& lhs, const structured_fv::Complex& rhs)
{
  return lhs > rhs ? lhs : rhs;
}

constexpr structured_fv::Complex min(const structured_fv::Complex& lhs, const structured_fv::Complex& rhs)
{
  return lhs < rhs ? lhs : rhs;
}

}


namespace structured_fv {

template <typename T,  UInt N>
constexpr FixedVec<T, N> operator+(const FixedVec<T, N>& a, const FixedVec<T, N>& b)
{
  FixedVec<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] + b[i];

  return c;
}

template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N> operator+(const FixedVec<T, N>& a, const T2& b)
{
  FixedVec<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] + b;

  return c;
}

template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N> operator+(const T2& a, const FixedVec<T, N>& b)
{
  return b + a;
}

template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N>& operator+=(FixedVec<T, N>& a, const FixedVec<T2, N>& b)
{
  for (int i=0; i < N; ++i)
    a[i] += b[i];

  return a;
}

template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N>& operator+=(FixedVec<T, N>& a, const T2& b)
{
  for (int i=0; i < N; ++i)
    a[i] += b;

  return a;
}

template <typename T,  UInt N>
constexpr FixedVec<T, N> operator-(const FixedVec<T, N>& a, const FixedVec<T, N>& b)
{
  FixedVec<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] - b[i];

  return c;
}

template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N> operator-(const FixedVec<T, N>& a, const T2& b)
{
  FixedVec<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] - b;

  return c;
}

template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N> operator-(const T2& a, const FixedVec<T, N>& b)
{
  FixedVec<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a - b[i];

  return c;
}

template <typename T,  UInt N>
constexpr FixedVec<T, N> operator-(const FixedVec<T, N>& a)
{
  FixedVec<T, N> b;
  for (int i=0; i < N; ++i)
    b[i] = -a[i];

  return b;
}

template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N>& operator-=(FixedVec<T, N>& a, const FixedVec<T2, N>& b)
{
  for (int i=0; i < N; ++i)
    a[i] -= b[i];

  return a;
}

template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N>& operator-=(FixedVec<T, N>& a, const T2& b)
{
  for (int i=0; i < N; ++i)
    a[i] -= b;

  return a;
}

template <typename T,  UInt N>
constexpr FixedVec<T, N> operator*(const FixedVec<T, N>& a, const FixedVec<T, N>& b)
{
  FixedVec<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] * b[i];

  return c;
}

template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N> operator*(const FixedVec<T, N>& a, const T2& b)
{
  FixedVec<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] * b;

  return c;
}

template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N> operator*(const T2& b, const FixedVec<T, N>& a)
{
  FixedVec<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = b * a[i];

  return c;
}

template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N>& operator*=(FixedVec<T, N>& a, const FixedVec<T2, N>& b)
{
  for (int i=0; i < N; ++i)
    a[i] *= b[i];

  return a;
}

template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N>& operator*=(FixedVec<T, N>& a, const T2& b)
{
  for (int i=0; i < N; ++i)
    a[i] *= b;

  return a;
}

template <typename T,  UInt N>
constexpr FixedVec<T, N> operator/(const FixedVec<T, N>& a, const FixedVec<T, N>& b)
{
  FixedVec<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] / b[i];

  return c;
}


template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N> operator/(const FixedVec<T, N>& a, const T2& b)
{
  FixedVec<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] / b;

  return c;
}

template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N> operator/(const T2& a, const FixedVec<T, N>& b)
{
  FixedVec<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a / b[i];

  return c;
}

template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N>& operator/=(FixedVec<T, N>& a, const FixedVec<T2, N>& b)
{
  for (int i=0; i < N; ++i)
    a[i] /= b[i];

  return a;
}

template <typename T,  UInt N, typename T2>
constexpr FixedVec<T, N>& operator/=(FixedVec<T, N>& a, const T2& b)
{
  for (int i=0; i < N; ++i)
    a[i] /= b;

  return a;
}

template <typename T,  UInt N>
constexpr T dot(const FixedVec<T, N>& a, const FixedVec<T, N>& b)
{
  T val = 0;
  for (int i=0; i < N; ++i)
    val += a[i] * b[i];

  return val;
}

/*
template <typename T,  UInt N>
std::ostream& operator<<(std::ostream& os, const FixedVec<T, N>& a)
{
  for (UInt i=0; i < N; ++i)
  {
    os << a[i];
    if (i < N - 1)
      os << ", ";
  }

  return os;
}
*/

template <typename T>
constexpr Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b)
{
  T c1 =   a[1] * b[2] - a[2] * b[1];
  T c2 = -(a[0] * b[2] - a[2] * b[0]);
  T c3 =   a[0] * b[1] - a[1] * b[0];


  return {c1, c2, c3};
}


template <typename T, typename T2, UInt N>
constexpr FixedVec<T, N> convert(const FixedVec<T2, N>& arr)
{
  FixedVec<T, N> arr2{};
  for (UInt i=0; i < N; ++i)
    arr2[i] = arr[i];

  return arr2;
}

constexpr Real degreesToRadians(Real degrees)
{
  return PI * (degrees / 180);
}

constexpr Real radiansToDegrees(Real radians)
{
  return radians * 180/ PI;
}

//TODO: finish this
template <typename T>
constexpr T smoothAbs(T x, Real delta=1e-6)
{
  //Real val1 = std::abs(x);
  T minus_x = -x;
  T abs_x = x > 0 ? x : minus_x;
  T val2 = x*x/delta;
  return abs_x > delta ? abs_x : val2;
}

template <typename T>
constexpr std::pair<T, T> smoothAbs_dot(T x, Real x_dot, Real delta=1e-6)
{
  T minus_x = -x;
  T minus_x_dot = -x_dot;

  T abs_x = x > 0 ? x : minus_x;
  T abs_x_dot = x > 0 ? x_dot : minus_x_dot;

  T val2 = x*x/delta;
  T val2_dot = 2*x*x_dot/delta;

  T v =  abs_x > delta ? abs_x : val2;
  T v_dot = abs_x > delta ? abs_x_dot : val2_dot;

  return {v, v_dot};
}

// coords should be in couter-clockwise order
constexpr Real computeQuadArea(const FixedVec<FixedVec<Real, 2>, 4>& coords)
{
  Vec3<Real> a1{0, 0, 0}, a2{0, 0, 0}, b1{0, 0, 0}, b2{0, 0, 0};

  for (int d=0; d < 2; ++d)
  {
    a1[d] = coords[1][d] - coords[0][d];
    a2[d] = coords[3][d] - coords[0][d];
    b1[d] = coords[3][d] - coords[2][d];
    b2[d] = coords[1][d] - coords[2][d];
  }

  return 0.5*(cross(a1, a2)[2] + cross(b1, b2)[2]);
}

template <typename T>
constexpr Int sgn(T val)
{
  return (T(0) < val) - (val < T(0));
}

inline Real imag(Real)
{
  return 0;
}

inline Real imag(Complex x)
{
  return x.imag();
}

inline Real real(Real x)
{
  return x;
}

inline Real real(Complex x)
{
  return x.real();
}

}


#endif