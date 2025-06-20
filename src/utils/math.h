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

constexpr Real smoothAbs(Real x, Real delta)
{
  Real val1 = std::abs(x);
  Real val2 = x*x/delta;
  return val1 > delta ? val1 : val2;
}

constexpr Real smoothAbsDeriv(Real x, Real delta)
{
  double val1 = std::abs(x);
  Real val1_deriv = x > 0 ? 1 : -1;

  //Real val2 = x*x/delta;
  Real val2_deriv = 2*x/delta;

  return val1 > delta ? val1_deriv : val2_deriv;
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

}


#endif