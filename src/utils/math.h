#ifndef UTILS_MATH_H
#define UTILS_MATH_H

#include <array>
#include <ostream>
#include "project_defs.h"

namespace {
  template <typename T>
  using Vec3 = std::array<T, 3>;

  template <typename T>
  using Vec2 = std::array<T, 2>;

  static constexpr double PI = 3.141592653589793238462643383279502884;
  using structured_fv::Real;
}

namespace std {

template <typename T,  size_t N>
constexpr std::array<T, N> operator+(const std::array<T, N>& a, const std::array<T, N>& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] + b[i];

  return c;
}

template <typename T,  size_t N, typename T2>
constexpr std::array<T, N> operator+(const std::array<T, N>& a, const T2& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] + b;

  return c;
}

template <typename T,  size_t N, typename T2>
constexpr std::array<T, N> operator+(const T2& a, const std::array<T, N>& b)
{
  return b + a;
}

template <typename T,  size_t N, typename T2>
constexpr std::array<T, N>& operator+=(std::array<T, N>& a, const std::array<T2, N>& b)
{
  for (int i=0; i < N; ++i)
    a[i] += b[i];

  return a;
}

template <typename T,  size_t N, typename T2>
constexpr std::array<T, N>& operator+=(std::array<T, N>& a, const T2& b)
{
  for (int i=0; i < N; ++i)
    a[i] += b;

  return a;
}

template <typename T,  size_t N>
constexpr std::array<T, N> operator-(const std::array<T, N>& a, const std::array<T, N>& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] - b[i];

  return c;
}

template <typename T,  size_t N, typename T2>
constexpr std::array<T, N> operator-(const std::array<T, N>& a, const T2& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] - b;

  return c;
}

template <typename T,  size_t N, typename T2>
constexpr std::array<T, N> operator-(const T2& a, const std::array<T, N>& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a - b[i];

  return c;
}

template <typename T,  size_t N, typename T2>
constexpr std::array<T, N>& operator-=(std::array<T, N>& a, const std::array<T2, N>& b)
{
  for (int i=0; i < N; ++i)
    a[i] -= b[i];

  return a;
}

template <typename T,  size_t N, typename T2>
constexpr std::array<T, N>& operator-=(std::array<T, N>& a, const T2& b)
{
  for (int i=0; i < N; ++i)
    a[i] -= b;

  return a;
}

template <typename T,  size_t N>
constexpr std::array<T, N> operator*(const std::array<T, N>& a, const std::array<T, N>& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] * b[i];

  return c;
}

template <typename T,  size_t N, typename T2>
constexpr std::array<T, N> operator*(const std::array<T, N>& a, const T2& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] * b;

  return c;
}

template <typename T,  size_t N, typename T2>
constexpr std::array<T, N> operator*(const T2& b, const std::array<T, N>& a)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = b * a[i];

  return c;
}

template <typename T,  size_t N, typename T2>
constexpr std::array<T, N>& operator*=(std::array<T, N>& a, const std::array<T2, N>& b)
{
  for (int i=0; i < N; ++i)
    a[i] *= b[i];

  return a;
}

template <typename T,  size_t N, typename T2>
constexpr std::array<T, N>& operator*=(std::array<T, N>& a, const T2& b)
{
  for (int i=0; i < N; ++i)
    a[i] *= b;

  return a;
}

template <typename T,  size_t N>
constexpr std::array<T, N> operator/(const std::array<T, N>& a, const std::array<T, N>& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] / b[i];

  return c;
}


template <typename T,  size_t N, typename T2>
constexpr std::array<T, N> operator/(const std::array<T, N>& a, const T2& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a[i] / b;

  return c;
}

template <typename T,  size_t N, typename T2>
constexpr std::array<T, N> operator/(const T2& a, const std::array<T, N>& b)
{
  std::array<T, N> c;
  for (int i=0; i < N; ++i)
    c[i] = a / b[i];

  return c;
}

template <typename T,  size_t N, typename T2>
constexpr std::array<T, N>& operator/=(std::array<T, N>& a, const std::array<T2, N>& b)
{
  for (int i=0; i < N; ++i)
    a[i] /= b[i];

  return a;
}

template <typename T,  size_t N, typename T2>
constexpr std::array<T, N>& operator/=(std::array<T, N>& a, const T2& b)
{
  for (int i=0; i < N; ++i)
    a[i] /= b;

  return a;
}

template <typename T,  size_t N>
constexpr T dot(const std::array<T, N>& a, const std::array<T, N>& b)
{
  T val = 0;
  for (int i=0; i < N; ++i)
    val += a[i] * b[i];

  return val;
}

template <typename T,  size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& a)
{
  for (size_t i=0; i < N; ++i)
  {
    os << a[i];
    if (i < N - 1)
      os << ", ";
  }

  return os;
}

template <typename T>
constexpr Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b)
{
  T c1 =   a[1] * b[2] - a[2] * b[1];
  T c2 = -(a[0] * b[2] - a[2] * b[0]);
  T c3 =   a[0] * b[1] - a[1] * b[0];


  return {c1, c2, c3};
}

}

template <typename T, typename T2, size_t N>
constexpr std::array<T, N> convert(const std::array<T2, N>& arr)
{
  std::array<T, N> arr2{};
  for (size_t i=0; i < N; ++i)
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
constexpr Real computeQuadArea(const std::array<std::array<Real, 2>, 4>& coords)
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


#endif