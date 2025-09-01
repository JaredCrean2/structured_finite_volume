#ifndef STRUCTURED_FINITE_VOLUME_UTILS_VEC_H
#define STRUCTURED_FINITE_VOLUME_UTILS_VEC_H

#include "project_defs.h"
#include "utils/error_handling.h"
#include <algorithm>
#include <iostream>

namespace structured_fv {

// this is a replacement for std::array that will work on device eventually
template <typename T, UInt N>
struct FixedVec
{
  using value_type = T;
  using size_type = UInt;
  using difference_type = std::ptrdiff_t;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;
  using iterator = T*;
  using const_iterator = const T*;

  // constructors use aggregate initialization
  constexpr reference at(UInt pos)
  {
    assertAlways(pos < N, "bounds error");
    return _m_vals[pos];
  }

  constexpr const_reference at(UInt pos) const
  {
    assertAlways(pos < N, "bounds error");
    return _m_vals[pos];
  }

  constexpr reference operator[](UInt pos) 
  {
    assert(pos < N);
    return _m_vals[pos];
  }

  constexpr const_reference operator[](UInt pos) const
  {
    assert(pos < N);
    return _m_vals[pos];
  }

  constexpr reference front()
  {
    static_assert(N > 0);
    return _m_vals[0];
  }

  constexpr const_reference front() const
  {
    static_assert(N > 0);
    return _m_vals[0];
  }

  constexpr reference back()
  {
    static_assert(N > 0);
    return _m_vals[N-1];
  }

  constexpr const_reference back() const
  {
    static_assert(N > 0);
    return _m_vals[N-1];
  }

  constexpr pointer data() noexcept { return _m_vals; }

  constexpr const_pointer data() const noexcept { return _m_vals; }

  constexpr iterator begin() noexcept { return _m_vals; }

  constexpr const_iterator begin() const noexcept { return _m_vals; }

  constexpr const_iterator cbegin() const noexcept { return _m_vals; }

  constexpr iterator end() noexcept { return begin() + N; }

  constexpr const_iterator end() const noexcept { return begin() + N; }

  constexpr const_iterator cend() const noexcept { return begin() + N; }

  constexpr bool empty() const noexcept { return N == 0; }

  constexpr size_type size() const noexcept { return N; }

  constexpr size_type max_size() const noexcept { return N; }

  constexpr void fill(const T& val) noexcept
  {
    for (UInt i=0; i < N; ++i)
      _m_vals[i] = val;
  }

  constexpr void swap(FixedVec<T, N>& other) noexcept
  {
    for (UInt i=0; i < N; ++i)
    {
      T tmp = _m_vals[i];
      _m_vals[i] = other[i];
      other[i] = tmp;
    }
  }

  T _m_vals[N];  // this can't be private if the class is an aggregate, but no one
                 // should access it directly
};

template <typename T>
using Vec2 = FixedVec<T, 2>;

template <typename T>
using Vec3 = FixedVec<T, 3>;

template <typename T>
using Vec4 = FixedVec<T, 4>;

template <typename T>
using Vec5 = FixedVec<T, 5>;

// std::array is only required to implement comparisons for arrays of the same T and
// N, which is a little weird
template <typename T, UInt N>
bool operator==(const FixedVec<T, N>& lhs, const FixedVec<T, N>& rhs)
{
  return std::equal(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

template <typename T, UInt N>
bool operator!=(const FixedVec<T, N>& lhs, const FixedVec<T, N>& rhs)
{
  return !(lhs == rhs);
}

template <typename T, UInt N>
bool operator<(const FixedVec<T, N>& lhs, const FixedVec<T, N>& rhs)
{
  return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
}

template <typename T, UInt N>
bool operator<=(const FixedVec<T, N>& lhs, const FixedVec<T, N>& rhs)
{
  return (lhs < rhs) || (lhs == rhs);
}

template <typename T, UInt N>
bool operator>(const FixedVec<T, N>& lhs, const FixedVec<T, N>& rhs)
{
  return !(lhs <= rhs);
}

template <typename T, UInt N>
bool operator>=(const FixedVec<T, N>& lhs, const FixedVec<T, N>& rhs)
{
  return !(lhs < rhs);
}

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

template <size_t I, typename T, UInt N>
constexpr auto&& get(const FixedVec<T, N>& arr)
{
  return arr[I];
}

template <size_t I, typename T, UInt N>
constexpr auto&& get(FixedVec<T, N>& arr)
{
  return arr[I];
}

template <size_t I, typename T, UInt N>
constexpr auto&& get(FixedVec<T, N>&& arr)
{
  return arr[I];
}

template <size_t I, typename T, UInt N>
constexpr auto&& get(const FixedVec<T, N>&& arr)
{
  return arr[I];
}


template <typename T, UInt N>
std::ostream& operator<<(std::ostream& os, const FixedVec<T, N>& arr)
{
  for (UInt i=0; i < N; ++i)
  {
    os << arr[i];
    if (i < N-1)
      os << ", ";
  }

  return os;
}

template <typename Tout, typename Tin, UInt N>
FixedVec<Tout, N> copy(const FixedVec<Tin, N>& vin)
{
  FixedVec<Tout, N> vout;
  for (UInt i=0; i < N; ++i)
    vout[i] = vin[i];

  return vout;
}

}

namespace std {
  
template<typename T, structured_fv::UInt N>
struct std::tuple_size<structured_fv::FixedVec<T, N>>
{
  static constexpr std::size_t value = N;
};

template <size_t M, typename T, structured_fv::UInt N>
struct std::tuple_element<M, structured_fv::FixedVec<T, N>>
{
  using type = T;
};

}

#endif