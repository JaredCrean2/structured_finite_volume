#include "vec.h"
#include "project_defs.h"
#include "math.h"

namespace structured_fv {

// a Dual number class that implements forward mode algorithmic differentiation with
// N dual part
template <typename T, UInt N>
class Dual
{
  public:
    using value_type = T;
    static constexpr UInt size = N;

    constexpr Dual() :
      m_dx{},
      m_val{}
    {}

    constexpr Dual(NoInit)
    {}

    constexpr Dual(const T& val) :
      m_dx{},
      m_val(val)
    {}

    constexpr Dual(const T& val, const FixedVec<T, N>& dx) :
      m_dx(dx),
      m_val(val)
    {}

    template <typename T2>
    constexpr Dual(const Dual<T2, N>& b) :
      m_val(b.m_val)
    {
      for (UInt i=0; i < N; ++i)
        m_dx[i] = b.m_dx[i];
    }

    constexpr const T& get_value() const { return m_val; }

    constexpr T& get_value() { return m_val; }

    constexpr const T& get_deriv(UInt i) const { return m_dx[i]; }

    constexpr T& get_deriv(UInt i) { return m_dx[i]; }

    constexpr const T& operator()() const { return get_value(); }

    constexpr T& operator()() { return get_value(); }

    constexpr const T& operator()(UInt i) const { return get_deriv(i); }

    constexpr T& operator()(UInt i) { return get_deriv(i); }

    constexpr Dual operator-() const
    {
      Dual<T, N> b(NoInit{});
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = -m_dx[i];
      b.m_val = -m_val;

      return b;
    }    

    // dual-dual operations
    template <typename T2>
    constexpr Dual& operator+=(const Dual<T2, N>& b)
    {
      m_dx  += b.m_dx;
      m_val += b.m_val;
      return *this;
    }   

    template <typename T2>
    constexpr Dual& operator-=(const Dual<T2, N>& b)
    {
      m_dx  -= b.m_dx;
      m_val -= b.m_val;
      return *this;
    }
      

    template <typename T2>
    constexpr Dual& operator*=(const Dual<T2, N>& b)
    {
      for (UInt i=0; i < N; ++i)
        m_dx[i] = m_dx[i]*b.m_val + m_val*b.m_dx[i];

      m_val *= b.m_val;
      return *this;
    }
      

    template <typename T2>
    constexpr Dual& operator/=(const Dual<T2, N>& b)
    {
      for (UInt i=0; i < N; ++i)
        m_dx[i] = m_dx[i]/b.m_val - (m_val/(b.m_val*b.m_val))*b.m_dx[i];

      m_val /= b.m_val;
      return *this;
    }

    // dual-non-dual operations
    template <typename T2>
    constexpr Dual& operator+=(const T2& b)
    {
      m_val += b;
      return *this;
    } 
    
    template <typename T2>
    constexpr Dual& operator-=(const T2& b)
    {
      m_val -= b;
      return *this;
    }
    
    template <typename T2>
    constexpr Dual& operator*=(const T2& b)
    {
      m_dx *= b;
      m_val *= b;
      return *this;
    }     

    template <typename T2>
    constexpr Dual& operator/=(const T2& b)
    {
      m_dx /= b;
      m_val /= b;
      return *this;
    }
    
    // overloads for functions in cmath.h
    
    Dual sqrt() const
    {
      Dual<T, N> b(NoInit{});
      T sqrt_v = std::sqrt(m_val);
      T deriv_fac = 1.0/(2*sqrt_v);
      b.m_val = sqrt_v;
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = deriv_fac*m_dx[i];

      return b;
    }

    Dual cbrt() const
    {
      Dual<T, N> b(NoInit{});
      T sqrt_v = std::cbrt(m_val);
      T deriv_fac = std::pow(m_val, -2.0/3)/3;
      b.m_val = sqrt_v;
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = deriv_fac*m_dx[i];

      return b;
    }

    Dual exp() const
    {
      Dual<T, N> b(NoInit{});
      T exp_v = std::exp(m_val);
      b.m_val = exp_v;
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = exp_v*m_dx[i];

      return b;
    }
  
    Dual exp2() const
    {
      Dual<T, N> b(NoInit{});
      T exp_v = std::exp2(m_val);
      T log2_v = std::log(2);  //TODO: can be constexpr in c++26
      T deriv_fac = exp_v * log2_v;
      b.m_val = exp_v;
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = deriv_fac*m_dx[i];

      return b;
    }

    Dual expm1() const
    {
      Dual<T, N> b(NoInit{});
      T exp_v = std::expm1(m_val);
      T deriv_fac = std::exp(m_val);
      b.m_val = exp_v;
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = deriv_fac*m_dx[i];

      return b;
    }    

    Dual log() const
    {
      Dual<T, N> b(NoInit{});
      b.m_val = std::log(m_val);
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = m_dx[i]/m_val;

      return b;
    }
    
    Dual log10() const
    {
      Dual<T, N> b(NoInit{});
      b.m_val = std::log10(m_val);
      T deriv_fac = std::log(10)*m_val;  //TODO: can be constexpr in c++26
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = m_dx[i]/deriv_fac;

      return b;
    } 
    
    Dual log2() const
    {
      Dual<T, N> b(NoInit{});
      b.m_val = std::log2(m_val);
      T deriv_fac = std::log(2)*m_val;  //TODO: can be constexpr in c++26
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = m_dx[i]/deriv_fac;

      return b;
    }
    
    Dual log1p() const
    {
      Dual<T, N> b(NoInit{});
      b.m_val = std::log1p(m_val);
      T deriv_fac = 1 + m_val;
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = m_dx[i]/deriv_fac;

      return b;
    } 
   
    template <typename T2>
    Dual pow(const Dual<T2, N>& b) const
    {
      Dual<T, N> c(NoInit{});
      c.m_val = std::pow(m_val, b.m_val);
      T deriv_fac1 = b.m_val + std::pow(m_val, b()-1);
      T deriv_fac2 = c.m_val * std::log(m_val);
      for (UInt i=0; i < N; ++i)
        c.m_dx[i] = deriv_fac1*m_dx[i] + deriv_fac2*b.m_dy[i];

      return c;
    }

    Dual sin() const
    {
      Dual<T, N> b(NoInit{});
      b.m_val = std::sin(m_val);
      T deriv_fac = std::cos(m_val);
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = deriv_fac*m_dx[i];

      return b;
    } 
    
    Dual cos() const
    {
      Dual<T, N> b(NoInit{});
      b.m_val = std::cos(m_val);
      T deriv_fac = -std::sin(m_val);
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = deriv_fac*m_dx[i];

      return b;
    }
    
    Dual tan() const
    {
      Dual<T, N> b(NoInit{});
      b.m_val = std::tan(m_val);
      T deriv_fac = std::cos(m_val);
      deriv_fac = deriv_fac*deriv_fac;
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = m_dx[i]/deriv_fac;

      return b;
    } 
    
    Dual asin() const
    {
      Dual<T, N> b(NoInit{});
      b.m_val = std::asin(m_val);
      T deriv_fac = std::sqrt(1 - m_val*m_val);
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = m_dx[i]/deriv_fac;

      return b;
    }

    Dual acos() const
    {
      Dual<T, N> b(NoInit{});
      b.m_val = std::acos(m_val);
      T deriv_fac = -std::sqrt(1 - m_val*m_val);
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = m_dx[i]/deriv_fac;

      return b;
    }
    
    Dual atan() const
    {
      Dual<T, N> b(NoInit{});
      b.m_val = std::atan(m_val);
      T deriv_fac = 1 + m_val*m_val;
      for (UInt i=0; i < N; ++i)
        b.m_dx[i] = m_dx[i]/deriv_fac;

      return b;
    }    

  private:
    FixedVec<T, N> m_dx;
    T m_val;

    // Only need to friend Dual<T2, N>, but the standard
    // doesnt allow that
    template <typename T2, UInt N2>
    friend class Dual;
};

namespace impl {

// get common element type, where T1 and T2 can be either regular numbers
// or dual numbers
template <typename T1, typename T2>
struct common_type
{
  using type = std::common_type_t<T1, T2>;
};

template <typename T1, typename T2, UInt N>
struct common_type<T1, Dual<T2, N>>
{
  using type = std::common_type_t<T1, typename Dual<T2, N>::value_type>;
};

template <typename T1, typename T2, UInt N>
struct common_type<Dual<T1, N>, T2>
{
  using type = std::common_type_t<T1, typename Dual<T2, N>::value_type>;
};

template <typename T1, typename T2, UInt N>
struct common_type<Dual<T1, N>, Dual<T2, N>>
{
  using type = std::common_type_t<typename Dual<T2, N>::value_type, typename Dual<T2, N>::value_type>;;
};

template <typename T1, typename T2>
using common_type_t = typename common_type<T1, T2>::type;


// if a type T is a dual number
template <typename T>
struct is_dual : std::false_type
{};

template <typename T, UInt N>
struct is_dual<Dual<T, N>> : std::true_type
{};

template <typename T>
constexpr bool is_dual_v = is_dual<T>::value;

template <typename T1, typename T2>
using EnableIfEitherIsDual = std::enable_if_t<is_dual_v<T1> || is_dual_v<T2>, bool>;

template <typename T1>
using EnableIfNotDual = std::enable_if_t<!is_dual_v<T1>, bool>;

template <typename T1, typename T2>
constexpr UInt get_dual_length()
{
  if constexpr (is_dual_v<T1> && is_dual_v<T2>)
  {
    static_assert(T1::size == T2::size, "dual numbers must have same length N");
    return T1::size;
  }  else if constexpr (is_dual_v<T1>)
  {
    return T1::size;
  } else if constexpr (is_dual_v<T2>)
  {
    return T2::size;
  } else
  {
    static_assert(false, "at least one of T1 and T2 must be a Dual<T, N>");
  }
}

}

template <typename T, UInt N>
T real(const Dual<T, N>& a)
{
  return a();
}

// non-member operations where at least one of the operands is a Dual
template <typename T1, typename T2, impl::EnableIfEitherIsDual<T1, T2> = true>
constexpr auto operator+(const T1& a, const T2& b)
{
  using T = impl::common_type_t<T1, T2>;
  constexpr UInt N = impl::get_dual_length<T1, T2>();
  Dual<T, N> c(a);
  c += b;
  return c;
}

template <typename T1, typename T2, impl::EnableIfEitherIsDual<T1, T2> = true>
constexpr auto operator-(const T1& a, const T2& b)
{
  using T = impl::common_type_t<T1, T2>;
  constexpr UInt N = impl::get_dual_length<T1, T2>();
  Dual<T, N> c(a);
  c -= b;
  return c;
}

template <typename T1, typename T2, impl::EnableIfEitherIsDual<T1, T2> = true>
constexpr auto operator*(const T1& a, const T2& b)
{
  using T = impl::common_type_t<T1, T2>;
  constexpr UInt N = impl::get_dual_length<T1, T2>();
  Dual<T, N> c(a);
  c *= b;
  return c;
}

template <typename T1, typename T2, impl::EnableIfEitherIsDual<T1, T2> = true>
constexpr auto operator/(const T1& a, const T2& b)
{
  using T = impl::common_type_t<T1, T2>;
  constexpr UInt N = impl::get_dual_length<T1, T2>();
  Dual<T, N> c(a);
  c /= b;
  return c;
}

template <typename T1, typename T2, UInt N>
constexpr auto pow(const Dual<T1, N>& a, const Dual<T2, N>& b)
{
  using T = impl::common_type_t<T1, T2>;
  Dual<T, N> c(NoInit{});

  c() = std::pow(a(), b());
  T deriv_fac1 = b() * std::pow(a(), b()-1);
  T deriv_fac2 = c() * std::log(a());
  for (UInt i=0; i < N; ++i)
    c(i) = deriv_fac1*a(i) + deriv_fac2*b(i);

  return c;
}

template <typename T1, typename T2, UInt N, impl::EnableIfNotDual<T2> = true>
constexpr auto pow(const Dual<T1, N>& a, const T2& b)
{
  using T = impl::common_type_t<T1, T2>;
  Dual<T, N> c(NoInit{});

  c() = std::pow(a(), b);
  T deriv_fac1 = b * std::pow(a(), b-1);
  for (UInt i=0; i < N; ++i)
    c(i) = deriv_fac1*a(i);

  return c;
}

template <typename T1, typename T2, UInt N, impl::EnableIfNotDual<T1> = true>
constexpr auto pow(const T1& a, const Dual<T2, N>& b)
{
  using T = impl::common_type_t<T1, T2>;
  Dual<T, N> c(NoInit{});

  c() = std::pow(a, b());
  T deriv_fac2 = c() * std::log(a);
  for (UInt i=0; i < N; ++i)
    c(i) = deriv_fac2*b(i);

  return c;
}

template <typename T1, typename T2, impl::EnableIfEitherIsDual<T1, T2> = true>
constexpr auto operator<(const T1& a, const T2& b)
{
  return real(a) < real(b);
}

template <typename T1, typename T2, impl::EnableIfEitherIsDual<T1, T2> = true>
constexpr auto operator<=(const T1& a, const T2& b)
{
  return real(a) <= real(b);
}

template <typename T1, typename T2, impl::EnableIfEitherIsDual<T1, T2> = true>
constexpr auto operator>(const T1& a, const T2& b)
{
  return real(a) > real(b);
}

template <typename T1, typename T2, impl::EnableIfEitherIsDual<T1, T2> = true>
constexpr auto operator>=(const T1& a, const T2& b)
{
  return real(a) >= real(b);
}


//TODO: make these nonmembers
template <typename T, UInt N>
Dual<T, N> sqrt(const Dual<T, N>& a)
{
  return a.sqrt();
}

template <typename T, UInt N>
Dual<T, N> cbrt(const Dual<T, N>& a)
{
  return a.cbrt();
}

template <typename T, UInt N>
Dual<T, N> exp(const Dual<T, N>& a)
{
  return a.exp();
}

template <typename T, UInt N>
Dual<T, N> exp2(const Dual<T, N>& a)
{
  return a.exp2();
}

template <typename T, UInt N>
Dual<T, N> expm1(const Dual<T, N>& a)
{
  return a.expm1();
}

template <typename T, UInt N>
Dual<T, N> log(const Dual<T, N>& a)
{
  return a.log();
}

template <typename T, UInt N>
Dual<T, N> log10(const Dual<T, N>& a)
{
  return a.log10();
}

template <typename T, UInt N>
Dual<T, N> log2(const Dual<T, N>& a)
{
  return a.log2();
}

template <typename T, UInt N>
Dual<T, N> log1p(const Dual<T, N>& a)
{
  return a.log1p();
}

template <typename T, UInt N>
Dual<T, N> sin(const Dual<T, N>& a)
{
  return a.sin();
}

template <typename T, UInt N>
Dual<T, N> cos(const Dual<T, N>& a)
{
  return a.cos();
}

template <typename T, UInt N>
Dual<T, N> tan(const Dual<T, N>& a)
{
  return a.tan();
}

template <typename T, UInt N>
Dual<T, N> asin(const Dual<T, N>& a)
{
  return a.asin();
}

template <typename T, UInt N>
Dual<T, N> acos(const Dual<T, N>& a)
{
  return a.acos();
}

template <typename T, UInt N>
Dual<T, N> atan(const Dual<T, N>& a)
{
  return a.atan();
}
}