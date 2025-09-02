#include "gtest/gtest.h"
#include "utils/dual_number.h"
#include "utils/project_defs.h"
#include <type_traits>

namespace {
using structured_fv::Real;
using structured_fv::UInt;


using structured_fv::FixedVec;
using structured_fv::Complex;
using structured_fv::Dual;
using RealDual2 = structured_fv::Dual<double, 2>;
using FloatDual2 = structured_fv::Dual<float, 2>;

template <typename T>
class DualNumberBinaryFunctionTester : public testing::Test
{};

template <typename T>
class DualNumberOutOfPlaceFunctionTester : public testing::Test
{};

template <typename T>
class DualNumberUnaryFunctionTester : public testing::Test
{};

struct Plus
{
  template <typename T, typename T2>
  auto operator()(T a, T2 b)
  {
    return a + b;
  }

  static constexpr const char* name = "Plus";
};

struct Minus
{
  template <typename T, typename T2>
  auto operator()(T a, T2 b)
  {
    return a - b;
  }

  static constexpr const char* name = "Minus";
};

struct Times
{
  template <typename T, typename T2>
  auto operator()(T a, T2 b)
  {
    return a * b;
  }

  static constexpr const char* name = "Times";
};

struct Div
{
  template <typename T, typename T2>
  auto operator()(T a, T2 b)
  {
    return a / b;
  }

  static constexpr const char* name = "Div";
};

struct Pow
{
  template <typename T, typename T2>
  auto operator()(T a, T2 b)
  {
    return structured_fv::pow(a, b);
  }

  static constexpr const char* name = "pow";
};

struct PlusEqual
{
  template <typename T, typename T2>
  auto operator()(T a, T2 b)
  {
    return a += b;
  }

  static constexpr const char* name = "PlusEqual";
};

struct MinusEqual
{
  template <typename T, typename T2>
  auto operator()(T a, T2 b)
  {
    return a -= b;
  }

  static constexpr const char* name = "MinusEqual";
};

struct TimesEqual
{
  template <typename T, typename T2>
  auto operator()(T a, T2 b)
  {
    return a *= b;
  }

  static constexpr const char* name = "TimesEqual";
};

struct DivEqual
{
  template <typename T, typename T2>
  auto operator()(T a, T2 b)
  {
    return a /= b;
  }

  static constexpr const char* name = "DivEqual";
};

// Unary functions

struct UnaryMinus
{
  template <typename T>
  T operator()(T a)
  {
    return -a;
  }

  static constexpr const char* name = "UnaryPlus";
  static constexpr bool HasComplexSupport = true;
};

struct Sqrt
{
  template <typename T>
  T operator()(T a)
  {
    return structured_fv::sqrt(a);
  }

  static constexpr const char* name = "Sqrt";
  static constexpr bool HasComplexSupport = true;
};

struct Cbrt
{
  template <typename T>
  T operator()(T a)
  {
    return structured_fv::cbrt(a);
  }

  static constexpr const char* name = "cbrt";
  static constexpr bool HasComplexSupport = false;
};

struct Exp
{
  template <typename T>
  T operator()(T a)
  {
    return structured_fv::exp(a);
  }

  static constexpr const char* name = "Exp";
  static constexpr bool HasComplexSupport = true;
};

struct Exp2
{
  template <typename T>
  T operator()(T a)
  {
    return structured_fv::exp2(a);
  }

  static constexpr const char* name = "Exp2";
  static constexpr bool HasComplexSupport = false;
};

struct Expm1
{
  template <typename T>
  T operator()(T a)
  {
    return structured_fv::expm1(a);
  }

  static constexpr const char* name = "Expm1";
  static constexpr bool HasComplexSupport = false;
};

struct NaturalLog
{
  template <typename T>
  T operator()(T a)
  {
    return structured_fv::log(a);
  }

  static constexpr const char* name = "Ln";
  static constexpr bool HasComplexSupport = true;
};

struct Log10
{
  template <typename T>
  T operator()(T a)
  {
    return structured_fv::log10(a);
  }

  static constexpr const char* name = "log10";
  static constexpr bool HasComplexSupport = true;
};

struct Log2
{
  template <typename T>
  T operator()(T a)
  {
    return structured_fv::log2(a);
  }

  static constexpr const char* name = "log2";
  static constexpr bool HasComplexSupport = false;
};

struct Log1p
{
  template <typename T>
  T operator()(T a)
  {
    return structured_fv::log1p(a);
  }

  static constexpr const char* name = "log1p";
  static constexpr bool HasComplexSupport = false;
};

struct Sin
{
  template <typename T>
  T operator()(T a)
  {
    return structured_fv::sin(a);
  }

  static constexpr const char* name = "sin";
  static constexpr bool HasComplexSupport = true;
};

struct Cos
{
  template <typename T>
  T operator()(T a)
  {
    return structured_fv::cos(a);
  }

  static constexpr const char* name = "cos";
  static constexpr bool HasComplexSupport = true;
};

struct Tan
{
  template <typename T>
  T operator()(T a)
  {
    return structured_fv::tan(a);
  }

  static constexpr const char* name = "tan";
  static constexpr bool HasComplexSupport = true;
};

struct Asin
{
  template <typename T>
  T operator()(T a)
  {
    return structured_fv::asin(a);
  }

  static constexpr const char* name = "asin";
  static constexpr bool HasComplexSupport = true;
};

struct Acos
{
  template <typename T>
  T operator()(T a)
  {
    return structured_fv::acos(a);
  }

  static constexpr const char* name = "acos";
  static constexpr bool HasComplexSupport = true;
};

struct Atan
{
  template <typename T>
  T operator()(T a)
  {
    return structured_fv::atan(a);
  }

  static constexpr const char* name = "atan";
  static constexpr bool HasComplexSupport = true;
};

template <typename Func>
constexpr bool HasComplexSupport = std::is_invocable_v<Func, std::complex<double>>;


using BinaryFunctions = ::testing::Types<Plus,
                                         Minus,
                                         Times,
                                         Div,
                                         PlusEqual,
                                         MinusEqual,
                                         TimesEqual,
                                         DivEqual>;

using OutOfPlaceFunctions = ::testing::Types<Plus,
                                             Minus,
                                             Times,
                                             Div,
                                             Pow>;

using UnaryFunctions = ::testing::Types<UnaryMinus,
                                        Sqrt,
                                        Cbrt,
                                        Exp,
                                        Exp2,
                                        Expm1,
                                        NaturalLog,
                                        Log10,
                                        Log2,
                                        Log1p,
                                        Sin,
                                        Cos,
                                        Tan,
                                        Asin,
                                        Acos,
                                        Atan>;
class NameGenerator
{
  public:
    template <typename T>
    static std::string GetName(int)
    {
      return T::name;     
    }
};

TYPED_TEST_SUITE(DualNumberBinaryFunctionTester, BinaryFunctions, NameGenerator);

TYPED_TEST_SUITE(DualNumberOutOfPlaceFunctionTester, OutOfPlaceFunctions, NameGenerator);

TYPED_TEST_SUITE(DualNumberUnaryFunctionTester, UnaryFunctions, NameGenerator);

template <typename Func>
constexpr bool IsComplexBinaryFunc = std::is_invocable_r_v<Complex, Func, Complex, Complex>;

template <typename Func>
constexpr bool IsComplexUnaryFunc = std::is_invocable_r_v<Complex, Func, Complex>;

template <typename Func>
FixedVec<Real, 2> compute_derivs_cs(const Dual<Real, 2>& a, const Dual<Real, 2>& b, Func func)
{
  static_assert(IsComplexBinaryFunc<Func>);

  Real h = 1e-40;
  Complex pert(0, h);

  FixedVec<Real, 2> derivs;
  for (UInt i=0; i < 2; ++i)
  {
    Complex ac(a(), h*a(i));
    Complex bc(b(), h*b(i));
    Complex cc = func(ac, bc);
    derivs[i] = cc.imag()/h;
  }

  return derivs;
}

template <typename Func>
FixedVec<Real, 2> compute_derivs_cs(const Dual<Real, 2>& a, int b, Func func)
{
  static_assert(IsComplexBinaryFunc<Func>);

  Real h = 1e-40;
  Complex pert(0, h); 

  FixedVec<Real, 2> derivs;
  for (UInt i=0; i < 2; ++i)
  {
    Complex ac(a(), h*a(i));
    Complex bc(b, 0);
    Complex cc = func(ac, bc);
    derivs[i] = cc.imag()/h;
  }

  return derivs;
}

template <typename Func>
FixedVec<Real, 2> compute_derivs_cs(int a, const Dual<Real, 2>& b, Func func)
{
  static_assert(IsComplexBinaryFunc<Func>);

  Real h = 1e-40;
  Complex pert(0, h);

  FixedVec<Real, 2> derivs;
  for (UInt i=0; i < 2; ++i)
  {
    Complex ac(a, 0);
    Complex bc(b(), h*b(i));
    Complex cc = func(ac, bc);
    derivs[i] = cc.imag()/h;
  }

  return derivs;
}

template <typename Func>
FixedVec<Real, 2> compute_derivs_cs(const Dual<Real, 2>& a, Func func)
{
  static_assert(IsComplexUnaryFunc<Func>);

  Real h = 1e-40;
  Complex pert(0, h);

  FixedVec<Real, 2> derivs;
  for (UInt i=0; i < 2; ++i)
  {
    if constexpr (Func::HasComplexSupport)
    {
      Complex ac(a(), h*a(i));
      Complex bc = func(ac);
      derivs[i] = bc.imag()/h;
    } else
    {
      Real h = 1e-7;

      Real f0 = func(a());
      Real f1 = func(a() + h);
      Real dfdx = (f1 - f0)/h;
      derivs[i] = dfdx * a(i);      
    }
  }

  return derivs;
}

}

TEST(DualNumber, Typedefs)
{
  static_assert(std::is_same_v<RealDual2::value_type, double>);
  static_assert(RealDual2::size == 2);
}

TEST(DualNumber, ConstructorZeroArg)
{
  RealDual2 val;
  EXPECT_EQ(val(), 0);
  EXPECT_EQ(val(0), 0);
  EXPECT_EQ(val(1), 0);
}

TEST(DualNumber, ConstructorOneArg)
{
  RealDual2 val(1);
  EXPECT_EQ(val(), 1);
  EXPECT_EQ(val(0), 0);
  EXPECT_EQ(val(1), 0);
}

TEST(DualNumber, ConstructorOneArgDual)
{
  FloatDual2 b(1, {2, 3});
  RealDual2 val(b);
  
  EXPECT_EQ(val(), 1);
  EXPECT_EQ(val(0), 2);
  EXPECT_EQ(val(1), 3);
}

TEST(DualNumber, ConstructorTwoArg)
{
  RealDual2 val(1, {2, 3});
  EXPECT_EQ(val(), 1);
  EXPECT_EQ(val(0), 2);
  EXPECT_EQ(val(1), 3);
}

TEST(DualNumber, Getters)
{
  RealDual2 val(1, {2, 3});
  EXPECT_EQ(val.get_value(), 1);
  EXPECT_EQ(val.get_deriv(0), 2);
  EXPECT_EQ(val.get_deriv(1), 3);

  EXPECT_EQ(val(), 1);
  EXPECT_EQ(val(0), 2);
  EXPECT_EQ(val(1), 3);  
}

TEST(DualNumber, UnaryMinus)
{
  RealDual2 val(1, {2, 3});

  RealDual2 val2 = -val;

  EXPECT_EQ(val2(),  -1);
  EXPECT_EQ(val2(0), -2);
  EXPECT_EQ(val2(1), -3);  
}

TEST(DualNumber, Plus)
{
  RealDual2 val(1, {2, 3});
  FloatDual2 val2(2, {4, 6});

  const RealDual2& val_ret = val + val2;

  EXPECT_NE(&val, &val_ret);
  EXPECT_EQ(val_ret(), 3);
  EXPECT_EQ(val_ret(0), 6);
  EXPECT_EQ(val_ret(1), 9);  
}

TEST(DualNumber, PlusEqual)
{
  RealDual2 val(1, {2, 3});
  FloatDual2 val2(2, {4, 6});

  RealDual2& val_ret = val += val2;

  EXPECT_EQ(&val, &val_ret);
  EXPECT_EQ(val(), 3);
  EXPECT_EQ(val(0), 6);
  EXPECT_EQ(val(1), 9);  
}

TEST(DualNumber, Minus)
{
  RealDual2 val(1, {2, 3});
  FloatDual2 val2(2, {4, 6});

  const RealDual2& val_ret = val - val2;
  
  EXPECT_NE(&val, &val_ret);
  EXPECT_EQ(val_ret(),  -1);
  EXPECT_EQ(val_ret(0), -2);
  EXPECT_EQ(val_ret(1), -3);  
}


TEST(DualNumber, MinusEqual)
{
  RealDual2 val(1, {2, 3});
  FloatDual2 val2(2, {4, 6});

  RealDual2& val_ret = val -= val2;
  
  EXPECT_EQ(&val, &val_ret);
  EXPECT_EQ(val(),  -1);
  EXPECT_EQ(val(0), -2);
  EXPECT_EQ(val(1), -3);  
}

TEST(DualNumber, TimesEqual)
{
  RealDual2 val(1, {2, 3});
  FloatDual2 val2(2, {4, 6});

  RealDual2& val_ret = val *= val2;
  
  EXPECT_EQ(&val, &val_ret);
  EXPECT_EQ(val(),  2);
  EXPECT_EQ(val(0), 4 + 4);
  EXPECT_EQ(val(1), 6 + 6);  
}


TEST(DualNumber, Times)
{
  RealDual2 val(1, {2, 3});
  FloatDual2 val2(2, {4, 6});

  const RealDual2& val_ret = val * val2;
  
  EXPECT_NE(&val, &val_ret);
  EXPECT_EQ(val_ret(),  2);
  EXPECT_EQ(val_ret(0), 4 + 4);
  EXPECT_EQ(val_ret(1), 6 + 6);  
}

TEST(DualNumber, DivEqual)
{
  RealDual2 val(1, {2, 3});
  FloatDual2 val2(2, {4, 6});

  FloatDual2& val2_ret = val2 /= val;
  
  EXPECT_EQ(&val2, &val2_ret);
  EXPECT_EQ(val2(),  2);
  EXPECT_EQ(val2(0), 0);
  EXPECT_EQ(val2(1), 0);
}

TEST(DualNumber, Div)
{
  RealDual2 val(1, {2, 3});
  FloatDual2 val2(2, {4, 6});

  const FloatDual2& val2_ret = val2 / val;
  
  EXPECT_NE(&val2, &val2_ret);
  EXPECT_EQ(val2_ret(),  2);
  EXPECT_EQ(val2_ret(0), 0);
  EXPECT_EQ(val2_ret(1), 0);  
}

TEST(DualNumber, ComparisonOperators)
{
  RealDual2 val(2, {2, 3});
  FloatDual2 val2(1, {4, 6});
  FloatDual2 val3(2, {4, 6});

  EXPECT_FALSE(val <  val2);
  EXPECT_FALSE(val <= val2);
  EXPECT_TRUE( val >  val2);
  EXPECT_TRUE( val >= val2);  
  
  EXPECT_FALSE(val <  val3);
  EXPECT_TRUE( val <= val3);
  EXPECT_FALSE(val >  val3);
  EXPECT_TRUE( val >= val3);
}

TEST(DualNumber, EqualityOperators)
{
  RealDual2 val1(2, {2, 3});
  RealDual2 val2(2, {3, 4});

  structured_fv::DualStrictlyEqual pred;

  EXPECT_EQ(val1, val2);
  EXPECT_FALSE(pred(val1, val2));

  EXPECT_FALSE(val1 != val2);
}


TYPED_TEST(DualNumberBinaryFunctionTester, DualDual)
{
  using T = TypeParam;
  T binary_func;
  RealDual2 val(1, {2, 3});
  RealDual2 val2(2, {4, 6});
  RealDual2 val3 = binary_func(val, val2);
  FixedVec<double, 2> derivs = compute_derivs_cs(val, val2, binary_func);
  
  EXPECT_EQ(val3(), binary_func(val(), val2()));
  EXPECT_EQ(val3(0), derivs[0]);
  EXPECT_EQ(val3(1), derivs[1]);
}

TYPED_TEST(DualNumberBinaryFunctionTester, DualInt)
{
  using T = TypeParam;
  T binary_func;
  RealDual2 val(1, {2, 3});
  int val2 = 2;
  RealDual2 val3 = binary_func(val, val2);
  FixedVec<double, 2> derivs = compute_derivs_cs(val, val2, binary_func);
  
  EXPECT_EQ(val3(), binary_func(val(), val2));
  EXPECT_EQ(val3(0), derivs[0]);
  EXPECT_EQ(val3(1), derivs[1]);
}


TYPED_TEST(DualNumberOutOfPlaceFunctionTester, IntDual)
{
  using T = TypeParam;
  T binary_func;
  int val = 2;
  RealDual2 val2(1, {2, 3});
  RealDual2 val3 = binary_func(val, val2);
  FixedVec<double, 2> derivs = compute_derivs_cs(val, val2, binary_func);
  
  EXPECT_EQ(val3(), binary_func(val, val2()));
  EXPECT_EQ(val3(0), derivs[0]);
  EXPECT_EQ(val3(1), derivs[1]);
}

TYPED_TEST(DualNumberOutOfPlaceFunctionTester, DualInt)
{
  using T = TypeParam;
  T binary_func;
  RealDual2 val(1, {2, 3});
  int val2 = 2;
  RealDual2 val3 = binary_func(val, val2);
  FixedVec<double, 2> derivs = compute_derivs_cs(val, val2, binary_func);
  
  EXPECT_EQ(val3(), binary_func(val(), val2));
  EXPECT_EQ(val3(0), derivs[0]);
  EXPECT_EQ(val3(1), derivs[1]);
}

TYPED_TEST(DualNumberOutOfPlaceFunctionTester, DualDual)
{
  using T = TypeParam;
  T binary_func;
  RealDual2 val(1, {2, 3});
  RealDual2 val2(4, {5, 6});
  RealDual2 val3 = binary_func(val, val2);
  FixedVec<double, 2> derivs = compute_derivs_cs(val, val2, binary_func);
  
  EXPECT_EQ(val3(), binary_func(val(), val2()));
  EXPECT_EQ(val3(0), derivs[0]);
  EXPECT_EQ(val3(1), derivs[1]);
}

TYPED_TEST(DualNumberUnaryFunctionTester, Dual)
{
  using T = TypeParam;
  T unary_func;
  RealDual2 val(0.5, {3, 4});
  RealDual2 val2 = unary_func(val);
  FixedVec<double, 2> derivs = compute_derivs_cs(val, unary_func);
  
  EXPECT_EQ(val2(), unary_func(val()));
  if constexpr (T::HasComplexSupport)
  {
    EXPECT_DOUBLE_EQ(val2(0), derivs[0]);
    EXPECT_DOUBLE_EQ(val2(1), derivs[1]);
  } else
  {
    EXPECT_NEAR(val2(0), derivs[0], 1e-5);
    EXPECT_NEAR(val2(1), derivs[1], 1e-5);    
  }
}
