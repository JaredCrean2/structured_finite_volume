
#include "gtest/gtest.h"
#include "linear_system/large_matrix.h"
#include "linear_system/large_matrix_dense.h"

namespace {
using namespace structured_fv;
using Mat3 = Matrix<Real, 3, 3>;
}

TEST(LargeMatrixDense, GeneralSolve)
{
  linear_system::LargeMatrixOpts opts;
  opts.is_structurally_symmetric = false;
  opts.is_value_symmetric        = false;
  opts.factor_in_place           = false;

  linear_system::LargeMatrixDense mat("A", 3, 3, opts);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  FixedVec<GlobalDof, 3> dofs{0, 1, 2};
  Mat3 jac({1, 2, 3,
            4, 5, 6,
            8, 8, 9});
  linear_system::VectorType x("x", 3);
  linear_system::VectorType b("b", 3);
  b(0) = 3;
  b(1) = 6;
  b(2) = 9;
  Vec3<Real> x_ex {0, 0, 1};

  mat.assembleValues(dofs, dofs, jac);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x(i), x_ex[i], 1e-13);
}

TEST(LargeMatrixDense, AssembleValuesAdditive)
{
  linear_system::LargeMatrixOpts opts;
  opts.is_structurally_symmetric = false;
  opts.is_value_symmetric        = false;
  opts.factor_in_place           = false;

  linear_system::LargeMatrixDense mat("A", 3, 3, opts);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  FixedVec<GlobalDof, 3> dofs{0, 1, 2};
  Mat3 vals({1, 2, 3,
             4, 5, 3,
             4, 4, 9});
  Mat3 vals2({0, 0, 0,
              0, 0, 3,
              4, 4, 0});
  linear_system::VectorType x("x", 3);
  linear_system::VectorType b("b", 3);
  b(0) = 3;
  b(1) = 6;
  b(2) = 9;
  Vec3<Real> x_ex = {0, 0, 1};

  mat.assembleValues(dofs, dofs, vals);
  mat.assembleValues(dofs, dofs, vals2);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x(i), x_ex[i], 1e-13);
}


TEST(LargeMatrixDense, AssembleValuesIgnore)
{
  linear_system::LargeMatrixOpts opts;
  opts.is_structurally_symmetric = false;
  opts.is_value_symmetric        = false;
  opts.factor_in_place           = false;

  linear_system::LargeMatrixDense mat("A", 3, 3, opts);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  FixedVec<GlobalDof, 4> dofs{0, 1, 2, -1};
  Matrix<Real, 4, 4> vals({1, 2, 3,       666,
                           4, 5, 6,       666,
                           8, 8, 9,       666,
                           666, 666, 666, 666});
  linear_system::VectorType x("x", 3);
  linear_system::VectorType b("b", 3);
  b(0) = 3;
  b(1) = 6;
  b(2) = 9;
  Vec3<Real> x_ex = {0, 0, 1};

  mat.assembleValues(dofs, dofs, vals);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);
}


TEST(LargeMatrixDense, ZeroMatrix)
{
  linear_system::LargeMatrixOpts opts;
  opts.is_structurally_symmetric = false;
  opts.is_value_symmetric        = false;
  opts.factor_in_place           = false;

  linear_system::LargeMatrixDense mat("A", 3, 3, opts);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  FixedVec<GlobalDof, 3> dofs{0, 1, 2};
  Mat3 vals({1, 2, 3,
             4, 5, 6,
             8, 8, 9});
  linear_system::VectorType x("x", 3);
  linear_system::VectorType b("b", 3);
  b(0) = 3;
  b(1) = 6;
  b(2) = 9;
  Vec3<Real> x_ex = {0, 0, 1};

  mat.assembleValues(dofs, dofs, vals);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x(i), x_ex[i], 1e-13);


  mat.zeroMatrix();
  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
      vals(i, j) *= 2;

  mat.assembleValues(dofs, dofs, vals);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x(i), 0.5 * x_ex[i], 1e-13);
}


TEST(LargeMatrixDense, FactorInPlace)
{
  linear_system::LargeMatrixOpts opts;
  opts.is_structurally_symmetric = false;
  opts.is_value_symmetric        = false;
  opts.factor_in_place           = true;

  linear_system::LargeMatrixDense mat("A", 3, 3, opts);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  FixedVec<GlobalDof, 3> dofs{0, 1, 2};
  Mat3 vals({1, 2, 3,
             4, 5, 6,
             8, 8, 9});
  linear_system::VectorType x("x", 3);
  linear_system::VectorType b("b", 3);
  b(0) = 3;
  b(1) = 6;
  b(2) = 9;
  Vec3<Real> x_ex = {0, 0, 1};

  mat.assembleValues(dofs, dofs, vals);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);
}


TEST(LargeMatrixDense, SPD)
{
  linear_system::LargeMatrixOpts opts;
  opts.is_structurally_symmetric = true;
  opts.is_value_symmetric        = true;
  opts.factor_in_place           = false;

  linear_system::LargeMatrixDense mat("A", 3, 3, opts);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  FixedVec<GlobalDof, 3> dofs{0, 1, 2};
  Mat3 vals({10, 2, 3,
             2, 12, 5,
             3, 5, 20});
  linear_system::VectorType x("x", 3);
  linear_system::VectorType b("b", 3);
  b(0) = 1;
  b(1) = 2;
  b(2) = 3;
  
  Vec3<Real> x_ex = {0.043027, 0.111276, 0.115727};

  mat.assembleValues(dofs, dofs, vals);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-5);
}


TEST(LargeMatrixDense, SPDInPlace)
{  
  linear_system::LargeMatrixOpts opts;
  opts.is_structurally_symmetric = true;
  opts.is_value_symmetric        = true;
  opts.factor_in_place           = true;

  linear_system::LargeMatrixDense mat("A", 3, 3, opts);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);
  
  FixedVec<GlobalDof, 3> dofs{0, 1, 2};
  Mat3 vals({10, 2, 3,
             2, 12, 5,
             3, 5, 20});
  linear_system::VectorType x("x", 3);
  linear_system::VectorType b("b", 3);
  b(0) = 1;
  b(1) = 2;
  b(2) = 3;
  
  Vec3<Real> x_ex = {0.043027, 0.111276, 0.115727};

  mat.assembleValues(dofs, dofs, vals);
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-5);
}

