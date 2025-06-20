#include "gtest/gtest.h"
#include "utils/array3.h"

using namespace structured_fv;

TEST(Array3, DefaultConstructor)
{
  Array3<Int, 3, 4, 5> mat;
  EXPECT_EQ(mat.extent0(), 3);
  EXPECT_EQ(mat.extent1(), 4);
  EXPECT_EQ(mat.extent2(), 5);

  EXPECT_EQ(mat.size()[0], 3);
  EXPECT_EQ(mat.size()[1], 4);
  EXPECT_EQ(mat.size()[2], 5);
  for (UInt i=0; i < 3; ++i)
    for (UInt j=0; j < 4; ++j)
      for (UInt k=0; k < 5; ++k)
        EXPECT_EQ(mat(i, j, k), 0);

  const Array3<Int, 3, 4, 5>& mat_const = mat;
  EXPECT_EQ(mat_const.extent0(), 3);
  EXPECT_EQ(mat_const.extent1(), 4);
  EXPECT_EQ(mat_const.extent2(), 5);
  EXPECT_EQ(mat_const.size()[0], 3);
  EXPECT_EQ(mat_const.size()[1], 4);
  EXPECT_EQ(mat_const.size()[2], 5);

  for (UInt i=0; i < 3; ++i)
    for (UInt j=0; j < 4; ++j)
      for (UInt k=0; k < 5; ++k)
        EXPECT_EQ(mat_const(i, j, k), 0);  
}

TEST(Array3, ValueConstructor)
{
  Array3<Int, 2, 3, 4> mat({1,  2,  3,  4,
                            5,  6,  7,  8,
                            9, 10, 11, 12,
                            
                            13, 14, 15, 16,
                            17, 18, 19, 20,
                            21, 22, 23, 24});
  
  Int val = 1;
  for (UInt i=0; i < 2; ++i)
    for (UInt j=0; j < 3; ++j)
      for (UInt k=0; k < 4; ++k)
        EXPECT_EQ(mat(i, j, k), val++);

  const Array3<Int, 2, 3, 4>& mat_const = mat;
  val = 1;
  for (UInt i=0; i < 2; ++i)
    for (UInt j=0; j < 3; ++j)
      for (UInt k=0; k < 4; ++k)
        EXPECT_EQ(mat_const(i, j, k), val++);
}

TEST(Array3, DefaultTemplateArgument)
{
  Array3<Int, 3> mat;
  EXPECT_EQ(mat.extent0(), 3);
  EXPECT_EQ(mat.extent1(), 3);
  EXPECT_EQ(mat.extent2(), 3);

  EXPECT_EQ(mat.size()[0], 3);
  EXPECT_EQ(mat.size()[1], 3);
  EXPECT_EQ(mat.size()[2], 3);
}