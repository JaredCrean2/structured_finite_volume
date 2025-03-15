#include "gtest/gtest.h"
#include "utils/bitwise.h"

using namespace structured_fv;

TEST(Bitwise, Masks)
{
  EXPECT_EQ(mask_if_equal(0, 0), UIntMax);
  EXPECT_EQ(mask_if_equal(1, 1), UIntMax);
  EXPECT_EQ(mask_if_equal(0, 1), 0U);

  EXPECT_EQ(mask_if_notequal(0, 0), 0U);
  EXPECT_EQ(mask_if_notequal(1, 1), 0U);
  EXPECT_EQ(mask_if_notequal(0, 1), UIntMax);
}

