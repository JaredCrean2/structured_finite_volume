#include "gtest/gtest.h"
#include "utils/range.h"
#include <vector>

using namespace structured_fv;

TEST(RangeIter, Dereference)
{
  RangeIter iter(2);
  EXPECT_EQ(*iter, 2);
}

TEST(RangeIter, PrefixIncrement)
{
  RangeIter iter(2);
  EXPECT_EQ(*(++iter), 3);
  EXPECT_EQ(*iter, 3);
}

TEST(RangeIter, PostfixIncrement)
{
  RangeIter iter(2);
  EXPECT_EQ(*(iter++), 2);
  EXPECT_EQ(*iter, 3);
}

TEST(RangeIter, PrefixDecrement)
{
  RangeIter iter(2);
  EXPECT_EQ(*(--iter), 1);
  EXPECT_EQ(*iter, 1);
}

TEST(RangeIter, PostfixDecrement)
{
  RangeIter iter(2);
  EXPECT_EQ(*(iter--), 2);
  EXPECT_EQ(*iter, 1);
}

TEST(RangeIter, PlusEqual)
{
  RangeIter iter(2);
  EXPECT_EQ(*(iter += 2), 4);
  EXPECT_EQ(*iter, 4);
}

TEST(RangeIter, MinusEqual)
{
  RangeIter iter(2);
  EXPECT_EQ(*(iter -= 2), 0);
  EXPECT_EQ(*iter, 0);
}

TEST(RangeIter, Plus)
{
  RangeIter iter(2);
  EXPECT_EQ(*(iter + 2), 4);
  EXPECT_EQ(*iter, 2);
}

TEST(RangeIter, Minus)
{
  RangeIter iter(2);
  EXPECT_EQ(*(iter - 2), 0);
  EXPECT_EQ(*iter, 2);
}

TEST(RangeIter, MinusIterator)
{
  RangeIter iter(2);
  EXPECT_EQ(iter - RangeIter(2), 0);
  EXPECT_EQ(*iter, 2);
}

TEST(RangeIter, Equality)
{
  EXPECT_TRUE(RangeIter(2) == RangeIter(2));
  EXPECT_FALSE(RangeIter(2) != RangeIter(2));

  EXPECT_FALSE(RangeIter(2) == RangeIter(3));
  EXPECT_TRUE(RangeIter(2) != RangeIter(3));
}

TEST(RangeIter, RelationalOperators)
{
  EXPECT_TRUE(RangeIter(2)  < RangeIter(3));
  EXPECT_TRUE(RangeIter(2)  <= RangeIter(3));
  EXPECT_FALSE(RangeIter(2) > RangeIter(3));
  EXPECT_FALSE(RangeIter(2) >= RangeIter(3));

  EXPECT_FALSE(RangeIter(3)  < RangeIter(2));
  EXPECT_FALSE(RangeIter(3)  <= RangeIter(2));
  EXPECT_TRUE(RangeIter(3) > RangeIter(2));
  EXPECT_TRUE(RangeIter(3) >= RangeIter(2));

  EXPECT_FALSE(RangeIter(2)  < RangeIter(2));
  EXPECT_TRUE(RangeIter(2)  <= RangeIter(2));
  EXPECT_FALSE(RangeIter(2) > RangeIter(2));
  EXPECT_TRUE(RangeIter(2) >= RangeIter(2));    
}

TEST(Range, CallOperator)
{
  Range range(2, 5);
  EXPECT_EQ(range(0), 2);
  EXPECT_EQ(range(1), 3);
  EXPECT_EQ(range(2), 4);
}

TEST(Range, RangeFor)
{
  Range range(2, 5);
  std::vector<UInt> vals;

  for (UInt i : range)
    vals.push_back(i);


  std::vector<UInt> vals_expected = {2, 3, 4};
  EXPECT_EQ(vals, vals_expected);
}

TEST(Range, RangeForEmptyRange)
{
  Range range(2, 2);
  std::vector<UInt> vals;

  for (UInt i : range)
    vals.push_back(i);

  EXPECT_EQ(vals.size(), 0);  
}

TEST(Range, RangeForNegativeRange)
{
  Range range(2, 1);
  std::vector<UInt> vals;

  for (UInt i : range)
    vals.push_back(i);

  EXPECT_EQ(vals.size(), 0);  
}

TEST(Range, Length)
{
  EXPECT_EQ(Range(2, 5).size(), 3);
}

TEST(Range, Equality)
{
  EXPECT_TRUE(Range(2, 4) == Range(2, 4));
  EXPECT_FALSE(Range(2, 4) == Range(2, 5));
  EXPECT_FALSE(Range(2, 5) == Range(2, 4));

  EXPECT_FALSE(Range(2, 4) != Range(2, 4));
  EXPECT_TRUE(Range(2, 4) != Range(2, 5));
  EXPECT_TRUE(Range(2, 5) != Range(2, 4));  
}

TEST(Range, In)
{
  EXPECT_FALSE(in(Range(2, 4), 1));
  EXPECT_TRUE( in(Range(2, 4), 2));
  EXPECT_TRUE( in(Range(2, 4), 3));
  EXPECT_FALSE(in(Range(2, 4), 4));  
}

TEST(Range2D, RangeFor)
{
  Range2D range(2, 5, 5, 7);
  std::vector<std::pair<UInt, UInt>> vals;

  for (UInt i : range.getXRange())
    for (UInt j : range.getYRange())
      vals.emplace_back(i, j);

  std::vector<std::pair<UInt, UInt>> vals_expected = { {2, 5}, {2, 6},
                                                       {3, 5}, {3, 6},
                                                       {4, 5}, {4, 6}};
  EXPECT_EQ(vals, vals_expected);
}

TEST(Range2D, Equality)
{
  EXPECT_TRUE(Range2D(0, 1, 0, 2) == Range2D(0, 1, 0, 2));
  EXPECT_FALSE(Range2D(0, 1, 0, 2) == Range2D(1, 1, 0, 2));
  EXPECT_FALSE(Range2D(0, 1, 0, 2) == Range2D(0, 2, 0, 2));
  EXPECT_FALSE(Range2D(0, 1, 0, 2) == Range2D(1, 1, 1, 2));
  EXPECT_FALSE(Range2D(0, 1, 0, 2) == Range2D(1, 1, 0, 3));

  EXPECT_FALSE(Range2D(0, 1, 0, 2) != Range2D(0, 1, 0, 2));
  EXPECT_TRUE(Range2D(0, 1, 0, 2)  != Range2D(1, 1, 0, 2));
  EXPECT_TRUE(Range2D(0, 1, 0, 2)  != Range2D(0, 2, 0, 2));
  EXPECT_TRUE(Range2D(0, 1, 0, 2)  != Range2D(1, 1, 1, 2));
  EXPECT_TRUE(Range2D(0, 1, 0, 2)  != Range2D(1, 1, 0, 3));  
}

TEST(Range2D, In)
{
  Range2D range(2, 4, 3, 5);
  EXPECT_FALSE(in(range, 1, 2));
  EXPECT_FALSE(in(range, 2, 2));
  EXPECT_FALSE(in(range, 2, 2));

  EXPECT_TRUE(in(range, 2, 3));
  EXPECT_TRUE(in(range, 3, 4));

  EXPECT_FALSE(in(range, 4, 4));  
  EXPECT_FALSE(in(range, 3, 5));
  EXPECT_FALSE(in(range, 4, 5));
}