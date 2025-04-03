#include "gtest/gtest.h"
#include "disc/disc_vector.h"

using namespace structured_fv;

TEST(DiscVector, Indexing)
{
  GlobalDof num_dofs = 10;
  disc::DiscVector<Real> vec(num_dofs, "foo");

  EXPECT_EQ(vec.size(), num_dofs);

  for (GlobalDof i=0; i < num_dofs; ++i)
    vec(i) = i;

  for (GlobalDof i=0; i < num_dofs; ++i)
    EXPECT_EQ(vec(i), i);
}

TEST(DiscVector, Set)
{
  GlobalDof num_dofs = 10;
  disc::DiscVector<Real> vec(num_dofs, "foo");

  vec.set(42);
  for (GlobalDof i=0; i < num_dofs; ++i)
    EXPECT_EQ(vec(i), 42);  
}