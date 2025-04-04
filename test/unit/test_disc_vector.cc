#include "disc/disc_vector.h"
#include "disc/discretization.h"
#include "mesh/block_spec.h"
#include "mesh/structured_mesh.h"
#include "gtest/gtest.h"

namespace {
using namespace structured_fv;

struct DiscVectorTester : public ::testing::Test {
  public:
    DiscVectorTester() {
      mesh::MeshSpec spec(1, 1);
      spec.blocks(0, 0) = mesh::MeshBlockSpec(2, 5);
      auto mesh = std::make_shared<mesh::StructuredMesh>(spec);
      m_disc = std::make_shared<disc::StructuredDisc>(mesh, 0, 1);
      m_num_dofs = m_disc->getNumDofs();
    }

    disc::StructuredDiscPtr m_disc;
    GlobalDof m_num_dofs = 0;
};

} // namespace

TEST_F(DiscVectorTester, Indexing)
{
  disc::DiscVector<Real> vec(m_disc, "foo");

  EXPECT_EQ(vec.size(), m_num_dofs);

  for (GlobalDof i = 0; i < m_num_dofs; ++i)
    vec(i) = i;

  for (GlobalDof i = 0; i < m_num_dofs; ++i)
    EXPECT_EQ(vec(i), i);
}

TEST_F(DiscVectorTester, Set)
{
  disc::DiscVector<Real> vec(m_disc, "foo");

  vec.set(42);
  for (GlobalDof i = 0; i < m_num_dofs; ++i)
    EXPECT_EQ(vec(i), 42);
}

