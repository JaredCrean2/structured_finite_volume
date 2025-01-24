#include "mesh/structured_mesh.h"
#include "gtest/gtest.h"

using namespace structured_fv;
using namespace structured_fv::mesh;

TEST(StructuredMesh, SingleBlock) {
  MeshSpec meshspec(1, 1, 2);
  meshspec.blocks(0, 0) = MeshBlockSpec(2, 3, 0);

  StructuredMesh mesh(meshspec);

  EXPECT_EQ(mesh.getNumBlocks(), 5);
  EXPECT_EQ(mesh.getNumRegularBlocks(), 1);
  EXPECT_EQ(mesh.getNumGhostBCBlocks(), 4);
  EXPECT_EQ(mesh.getNumBlockInterfaces(), 4);
  EXPECT_EQ(mesh.getNumRegularBlockInterfaces(), 0);
  EXPECT_EQ(mesh.getNumGhostBCBlockInterfaces(), 4);
}