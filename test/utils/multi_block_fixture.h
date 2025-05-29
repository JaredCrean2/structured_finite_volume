#ifndef STRUCTURE_FINITE_VOLUME_TEST_UTILS_MULTI_BLOCK_FIXTURE_H
#define STRUCTURE_FINITE_VOLUME_TEST_UTILS_MULTI_BLOCK_FIXTURE_H

#include "utils/project_defs.h"
#include  "mesh/structured_mesh.h"
#include "disc/discretization.h"
#include "disc/disc_vector.h"

using namespace structured_fv;

namespace test_utils {

using structured_fv::UInt;
using structured_fv::Real;

// sets up two meshes, the first one with a single block, the second one
// that uniformly subdivides that block into 4 blocks.
// Provides functions that map between indices of the two meshes, to check
// that physics models give the same results on both meshes
class MultiBlockFixture 
{
  public:
  MultiBlockFixture();

  void setup_disc(UInt dofs_per_cell); 


  // given the indices of a cell on the single block of mesh1, compute the
  // (block, i, j) for the same cell on mesh 2
  std::tuple<UInt, UInt, UInt> get_idx(UInt i, UInt j) const;

  void test_residuals_equal(std::shared_ptr<disc::DiscVector<Real>> residual1,
                            std::shared_ptr<disc::DiscVector<Real>> residual2) const;

  public:
    UInt m_dofs_per_cell = 1;
    int m_num_bc_ghost_cells = 2;
    mesh::MeshSpec spec1;
    std::shared_ptr<mesh::StructuredMesh> m_mesh1;
    std::shared_ptr<disc::StructuredDisc> m_disc1;
    disc::DiscVectorPtr<Real> m_solution1;
    disc::DiscVectorPtr<Real> m_residual1;

    mesh::MeshSpec spec2;
    std::shared_ptr<mesh::StructuredMesh> m_mesh2;
    std::shared_ptr<disc::StructuredDisc> m_disc2;
    disc::DiscVectorPtr<Real> m_solution2;
    disc::DiscVectorPtr<Real> m_residual2;    
};

}

#endif