#include "linear_system/sparsity_pattern_mesh.h"
#include "disc/discretization.h"
#include "disc/elem_field.h"
#include "utils/neighbor_direction.h"

namespace structured_fv {
namespace linear_system {

SparsityPatternDisc::SparsityPatternDisc(disc::StructuredDiscPtr disc, UInt stencil_size) :
  m_disc(disc),
  m_stencil_size(stencil_size)
{}

PetscInt SparsityPatternDisc::getNumOwnedDofs() const
{
  return m_disc->getNumDofs();
}

void SparsityPatternDisc::computePattern(bool symmetric)
{
  auto& onproc_dofs  = symmetric ? m_onproc_dofs_sym : m_onproc_dofs;
  auto& remote_dofs = symmetric ? m_remote_dofs_sym : m_remote_dofs;

  onproc_dofs.resize(m_disc->getNumDofs());
  remote_dofs.resize(m_disc->getNumDofs());

  for (UInt block_id : m_disc->getRegularBlocksIds())
  {
    std::cout << "block = " << block_id << std::endl;
    const disc::StructuredBlock& block = m_disc->getBlock(block_id);
    const auto dof_nums = m_disc->getDofNumbering()->getData(block_id);
    for (UInt i : block.getOwnedCells().getXRange())
      for (UInt j : block.getOwnedCells().getYRange())
        for (UInt k=0; k < m_disc->getNumDofsPerNode(); ++k)
        {
          std::cout << "i, j, k = " << i << ", " << j << ", " << k << std::endl;
          std::cout << "stencil size = " << m_stencil_size << std::endl;
          UInt num_connected_dofs = m_disc->getNumDofsPerNode();
          GlobalDof dof_k = dof_nums(i, j, k);
          for (UInt p=0; p < m_stencil_size; ++p)
          {
            std::cout << "p = " << p << std::endl;
            for (NeighborDirection dir : {NeighborDirection::North, NeighborDirection::East,
                                          NeighborDirection::South, NeighborDirection::West})
            {
              auto [i2, j2] = computeIndices(dir, p+1, i, j);
              for (UInt k2=0; k2 < m_disc->getNumDofsPerNode(); ++k2)
              {
                std::cout << "k2 = " << k2 << std::endl;
                if (symmetric && dof_nums(i2, j2, k2) > dof_k)
                {
                  num_connected_dofs++;
                } else
                {
                  num_connected_dofs++;
                }
              }
            }
          }

          onproc_dofs[dof_k] = num_connected_dofs;
          remote_dofs[dof_k] = 0;
        }
  }

  if (symmetric)
    m_computed_symmetric = true;
  else
    m_computed_nonsymmetric = true;
}


}  // namespace
}