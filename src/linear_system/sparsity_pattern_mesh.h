#ifndef STRUCTURED_FINITE_VOLUME_LINEAR_SYSTEM_SPARSITY_PATTERN_MESH_H
#define STRUCTURED_FINITE_VOLUME_LINEAR_SYSTEM_SPARSITY_PATTERN_MESH_H

#include "disc/discretization.h"
#include "linear_system/sparsity_pattern.h"
#include "utils/project_defs.h"
#include <memory>


namespace structured_fv {
namespace linear_system {

class SparsityPatternDisc : public SparsityPattern
{
  public:
    explicit SparsityPatternDisc(disc::StructuredDiscPtr disc, UInt stencil_size);

    PetscInt getNumOwnedDofs() const override;

    const std::vector<PetscInt>& getDiagonalCounts() override
    { 
      if (!m_computed_nonsymmetric)
        computePattern(false);

      return m_onproc_dofs;
    }

    const std::vector<PetscInt>& getDiagonalCountsSym() override
    { 
      if (!m_computed_symmetric)
        computePattern(true);

      return m_onproc_dofs_sym;
    }

    const std::vector<PetscInt>& getOffProcCounts() override 
    {
      if (!m_computed_nonsymmetric)
        computePattern(false);

      return m_remote_dofs; 
    }

    const std::vector<PetscInt>& getOffProcCountsSym() override 
    {
      if (!m_computed_symmetric)
        computePattern(true);

      return m_remote_dofs_sym; 
    }

    const std::vector<PetscInt>& getGhostGlobalIndices() override
    {
      //return m_ghost_global_dofs;
      return m_empty_vector;
    }
/*
    const std::vector<PetscInt>& getGhostLocalIndices() override
    {
      return m_ghost_onproc_dofs;
    }

    const std::vector<PetscInt>& getOwnedToLocalInfo() override
    {
      return m_owned_dof_to_local;
    }

    const std::vector<PetscInt>& getLocalToGlobalDofs() override
    {
      return m_local_to_global_dofs;
    }
*/

  private:
    void computePattern(bool symmetric);

    bool m_computed_nonsymmetric = false;
    bool m_computed_symmetric    = false;
    disc::StructuredDiscPtr m_disc;
    UInt m_stencil_size;
    std::vector<PetscInt> m_onproc_dofs;
    std::vector<PetscInt> m_remote_dofs;
    std::vector<PetscInt> m_onproc_dofs_sym;
    std::vector<PetscInt> m_remote_dofs_sym;
    std::vector<PetscInt> m_ghost_global_dofs;
    std::vector<PetscInt> m_empty_vector;
    //std::vector<PetscInt> m_ghost_onproc_dofs;
    //std::vector<PetscInt> m_owned_dof_to_local;
    //std::vector<PetscInt> m_local_to_global_dofs;
};


} // namespace
}

#endif