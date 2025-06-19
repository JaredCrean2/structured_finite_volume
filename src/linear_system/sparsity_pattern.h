#ifndef STRUCTURED_FINITE_VOLUME_LINEAR_SYSTEM_SPARSITY_PATTERN_H
#define STRUCTURED_FINITE_VOLUME_LINEAR_SYSTEM_SPARSITY_PATTERN_H

#include <petscsystypes.h>
#include <vector>

namespace structured_fv {
namespace linear_system {

class SparsityPattern
{
  public:
    virtual ~SparsityPattern() = default;

    virtual PetscInt getNumOwnedDofs() const = 0;

    virtual PetscInt getNumLocalDofs() { return getNumOwnedDofs(); /*+ getGhostLocalIndices().size();*/ }

    // returns vector giving number of local dofs connected to each dof
    virtual const std::vector<PetscInt>& getDiagonalCounts() = 0;

    // similar to the above, but returns the count for only the matrix diagonal + upper triangle
    virtual const std::vector<PetscInt>& getDiagonalCountsSym() = 0;

    // returns vector giving number of remote dofs connected to each dof
    virtual const std::vector<PetscInt>& getOffProcCounts() = 0;
        
    // similar to the above, but returns the count for only the upper triangle
    virtual const std::vector<PetscInt>& getOffProcCountsSym() = 0;

    virtual const std::vector<PetscInt>& getGhostGlobalIndices() = 0;

    //virtual const std::vector<PetscInt>& getGhostLocalIndices() = 0;

    //virtual const std::vector<PetscInt>& getOwnedToLocalInfo() = 0;

    //virtual const std::vector<PetscInt>& getLocalToGlobalDofs() = 0;
};

} // namespace
}

#endif