#ifndef STRUCTURED_FINITE_VOLUME_LINEAR_SYSTEM_ASSEMBLER_BASE_H
#define STRUCTURED_FINITE_VOLUME_LINEAR_SYSTEM_ASSEMBLER_BASE_H

#include <memory>
#include "utils/project_defs.h"


namespace structured_fv {
namespace linear_system {

// the idea for Assembler and AssemblerBase is they provide the two different
// interfaces for Assemblers.  AssemblerBase is the interface for the linear
// solvers and nonlinear solvers, while Assembler is the interface for physics
// modules to put values into the matrix


class AssemblerBase
{
  public:    
  
    virtual ~AssemblerBase() = default;

    virtual void setAlpha(Real alpha) = 0;

    virtual Real getAlpha() const = 0;
};

using AssemblerBasePtr = std::shared_ptr<AssemblerBase>;


} // namespace
}

#endif
