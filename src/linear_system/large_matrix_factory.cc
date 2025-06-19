#include "linear_system/large_matrix_factory.h"
#include "linear_system/large_matrix_dense.h"
#include "linear_system/large_matrix_petsc.h"

namespace structured_fv {
namespace linear_system {

std::shared_ptr<LargeMatrix> largeMatrixFactory(LargeMatrixType type,
                                               std::shared_ptr<LargeMatrixOpts> opts, 
                                               std::shared_ptr<SparsityPattern> sparsity)
{
  switch(type)
  {
    case LargeMatrixType::Dense:
    {
      return std::make_shared<LargeMatrixDense>("A", sparsity->getNumOwnedDofs(), sparsity->getNumOwnedDofs(), *opts);
    }

    case LargeMatrixType::Petsc:
    {
      auto& opts_petsc = dynamic_cast<LargeMatrixOptsPetsc&>(*opts);
      return std::make_shared<LargeMatrixPetsc>("A", opts_petsc, sparsity);
    }

    default:
      throw std::runtime_error("unrecognized type");
  }
}

}
}