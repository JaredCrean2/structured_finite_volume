#include "linear_system/large_matrix.h"

namespace structured_fv {
namespace linear_system {

void solve(LargeMatrixPtr mat, const disc::DiscVectorPtr<Real> b, disc::DiscVectorPtr<Real> x)
{
  auto b_vec = b->getData();
  auto x_vec = x->getData();

  mat->solve(b_vec, x_vec);
}

} // namespace
}