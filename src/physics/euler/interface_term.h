#include "euler_model.h"
#include "typedefs.h"

namespace structured_fv {
namespace euler {

void evaluateInterfaceTermsEntryPoint(const EulerOpts& opts, const ElementFieldPtr<Real>& solution, Real t, const StructuredDiscPtr disc,
                                      ElementFieldPtr<Real> residual);

}
}