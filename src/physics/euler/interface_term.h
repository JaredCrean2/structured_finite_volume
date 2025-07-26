#include "euler_model.h"
#include "typedefs.h"

namespace structured_fv {
namespace euler {

template <typename T>
void evaluateInterfaceTerms(const EulerOpts& opts, Fields<T>& fields, Real t, const StructuredDiscPtr disc);

}
}