#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_INTERFACE_TERM_JAC_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_EULER_INTERFACE_TERM_JAC_H

#include "utils/project_defs.h"
#include "euler_model.h"
#include "linear_system/assembler_base.h"

namespace structured_fv {
namespace euler {


void evaluateInterfaceTermsJac(const EulerOpts& opts, Fields<Real>& fields, Real t, const StructuredDiscPtr disc,
                               linear_system::AssemblerBasePtr assembler);

}
}

#endif