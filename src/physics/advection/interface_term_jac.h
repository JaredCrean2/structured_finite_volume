#ifndef STRUCTURED_FINITE_VOLUME_PHYSICS_ADVECTION_INTERFACE_TERM_JAC_H
#define STRUCTURED_FINITE_VOLUME_PHYSICS_ADVECTION_INTERFACE_TERM_JAC_H

#include "utils/project_defs.h"
#include "advection_model.h"
#include "interface_term.h"
#include "linear_system/assembler.h"
#include "linear_system/large_matrix_dense.h"
#include "linear_system/large_matrix_petsc.h"

#include "utils/matrix.h"

namespace structured_fv {
namespace advection {


void evaluateInterfaceTermJac(const AdvectionOpts& opts, Fields<Real> fields, Real t, StructuredDiscPtr disc,
                              linear_system::AssemblerBasePtr assembler);

}
}

#endif