#include "linear_system/large_matrix_petsc.h"
#include "linear_system/sparsity_pattern.h"
#include <petscvec.h>

namespace structured_fv {
namespace linear_system {

LargeMatrixPetsc::LargeMatrixPetsc(const std::string& name, LargeMatrixOptsPetsc opts, std::shared_ptr<SparsityPattern> sparsity_pattern) :
  LargeMatrix(sparsity_pattern->getNumOwnedDofs(), sparsity_pattern->getNumOwnedDofs(), sparsity_pattern),
  m_opts(opts)
  //m_owned_dof_to_local(sparsity_pattern->getOwnedToLocalInfo()),
  //m_ghost_dofs_to_local(sparsity_pattern->getGhostLocalIndices())
{  
  bool is_symmetric_mattype = false;
  if (opts.petsc_opts.count("mat_type") == 0)
  {
    opts.petsc_opts["mat_type"] = opts.is_value_symmetric ? "sbaij" : "aij";
  } else
  {
    if (opts.petsc_opts["mat_type"] == "sbaij")
      is_symmetric_mattype = true;
  }  
  setPetscOptions(opts);


  const auto& ghost_global_dofs = sparsity_pattern->getGhostGlobalIndices();
    
  //VecCreate(PETSC_COMM_WORLD, &m_x);  
  VecCreateGhost(PETSC_COMM_WORLD, getNLocal(), PETSC_DECIDE, ghost_global_dofs.size(), ghost_global_dofs.data(), &m_x);
  PetscObjectSetName((PetscObject)m_x, (opts.opts_prefix + "_Solution").c_str());
  VecSetOptionsPrefix(m_x, opts.opts_prefix.c_str());
  //VecSetSizes(m_x, getNLocal(), PETSC_DECIDE);
  VecSetFromOptions(m_x);

  VecCreate(PETSC_COMM_WORLD, &m_b);
  //VecCreateGhost(PETSC_COMM_WORLD, getMLocal(), PETSC_DECIDE, ghost_global_dofs.size(), ghost_global_dofs.data(), &m_b);
  PetscObjectSetName((PetscObject)m_b, (opts.opts_prefix + "_rhs").c_str());
  VecSetOptionsPrefix(m_b, opts.opts_prefix.c_str());
  VecSetSizes(m_b, getMLocal(), PETSC_DECIDE);
  VecSetFromOptions(m_b);

  MatCreate(PETSC_COMM_WORLD, &m_A);
  PetscObjectSetName((PetscObject)m_A, (opts.opts_prefix + "_matrix").c_str());
  MatSetSizes(m_A, getMLocal(), getNLocal(), PETSC_DECIDE, PETSC_DECIDE);
  MatSetOptionsPrefix(m_A, opts.opts_prefix.c_str());
  MatSetFromOptions(m_A);

  const PetscInt *diagonal_counts = nullptr, *offproc_counts=nullptr;
  const PetscInt *diagonal_counts_sym = nullptr, *offproc_counts_sym=nullptr;
  if (is_symmetric_mattype)
  {
    diagonal_counts_sym = sparsity_pattern->getDiagonalCountsSym().data();
    offproc_counts_sym  = sparsity_pattern->getOffProcCountsSym().data();
  } else
  {
    diagonal_counts = sparsity_pattern->getDiagonalCounts().data();
    offproc_counts  = sparsity_pattern->getOffProcCounts().data();
  }

  MatXAIJSetPreallocation(m_A, 1, diagonal_counts, offproc_counts,
                                  diagonal_counts_sym, offproc_counts_sym);

  //TODO: what happens if MAT_SYMMETRIC is true but a non-symmetric matrix format is used?
  if (opts.is_value_symmetric)
    MatSetOption(m_A, MAT_SYMMETRIC, PETSC_TRUE);
  else if (opts.is_structurally_symmetric)
    MatSetOption(m_A, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE);

  if (opts.is_value_symmetric || opts.is_structurally_symmetric)
    MatSetOption(m_A, MAT_SYMMETRY_ETERNAL, PETSC_TRUE);

  MatSetOption(m_A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
  //MatSetOption(m_A, MAT_NEW_NONZERO_LOCATIONS,      PETSC_FALSE);
  //MatSetOption(m_A, MAT_NEW_NONZERO_LOCATION_ERR,   PETSC_TRUE);
  MatSetOption(m_A, MAT_IGNORE_OFF_PROC_ENTRIES, PETSC_TRUE);
  // TODO: MAT_SORTED_FULL for preallocation
  // TODO: MatSetOption(m_A, MAT_ROWS_SORTED) and similarly for MAT_COLUMNS_SORTED (see page 56 of manual)
  



  KSPCreate(PETSC_COMM_WORLD, &m_ksp);
  PetscObjectSetName((PetscObject)m_ksp, (opts.opts_prefix + "_ksp").c_str());
  KSPSetOperators(m_ksp, m_A, m_A);
  KSPSetOptionsPrefix(m_ksp, opts.opts_prefix.c_str());
  KSPSetFromOptions(m_ksp);

  KSPGetPC(m_ksp, &m_pc);
  PetscObjectSetName((PetscObject)m_pc, (opts.opts_prefix + "_pc").c_str());
  PCSetOptionsPrefix(m_pc, opts.opts_prefix.c_str());
  PCSetFromOptions(m_pc);

  if (opts.factor_in_place)
    PCFactorSetUseInPlace(m_pc, PETSC_TRUE);
}

LargeMatrixPetsc::LargeMatrixPetsc(std::shared_ptr<SparsityPattern> sparsity_pattern) :
  LargeMatrix(sparsity_pattern->getNumOwnedDofs(), sparsity_pattern->getNumOwnedDofs(), sparsity_pattern)
{}



LargeMatrixPetsc::~LargeMatrixPetsc()
{
  if (m_x)
    VecDestroy(&m_x);

  if (m_b)
    VecDestroy(&m_b);

  if (m_b)
    MatDestroy(&m_A);

  if (m_ksp)
    KSPDestroy(&m_ksp);
  // KSP manages the PC, no need to destroy it explicitly
}

std::shared_ptr<LargeMatrix> LargeMatrixPetsc::clone()
{
  auto copy = std::make_shared<LargeMatrixPetsc>(getSparsityPattern());
  MatDuplicate(m_A, MAT_SHARE_NONZERO_PATTERN, &(copy->m_A));

  return copy;
}

AssemblerPtr<LargeMatrixPetsc> LargeMatrixPetsc::getAssembler(disc::StructuredDiscPtr disc) const
{
  return std::make_shared<Assembler<LargeMatrixPetsc>>(disc, *this);
}

void LargeMatrixPetsc::finishMatrixAssembly_impl()
{
  MatAssemblyBegin(m_A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m_A, MAT_FINAL_ASSEMBLY);
  //MatView(m_A, PETSC_VIEWER_STDOUT_WORLD);
}


void LargeMatrixPetsc::factor_impl()
{
  if (m_new_matrix_values)
  {
    PCSetReusePreconditioner(m_pc, PETSC_FALSE);
    PCSetUp(m_pc);
    PCSetReusePreconditioner(m_pc, PETSC_TRUE);
  }
}


void LargeMatrixPetsc::solve_impl(ConstVectorType& b, VectorType& x)
{
  //assertAlways(b.shape()[0] == m_owned_dof_to_local.size() + m_ghost_dofs_to_local.size(), "vector must be a local (owned + ghost) vector");
  //assertAlways(x.shape()[0] == m_owned_dof_to_local.size() + m_ghost_dofs_to_local.size(), "vector must be a local (owned + ghost) vector");
  if (m_new_matrix_values)
  {
    KSPSetUp(m_ksp);
    m_new_matrix_values = false;
  }

  copyVec(b, m_b, false);
  copyVec(x, m_x, true);

  KSPSolve(m_ksp, m_b, m_x);

  
  KSPConvergedReason reason;
  KSPGetConvergedReason(m_ksp, &reason);

  if (reason < 0)
  {
    const char* reason_str;
    KSPGetConvergedReasonString(m_ksp, &reason_str);
    throw std::runtime_error(std::string("KSP failed to converge for reason: ") + reason_str);
  }
  
  copyVec(m_x, x, true);
}

void LargeMatrixPetsc::matVec_impl(ConstVectorType& x, VectorType& b)
{
  //assertAlways(b.shape()[0] == m_owned_dof_to_local.size() + m_ghost_dofs_to_local.size(), "vector must be a local (owned + ghost) vector");
  //assertAlways(x.shape()[0] == m_owned_dof_to_local.size() + m_ghost_dofs_to_local.size(), "vector must be a local (owned + ghost) vector");
  
  copyVec(x, m_x, true);
  MatMult(m_A, m_x, m_b);
  copyVec(m_b, b, false);
}

void LargeMatrixPetsc::axpy_impl(Real alpha, std::shared_ptr<LargeMatrix> x)
{
  auto x_petsc = std::dynamic_pointer_cast<LargeMatrixPetsc>(x);
  assert(x_petsc);
  MatAXPY(m_A, alpha, x_petsc->m_A, SAME_NONZERO_PATTERN);
}

AssemblerBasePtr LargeMatrixPetsc::getAssembler_impl(disc::StructuredDiscPtr disc) const
{
  return getAssembler(disc);
}


void LargeMatrixPetsc::copyVec(ConstVectorType& vec_in, Vec vec_out, bool set_ghosts)
{
  Vec vec_out_local;
  if (set_ghosts)
    VecGhostGetLocalForm(vec_out, &vec_out_local);
  else
    vec_out_local = vec_out;


  PetscScalar* vec_p;
  VecGetArray(vec_out_local, &vec_p);
  for (GlobalDof i=0; i < vec_in.extent(0); ++i)
    vec_p[i] = vec_in(i);
  /*
  for (size_t owned_dof=0; owned_dof < m_owned_dof_to_local.size(); ++owned_dof)
  {
    auto local_dof = m_owned_dof_to_local[owned_dof];
    vec_p[owned_dof] = vec_in[local_dof];
  }

  if (set_ghosts)
  {
    for (size_t ghost_idx=0; ghost_idx < m_ghost_dofs_to_local.size(); ++ghost_idx)
    {
      auto local_dof = m_ghost_dofs_to_local[ghost_idx];
      vec_p[ghost_idx + m_owned_dof_to_local.size()] = vec_in[local_dof];
    }
  }
  */

  VecRestoreArray(vec_out_local, &vec_p);
  if (set_ghosts)
    VecGhostRestoreLocalForm(vec_out, &vec_out_local);
}


void LargeMatrixPetsc::copyVec(const Vec vec_in, VectorType& vec_out, bool set_ghosts)
{
  Vec vec_in_local;
  if (set_ghosts)
  {
    VecGhostGetLocalForm(vec_in, &vec_in_local);
    //VecGhostUpdateBegin(vec_in, INSERT_VALUES, SCATTER_FORWARD);
  } else
    vec_in_local = vec_in;
  
  const PetscScalar* vec_p;
  //int len = vec_out.shape()[0];
  VecGetArrayRead(vec_in_local, &vec_p);
  //std::copy(vec_p, vec_p + len, &(vec_out[0]));
  for (GlobalDof i=0; i < vec_out.extent(0); ++i)
    vec_out(i) = vec_p[i];

/*
  for (size_t owned_dof=0; owned_dof < m_owned_dof_to_local.size(); ++owned_dof)
  {
    auto local_dof = m_owned_dof_to_local[owned_dof];
    vec_out[local_dof] = vec_p[owned_dof];
  }

  if (set_ghosts)
  {
    VecGhostUpdateEnd(vec_in, INSERT_VALUES, SCATTER_FORWARD);
    for (size_t ghost_idx=0; ghost_idx < m_ghost_dofs_to_local.size(); ++ghost_idx)
    {
      auto local_dof = m_ghost_dofs_to_local[ghost_idx];
      vec_out[local_dof] = vec_p[ghost_idx + m_owned_dof_to_local.size()];
    }
  }
*/
  VecRestoreArrayRead(vec_in_local, &vec_p);
  if (set_ghosts)
    VecGhostRestoreLocalForm(vec_in, &vec_in_local);
}



void setPetscOptions(const LargeMatrixOptsPetsc& opts)
{
  for (auto& p : opts.petsc_opts)
  {
    std::string full_key = std::string("-") + opts.opts_prefix + p.first;
    PetscOptionsSetValue(NULL, full_key.c_str(), p.second.c_str());
  }
}

void setPetscGlobalOption(const std::string& key, const std::string& val)
{
  std::string full_key = std::string("-") + key;
  PetscOptionsSetValue(NULL, full_key.c_str(), val.c_str());
}

}
}