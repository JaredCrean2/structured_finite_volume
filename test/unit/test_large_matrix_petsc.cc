#include "gtest/gtest.h"
#include "linear_system/large_matrix.h"
#include "linear_system/large_matrix_petsc.h"
#include "linear_system/sparsity_pattern.h"
#include "petscsys.h"
#include "utils/project_defs.h"

namespace {

using namespace structured_fv;

using Mat3 = Matrix<Real, 3, 3>;

class SparsityPatternTest : public linear_system::SparsityPattern
{
  public:
    explicit SparsityPatternTest(int size) :
      m_local(size, size),
      m_offproc(size, size)
      //m_owned_to_local_dofs(size),
      //m_local_dof_to_global(size)
    {
      //for (int i=0; i < size; ++i)
      //  m_owned_to_local_dofs[i] = i;
    }

    PetscInt getNumOwnedDofs() const override { return m_local.size(); }

    const std::vector<PetscInt>& getDiagonalCounts() override { return m_local; }

    const std::vector<PetscInt>& getDiagonalCountsSym() override { return m_local; }

    const std::vector<PetscInt>& getOffProcCounts() override { return m_offproc; }
        
    const std::vector<PetscInt>& getOffProcCountsSym() override { return m_offproc; }

    const std::vector<PetscInt>& getGhostGlobalIndices() override { return m_ghost_dofs; }

    //const std::vector<PetscInt>& getGhostLocalIndices() override { return m_ghost_dofs; }

    //const std::vector<PetscInt>& getOwnedToLocalInfo() override { return m_owned_to_local_dofs; };

    //const std::vector<PetscInt>& getLocalToGlobalDofs() override { return m_local_dof_to_global; }


  private:
    std::vector<PetscInt> m_local;
    std::vector<PetscInt> m_offproc;
    std::vector<PetscInt> m_ghost_dofs;
    //std::vector<PetscInt> m_owned_to_local_dofs;
    //std::vector<PetscInt> m_local_dof_to_global;
};

linear_system::LargeMatrixOptsPetsc get_options()
{
  linear_system::LargeMatrixOptsPetsc opts;
  opts.is_structurally_symmetric = false;
  opts.is_value_symmetric        = false;
  opts.factor_in_place           = false;
  opts.petsc_opts["ksp_atol"] = "1e-15";
  opts.petsc_opts["ksp_rtol"] = "1e-50";
  opts.petsc_opts["ksp_monitor"] = "";

  return opts;
}

}


TEST(LargeMatrixPetsc, GeneralSolve)
{
  //linear_system::setPetscGlobalOption("on_error_abort", "true");
  //PetscPushErrorHandler(&PetscAbortErrorHandler, nullptr);
  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  auto opts = get_options();
  auto sparsity_pattern = std::make_shared<SparsityPatternTest>(3);
  linear_system::LargeMatrixPetsc mat("A", opts, sparsity_pattern);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::array<GlobalDof, 3> dofs{0 + 3*comm_rank, 1 + 3*comm_rank, 2 + 3*comm_rank};
  Mat3 jac({1, 2, 3,
            4, 5, 6,
            8, 8, 9});
  linear_system::VectorType x("x", 3);
  linear_system::VectorType b("b", 3);
  b(0) = 3;
  b(1) = 6;
  b(2) = 9;
  Vec3<Real> x_ex {0, 0, 1};

  mat.assembleValues(dofs, dofs, jac);
  mat.finishMatrixAssembly();
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);
}


TEST(LargeMatrixPetsc, AssembleValuesAdditive)
{
  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  auto opts = get_options();
  auto sparsity_pattern = std::make_shared<SparsityPatternTest>(3);
  linear_system::LargeMatrixPetsc mat("A", opts, sparsity_pattern);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::array<GlobalDof, 3> dofs{0 + 3*comm_rank, 1 + 3*comm_rank, 2 + 3*comm_rank};
  Mat3 vals({1, 2, 3,
             4, 5, 3,
             4, 4, 9});
  Mat3 vals2({0, 0, 0,
              0, 0, 3,
              4, 4, 0});
  linear_system::VectorType x("x", 3);
  linear_system::VectorType b("b", 3);
  b(0) = 3;
  b(1) = 6;
  b(2) = 9;
  Vec3<Real> x_ex = {0, 0, 1};

  mat.assembleValues(dofs, dofs, vals);
  mat.assembleValues(dofs, dofs, vals2);
  mat.finishMatrixAssembly();
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);
}


TEST(LargeMatrixPetsc, AssembleValuesIgnore)
{
  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  auto opts = get_options();
  auto sparsity_pattern = std::make_shared<SparsityPatternTest>(3);
  linear_system::LargeMatrixPetsc mat("A", opts, sparsity_pattern);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::array<GlobalDof, 4> dofs{0 + 3*comm_rank, 1 + 3*comm_rank, 2 + 3*comm_rank, -1};
  Matrix<Real, 4> vals({1, 2, 3,       666,
                        4, 5, 6,       666,
                        8, 8, 9,       666,
                        666, 666, 666, 666});
  linear_system::VectorType x("x", 3);
  linear_system::VectorType b("b", 3);
  b(0) = 3;
  b(1) = 6;
  b(2) = 9;
  Vec3<Real> x_ex = {0, 0, 1};

  mat.assembleValues(dofs, dofs, vals);
  mat.finishMatrixAssembly();
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);
}


TEST(LargeMatrixPetsc, ZeroMatrix)
{
  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  auto opts = get_options();
  auto sparsity_pattern = std::make_shared<SparsityPatternTest>(3);
  linear_system::LargeMatrixPetsc mat("A", opts, sparsity_pattern);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::array<GlobalDof, 3> dofs{0 + 3*comm_rank, 1 + 3*comm_rank, 2 + 3*comm_rank};
  Mat3 vals({1, 2, 3,
             4, 5, 6,
             8, 8, 9});
  linear_system::VectorType x("x", 3);
  linear_system::VectorType b("b", 3);
  b(0) = 3;
  b(1) = 6;
  b(2) = 9;
  Vec3<Real> x_ex = {0, 0, 1};


  mat.assembleValues(dofs, dofs, vals);
  mat.finishMatrixAssembly();
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);


  mat.zeroMatrix();
  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
      vals(i, j) *= 2;

  mat.assembleValues(dofs, dofs, vals);
  mat.finishMatrixAssembly();
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], 0.5 * x_ex[i], 1e-13);
}


TEST(LargeMatrixPetsc, FactorInPlace)
{
  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  auto opts = get_options();
  opts.factor_in_place  = true;
  opts.petsc_opts["ksp_type"] = "preonly";
  opts.petsc_opts.erase("ksp_monitor");
  opts.petsc_opts.erase("ksp_atol");
  opts.petsc_opts.erase("ksp_rtol");
  auto sparsity_pattern = std::make_shared<SparsityPatternTest>(3);
  linear_system::LargeMatrixPetsc mat("A", opts, sparsity_pattern);


  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::array<GlobalDof, 3> dofs{0 + 3*comm_rank, 1 + 3*comm_rank, 2 + 3*comm_rank};
  Mat3 vals({1, 2, 3,
             4, 5, 6,
             8, 8, 9});
  linear_system::VectorType x("x", 3);
  linear_system::VectorType b("b", 3);
  b(0) = 3;
  b(1) = 6;
  b(2) = 9;
  Vec3<Real> x_ex = {0, 0, 1};

  mat.assembleValues(dofs, dofs, vals);
  mat.finishMatrixAssembly();
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-13);
}


TEST(LargeMatrixPetsc, SPD)
{
  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  
  auto opts = get_options();
  opts.is_structurally_symmetric = true;
  opts.is_value_symmetric        = true;
  opts.petsc_opts["pc_type"] = "bjacobi";
  opts.petsc_opts["pc_sub_type"] = "icc";
  opts.petsc_opts["mat_type"] = "aij";  //TODO: need to update assembly procedure to filter
                                        //      out below-diagonal entries
  auto sparsity_pattern = std::make_shared<SparsityPatternTest>(3);
  linear_system::LargeMatrixPetsc mat("A", opts, sparsity_pattern);

  EXPECT_EQ(mat.getMLocal(), 3);
  EXPECT_EQ(mat.getNLocal(), 3);

  std::array<GlobalDof, 3> dofs{0 + 3*comm_rank, 1 + 3*comm_rank, 2 + 3*comm_rank};
  Mat3 vals({10, 2, 3,
             2, 12, 5,
             3, 5, 20});
  linear_system::VectorType x("x", 3);
  linear_system::VectorType b("b", 3);
  b(0) = 1;
  b(1) = 2;
  b(2) = 3;
  
  Vec3<Real> x_ex = {0.043027, 0.111276, 0.115727};

  mat.assembleValues(dofs, dofs, vals);
  mat.finishMatrixAssembly();
  mat.factor();
  mat.solve(b, x);

  for (int i=0; i < 3; ++i)
    EXPECT_NEAR(x[i], x_ex[i], 1e-5);
}