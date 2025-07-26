#include "gtest/gtest.h"
#include "disc/disc_block.h"
#include "disc/disc_interface.h"
#include "physics/common/slope_limiters.h"
#include "physics/euler/euler_model.h"
#include "nonlinear_solvers/explicit_euler.h"
#include "mesh/structured_mesh.h"
#include "physics/euler/flux_function_enums.h"
#include "physics/euler/typedefs.h"
#include "linear_system/large_matrix_dense.h"
#include "linear_system/assembler.h"
#include "utils/math.h"
#include "multi_block_fixture.h"
#include "jacobian.h"

namespace {

using namespace structured_fv;
using MatrixType = linear_system::LargeMatrixDense;
using AssemblerType = linear_system::Assembler<MatrixType>;

class EulerTester : public ::testing::Test
{
  public:
    EulerTester() :
      spec(1, 1, m_num_bc_ghost_cells)
    {
      spec.blocks(0, 0) = mesh::MeshBlockSpec(3, 4, 0, [](Real x, Real y) { return FixedVec<Real, 2>{x, y}; });
      m_mesh = std::make_shared<mesh::StructuredMesh>(spec);
      m_disc = std::make_shared<disc::StructuredDisc>(m_mesh, m_num_bc_ghost_cells, m_dofs_per_cell);

      linear_system::LargeMatrixOpts opts;
      m_matrix = std::make_shared<MatrixType>("A", m_disc->getNumDofs(), m_disc->getNumDofs(), opts);
      m_assembler = std::make_shared<AssemblerType>(m_disc, *m_matrix);        
    }

    void setup(const euler::EulerOpts& opts, euler::Fxyt bc_func,
               euler::Fxyt source_func)
    {
      std::vector<euler::Fxyt> bc_funcs{bc_func, bc_func, bc_func, bc_func};
      m_euler_model = std::make_shared<euler::EulerModel>(opts, m_disc, bc_funcs, source_func);
    }


    UInt m_dofs_per_cell = 4;
    int m_num_bc_ghost_cells = 2;
    mesh::MeshSpec spec;
    std::shared_ptr<mesh::StructuredMesh> m_mesh;
    std::shared_ptr<disc::StructuredDisc> m_disc;
    std::shared_ptr<MatrixType> m_matrix;
    std::shared_ptr<AssemblerType> m_assembler;    
    std::shared_ptr<euler::EulerModel> m_euler_model;
};

class EulerTesterMultiBlock : public test_utils::MultiBlockFixture,
                              public ::testing::Test
{
  public:
    void setup(const euler::EulerOpts &opts, euler::Fxyt bc_func,
               euler::Fxyt source_func)
    {
      setup_disc(4);
      std::vector<euler::Fxyt> bc_funcs{bc_func, bc_func, bc_func, bc_func};
      m_euler_model1 = std::make_shared<euler::EulerModel>(
          opts, m_disc1, bc_funcs, source_func);

      linear_system::LargeMatrixOpts matrix_opts;
      m_matrix1 = std::make_shared<MatrixType>("A", m_disc1->getNumDofs(), m_disc1->getNumDofs(), matrix_opts);
      m_assembler1 = std::make_shared<AssemblerType>(m_disc1, *m_matrix1);           

      std::vector<euler::Fxyt> bc_funcs2{bc_func, bc_func, bc_func, bc_func,
                                         bc_func, bc_func, bc_func, bc_func};

      m_euler_model2 = std::make_shared<euler::EulerModel>(
          opts, m_disc2, bc_funcs2, source_func);

      m_matrix2 = std::make_shared<MatrixType>("A", m_disc2->getNumDofs(), m_disc2->getNumDofs(), matrix_opts);
      m_assembler2 = std::make_shared<AssemblerType>(m_disc2, *m_matrix2);
    }

    std::shared_ptr<MatrixType> m_matrix1;
    std::shared_ptr<AssemblerType> m_assembler1;
    std::shared_ptr<euler::EulerModel> m_euler_model1;

    std::shared_ptr<MatrixType> m_matrix2;
    std::shared_ptr<AssemblerType> m_assembler2;
    std::shared_ptr<euler::EulerModel> m_euler_model2;
};
}

TEST_F(EulerTester, SourceTerm)
{
  auto solution = std::make_shared<disc::DiscVector<Real>>(m_disc, "solution");
  auto residual = std::make_shared<disc::DiscVector<Real>>(m_disc, "residual");
  solution->set({1, 0, 0, 10000});
  residual->set(0);

  Real t = 1.0;
  euler::EulerOpts opts;
  auto bc_func = [](Real x, Real y, Real t){ return FixedVec<Real, 4>{1, 0, 0, 10000}; };
  auto source_func = [](Real x, Real y, Real t){ return FixedVec<Real, 4>{x+y+t, x+y+t+1, x+y+t+2, x+y+t+3}; };
  setup(opts, bc_func, source_func);
  m_euler_model->evaluateRhs(solution, t, residual);

  const disc::StructuredBlock& block = m_disc->getBlock(0);
  const auto& coords = m_disc->getCoordField()->getData(0);
  const auto& dof_nums = m_disc->getDofNumbering()->getData(0);
  for (UInt i : block.getOwnedCells().getXRange())
    for (UInt j : block.getOwnedCells().getYRange())
    {
      Vec2<Real> x = disc::computeCellCentroid(coords, i, j);
      for (UInt k=0; k < 4; ++k)
      {
        GlobalDof dof = dof_nums(i, j, k);
        EXPECT_NEAR((*residual)(dof), source_func(x[0], x[1], t)[k], 1e-13);
      }
    }
}

//TODO: write tests where qL and qR can be connected by a single wave
//      both HLLE and Roe fluxes should be exact for this case

//TODO: check for nans
TEST_F(EulerTester, BoundaryTerm)
{
  auto solution = std::make_shared<disc::DiscVector<Real>>(m_disc, "solution");
  auto residual = std::make_shared<disc::DiscVector<Real>>(m_disc, "residual");
  solution->set({1, 2, 3, 40000});
  residual->set(0);

  Real t = 1.0;
  euler::EulerOpts opts;
  auto bc_func = [](Real x, Real y, Real t){ return FixedVec<Real, 4>{1, 2, 3, 40000}; };
  auto source_func = [](Real x, Real y, Real t){ return FixedVec<Real, 4>{0, 0, 0, 0}; };
  setup(opts, bc_func, source_func);
  m_euler_model->evaluateRhs(solution, t, residual);

  const disc::StructuredBlockInterface& left_iface = m_disc->getBlockInterface(3);
  const auto& dof_nums = m_disc->getDofNumbering()->getData(left_iface.getBlockIdL());
  for (UInt i : left_iface.getOwnedBoundaryCellsL().getXRange())
    for (UInt j : left_iface.getOwnedBoundaryCellsL().getYRange())
      for (UInt k=0; k < 4; ++k)
        EXPECT_NEAR((*residual)(dof_nums(i, j, k)), 0.0, 1e-13);

  for (UInt iface_id=0; iface_id < 3; ++iface_id)
  {
    const disc::StructuredBlockInterface& iface = m_disc->getBlockInterface(iface_id);
    const auto& dof_nums = m_disc->getDofNumbering()->getData(iface.getBlockIdL());

    for (UInt i : iface.getOwnedBoundaryCellsL().getXRange())
      for (UInt j : iface.getOwnedBoundaryCellsL().getYRange())
      {
        if (!in(left_iface.getOwnedBoundaryCellsL(), i, j))
        {
          EXPECT_NEAR((*residual)(dof_nums(i, j, 0)), 0, 1e-13);
        }
      }
  }
}

TEST_F(EulerTesterMultiBlock, ConsistencyFirstOrder) {
  Real t = 1.0;
  euler::EulerOpts opts;

  Real u = 40;
  Real T = 298;
  auto u_ex     = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };
  auto bc_func  = [&](Real x, Real y, Real t) { return u_ex(x, y, t); };
  auto u0       = [&](Real x, Real y) { return u_ex(x, y, 0); };

  setup(opts, bc_func, src_term);

  m_solution1->set(u0);
  m_solution2->set(u0);

  m_euler_model1->evaluateRhs(m_solution1, t, m_residual1);
  m_euler_model2->evaluateRhs(m_solution2, t, m_residual2);

  test_residuals_equal(m_residual1, m_residual2);
}

TEST_F(EulerTesterMultiBlock, ConsistencySecondOrder) {
  Real t = 1.0;
  euler::EulerOpts opts;
  opts.limiter = common::SlopeLimiter::VanLeer;

  Real u = 40;
  Real T = 298;
  auto u_ex     = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };
  auto bc_func  = [&](Real x, Real y, Real t) { return u_ex(x, y, t); };
  auto u0       = [&](Real x, Real y) { return u_ex(x, y, 0); };

  setup(opts, bc_func, src_term);

  m_solution1->set(u0);
  m_solution2->set(u0);

  m_euler_model1->evaluateRhs(m_solution1, t, m_residual1);
  m_euler_model2->evaluateRhs(m_solution2, t, m_residual2);

  test_residuals_equal(m_residual1, m_residual2);
}


TEST_F(EulerTester, JacobianFirstOrder)
{
  Real t = 1.0;
  euler::EulerOpts opts;
  opts.flux = euler::FluxFunction::LLF;
  opts.limiter = common::SlopeLimiter::FirstOrder; //common::SlopeLimiter::VanAlba;

  Real u = 40;
  Real T = 298;
  auto u_ex     = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };
  auto bc_func  = [&](Real x, Real y, Real t) { return u_ex(x, y, t); };
  auto u0       = [&](Real x, Real y) { return u_ex(x, y, 0); };

  setup(opts, bc_func, src_term);

  auto solution = std::make_shared<disc::DiscVector<Real>>(m_disc, "solution");
  auto v = std::make_shared<disc::DiscVector<Real>>(m_disc, "v");
  auto residual1 = std::make_shared<disc::DiscVector<Real>>(m_disc, "residual1");
  auto dRdq_v = std::make_shared<disc::DiscVector<Real>>(m_disc, "dR/dq v");

  solution->set(u0);
  v->set(0);
  m_euler_model->evaluateJacobian(solution, t, residual1, m_assembler);
  for (GlobalDof dof=0; dof < m_disc->getNumDofs(); ++dof)
  {
    (*v)(dof) = 1;
    m_euler_model->computeJacVecProduct(solution, t, v, dRdq_v);
    (*v)(dof) = 0;

    for (GlobalDof dof2=0; dof2 < m_disc->getNumDofs(); ++dof2)
    {
      std::cout << "diff = " << (*m_matrix)(dof2, dof) - (*dRdq_v)(dof2) << std::endl;
      EXPECT_DOUBLE_EQ_CUSTOM((*m_matrix)(dof2, dof), (*dRdq_v)(dof2), 300);
    }
  }
}

TEST_F(EulerTester, JacobianSecondOrder)
{
  Real t = 1.0;
  euler::EulerOpts opts;
  opts.flux = euler::FluxFunction::LLF;
  opts.limiter = common::SlopeLimiter::VanAlba;

  Real u = 40;
  Real T = 298;
  auto u_ex     = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };
  auto bc_func  = [&](Real x, Real y, Real t) { return u_ex(x, y, t); };
  auto u0       = [&](Real x, Real y) { return u_ex(x, y, 0); };

  setup(opts, bc_func, src_term);

  auto solution = std::make_shared<disc::DiscVector<Real>>(m_disc, "solution");
  auto v = std::make_shared<disc::DiscVector<Real>>(m_disc, "v");
  auto residual1 = std::make_shared<disc::DiscVector<Real>>(m_disc, "residual1");
  auto dRdq_v = std::make_shared<disc::DiscVector<Real>>(m_disc, "dR/dq v");

  solution->set(u0);
  v->set(0);
  m_euler_model->evaluateJacobian(solution, t, residual1, m_assembler);
  for (GlobalDof dof=0; dof < m_disc->getNumDofs(); ++dof)
  {
    (*v)(dof) = 1;
    m_euler_model->computeJacVecProduct(solution, t, v, dRdq_v);
    (*v)(dof) = 0;

    for (GlobalDof dof2=0; dof2 < m_disc->getNumDofs(); ++dof2)
    {
      //std::cout << "\nentry " << dof2 << ", " << dof << std::endl;
      //std::cout << "diff = " << (*m_matrix)(dof2, dof) - (*dRdq_v)(dof2) << std::endl;
      EXPECT_DOUBLE_EQ_CUSTOM((*m_matrix)(dof2, dof), (*dRdq_v)(dof2), 300);
    }
  }
}

TEST_F(EulerTesterMultiBlock, JacobianFirstOrder)
{
  Real t = 1.0;
  euler::EulerOpts opts;
  opts.flux = euler::FluxFunction::LLF;
  opts.limiter = common::SlopeLimiter::FirstOrder; //common::SlopeLimiter::VanAlba;

  Real u = 40;
  Real T = 298;
  auto u_ex     = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };
  auto bc_func  = [&](Real x, Real y, Real t) { return u_ex(x, y, t); };
  auto u0       = [&](Real x, Real y) { return u_ex(x, y, 0); };

  setup(opts, bc_func, src_term);

  auto solution = std::make_shared<disc::DiscVector<Real>>(m_disc2, "solution");
  auto v = std::make_shared<disc::DiscVector<Real>>(m_disc2, "v");
  auto residual1 = std::make_shared<disc::DiscVector<Real>>(m_disc2, "residual1");
  auto dRdq_v = std::make_shared<disc::DiscVector<Real>>(m_disc2, "dR/dq v");

  solution->set(u0);
  v->set(0);
  m_euler_model2->evaluateJacobian(solution, t, residual1, m_assembler2);
  for (GlobalDof dof=0; dof < m_disc2->getNumDofs(); ++dof)
  {
    (*v)(dof) = 1;
    m_euler_model2->computeJacVecProduct(solution, t, v, dRdq_v);
    (*v)(dof) = 0;

    for (GlobalDof dof2=0; dof2 < m_disc2->getNumDofs(); ++dof2)
    {
      std::cout << "entry " << dof2 << ", " << dof << std::endl;
      std::cout << "diff = " << (*m_matrix2)(dof2, dof) - (*dRdq_v)(dof2) << std::endl;
      EXPECT_DOUBLE_EQ_CUSTOM((*m_matrix2)(dof2, dof), (*dRdq_v)(dof2), 300);
    }
  }
}

TEST_F(EulerTesterMultiBlock, JacobianSecondOrder)
{
  Real t = 1.0;
  euler::EulerOpts opts;
  opts.flux = euler::FluxFunction::LLF;
  opts.limiter = common::SlopeLimiter::VanAlba;

  Real u = 40;
  Real T = 298;
  auto u_ex     = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };
  auto bc_func  = [&](Real x, Real y, Real t) { return u_ex(x, y, t); };
  auto u0       = [&](Real x, Real y) { return u_ex(x, y, 0); };

  setup(opts, bc_func, src_term);

  auto solution = std::make_shared<disc::DiscVector<Real>>(m_disc2, "solution");
  auto v = std::make_shared<disc::DiscVector<Real>>(m_disc2, "v");
  auto residual1 = std::make_shared<disc::DiscVector<Real>>(m_disc2, "residual1");
  auto dRdq_v = std::make_shared<disc::DiscVector<Real>>(m_disc2, "dR/dq v");

  solution->set(u0);
  v->set(0);
  m_euler_model2->evaluateJacobian(solution, t, residual1, m_assembler2);
  for (GlobalDof dof=0; dof < m_disc2->getNumDofs(); ++dof)
  {
    (*v)(dof) = 1;
    m_euler_model2->computeJacVecProduct(solution, t, v, dRdq_v);
    (*v)(dof) = 0;

    for (GlobalDof dof2=0; dof2 < m_disc2->getNumDofs(); ++dof2)
    {
      std::cout << "entry " << dof2 << ", " << dof << std::endl;
      std::cout << "diff = " << (*m_matrix2)(dof2, dof) - (*dRdq_v)(dof2) << std::endl;
      EXPECT_DOUBLE_EQ_CUSTOM((*m_matrix2)(dof2, dof), (*dRdq_v)(dof2), 300);
    }
  }
}