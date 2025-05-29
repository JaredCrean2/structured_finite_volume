#include "gtest/gtest.h"
#include "disc/disc_block.h"
#include "disc/disc_interface.h"
#include "physics/common/slope_limiters.h"
#include "physics/euler/euler_model.h"
#include "nonlinear_solvers/explicit_euler.h"
#include "mesh/structured_mesh.h"
#include "physics/euler/typedefs.h"
#include "utils/math.h"
#include "multi_block_fixture.h"

namespace {

using namespace structured_fv;

class EulerTester : public ::testing::Test
{
  public:
    EulerTester() :
      spec(1, 1, m_num_bc_ghost_cells)
    {
      spec.blocks(0, 0) = mesh::MeshBlockSpec(3, 4, 0, [](Real x, Real y) { return std::array<Real, 2>{x, y}; });
      m_mesh = std::make_shared<mesh::StructuredMesh>(spec);
      m_disc = std::make_shared<disc::StructuredDisc>(m_mesh, m_num_bc_ghost_cells, m_dofs_per_cell);
    }

    void setup(const euler::EulerOpts& opts, euler::Fxyt bc_func,
               euler::Fxyt source_func)
    {
      std::vector<euler::Fxyt> bc_funcs{bc_func, bc_func, bc_func, bc_func};
      m_euler_model = std::make_shared<euler::EulerModel>(opts, m_disc, bc_funcs, source_func);
    }



    UInt m_dofs_per_cell = 4;
    int m_num_bc_ghost_cells = 1;
    mesh::MeshSpec spec;
    std::shared_ptr<mesh::StructuredMesh> m_mesh;
    std::shared_ptr<disc::StructuredDisc> m_disc;
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
      m_adv_model1 = std::make_shared<euler::EulerModel>(
          opts, m_disc1, bc_funcs, source_func);

      std::vector<euler::Fxyt> bc_funcs2{bc_func, bc_func, bc_func, bc_func,
                                         bc_func, bc_func, bc_func, bc_func};

      m_adv_model2 = std::make_shared<euler::EulerModel>(
          opts, m_disc2, bc_funcs2, source_func);
    }

    std::shared_ptr<euler::EulerModel> m_adv_model1;
    std::shared_ptr<euler::EulerModel> m_adv_model2;
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
  auto bc_func = [](Real x, Real y, Real t){ return std::array<Real, 4>{1, 0, 0, 10000}; };
  auto source_func = [](Real x, Real y, Real t){ return std::array<Real, 4>{x+y+t, x+y+t+1, x+y+t+2, x+y+t+3}; };
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
  auto bc_func = [](Real x, Real y, Real t){ return std::array<Real, 4>{1, 2, 3, 40000}; };
  auto source_func = [](Real x, Real y, Real t){ return std::array<Real, 4>{0, 0, 0, 0}; };
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
  auto u_ex     = [&](Real x, Real y, Real t) { return std::array<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return std::array<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };
  auto bc_func  = [&](Real x, Real y, Real t) { return u_ex(x, y, t); };
  auto u0       = [&](Real x, Real y) { return u_ex(x, y, 0); };

  setup(opts, bc_func, src_term);

  m_solution1->set(u0);
  m_solution2->set(u0);

  m_adv_model1->evaluateRhs(m_solution1, t, m_residual1);
  m_adv_model2->evaluateRhs(m_solution2, t, m_residual2);

  test_residuals_equal(m_residual1, m_residual2);
}

TEST_F(EulerTesterMultiBlock, ConsistencySecondOrder) {
  Real t = 1.0;
  euler::EulerOpts opts;
  opts.limiter = common::SlopeLimiter::VanLeer;

  Real u = 40;
  Real T = 298;
  auto u_ex     = [&](Real x, Real y, Real t) { return std::array<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return std::array<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };
  auto bc_func  = [&](Real x, Real y, Real t) { return u_ex(x, y, t); };
  auto u0       = [&](Real x, Real y) { return u_ex(x, y, 0); };

  setup(opts, bc_func, src_term);

  m_solution1->set(u0);
  m_solution2->set(u0);

  m_adv_model1->evaluateRhs(m_solution1, t, m_residual1);
  m_adv_model2->evaluateRhs(m_solution2, t, m_residual2);

  test_residuals_equal(m_residual1, m_residual2);
}


