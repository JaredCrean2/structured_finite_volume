#include "disc/disc_block.h"
#include "disc/disc_interface.h"
#include "disc/disc_vector.h"
#include "mesh/structured_mesh.h"
#include "nonlinear_solvers/explicit_euler.h"
#include "physics/advection/advection_model.h"
#include "physics/common/slope_limiters.h"
#include "utils/math.h"
#include "rotation_utils.h"
#include "gtest/gtest.h"
#include "multi_block_fixture.h"

namespace {

using namespace structured_fv;

class AdvectionTester : public ::testing::Test
{
  public:
    AdvectionTester() : 
      spec(1, 1, m_num_bc_ghost_cells) 
    {
      spec.blocks(0, 0) = mesh::MeshBlockSpec(3, 4, 0, 
                           [](Real x, Real y) { return FixedVec<Real, 2>{x, y}; });
      m_mesh = std::make_shared<mesh::StructuredMesh>(spec);
      m_disc = std::make_shared<disc::StructuredDisc>(m_mesh, m_num_bc_ghost_cells, m_dofs_per_cell);
    }

    void setup(const advection::AdvectionOpts &opts, advection::Fxyt bc_func,
               advection::Fxyt source_func)
    {
      std::vector<advection::Fxyt> bc_funcs{bc_func, bc_func, bc_func, bc_func};
      m_adv_model = std::make_shared<advection::AdvectionModel>(
          opts, m_disc, bc_funcs, source_func);
    }

    UInt m_dofs_per_cell = 1;
    int m_num_bc_ghost_cells = 2;
    mesh::MeshSpec spec;
    std::shared_ptr<mesh::StructuredMesh> m_mesh;
    std::shared_ptr<disc::StructuredDisc> m_disc;
    std::shared_ptr<advection::AdvectionModel> m_adv_model;
};

class AdvectionTesterMultiBlock : public test_utils::MultiBlockFixture,
                                  public ::testing::Test
{
  public:
    void setup(const advection::AdvectionOpts &opts, advection::Fxyt bc_func,
              advection::Fxyt source_func)
    {
      setup_disc(1);
      std::vector<advection::Fxyt> bc_funcs{bc_func, bc_func, bc_func, bc_func};
      m_adv_model1 = std::make_shared<advection::AdvectionModel>(
          opts, m_disc1, bc_funcs, source_func);

      std::vector<advection::Fxyt> bc_funcs2{bc_func, bc_func, bc_func, bc_func,
                                            bc_func, bc_func, bc_func, bc_func};

      m_adv_model2 = std::make_shared<advection::AdvectionModel>(
          opts, m_disc2, bc_funcs2, source_func);
    }

    std::shared_ptr<advection::AdvectionModel> m_adv_model1;
    std::shared_ptr<advection::AdvectionModel> m_adv_model2;
};
} // namespace

TEST_F(AdvectionTester, SourceTerm)
{
  auto solution = std::make_shared<disc::DiscVector<Real>>(m_disc, "solution");
  auto residual = std::make_shared<disc::DiscVector<Real>>(m_disc, "residual");
  solution->set(0);
  residual->set(0);

  Real t = 1.0;
  advection::AdvectionOpts opts;
  opts.adv_velocity = {1, 0};
  auto bc_func = [](Real x, Real y, Real t) { return 0; };
  auto source_func = [](Real x, Real y, Real t) { return x + y + t; };
  setup(opts, bc_func, source_func);
  m_adv_model->evaluateRhs(solution, t, residual);

  const disc::StructuredBlock &block = m_disc->getBlock(0);
  const auto &coords = m_disc->getCoordField()->getData(0);
  const auto &dof_nums = m_disc->getDofNumbering()->getData(0);
  for (UInt i : block.getOwnedCells().getXRange())
    for (UInt j : block.getOwnedCells().getYRange()) {
      Vec2<Real> x = disc::computeCellCentroid(coords, i, j);
      GlobalDof dof = dof_nums(i, j, 0);
      EXPECT_NEAR((*residual)(dof), source_func(x[0], x[1], t), 1e-13);
    }
}

TEST_F(AdvectionTester, BoundaryTerm)
{
  auto solution = std::make_shared<disc::DiscVector<Real>>(m_disc, "solution");
  auto residual = std::make_shared<disc::DiscVector<Real>>(m_disc, "residual");
  solution->set(0);
  residual->set(0);

  Real t = 1.0;
  advection::AdvectionOpts opts;
  opts.adv_velocity = {1, 0};
  auto bc_func = [](Real x, Real y, Real t) { return 1; };
  auto source_func = [](Real x, Real y, Real t) { return 0; };
  setup(opts, bc_func, source_func);
  m_adv_model->evaluateRhs(solution, t, residual);

  const disc::StructuredBlockInterface &left_iface =
      m_disc->getBlockInterface(3);
  const auto &coords =
      m_disc->getCoordField()->getData(left_iface.getBlockIdL());
  const auto &dof_nums =
      m_disc->getDofNumbering()->getData(left_iface.getBlockIdL());
  Real dx = 1.0 / 3;
  Real dy = 1.0 / 4;
  Real cell_area = dx * dy;
  for (UInt i : left_iface.getOwnedBoundaryCellsL().getXRange())
    for (UInt j : left_iface.getOwnedBoundaryCellsL().getYRange()) {
      Vec2<Real> x = disc::computeCellCentroid(coords, i, j);
      Real flux = bc_func(x[0], x[1], t) * dy;
      EXPECT_NEAR((*residual)(dof_nums(i, j, 0)), flux / cell_area, 1e-13);
    }

  for (UInt iface_id = 0; iface_id < 3; ++iface_id) {
    const disc::StructuredBlockInterface &iface =
        m_disc->getBlockInterface(iface_id);
    const auto &dof_nums =
        m_disc->getDofNumbering()->getData(iface.getBlockIdL());

    for (UInt i : iface.getOwnedBoundaryCellsL().getXRange())
      for (UInt j : iface.getOwnedBoundaryCellsL().getYRange()) {
        if (!in(left_iface.getOwnedBoundaryCellsL(), i, j)) {
          EXPECT_NEAR((*residual)(dof_nums(i, j, 0)), 0, 1e-13);
        }
      }
  }
}


TEST_F(AdvectionTesterMultiBlock, ConsistencyFirstOrder)
{
  Real t = 1.0;
  advection::AdvectionOpts opts;
  opts.adv_velocity = {1, 2};
  auto u           = [] (Real x, Real y, Real t) { return 2*x*x + 3*y*y + t; };
  auto u0          = [&] (Real x, Real y) { return FixedVec<Real, 1>{u(x, y, 0)}; };
  auto bc_func     = [&](Real x, Real y, Real t) { return u(x, y, t); };
  auto source_func = [] (Real x, Real y, Real t) { return x*x + y*x + t*t; };

  setup(opts, bc_func, source_func);

  m_solution1->set(u0);
  m_solution2->set(u0);

  m_adv_model1->evaluateRhs(m_solution1, t, m_residual1);
  m_adv_model2->evaluateRhs(m_solution2, t, m_residual2);

  test_residuals_equal(m_residual1, m_residual2);
}

TEST_F(AdvectionTesterMultiBlock, ConsistencySecondOrder)
{
  Real t = 1.0;
  advection::AdvectionOpts opts;
  opts.limiter = common::SlopeLimiter::VanLeer;
  opts.adv_velocity = {1, 2};
  auto u           = [] (Real x, Real y, Real t) { return 2*x*x + 3*y*y + t; };
  auto u0          = [&] (Real x, Real y) { return FixedVec<Real, 1>{u(x, y, 0)}; };
  auto bc_func     = [&](Real x, Real y, Real t) { return u(x, y, t); };
  auto source_func = [] (Real x, Real y, Real t) { return x*x + y*x + t*t; };

  setup(opts, bc_func, source_func);

  m_solution1->set(u0);
  m_solution2->set(u0);

  m_adv_model1->evaluateRhs(m_solution1, t, m_residual1);
  m_adv_model2->evaluateRhs(m_solution2, t, m_residual2);

  test_residuals_equal(m_residual1, m_residual2);
}