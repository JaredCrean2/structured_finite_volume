#include "gtest/gtest.h"
#include "disc/disc_block.h"
#include "disc/disc_interface.h"
#include "nonlinear_solvers/explicit_euler.h"
#include "physics/advection/advection_model.h"
#include "nonlinear_solvers/explicit_euler.h"
#include "mesh/structured_mesh.h"
#include "utils/math.h"

namespace {

using namespace structured_fv;

class AdvectionTester : public ::testing::Test
{
  public:
    AdvectionTester() :
      spec(1, 1, m_num_bc_ghost_cells)
    {
      spec.blocks(0, 0) = mesh::MeshBlockSpec(3, 4, 0, [](Real x, Real y) { return std::array<Real, 2>{x, y}; });
      m_mesh = std::make_shared<mesh::StructuredMesh>(spec);
      m_disc = std::make_shared<disc::StructuredDisc>(m_mesh, m_num_bc_ghost_cells, m_dofs_per_cell);
    }

    void setup(const advection::AdvectionOpts& opts, advection::Fxyt bc_func,
               advection::Fxyt source_func)
    {
      std::vector<advection::Fxyt> bc_funcs{bc_func, bc_func, bc_func, bc_func};
      m_adv_model = std::make_shared<advection::AdvectionModel>(opts, m_disc, bc_funcs, source_func);
    }



    UInt m_dofs_per_cell = 1;
    int m_num_bc_ghost_cells = 1;
    mesh::MeshSpec spec;
    std::shared_ptr<mesh::StructuredMesh> m_mesh;
    std::shared_ptr<disc::StructuredDisc> m_disc;
    std::shared_ptr<advection::AdvectionModel> m_adv_model;
};
}

TEST_F(AdvectionTester, SourceTerm)
{
  auto solution = std::make_shared<disc::DiscVector<Real>>(m_disc, "solution");
  auto residual = std::make_shared<disc::DiscVector<Real>>(m_disc, "residual");
  solution->set(0);
  residual->set(0);

  Real t = 1.0;
  advection::AdvectionOpts opts;
  opts.adv_velocity = {1, 0};
  auto bc_func = [](Real x, Real y, Real t){ return 0; };
  auto source_func = [](Real x, Real y, Real t){ return x + y + t; };
  setup(opts, bc_func, source_func);
  m_adv_model->evaluateRhs(solution, t, residual);

  const disc::StructuredBlock& block = m_disc->getBlock(0);
  const auto& coords = m_disc->getCoordField()->getData(0);
  const auto& dof_nums = m_disc->getDofNumbering()->getData(0);
  for (UInt i : block.getOwnedCells().getXRange())
    for (UInt j : block.getOwnedCells().getYRange())
    {
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
  auto bc_func = [](Real x, Real y, Real t){ return 1; };
  auto source_func = [](Real x, Real y, Real t){ return 0; };
  setup(opts, bc_func, source_func);
  m_adv_model->evaluateRhs(solution, t, residual);

  const disc::StructuredBlockInterface& left_iface = m_disc->getBlockInterface(3);
  const auto& coords = m_disc->getCoordField()->getData(left_iface.getBlockIdL());
  const auto& dof_nums = m_disc->getDofNumbering()->getData(left_iface.getBlockIdL());
  Real dx = 1.0/3;
  Real dy = 1.0/4;
  Real cell_area = dx*dy;
  for (UInt i : left_iface.getOwnedBoundaryCellsL().getXRange())
    for (UInt j : left_iface.getOwnedBoundaryCellsL().getYRange())
    {
      Vec2<Real> x = disc::computeCellCentroid(coords, i, j);
      Real flux = bc_func(x[0], x[1], t) * dy;
      EXPECT_NEAR((*residual)(dof_nums(i, j, 0)), flux/cell_area, 1e-13);
    }

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

TEST(Advection, XSolutionConvergence)
{
  std::vector<UInt> ncells = {3, 6, 12};
  std::vector<Real> errors, ratios;

  for (UInt ncell : ncells)
  {
    UInt num_bc_ghost_cells = 1;
    UInt dofs_per_cell = 1;
    mesh::MeshSpec spec(1, 1);
    spec.blocks(0, 0) = mesh::MeshBlockSpec(ncell, ncell, 0, [](Real x, Real y) { return std::array<Real, 2>{x, y}; });
    auto mesh = std::make_shared<mesh::StructuredMesh>(spec);
    auto disc = std::make_shared<disc::StructuredDisc>(mesh, num_bc_ghost_cells, dofs_per_cell);

    advection::AdvectionOpts opts;
    opts.adv_velocity = {1, 0};

    auto u_ex = [&](Real x, Real y, Real t) { return std::exp(x); };
    auto src_term = [&](Real x, Real y, Real t) { return opts.adv_velocity[0]*std::exp(x); };
    auto u0 = [&](Real x, Real y) { return std::array<Real, 1>{u_ex(x, y, 0)}; };  
    std::vector<advection::Fxyt> bc_funcs{u_ex, u_ex, u_ex, u_ex};

    auto advection_model = std::make_shared<advection::AdvectionModel>(opts, disc, bc_funcs, src_term);

    nlsolvers::ExplicitEulerOpts time_solver_opts;
    time_solver_opts.t_end = 9999;
    time_solver_opts.residual_tol = 1e-10;
    Real delta_x = 1.0/spec.blocks(0, 0).num_cells_x;
    time_solver_opts.delta_t = 0.9 * delta_x/opts.adv_velocity[0];
    time_solver_opts.itermax = 500;
    auto sol = std::make_shared<disc::DiscVector<Real>>(disc, "solution");
    sol->set(u0);

    nlsolvers::explicitEuler(time_solver_opts, advection_model, sol);

    Real error = 0.0;
    for (UInt block_id : disc->getRegularBlocksIds())
    {
      const disc::StructuredBlock& block = disc->getBlock(block_id);
      const auto& dof_nums = disc->getDofNumbering()->getData(block_id);
      const auto& coords = disc->getCoordField()->getData(block_id);

      for (UInt i : block.getOwnedCells().getXRange())
        for (UInt j : block.getOwnedCells().getYRange())
        {
          Vec2<Real> x = disc::computeCellCentroid(coords, i, j);
          Real u_ex_val = u_ex(x[0], x[1], time_solver_opts.t_end);
          error = std::max(error, std::abs(u_ex_val - (*sol)(dof_nums(i, j, 0))));
        }
    }

    errors.push_back(error);
    ratios.push_back(errors.size() > 1 ? errors[errors.size()-2]/errors[errors.size()-1] : 0);
  }

  std::cout << "errors:" << std::endl;
  for (UInt i=0; i < ncells.size(); ++i)
  {
    std::cout << ncells[i] << ": " << errors[i] << ", ratios = " << ratios[i] << std::endl;
  }

  EXPECT_NEAR(ratios.back(), 2, 0.1);
}

