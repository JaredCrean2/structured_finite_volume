#include "gtest/gtest.h"
#include "disc/disc_block.h"
#include "disc/disc_interface.h"
#include "nonlinear_solvers/explicit_euler.h"
#include "physics/euler/euler_model.h"
#include "nonlinear_solvers/explicit_euler.h"
#include "mesh/structured_mesh.h"
#include "physics/euler/typedefs.h"
#include "utils/math.h"

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
  //const auto& coords = m_disc->getCoordField()->getData(left_iface.getBlockIdL());
  const auto& dof_nums = m_disc->getDofNumbering()->getData(left_iface.getBlockIdL());
  //Real dx = 1.0/3;
  //Real dy = 1.0/4;
  //Real cell_area = dx*dy;
  for (UInt i : left_iface.getOwnedBoundaryCellsL().getXRange())
    for (UInt j : left_iface.getOwnedBoundaryCellsL().getYRange())
    {
      //Vec2<Real> x = disc::computeCellCentroid(coords, i, j);
      //auto flux = bc_func(x[0], x[1], t) * dy;
      for (UInt k=0; k < 4; ++k)
        EXPECT_NEAR((*residual)(dof_nums(i, j, k)), 0.0, 1e-13);
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


TEST(Euler, XSolutionConvergence)
{
  std::vector<UInt> ncells = {3, 6, 12, 24, 48, 96}; // {3, 6, 12};
  std::vector<Real> errors, ratios;

  for (UInt ncell : ncells)
  {
    UInt num_bc_ghost_cells = 1;
    UInt dofs_per_cell = 4;
    mesh::MeshSpec spec(1, 1, num_bc_ghost_cells);
    spec.blocks(0, 0) = mesh::MeshBlockSpec(ncell, ncell, 0, [](Real x, Real y) { return std::array<Real, 2>{x, y}; });
    auto mesh = std::make_shared<mesh::StructuredMesh>(spec);
    auto disc = std::make_shared<disc::StructuredDisc>(mesh, num_bc_ghost_cells, dofs_per_cell);

    euler::EulerOpts opts;

    auto u_ex = [&](Real x, Real y, Real t) { return std::array<Real, 4>{x+1, x+1, 0, 100}; };
    auto src_term = [&](Real x, Real y, Real t) { return std::array<Real, 4>{1, 1 - euler::Gamma_m1/2, 0, -euler::Gamma_m1/2}; };
    auto u0 = [&](Real x, Real y) { return u_ex(x, y, 0); };  
    std::vector<euler::Fxyt> bc_funcs{u_ex, u_ex, u_ex, u_ex};

    auto euler_model = std::make_shared<euler::EulerModel>(opts, disc, bc_funcs, src_term);

    nlsolvers::ExplicitEulerOpts time_solver_opts;
    time_solver_opts.t_end = 9999;
    time_solver_opts.residual_tol = 1e-5;
    Real delta_x = 1.0/spec.blocks(0, 0).num_cells_x;
    Real delta_y = 1.0/spec.blocks(0, 0).num_cells_y;
    time_solver_opts.delta_t = 0.25 * std::min(delta_x, delta_y)/(1 + 6.1);
    time_solver_opts.itermax = 50000;
    auto sol = std::make_shared<disc::DiscVector<Real>>(disc, "solution");
    sol->set(u0);

    nlsolvers::explicitEuler(time_solver_opts, euler_model, sol);

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
          auto  u_ex_val = u_ex(x[0], x[1], time_solver_opts.t_end);
          for (UInt k=0; k < 4; ++k)
          {
            Real error_k = std::abs(u_ex_val[k] - (*sol)(dof_nums(i, j, k)));
            std::cout << i << ", " << j << ", " << k << " error = " << error_k << std::endl;
            if (error_k > error)
            {
              std::cout << "found new max error" << std::endl;
            }
            error = std::max(error, error_k);
          }
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

  EXPECT_NEAR(ratios.back(), 2, 0.25);
}


