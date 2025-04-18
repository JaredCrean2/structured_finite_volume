#include "gtest/gtest.h"
#include "physics/advection/advection_model.h"
#include "nonlinear_solvers/explicit_euler.h"
#include "mesh/structured_mesh.h"
#include "utils/math.h"

using namespace structured_fv;

namespace {

template <typename SolEx, typename Src>
void runConvergenceStudy(const advection::AdvectionOpts opts, const std::vector<UInt>& ncells, Real expected_ratio,
                         SolEx u_ex, Src src_term)
{
  std::vector<Real> errors, ratios;
  for (UInt ncell : ncells)
  {
    std::cout << "\nncells = " << ncell << std::endl;
    UInt num_bc_ghost_cells = 1;
    UInt dofs_per_cell = 1;
    mesh::MeshSpec spec(1, 1, num_bc_ghost_cells);
    spec.blocks(0, 0) = mesh::MeshBlockSpec(ncell, ncell, 0, [](Real x, Real y) { return std::array<Real, 2>{x, y}; });
    auto mesh = std::make_shared<mesh::StructuredMesh>(spec);
    auto disc = std::make_shared<disc::StructuredDisc>(mesh, num_bc_ghost_cells, dofs_per_cell);

    auto u0 = [&](Real x, Real y) { return std::array<Real, 1>{u_ex(x, y, 0)}; };  
    std::vector<advection::Fxyt> bc_funcs{u_ex, u_ex, u_ex, u_ex};

    auto advection_model = std::make_shared<advection::AdvectionModel>(opts, disc, bc_funcs, src_term);

    nlsolvers::ExplicitEulerOpts time_solver_opts;
    time_solver_opts.t_end = 9999;
    time_solver_opts.residual_tol = 1e-14;
    Real delta_x = 1.0/spec.blocks(0, 0).num_cells_x;
    time_solver_opts.delta_t = 0.5 * delta_x/std::sqrt(opts.adv_velocity[0]*opts.adv_velocity[0] + opts.adv_velocity[1]*opts.adv_velocity[1]);  //TODO: include delta_y
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

  EXPECT_NEAR(ratios.back(), expected_ratio, 0.1);

}

}

TEST(Advection, XSolutionConvergence)
{
  advection::AdvectionOpts opts;
  opts.adv_velocity = {2, 0};

  auto u_ex = [&](Real x, Real y, Real t) { return std::exp(x); };
  auto src_term = [&](Real x, Real y, Real t) { return opts.adv_velocity[0]*std::exp(x); };

  std::vector<UInt> ncells = {3, 6, 12, 24};
  runConvergenceStudy(opts, ncells, 2, u_ex, src_term);
}

TEST(Advection, YSolutionConvergence)
{
  advection::AdvectionOpts opts;
  opts.adv_velocity = {0, 2};

  auto u_ex = [&](Real x, Real y, Real t) { return std::exp(y); };
  auto src_term = [&](Real x, Real y, Real t) { return opts.adv_velocity[1]*std::exp(y); };

  std::vector<UInt> ncells = {3, 6, 12, 24};
  runConvergenceStudy(opts, ncells, 2, u_ex, src_term);
}

TEST(Advection, XYSolutionConvergence)
{
  advection::AdvectionOpts opts;
  opts.adv_velocity = {2, 3};

  auto u_ex = [&](Real x, Real y, Real t) { return std::exp(x+y); };
  auto src_term = [&](Real x, Real y, Real t) { return opts.adv_velocity[0]*std::exp(x+y) + opts.adv_velocity[1]*std::exp(x+y); };

  std::vector<UInt> ncells = {3, 6, 12, 24};
  runConvergenceStudy(opts, ncells, 2, u_ex, src_term);
}