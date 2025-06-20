#include "gtest/gtest.h"
#include "nonlinear_solvers/explicit_timestepper_opts.h"
#include "physics/advection/advection_model.h"
#include "nonlinear_solvers/explicit_euler.h"
#include "nonlinear_solvers/rk2_ssp.h"
#include "mesh/structured_mesh.h"
#include "physics/common/slope_limiters.h"
#include "utils/math.h"

using namespace structured_fv;

namespace {

class ConvergenceTester
{
  public:

    template <typename SolEx, typename Src>
    void runConvergenceStudy(const advection::AdvectionOpts opts, const std::vector<UInt>& ncells,
                            SolEx u_ex, Src src_term, bool unsteady = false)
    {
      int expected_ratio = opts.limiter == common::SlopeLimiter::FirstOrder ? 2 : 4;
      std::vector<Real> errors, ratios;
      for (UInt ncell : ncells)
      {
        std::cout << "\nncells = " << ncell << std::endl;
        UInt num_bc_ghost_cells = 2;
        UInt dofs_per_cell = 1;
        mesh::MeshSpec spec(1, 1, num_bc_ghost_cells);
        spec.blocks(0, 0) = mesh::MeshBlockSpec(ncell, ncell, 0, [](Real x, Real y) { return FixedVec<Real, 2>{x, y}; });
        auto mesh = std::make_shared<mesh::StructuredMesh>(spec);
        auto disc = std::make_shared<disc::StructuredDisc>(mesh, num_bc_ghost_cells, dofs_per_cell);

        auto u0 = [&](Real x, Real y) { return FixedVec<Real, 1>{u_ex(x, y, 0)}; };  
        std::vector<advection::Fxyt> bc_funcs{u_ex, u_ex, u_ex, u_ex};

        auto advection_model = std::make_shared<advection::AdvectionModel>(opts, disc, bc_funcs, src_term);
        auto sol = std::make_shared<disc::DiscVector<Real>>(disc, "solution");
        sol->set(u0);

        nlsolvers::ExplicitTimestepperOpts time_solver_opts_base;
        Real delta_x = 1.0/spec.blocks(0, 0).num_cells_x;
        Real delta_t = 0.375 * delta_x/std::sqrt(opts.adv_velocity[0]*opts.adv_velocity[0] + opts.adv_velocity[1]*opts.adv_velocity[1]);  //TODO: include delta_y
        time_solver_opts_base.t_end = unsteady ? 0.5/*1.0*/ : 9999;
        time_solver_opts_base.residual_tol = unsteady ? -1 : 1e-12;
        time_solver_opts_base.delta_t = delta_t;
        time_solver_opts_base.itermax = 9999;
        if (!unsteady || opts.limiter == common::SlopeLimiter::FirstOrder)
        {
          nlsolvers::ExplicitEulerOpts time_solver_opts;
          time_solver_opts = time_solver_opts_base;
          nlsolvers::explicitEuler(time_solver_opts, advection_model, sol);
        } else
        {
          nlsolvers::RK2SSPOpts time_solver_opts;
          time_solver_opts = time_solver_opts_base;
          nlsolvers::rk2ssp(time_solver_opts, advection_model, sol);
        }

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
              Real u_ex_val = u_ex(x[0], x[1], time_solver_opts_base.t_end);
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

      EXPECT_NEAR(ratios.back(), expected_ratio, 0.3);

    }

};

class Steady : public ConvergenceTester,
               public testing::TestWithParam<common::SlopeLimiter>
{};

class Unsteady : public ConvergenceTester,
                 public testing::TestWithParam<common::SlopeLimiter>
{};

class NameGenerator
{
  public:
    std::string operator()(const testing::TestParamInfo<common::SlopeLimiter>& info)
    {
      return common::get_name(info.param);
    }

};

INSTANTIATE_TEST_SUITE_P(AdvectionTests, Unsteady,
                         testing::Values(common::SlopeLimiter::FirstOrder,
                                         common::SlopeLimiter::MinMod,
                                         common::SlopeLimiter::SuperBee,
                                         common::SlopeLimiter::VanAlba,
                                         common::SlopeLimiter::VanLeer),
                         NameGenerator());

// the non-smooth limiters dont converge for several of these problems, even
// for linear advection
INSTANTIATE_TEST_SUITE_P(AdvectionTests, Steady,
                         testing::Values(common::SlopeLimiter::FirstOrder,
                                         common::SlopeLimiter::VanAlba,
                                         common::SlopeLimiter::VanLeer),
                         NameGenerator());
}

TEST_P(Steady, XSolution)
{
  advection::AdvectionOpts opts;
  opts.limiter = GetParam();
  opts.adv_velocity = {2, 0};

  auto u_ex = [&](Real x, Real y, Real t) { return std::exp(x); };
  auto src_term = [&](Real x, Real y, Real t) { return opts.adv_velocity[0]*std::exp(x); };

  std::vector<UInt> ncells = {3, 6, 12, 24};
  runConvergenceStudy(opts, ncells, u_ex, src_term);
}

TEST_P(Unsteady, XSolutionUnsteady)
{
  advection::AdvectionOpts opts;
  opts.limiter = GetParam();
  opts.adv_velocity = {2, 0};

  auto u_ex = [&](Real x, Real y, Real t) { return std::exp(x + t); };
  auto src_term = [&](Real x, Real y, Real t) { return std::exp(x+t) + opts.adv_velocity[0]*std::exp(x+t); };

  //auto u_ex = [&](Real x, Real y, Real t) { return std::sin(x - opts.adv_velocity[0]*t); };
  //auto src_term = [&](Real x, Real y, Real t) { return 0; };


  std::vector<UInt> ncells = {3, 6, 12, 24};
  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

TEST_P(Steady, XSolutionNegativeVelocity)
{
  advection::AdvectionOpts opts;
  opts.limiter = GetParam();
  opts.adv_velocity = {-2, 0};

  auto u_ex = [&](Real x, Real y, Real t) { return std::exp(x); };
  auto src_term = [&](Real x, Real y, Real t) { return opts.adv_velocity[0]*std::exp(x); };

  std::vector<UInt> ncells = {3, 6, 12, 24};
  runConvergenceStudy(opts, ncells, u_ex, src_term);
}

TEST_P(Unsteady, XSolutionNegativeVelocityUnsteady)
{
  advection::AdvectionOpts opts;
  opts.limiter = GetParam();
  opts.adv_velocity = {-2, 0};

  auto u_ex = [&](Real x, Real y, Real t) { return std::exp(x+t); };
  auto src_term = [&](Real x, Real y, Real t) { return std::exp(x+t) + opts.adv_velocity[0]*std::exp(x+t); };

  std::vector<UInt> ncells = {3, 6, 12, 24};
  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

TEST_P(Steady, XSolutionDecreasing)
{
  advection::AdvectionOpts opts;
  opts.limiter = GetParam();
  opts.adv_velocity = {2, 0};

  auto u_ex = [&](Real x, Real y, Real t) { return std::exp(1.0) - std::exp(x); };
  auto src_term = [&](Real x, Real y, Real t) { return -opts.adv_velocity[0]*std::exp(x); };

  std::vector<UInt> ncells = {3, 6, 12, 24};
  runConvergenceStudy(opts, ncells, u_ex, src_term);
}

TEST_P(Unsteady, XSolutionDecreasingUnsteady)
{
  advection::AdvectionOpts opts;
  opts.limiter = GetParam();
  opts.adv_velocity = {2, 0};

  auto u_ex = [&](Real x, Real y, Real t) { return std::exp(1.0) - std::exp(x+t); };
  auto src_term = [&](Real x, Real y, Real t) { return -std::exp(x+t) -opts.adv_velocity[0]*std::exp(x+t); };

  std::vector<UInt> ncells = {3, 6, 12, 24};
  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

TEST_P(Steady, YSolution)
{
  advection::AdvectionOpts opts;
  opts.limiter = GetParam();
  opts.adv_velocity = {0, 2};

  auto u_ex = [&](Real x, Real y, Real t) { return std::exp(y); };
  auto src_term = [&](Real x, Real y, Real t) { return opts.adv_velocity[1]*std::exp(y); };

  std::vector<UInt> ncells = {3, 6, 12, 24};
  runConvergenceStudy(opts, ncells, u_ex, src_term);
}

TEST_P(Unsteady, YSolutionUnsteady)
{
  advection::AdvectionOpts opts;
  opts.limiter = GetParam();
  opts.adv_velocity = {0, 2};

  auto u_ex = [&](Real x, Real y, Real t) { return std::exp(y+t); };
  auto src_term = [&](Real x, Real y, Real t) { return std::exp(y+t) + opts.adv_velocity[1]*std::exp(y+t); };

  std::vector<UInt> ncells = {3, 6, 12, 24};
  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

TEST_P(Steady, XYSolution)
{
  advection::AdvectionOpts opts;
  opts.limiter = GetParam();
  opts.adv_velocity = {2, 3};

  auto u_ex = [&](Real x, Real y, Real t) { return std::exp(x+y); };
  auto src_term = [&](Real x, Real y, Real t) { return opts.adv_velocity[0]*std::exp(x+y) + opts.adv_velocity[1]*std::exp(x+y); };

  std::vector<UInt> ncells = {3, 6, 12, 24};
  runConvergenceStudy(opts, ncells, u_ex, src_term);
}

TEST_P(Unsteady, XYSolutionUnsteady)
{
  advection::AdvectionOpts opts;
  opts.limiter = GetParam();
  opts.adv_velocity = {2, 3};

  auto u_ex = [&](Real x, Real y, Real t) { return std::exp(x+y+t); };
  auto src_term = [&](Real x, Real y, Real t) { return std::exp(x+y+t) + opts.adv_velocity[0]*std::exp(x+y+t) + opts.adv_velocity[1]*std::exp(x+y+t); };

  std::vector<UInt> ncells = {3, 6, 12, 24};
  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

