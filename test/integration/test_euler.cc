#include "gtest/gtest.h"
#include "linear_system/large_matrix_factory.h"
#include "physics/common/slope_limiters.h"
#include "physics/euler/euler_flux.h"
#include "physics/euler/euler_model.h"
#include "nonlinear_solvers/explicit_euler.h"
#include "nonlinear_solvers/newton.h"
#include "linear_system/large_matrix_petsc.h"

using namespace structured_fv;

namespace {

class EulerConvergenceTest : public testing::TestWithParam<euler::FluxFunction>
{};

template <typename Tsol, typename Tsrc>
void runConvergenceStudy(const euler::EulerOpts opts, const std::vector<UInt>& ncells,
                         Tsol u_ex, Tsrc src_term, bool use_newton=false)
{
  std::vector<Real> errors, ratios;
  int expected_ratio = opts.limiter == common::SlopeLimiter::FirstOrder ? 2 : 4;

  for (UInt ncell : ncells)
  {
    std::cout << "\nncells = " << ncell << std::endl;
    UInt num_bc_ghost_cells = 2;
    UInt dofs_per_cell = 4;
    mesh::MeshSpec spec(1, 1, num_bc_ghost_cells);
    spec.blocks(0, 0) = mesh::MeshBlockSpec(ncell, ncell, 0, [](Real x, Real y) { return FixedVec<Real, 2>{x, y}; });
    auto mesh = std::make_shared<mesh::StructuredMesh>(spec);
    auto disc = std::make_shared<disc::StructuredDisc>(mesh, num_bc_ghost_cells, dofs_per_cell);

    //euler::EulerOpts opts;

    //Real u = 40;
    //Real T = 298;
    //auto u_ex = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
    //auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{1, 1 - euler::Gamma_m1/2, 0, -euler::Gamma_m1/2}; };
    //auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };

    auto u0 = [&](Real x, Real y) { return u_ex(x, y, 0); };  
    std::vector<euler::Fxyt> bc_funcs{u_ex, u_ex, u_ex, u_ex};

    auto euler_model = std::make_shared<euler::EulerModel>(opts, disc, bc_funcs, src_term);
    auto sol = std::make_shared<disc::DiscVector<Real>>(disc, "solution");
    sol->set(u0);

    Real t_end = 9999;
    if (use_newton)
    {
      std::cout << "\nsolving with newton" << std::endl;
      nlsolvers::NewtonOpts newton_opts;
      newton_opts.mat_type = linear_system::LargeMatrixType::Petsc;
      auto mat_opts = std::make_shared<linear_system::LargeMatrixOptsPetsc>();
      mat_opts->petsc_opts["ksp_atol"] = "1e-15";
      mat_opts->petsc_opts["ksp_rtol"] = "1e-50";
      mat_opts->petsc_opts["ksp_monitor"] = "";
      mat_opts->petsc_opts["on_error_abort"] = "";
      newton_opts.mat_opts = mat_opts;
      newton_opts.mat_opts = std::make_shared<linear_system::LargeMatrixOptsPetsc>();     
      nlsolvers::NewtonsMethod newton(newton_opts, euler_model, sol);
      newton.solve();
    } else
    {
      std::cout << "\n solving using explicit euler" << std::endl;
      nlsolvers::ExplicitEulerOpts time_solver_opts;
      time_solver_opts.t_end = t_end;
      time_solver_opts.residual_tol = 2e-7;
      Real delta_x = 1.0/spec.blocks(0, 0).num_cells_x;
      Real delta_y = 1.0/spec.blocks(0, 0).num_cells_y;
      time_solver_opts.delta_t = 0.25 * std::min(delta_x, delta_y)/(746.0);
      time_solver_opts.itermax = 100000;

      //TODO: use rk2ssp for 2nd order
      nlsolvers::explicitEuler(time_solver_opts, euler_model, sol);
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
          auto  u_ex_val = u_ex(x[0], x[1], t_end);
          for (UInt k=0; k < 4; ++k)
          {
            Real error_k = std::abs(u_ex_val[k] - (*sol)(dof_nums(i, j, k)));
            //std::cout << i << ", " << j << ", " << k << " error = " << error_k << std::endl;
            //if (error_k > error)
            //{
            //  std::cout << "found new max error" << std::endl;
            //}
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

  EXPECT_NEAR(ratios.back(), expected_ratio, 0.25);
}
}

TEST_P(EulerConvergenceTest, XSolutionSupersonic)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real u = 400;
  Real T = 298;
  auto u_ex = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };

  euler::EulerOpts opts;
  opts.flux = GetParam();

  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

TEST_P(EulerConvergenceTest, YSolutionSupersonic)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real v = 400;
  Real T = 298;
  auto u_ex = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{y*y+1, 0, v*(y*y+1), (y*y+1)*(euler::Cv*T + 0.5*v*v)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*y*v, 0, 2*y*v*v + 2*y*euler::R*T, 2*y*(euler::Cv*T + 0.5*v*v + euler::R*T)*v}; };

  euler::EulerOpts opts;
  opts.flux = GetParam();

  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

TEST_P(EulerConvergenceTest, XSolutionSubsonic)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real u = 40;
  Real T = 298;
  auto u_ex = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };

  euler::EulerOpts opts;
  opts.flux = GetParam();

  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

TEST_P(EulerConvergenceTest, XSolutionSubsonicVanAlba)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real u = 40;
  Real T = 298;
  auto u_ex = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };

  euler::EulerOpts opts;
  opts.flux = GetParam();
  opts.limiter = common::SlopeLimiter::VanAlba;

  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

TEST_P(EulerConvergenceTest, XSolutionSubsonicVanLeer)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real u = 40;
  Real T = 298;
  auto u_ex = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };

  euler::EulerOpts opts;
  opts.flux = GetParam();
  opts.limiter = common::SlopeLimiter::VanLeer;

  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

TEST_P(EulerConvergenceTest, YSolutionSubsonic)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real v = 40;
  Real T = 298;
  auto u_ex = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{y*y+1, 0, v*(y*y+1), (y*y+1)*(euler::Cv*T + 0.5*v*v)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*y*v, 0, 2*y*v*v + 2*y*euler::R*T, 2*y*(euler::Cv*T + 0.5*v*v + euler::R*T)*v}; };

  euler::EulerOpts opts;
  opts.flux = GetParam();

  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

TEST_P(EulerConvergenceTest, XSolutionSupersonicNegative)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real u = -400;
  Real T = 298;
  auto u_ex = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };

  euler::EulerOpts opts;
  opts.flux = GetParam();


  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

TEST_P(EulerConvergenceTest, YSolutionSupersonicNegative)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real v = -400;
  Real T = 298;
  auto u_ex = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{y*y+1, 0, v*(y*y+1), (y*y+1)*(euler::Cv*T + 0.5*v*v)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*y*v, 0, 2*y*v*v + 2*y*euler::R*T, 2*y*(euler::Cv*T + 0.5*v*v + euler::R*T)*v}; };

  euler::EulerOpts opts;
  opts.flux = GetParam();


  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

TEST_P(EulerConvergenceTest, XSolutionSubsonicNegative)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real u = -40;
  Real T = 298;
  auto u_ex = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };

  euler::EulerOpts opts;
  opts.flux = GetParam();


  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

TEST_P(EulerConvergenceTest, YSolutionSubsonicNegative)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real v = -40;
  Real T = 298;
  auto u_ex = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{y*y+1, 0, v*(y*y+1), (y*y+1)*(euler::Cv*T + 0.5*v*v)}; };
  auto src_term = [&](Real x, Real y, Real t) { return FixedVec<Real, 4>{2*y*v, 0, 2*y*v*v + 2*y*euler::R*T, 2*y*(euler::Cv*T + 0.5*v*v + euler::R*T)*v}; };

  euler::EulerOpts opts;
  opts.flux = GetParam();


  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}


TEST_P(EulerConvergenceTest, XSolutionSubsonicVelocity)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real v = 0;
  Real T = 298;

  auto rho   = [](auto x, auto y, auto t) { return (x*x + 1.0); };
  auto u     = [](auto x, auto y, auto t) { return 20.0*x + 200.0; };

  auto u_ex = [&](auto x, auto y, auto t) { return euler::compute_conservative_variables(Vec4<decltype(x)>{rho(x, y, t), u(x, y, t), v, T}, euler::PrimitiveVarTag()); };
  auto src_term = [&](Real x, Real y, Real t)
  {
    Real h = 1e-20;
    std::complex<Real> x_dot(x, h);
    FixedVec<std::complex<Real>, 4> u_dot = u_ex(x_dot, y, t);
    FixedVec<std::complex<Real>, 4> flux_dot = euler::compute_euler_flux(u_dot, {1, 0});
    FixedVec<Real, 4> src;
    for (UInt i=0; i < 4; ++i)
      src[i] = flux_dot[i].imag()/h;

    return src;
  };

  euler::EulerOpts opts;
  opts.flux = GetParam();

  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

TEST_P(EulerConvergenceTest, YSolutionSubsonicVelocity)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real u = 0;
  Real T = 298;

  auto rho   = [](auto x, auto y, auto t) { return (y*y + 1.0); };
  auto v     = [](auto x, auto y, auto t) { return 20.0*y + 200.0; };

  auto u_ex = [&](auto x, auto y, auto t) { return euler::compute_conservative_variables(Vec4<decltype(y)>{rho(x, y, t), u, v(x, y, t), T}, euler::PrimitiveVarTag()); };
  auto src_term = [&](Real x, Real y, Real t)
  {
    Real h = 1e-20;
    std::complex<Real> y_dot(y, h);
    FixedVec<std::complex<Real>, 4> u_dot = u_ex(x, y_dot, t);
    FixedVec<std::complex<Real>, 4> flux_dot = euler::compute_euler_flux(u_dot, {0, 1});
    FixedVec<Real, 4> src;
    for (UInt i=0; i < 4; ++i)
      src[i] = flux_dot[i].imag()/h;

    return src;
  };

  euler::EulerOpts opts;
  opts.flux = GetParam();

  runConvergenceStudy(opts, ncells, u_ex, src_term, true);
}

INSTANTIATE_TEST_SUITE_P(, EulerConvergenceTest, 
                         testing::Values(//euler::FluxFunction::Roe,
                                         //euler::FluxFunction::RoeHH,
                                         //euler::FluxFunction::HLLE,
                                         euler::FluxFunction::LLF
                                         /*euler::FluxFunction::HLLC*/),
                        [](const testing::TestParamInfo<euler::FluxFunction>& info)
                        { return euler::get_name(info.param);
                        });