#include "gtest/gtest.h"
#include "physics/euler/euler_flux.h"
#include "physics/euler/euler_model.h"
#include "nonlinear_solvers/explicit_euler.h"

using namespace structured_fv;

namespace {

class EulerConvergenceTest : public testing::TestWithParam<euler::FluxFunction>
{};

template <typename Tsol, typename Tsrc>
void runConvergenceStudy(const euler::EulerOpts opts, const std::vector<UInt>& ncells, Real expected_ratio,
                         Tsol u_ex, Tsrc src_term)
{
  std::vector<Real> errors, ratios;

  for (UInt ncell : ncells)
  {
    std::cout << "\nncells = " << ncell << std::endl;
    UInt num_bc_ghost_cells = 1;
    UInt dofs_per_cell = 4;
    mesh::MeshSpec spec(1, 1, num_bc_ghost_cells);
    spec.blocks(0, 0) = mesh::MeshBlockSpec(ncell, ncell, 0, [](Real x, Real y) { return std::array<Real, 2>{x, y}; });
    auto mesh = std::make_shared<mesh::StructuredMesh>(spec);
    auto disc = std::make_shared<disc::StructuredDisc>(mesh, num_bc_ghost_cells, dofs_per_cell);

    //euler::EulerOpts opts;

    //Real u = 40;
    //Real T = 298;
    //auto u_ex = [&](Real x, Real y, Real t) { return std::array<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
    //auto src_term = [&](Real x, Real y, Real t) { return std::array<Real, 4>{1, 1 - euler::Gamma_m1/2, 0, -euler::Gamma_m1/2}; };
    //auto src_term = [&](Real x, Real y, Real t) { return std::array<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };

    auto u0 = [&](Real x, Real y) { return u_ex(x, y, 0); };  
    std::vector<euler::Fxyt> bc_funcs{u_ex, u_ex, u_ex, u_ex};

    auto euler_model = std::make_shared<euler::EulerModel>(opts, disc, bc_funcs, src_term);

    nlsolvers::ExplicitEulerOpts time_solver_opts;
    time_solver_opts.t_end = 9999;
    time_solver_opts.residual_tol = 1e-7;
    Real delta_x = 1.0/spec.blocks(0, 0).num_cells_x;
    Real delta_y = 1.0/spec.blocks(0, 0).num_cells_y;
    time_solver_opts.delta_t = 0.25 * std::min(delta_x, delta_y)/(746.0);
    time_solver_opts.itermax = 100000;
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

  EXPECT_NEAR(ratios.back(), 2, 0.25);
}
}

TEST_P(EulerConvergenceTest, XSolutionSupersonic)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real u = 400;
  Real T = 298;
  auto u_ex = [&](Real x, Real y, Real t) { return std::array<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return std::array<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };

  euler::EulerOpts opts;
  opts.flux = GetParam();

  runConvergenceStudy(opts, ncells, 2, u_ex, src_term);
}

TEST_P(EulerConvergenceTest, XSolutionSubsonic)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real u = 40;
  Real T = 298;
  auto u_ex = [&](Real x, Real y, Real t) { return std::array<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return std::array<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };

  euler::EulerOpts opts;
  opts.flux = GetParam();

  runConvergenceStudy(opts, ncells, 2, u_ex, src_term);
}

TEST_P(EulerConvergenceTest, XSolutionSupersonicNegative)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real u = -400;
  Real T = 298;
  auto u_ex = [&](Real x, Real y, Real t) { return std::array<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return std::array<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };

  euler::EulerOpts opts;
  opts.flux = GetParam();


  runConvergenceStudy(opts, ncells, 2, u_ex, src_term);
}

TEST_P(EulerConvergenceTest, XSolutionSubsonicNegative)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real u = -40;
  Real T = 298;
  auto u_ex = [&](Real x, Real y, Real t) { return std::array<Real, 4>{x*x+1, u*(x*x+1), 0, (x*x+1)*(euler::Cv*T + 0.5*u*u)}; };
  auto src_term = [&](Real x, Real y, Real t) { return std::array<Real, 4>{2*x*u, 2*x*u*u + 2*x*euler::R*T, 0, 2*x*(euler::Cv*T + 0.5*u*u + euler::R*T)*u}; };

  euler::EulerOpts opts;
  opts.flux = GetParam();


  runConvergenceStudy(opts, ncells, 2, u_ex, src_term);
}


TEST_P(EulerConvergenceTest, XSolutionSubsonicVelocity)
{
  std::vector<UInt> ncells = {3, 6, 12, 24}; // {3, 6, 12};
  Real v = 0;
  Real T = 298;
  /*
  auto rho   = [](Real x, Real y, Real t) { return (x*x + 1); };
  auto rho_x = [](Real x, Real y, Real t) { return 2*x; };
  auto u     = [](Real x, Real y, Real t) { return 20*x + 400; };
  auto u_x   = [](Real x, Real y, Real t) { return 20; };

  auto p = [&](Real x, Real y, Real t) { return rho(x, y, t) * euler::R * T; };
  auto p_x = [&](Real x, Real y, Real t) { return rho_x(x, y, t) * euler::R * T; };

  auto E = [&](Real x, Real y, Real t) { return rho(x, y, t)*euler::Cv*T + 0.5*rho(x, y, t)*(u(x, y, t)*u(x, y, t) + v*v); };
  auto E_x = [&](Real x, Real y, Real t)
  {
    Real velocity = u(x, y, t)*u(x, y, t) + v*v;
    Real velocity_x = 2*u(x, y, t)*u_x(x, y, t); 
    return rho_x(x, y, t)*euler::Cv*T + 0.5*(rho_x(x, y, t)*velocity + rho(x, y, t)*velocity_x);
  };
  auto u_ex = [&](Real x, Real y, Real t) { return std::array<Real, 4>{rho(x, y, t), rho(x, y, t)*u(x, y, t), rho(x, y, t)*v, E(x, y, t)}; };
  auto src_term = [&](Real x, Real y, Real t) { return std::array<Real, 4>{rho(x, y, t)*u_x(x, y, t) + rho_x(x, y, t)*u(x, y, t),
                                                                           rho_x(x, y, t)*u(x, y, t)*u(x, y, t) + rho(x, y, t)*u(x, y, t)*u_x(x, y, t) + p_x(x, y, t),
                                                                           rho_x(x, y, t)*u(x, y, t)*v + rho(x, y, t)*u_x(x, y, t)*v,
                                                                           (E(x, y, t) + p(x, y, t))*u_x(x, y, t) + (E_x(x, y, t) + p_x(x, y, t))*u(x, y, t)}; };
*/
  auto rho   = [](auto x, auto y, auto t) { return (x*x + 1.0); };
  auto u     = [](auto x, auto y, auto t) { return 20.0*x + 200.0; };
  //auto p = [&](auto x, auto y, auto t) { return rho(x, y, t) * euler::R * T; };
  //auto E = [&](auto x, auto y, auto t) { return rho(x, y, t)*euler::Cv*T + 0.5*rho(x, y, t)*(u(x, y, t)*u(x, y, t) + v*v); };
  //auto u_ex = [&](auto x, auto y, auto t) { return euler::Vec4<Real>{rho(x, y, t), rho(x, y, t)*u(x, y, t), rho(x, y, t)*v, E(x, y, t)}; };
  auto u_ex = [&](auto x, auto y, auto t) { return euler::compute_conservative_variables(euler::Vec4<decltype(x)>{rho(x, y, t), u(x, y, t), v, T}, euler::PrimitiveVarTag()); };

  auto src_term = [&](Real x, Real y, Real t)
  {
    Real h = 1e-20;
    std::complex<Real> x_dot(x, h);
    std::array<std::complex<Real>, 4> u_dot = u_ex(x_dot, y, t);
    std::array<std::complex<Real>, 4> flux_dot = euler::compute_euler_flux(u_dot, {1, 0});
    std::array<Real, 4> src;
    for (UInt i=0; i < 4; ++i)
      src[i] = flux_dot[i].imag()/h;

    return src;
  };

  euler::EulerOpts opts;
  opts.flux = GetParam();

  runConvergenceStudy(opts, ncells, 2, u_ex, src_term);
}

INSTANTIATE_TEST_SUITE_P(, EulerConvergenceTest, 
                         testing::Values(euler::FluxFunction::Roe,
                                         euler::FluxFunction::HLLE,
                                         euler::FluxFunction::LLF),
                        [](const testing::TestParamInfo<euler::FluxFunction>& info)
                        { return euler::get_name(info.param);
                        });