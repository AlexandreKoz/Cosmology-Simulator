// Multirate outer PM level: driver-level convergence suites (2026-10-26
// round 3; AUTHORED, NOT RUN — owner validation required).
//
// Analytic split forces are driven through the production
// StepOrchestrator::executeHierarchicalBlockWithDispatcher with a
// HierarchicalPmCyclePlan, exactly as the time coordinator drives it, so the
// executed operator composition
//   K_PM[t0,t0+T/2] ∘ B_1 ∘ … ∘ B_m ∘ K_PM[t0+T/2,t0+T]
// (each B_j a power-of-two tree KDK block) is measured against exact
// solutions. Non-cosmological (no background): every kick/drift factor is a
// plain time interval, so the dynamics are Newtonian and exact solutions are
// available. Suites:
//   * PM-only    (tree force zero): equals leapfrog with step T; order 2 in T.
//   * tree-only  (PM force zero): identical to the reference scheme; order 2.
//   * combined   (fast tree + slow PM oscillator): r-RESPA order 2 in T at
//                fixed blocks per cycle; splitting error vs reference -> 0.
//   * variable   non-uniform block intervals inside each cycle: order 2.
//   * adaptive   per-particle rungs that change between blocks: order 2.
//   * two-body   circular orbit under a smooth Tree/PM force split: order 2
//                in phase; linear and angular momentum conserved to roundoff.
// Gates are observed orders from three successive halvings; they are the
// owner's validation targets for the claims in docs/time_integration.md.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "../support/test_require.hpp"
#include "cosmosim/core/cosmology.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/core/time_integration.hpp"

namespace {
using namespace cosmosim::core;
using Vec3 = std::array<double, 3>;

constexpr std::uint8_t k_max_bin = 2U;  // 4 fine ticks per block

struct Problem {
  std::size_t particle_count = 1U;
  // Initial positions/velocities (state lanes are zero-initialized).
  std::function<void(SimulationState&)> initialize;
  // Tree (short-range) and PM (long-range) acceleration of one row from the
  // current positions of all rows.
  std::function<void(const SimulationState&, std::size_t, Vec3&, Vec3&)> force;
  // Rung of a row for a block (applied at the block's synchronization point).
  std::function<std::uint8_t(std::size_t, std::uint64_t)> rung =
      [](std::size_t, std::uint64_t) { return std::uint8_t{0U}; };
};

struct Outcome {
  std::vector<Vec3> position;
  std::vector<Vec3> velocity;
  double time = 0.0;
  double scale_factor = 1.0;
};

// Integrates `cycles` PM cycles of length T. Each cycle is tiled by blocks
// whose intervals are fractions*T. multirate=false runs the reference scheme
// (PM surfaces at both endpoints of every block) on the same blocks.
// Optional cosmological background: code time in units of 1/H0, the outer
// PM plan factors computed exactly as TimeCoordinator computes them.
struct Cosmology {
  const LambdaCdmBackground* background = nullptr;
  double a_begin = 1.0;
  double t_begin = 0.0;
};

Outcome integrate(const Problem& problem, double cycle_length, std::size_t cycles,
                  const std::vector<double>& fractions, bool multirate, const Cosmology& cosmology = {}) {
  const std::size_t n = problem.particle_count;
  SimulationState state;
  state.resizeParticles(n);
  for (std::size_t row = 0U; row < n; ++row) {
    state.particle_sidecar.particle_id[row] = row + 1U;
    state.particle_sidecar.species_tag[row] = static_cast<std::uint32_t>(ParticleSpecies::kDarkMatter);
    state.particle_sidecar.owning_rank[row] = 0U;
    state.particles.mass_code[row] = 1.0;
  }
  problem.initialize(state);
  state.rebuildSpeciesIndex();
  state.updateAllParticleDriftEpoch(cosmology.t_begin, cosmology.a_begin);
  CHUI_REQUIRE(state.compactHomogeneousDmoMetadata(0U));
  IntegratorState integrator;
  integrator.current_time_code = cosmology.t_begin;
  integrator.current_scale_factor = cosmology.a_begin;
  if (cosmology.background != nullptr) integrator.time_si_per_code = 1.0 / cosmology.background->hubble0Si();
  HierarchicalTimeBinScheduler particles(k_max_bin), cells(k_max_bin);
  particles.reset(static_cast<std::uint32_t>(n), 0U, 0U);
  cells.reset(0U, 0U, 0U);
  const double periods = static_cast<double>(particles.binPeriodTicks(k_max_bin));

  std::vector<Vec3> tree(n), pm(n);
  for (std::size_t row = 0U; row < n; ++row) problem.force(state, row, tree[row], pm[row]);  // bootstrap
  const auto lanes = [&state](std::size_t row) {
    return std::array<double*, 6>{&state.particles.position_x_comoving[row], &state.particles.position_y_comoving[row],
                                  &state.particles.position_z_comoving[row], &state.particles.velocity_x_peculiar[row],
                                  &state.particles.velocity_y_peculiar[row], &state.particles.velocity_z_peculiar[row]};
  };
  const StageDispatchFunction dispatcher = [&](StepContext& context, bool) {
    const auto& kdk = context.hierarchical_kdk;
    if (context.stage == IntegrationStage::kGravityKickPre || context.stage == IntegrationStage::kGravityKickPost) {
      for (const auto row : context.active_set.particle_indices) {
        const double tree_factor = kdk.tree_kick_factor_code[particles.binIndex(row)];
        const auto lane = lanes(row);
        for (std::size_t axis = 0U; axis < 3U; ++axis) {
          *lane[3U + axis] += tree_factor * tree[row][axis] + kdk.pm_kick_factor_code * pm[row][axis];
        }
      }
    } else if (context.stage == IntegrationStage::kDrift) {
      for (const auto row : context.active_set.particle_indices) {
        const auto lane = lanes(row);
        for (std::size_t axis = 0U; axis < 3U; ++axis) {
          *lane[axis] += context.timeline_step.drift_factor_code * *lane[3U + axis];
          *lane[3U + axis] *= context.timeline_step.hubble_drag_factor;
        }
      }
    } else if (context.stage == IntegrationStage::kForceRefresh) {
      // Tree component for the rows closing a half step; the PM component is
      // refreshed (for every row) only on a closing PM surface.
      Vec3 unused{};
      for (const auto row : context.active_set.particle_indices) problem.force(state, row, tree[row], unused);
      if (kdk.include_long_range_force) {
        for (std::size_t row = 0U; row < n; ++row) problem.force(state, row, unused, pm[row]);
      }
    }
  };

  StepOrchestrator orchestrator;
  TransientStepWorkspace workspace;
  std::uint64_t block = 0U;
  for (std::size_t cycle = 0U; cycle < cycles; ++cycle) {
    for (std::size_t j = 0U; j < fractions.size(); ++j, ++block) {
      for (std::size_t row = 0U; row < n; ++row) {
        particles.setElementBin(static_cast<std::uint32_t>(row), problem.rung(row, block), particles.currentTick());
      }
      workspace.prepareGravityParticleIndexScratch(n);
      std::iota(workspace.gravity_particle_index_scratch.begin(), workspace.gravity_particle_index_scratch.end(), 0U);
      integrator.dt_time_code = fractions[j] * cycle_length / periods;
      HierarchicalPmCyclePlan plan{};
      if (multirate) {
        plan.enabled = true;
        plan.opens_cycle = j == 0U;
        plan.closes_cycle = j + 1U == fractions.size();
        plan.opening_pm_kick_factor_code = 0.5 * cycle_length;
        plan.closing_half_interval_time_code = 0.5 * cycle_length;
        if (cosmology.background != nullptr) {
          const double a0 = integrator.current_scale_factor;
          plan.cycle_mid_scale_factor = advanceScaleFactorByCosmicTime(
              *cosmology.background, a0, 0.5 * cycle_length * integrator.time_si_per_code);
          plan.opening_pm_kick_factor_code = computeComovingDriftFactor(
              *cosmology.background, a0, plan.cycle_mid_scale_factor, 64U) / (integrator.time_si_per_code * a0);
        }
      }
      orchestrator.executeHierarchicalBlockWithDispatcher(state, integrator, particles, cells, dispatcher,
                                                          cosmology.background, workspace, nullptr, nullptr, {}, plan);
      CHUI_REQUIRE(!integrator.inside_kdk_step);
      CHUI_REQUIRE(integrator.last_completed_restart_safe == (!multirate || j + 1U == fractions.size()));
    }
  }
  Outcome outcome;
  outcome.time = integrator.current_time_code;
  outcome.scale_factor = integrator.current_scale_factor;
  for (std::size_t row = 0U; row < n; ++row) {
    const auto lane = lanes(row);
    outcome.position.push_back({*lane[0], *lane[1], *lane[2]});
    outcome.velocity.push_back({*lane[3], *lane[4], *lane[5]});
  }
  return outcome;
}

double difference(const Outcome& a, const Outcome& b) {
  double e = 0.0;
  for (std::size_t row = 0U; row < a.position.size(); ++row) {
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      e = std::max(e, std::abs(a.position[row][axis] - b.position[row][axis]));
      e = std::max(e, std::abs(a.velocity[row][axis] - b.velocity[row][axis]));
    }
  }
  return e;
}

// Observed orders from successive halvings of T (errors[k] at T/2^k).
void requireOrder(const std::string& label, const std::vector<double>& errors, double low, double high) {
  for (std::size_t k = 1U; k < errors.size(); ++k) {
    CHUI_REQUIRE(errors[k] > 0.0 && std::isfinite(errors[k]));
    const double order = std::log2(errors[k - 1U] / errors[k]);
    std::cout << "multirate_convergence " << label << " halving=" << k << " error=" << errors[k]
              << " observed_order=" << order << '\n';
    CHUI_REQUIRE(order >= low && order <= high);
  }
}

// Independent 1D oscillators along x: A = -(ws^2 + wl^2) x, split into a
// tree part -ws^2 x and a PM part -wl^2 x. Exact: x = cos(Wt), u = -W sin(Wt).
Problem oscillators(std::vector<double> ws, std::vector<double> wl) {
  Problem problem;
  problem.particle_count = ws.size();
  problem.initialize = [](SimulationState& state) {
    for (std::size_t row = 0U; row < state.particles.size(); ++row) state.particles.position_x_comoving[row] = 1.0;
  };
  problem.force = [ws, wl](const SimulationState& state, std::size_t row, Vec3& tree, Vec3& pm) {
    const double x = state.particles.position_x_comoving[row];
    tree = {-ws[row] * ws[row] * x, 0.0, 0.0};
    pm = {-wl[row] * wl[row] * x, 0.0, 0.0};
  };
  return problem;
}

Outcome oscillatorExact(const std::vector<double>& ws, const std::vector<double>& wl, double t) {
  Outcome exact;
  exact.time = t;
  for (std::size_t row = 0U; row < ws.size(); ++row) {
    const double w = std::sqrt(ws[row] * ws[row] + wl[row] * wl[row]);
    exact.position.push_back({std::cos(w * t), 0.0, 0.0});
    exact.velocity.push_back({-w * std::sin(w * t), 0.0, 0.0});
  }
  return exact;
}

const std::vector<double> k_uniform4{0.25, 0.25, 0.25, 0.25};

void testPmOnly() {
  const std::vector<double> ws{0.0}, wl{1.0};
  const auto problem = oscillators(ws, wl);
  std::vector<double> errors;
  for (const std::size_t cycles : {20U, 40U, 80U}) {
    const double t_final = 2.0, cycle_length = t_final / static_cast<double>(cycles);
    const auto multirate = integrate(problem, cycle_length, cycles, k_uniform4, true);
    // With no tree force the composition is K(T/2) D(T) K(T/2): leapfrog at T.
    const auto leapfrog = integrate(problem, cycle_length, cycles, {1.0}, false);
    CHUI_REQUIRE(difference(multirate, leapfrog) < 1e-12);
    CHUI_REQUIRE(std::abs(multirate.time - t_final) < 1e-12);
    errors.push_back(difference(multirate, oscillatorExact(ws, wl, multirate.time)));
  }
  requireOrder("pm_only", errors, 1.9, 2.1);
}

void testTreeOnly() {
  const std::vector<double> ws{3.0}, wl{0.0};
  const auto problem = oscillators(ws, wl);
  std::vector<double> errors;
  for (const std::size_t cycles : {20U, 40U, 80U}) {
    const double cycle_length = 2.0 / static_cast<double>(cycles);
    const auto multirate = integrate(problem, cycle_length, cycles, k_uniform4, true);
    // No PM force: the PM surfaces carry zero and the schemes coincide.
    const auto reference = integrate(problem, cycle_length, cycles, k_uniform4, false);
    CHUI_REQUIRE(difference(multirate, reference) <= 1e-15);
    errors.push_back(difference(multirate, oscillatorExact(ws, wl, multirate.time)));
  }
  requireOrder("tree_only", errors, 1.9, 2.1);
}

void testCombinedAndVariableBlocks() {
  const std::vector<double> ws{4.0}, wl{1.0};
  const auto problem = oscillators(ws, wl);
  for (const auto& [label, fractions] : std::vector<std::pair<std::string, std::vector<double>>>{
           {"combined_uniform", k_uniform4}, {"combined_variable", {0.1, 0.4, 0.2, 0.3}}}) {
    std::vector<double> errors, splitting;
    for (const std::size_t cycles : {20U, 40U, 80U}) {
      const double cycle_length = 2.0 / static_cast<double>(cycles);
      const auto multirate = integrate(problem, cycle_length, cycles, fractions, true);
      const auto reference = integrate(problem, cycle_length, cycles, fractions, false);
      errors.push_back(difference(multirate, oscillatorExact(ws, wl, multirate.time)));
      splitting.push_back(difference(multirate, reference));
    }
    requireOrder(label, errors, 1.85, 2.15);
    requireOrder(label + "_splitting", splitting, 1.8, 2.2);
  }
}

void testAdaptiveRungs() {
  // Row 0: fast tree oscillator on the finest rung. Row 1: slow, alternating
  // between rungs 1 and 2 at successive block synchronization points.
  const std::vector<double> ws{6.0, 1.0}, wl{0.5, 0.5};
  auto problem = oscillators(ws, wl);
  problem.rung = [](std::size_t row, std::uint64_t block) {
    return row == 0U ? std::uint8_t{0U} : static_cast<std::uint8_t>(block % 2U == 0U ? 2U : 1U);
  };
  std::vector<double> errors;
  for (const std::size_t cycles : {20U, 40U, 80U}) {
    const double cycle_length = 2.0 / static_cast<double>(cycles);
    const auto multirate = integrate(problem, cycle_length, cycles, k_uniform4, true);
    errors.push_back(difference(multirate, oscillatorExact(ws, wl, multirate.time)));
  }
  requireOrder("adaptive_rungs", errors, 1.8, 2.2);
}

void testTwoBodyCircularOrbit() {
  // G = 1, two bodies of mass 1/2 on a circular orbit of separation 1
  // (angular frequency 1). Pair force F = -(1/2) r/|r|^3 on body 0, split by
  // S(r) = exp(-(r/r_s)^2): tree = S F, PM = (1 - S) F.
  constexpr double r_s = 0.8;
  Problem problem;
  problem.particle_count = 2U;
  problem.initialize = [](SimulationState& state) {
    state.particles.position_x_comoving[0] = 0.5;
    state.particles.position_x_comoving[1] = -0.5;
    state.particles.velocity_y_peculiar[0] = 0.5;
    state.particles.velocity_y_peculiar[1] = -0.5;
  };
  problem.force = [](const SimulationState& state, std::size_t row, Vec3& tree, Vec3& pm) {
    const std::size_t other = 1U - row;
    const Vec3 d{state.particles.position_x_comoving[row] - state.particles.position_x_comoving[other],
                 state.particles.position_y_comoving[row] - state.particles.position_y_comoving[other],
                 state.particles.position_z_comoving[row] - state.particles.position_z_comoving[other]};
    const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    const double split = std::exp(-(r / r_s) * (r / r_s));
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      const double f = -0.5 * d[axis] / (r * r * r);
      tree[axis] = split * f;
      pm[axis] = (1.0 - split) * f;
    }
  };
  const double pi = std::acos(-1.0);
  std::vector<double> errors;
  for (const std::size_t cycles : {50U, 100U, 200U}) {
    const double cycle_length = 2.0 * pi / static_cast<double>(cycles);
    const auto orbit = integrate(problem, cycle_length, cycles, k_uniform4, true);
    Outcome exact;
    const double t = orbit.time;
    exact.position = {{0.5 * std::cos(t), 0.5 * std::sin(t), 0.0}, {-0.5 * std::cos(t), -0.5 * std::sin(t), 0.0}};
    exact.velocity = {{-0.5 * std::sin(t), 0.5 * std::cos(t), 0.0}, {0.5 * std::sin(t), -0.5 * std::cos(t), 0.0}};
    errors.push_back(difference(orbit, exact));
    // Equal-mass central pair force on a common rung: every kick is
    // antisymmetric along the separation, so linear and angular momentum are
    // conserved to roundoff by every operator of the composition.
    for (std::size_t axis = 0U; axis < 3U; ++axis) {
      CHUI_REQUIRE(std::abs(orbit.velocity[0][axis] + orbit.velocity[1][axis]) < 1e-13);
    }
    double angular_momentum = 0.0;
    for (std::size_t row = 0U; row < 2U; ++row) {
      angular_momentum += 0.5 * (orbit.position[row][0] * orbit.velocity[row][1] -
                                 orbit.position[row][1] * orbit.velocity[row][0]);
    }
    CHUI_REQUIRE(std::abs(angular_momentum - 0.25) < 1e-12);
  }
  requireOrder("two_body_circular", errors, 1.8, 2.2);
}


// Einstein-de Sitter background, code time in 1/H0:
// a(t)^{3/2} = a0^{3/2} + 1.5 (t - t0), integral(dt/a) = 2 (sqrt(a1) - sqrt(a0)).
LambdaCdmBackground edsBackground() {
  CosmologyBackgroundConfig config;
  config.omega_matter = 1.0;
  config.omega_lambda = 0.0;
  config.omega_radiation = 0.0;
  config.omega_curvature = 0.0;
  return LambdaCdmBackground(config);
}

// Cosmological kick coefficients and long-range kick ownership. With constant
// comoving accelerations the canonical momentum p = a u is exact:
//   p(a1) = p0 + (A_tree + A_PM) * 2 (sqrt(a1) - sqrt(a0)),
// for any block tiling, rung history and PM cycle length. A duplicated or
// missing PM half-kick, a kick factor at the wrong epoch, or a gap/overlap
// between the opening [a0, a_mid] and closing [a_mid, a_end] PM integrals
// changes p by O(A * integral) and fails this test. Positions converge at
// second order to the analytic integral of p a^{-3/2} da.
void testCosmologicalKickCoefficients() {
  const auto background = edsBackground();
  const double a0 = 0.5, s0 = std::sqrt(a0), u0_z = 0.3;
  const Vec3 tree_acceleration{1.0, 0.0, 0.0}, pm_acceleration{0.0, 2.0, 0.0};
  Problem problem;
  problem.initialize = [&](SimulationState& state) { state.particles.velocity_z_peculiar[0] = u0_z; };
  problem.force = [&](const SimulationState&, std::size_t, Vec3& tree, Vec3& pm) {
    tree = tree_acceleration;
    pm = pm_acceleration;
  };
  problem.rung = [](std::size_t, std::uint64_t block) { return static_cast<std::uint8_t>(block % 3U); };
  const Cosmology cosmology{&background, a0, 0.0};
  const double span = 0.08;
  for (const auto& fractions : std::vector<std::vector<double>>{k_uniform4, {0.1, 0.4, 0.2, 0.3}}) {
    for (const bool multirate : {true, false}) {
      std::vector<double> position_errors;
      for (const std::size_t cycles : {8U, 16U, 32U}) {
        const auto outcome = integrate(problem, span / static_cast<double>(cycles), cycles, fractions, multirate, cosmology);
        const double a1 = outcome.scale_factor, s1 = std::sqrt(a1);
        CHUI_REQUIRE(std::abs(outcome.time - span) < 1e-12);
        CHUI_REQUIRE(std::abs(a1 - std::pow(std::pow(a0, 1.5) + 1.5 * span, 2.0 / 3.0)) < 1e-7);
        const double integral = 2.0 * (s1 - s0);
        const double p0[3] = {0.0, 0.0, a0 * u0_z};
        double position_error = 0.0;
        for (std::size_t axis = 0U; axis < 3U; ++axis) {
          const double acceleration = tree_acceleration[axis] + pm_acceleration[axis];
          const double p_expected = p0[axis] + acceleration * integral;
          const double p = a1 * outcome.velocity[0][axis];
          CHUI_REQUIRE(std::abs(p - p_expected) <= 1e-8 * (std::abs(acceleration) * integral + std::abs(p0[axis])) + 1e-15);
          const double x_expected = 2.0 * (p0[axis] - 2.0 * acceleration * s0) * (1.0 / s0 - 1.0 / s1) +
                                    4.0 * acceleration * std::log(s1 / s0);
          position_error = std::max(position_error, std::abs(outcome.position[0][axis] - x_expected));
        }
        position_errors.push_back(position_error);
      }
      if (position_errors.front() > 1e-12) {
        requireOrder(multirate ? "cosmological_positions_multirate" : "cosmological_positions_reference",
                     position_errors, 1.8, 2.2);
      } else {
        for (const double error : position_errors) CHUI_REQUIRE(error < 1e-10);
      }
    }
  }
}

// Cosmological split oscillator (comoving): self-convergence of the multirate
// composition, with rung changes, against a 32x finer reference-scheme run.
void testCosmologicalOscillatorConvergence() {
  const auto background = edsBackground();
  const std::vector<double> ws{6.0, 2.0}, wl{1.0, 1.0};
  auto problem = oscillators(ws, wl);
  problem.rung = [](std::size_t row, std::uint64_t block) {
    return row == 0U ? std::uint8_t{0U} : static_cast<std::uint8_t>(block % 2U == 0U ? 2U : 1U);
  };
  const Cosmology cosmology{&background, 0.5, 0.0};
  const double span = 0.2;
  const auto reference = integrate(problem, span / 2560.0, 2560U, {1.0}, false, cosmology);
  std::vector<double> errors;
  for (const std::size_t cycles : {20U, 40U, 80U}) {
    const auto multirate = integrate(problem, span / static_cast<double>(cycles), cycles, k_uniform4, true, cosmology);
    CHUI_REQUIRE(std::abs(multirate.scale_factor - reference.scale_factor) < 1e-9);
    errors.push_back(difference(multirate, reference));
  }
  requireOrder("cosmological_oscillator", errors, 1.8, 2.2);
}
}  // namespace

int main() {
  try {
    testPmOnly();
    testTreeOnly();
    testCombinedAndVariableBlocks();
    testAdaptiveRungs();
    testTwoBodyCircularOrbit();
    testCosmologicalKickCoefficients();
    testCosmologicalOscillatorConvergence();
  } catch (const std::exception& error) {
    std::cerr << "hierarchical multirate convergence: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
