#include "cosmosim/physics/metal_diffusion.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace cosmosim::physics {
namespace {

constexpr double k_small = 1.0e-30;

[[nodiscard]] bool finiteNonnegative(double value) noexcept {
  return std::isfinite(value) && value >= 0.0;
}

[[nodiscard]] double totalMetalMass(std::span<const double> metal_mass) {
  return std::accumulate(metal_mass.begin(), metal_mass.end(), 0.0);
}

[[nodiscard]] double harmonicMean(double left, double right) noexcept {
  if (!(left > 0.0) || !(right > 0.0)) {
    return 0.0;
  }
  return 2.0 * left * right / (left + right);
}


[[nodiscard]] double faceConductanceFromRhoKappa(
    std::span<const double> rho_kappa_code,
    const MetalDiffusionFace& face) noexcept {
  if (face.boundary_kind == MetalDiffusionBoundaryKind::kReflective ||
      face.boundary_kind == MetalDiffusionBoundaryKind::kOpen ||
      face.left_cell >= rho_kappa_code.size() ||
      face.right_cell >= rho_kappa_code.size() ||
      !(face.area_code > 0.0) || !(face.center_distance_code > 0.0)) {
    return 0.0;
  }
  const double rho_kappa_face = harmonicMean(
      rho_kappa_code[face.left_cell], rho_kappa_code[face.right_cell]);
  return rho_kappa_face * face.area_code / face.center_distance_code;
}

void validateCompactInputs(
    const MetalDiffusionConfig& config,
    MetalDiffusionFieldView view,
    std::span<const MetalDiffusionFace> faces,
    double dt_code) {
  if (!(dt_code >= 0.0) || !std::isfinite(dt_code)) {
    throw std::invalid_argument(
        "MetalDiffusionModel: dt_code must be finite and nonnegative");
  }
  if (view.gas_mass_code.size() != view.metal_mass_code.size() ||
      view.rho_kappa_code.size() != view.gas_mass_code.size()) {
    throw std::invalid_argument(
        "MetalDiffusionModel: compact field extents disagree");
  }
  for (std::size_t i = 0U; i < view.gas_mass_code.size(); ++i) {
    if (!finiteNonnegative(view.gas_mass_code[i]) ||
        !finiteNonnegative(view.metal_mass_code[i]) ||
        !finiteNonnegative(view.rho_kappa_code[i])) {
      throw std::invalid_argument(
          "MetalDiffusionModel: compact field contains invalid state");
    }
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, view.gas_mass_code[i]);
    if (view.metal_mass_code[i] > view.gas_mass_code[i] + tolerance) {
      throw std::invalid_argument(
          "MetalDiffusionModel: metal mass exceeds gas mass");
    }
  }
  for (const auto& face : faces) {
    if (face.left_cell >= view.gas_mass_code.size() ||
        face.right_cell >= view.gas_mass_code.size() ||
        !finiteNonnegative(face.area_code) ||
        !finiteNonnegative(face.center_distance_code)) {
      throw std::invalid_argument(
          "MetalDiffusionModel: invalid compact face geometry or index");
    }
  }
  if (config.enabled && config.model == core::MetalDiffusionModel::kNone) {
    throw std::invalid_argument(
        "MetalDiffusionModel: enabled diffusion requires a model");
  }
}

[[nodiscard]] std::uint64_t checkedWorkspaceBytes(
    std::size_t count, std::size_t buffer_count) {
  if (count != 0U && sizeof(double) >
      std::numeric_limits<std::uint64_t>::max() / count) {
    throw std::overflow_error("MetalDiffusionModel workspace byte count overflow");
  }
  const std::uint64_t one = static_cast<std::uint64_t>(count) * sizeof(double);
  if (buffer_count != 0U && one >
      std::numeric_limits<std::uint64_t>::max() / buffer_count) {
    throw std::overflow_error("MetalDiffusionModel workspace byte count overflow");
  }
  return one * static_cast<std::uint64_t>(buffer_count);
}

[[nodiscard]] double faceConductance(
    const MetalDiffusionConfig& config,
    std::span<const MetalDiffusionCell> cells,
    const MetalDiffusionFace& face) noexcept {
  if (face.boundary_kind == MetalDiffusionBoundaryKind::kReflective ||
      face.boundary_kind == MetalDiffusionBoundaryKind::kOpen ||
      face.left_cell >= cells.size() || face.right_cell >= cells.size() ||
      !(face.area_code > 0.0) || !(face.center_distance_code > 0.0)) {
    return 0.0;
  }
  const MetalDiffusionCell& left = cells[face.left_cell];
  const MetalDiffusionCell& right = cells[face.right_cell];
  if (!left.is_owned_leaf || !right.is_owned_leaf) {
    return 0.0;
  }
  const double left_rho_kappa = left.density_code *
      smagorinskyMetalDiffusivityCode(config, left);
  const double right_rho_kappa = right.density_code *
      smagorinskyMetalDiffusivityCode(config, right);
  const double rho_kappa_face = harmonicMean(left_rho_kappa, right_rho_kappa);
  return rho_kappa_face * face.area_code / face.center_distance_code;
}

void validateInputs(
    const MetalDiffusionConfig& config,
    std::span<const MetalDiffusionCell> cells,
    std::span<const MetalDiffusionFace> faces,
    double dt_code) {
  if (!(dt_code >= 0.0) || !std::isfinite(dt_code)) {
    throw std::invalid_argument("MetalDiffusionModel: dt_code must be finite and nonnegative");
  }
  for (std::size_t i = 0; i < cells.size(); ++i) {
    const auto& cell = cells[i];
    if (!finiteNonnegative(cell.gas_mass_code) ||
        !finiteNonnegative(cell.metal_mass_code) ||
        !finiteNonnegative(cell.density_code) ||
        !finiteNonnegative(cell.volume_code) ||
        !finiteNonnegative(cell.filter_length_code)) {
      throw std::invalid_argument("MetalDiffusionModel: non-finite or negative cell state");
    }
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
        std::max(1.0, cell.gas_mass_code);
    if (cell.metal_mass_code > cell.gas_mass_code + tolerance) {
      throw std::invalid_argument("MetalDiffusionModel: metal mass exceeds gas mass");
    }
  }
  for (const auto& face : faces) {
    if (face.left_cell >= cells.size() || face.right_cell >= cells.size() ||
        !finiteNonnegative(face.area_code) ||
        !finiteNonnegative(face.center_distance_code)) {
      throw std::invalid_argument("MetalDiffusionModel: invalid face geometry or index");
    }
  }
  if (config.enabled && config.model == core::MetalDiffusionModel::kNone) {
    throw std::invalid_argument("MetalDiffusionModel: enabled diffusion requires a model");
  }
}

[[nodiscard]] std::uint32_t rklStageCount(double dt, double explicit_dt, std::uint32_t maximum) {
  if (!(dt > explicit_dt) || !(explicit_dt > 0.0)) {
    return 1U;
  }
  const double ratio = dt / explicit_dt;
  // dt / dt_explicit <= (s^2+s-2)/4 for RKL2.
  const double root = 0.5 * (-1.0 + std::sqrt(9.0 + 16.0 * ratio));
  const std::uint32_t stages = static_cast<std::uint32_t>(std::ceil(root));
  return std::clamp(stages, 2U, maximum);
}

}  // namespace

double traceFreeStrainMagnitude(
    const MetalDiffusionVelocityGradient& gradient) noexcept {
  double symmetric[3][3]{};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      symmetric[i][j] = 0.5 * (gradient.grad[i][j] + gradient.grad[j][i]);
    }
  }
  const double trace_third =
      (symmetric[0][0] + symmetric[1][1] + symmetric[2][2]) / 3.0;
  for (int i = 0; i < 3; ++i) {
    symmetric[i][i] -= trace_third;
  }
  double contraction = 0.0;
  for (const auto& row : symmetric) {
    for (const double value : row) {
      contraction += value * value;
    }
  }
  return std::sqrt(std::max(0.0, 2.0 * contraction));
}

double smagorinskyMetalDiffusivityCode(
    const MetalDiffusionConfig& config,
    const MetalDiffusionCell& cell) noexcept {
  if (!config.enabled || config.model != core::MetalDiffusionModel::kSmagorinsky ||
      !(cell.filter_length_code > 0.0)) {
    return 0.0;
  }
  const double strain = traceFreeStrainMagnitude(cell.velocity_gradient);
  if (!std::isfinite(strain)) {
    return 0.0;
  }
  const double raw = config.smagorinsky_coefficient *
      cell.filter_length_code * cell.filter_length_code * strain;
  return std::clamp(
      raw, config.diffusivity_floor_code, config.diffusivity_ceiling_code);
}


double smagorinskyRhoKappaCode(
    const MetalDiffusionConfig& config,
    double density_code,
    double filter_length_code,
    double trace_free_strain_magnitude_code) noexcept {
  if (!config.enabled || config.model != core::MetalDiffusionModel::kSmagorinsky ||
      !(density_code > 0.0) || !(filter_length_code > 0.0) ||
      !std::isfinite(trace_free_strain_magnitude_code)) {
    return 0.0;
  }
  const double raw = config.smagorinsky_coefficient *
      filter_length_code * filter_length_code *
      trace_free_strain_magnitude_code;
  const double kappa = std::clamp(
      raw, config.diffusivity_floor_code, config.diffusivity_ceiling_code);
  return density_code * kappa;
}

std::uint64_t MetalDiffusionWorkspace::ownedCapacityBytes() const noexcept {
  std::uint64_t total = 0U;
  for (const auto& buffer : m_buffers) {
    const std::uint64_t bytes =
        static_cast<std::uint64_t>(buffer.capacity()) * sizeof(double);
    if (bytes > std::numeric_limits<std::uint64_t>::max() - total) {
      return std::numeric_limits<std::uint64_t>::max();
    }
    total += bytes;
  }
  return total;
}

std::uint64_t MetalDiffusionWorkspace::highWaterBytes() const noexcept {
  return m_high_water_bytes;
}

MetalDiffusionModel::MetalDiffusionModel(MetalDiffusionConfig config)
    : m_config(std::move(config)) {
  if (m_config.smagorinsky_coefficient < 0.0 ||
      !std::isfinite(m_config.smagorinsky_coefficient) ||
      !(m_config.parabolic_cfl > 0.0) || m_config.parabolic_cfl > 0.5 ||
      m_config.max_subcycles == 0U || m_config.max_rkl_stages < 2U ||
      m_config.diffusivity_floor_code < 0.0 ||
      !(m_config.diffusivity_ceiling_code > m_config.diffusivity_floor_code)) {
    throw std::invalid_argument("MetalDiffusionModel: invalid configuration");
  }
}

const MetalDiffusionConfig& MetalDiffusionModel::config() const noexcept {
  return m_config;
}

double MetalDiffusionModel::stableExplicitTimestepCode(
    std::span<const MetalDiffusionCell> cells,
    std::span<const MetalDiffusionFace> faces) const {
  if (!m_config.enabled || cells.empty() || faces.empty()) {
    return std::numeric_limits<double>::infinity();
  }
  std::vector<double> coefficient_sum(cells.size(), 0.0);
  for (const auto& face : faces) {
    const double conductance = faceConductance(m_config, cells, face);
    if (!(conductance > 0.0)) {
      continue;
    }
    coefficient_sum[face.left_cell] += conductance;
    coefficient_sum[face.right_cell] += conductance;
  }
  double stable = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < cells.size(); ++i) {
    if (coefficient_sum[i] > 0.0 && cells[i].gas_mass_code > 0.0) {
      stable = std::min(
          stable,
          m_config.parabolic_cfl * cells[i].gas_mass_code / coefficient_sum[i]);
    }
  }
  return stable;
}

std::vector<double> MetalDiffusionModel::evaluateDerivative(
    std::span<const MetalDiffusionCell> cells,
    std::span<const MetalDiffusionFace> faces,
    std::span<const double> metal_mass_code,
    std::uint64_t* faces_evaluated) const {
  std::vector<double> derivative(cells.size(), 0.0);
  for (const auto& face : faces) {
    const double conductance = faceConductance(m_config, cells, face);
    if (!(conductance > 0.0)) {
      continue;
    }
    const auto left = static_cast<std::size_t>(face.left_cell);
    const auto right = static_cast<std::size_t>(face.right_cell);
    const double z_left = cells[left].gas_mass_code > 0.0
        ? metal_mass_code[left] / cells[left].gas_mass_code : 0.0;
    const double z_right = cells[right].gas_mass_code > 0.0
        ? metal_mass_code[right] / cells[right].gas_mass_code : 0.0;
    const double transfer_rate_left_to_right = conductance * (z_left - z_right);
    derivative[left] -= transfer_rate_left_to_right;
    derivative[right] += transfer_rate_left_to_right;
    if (faces_evaluated != nullptr) {
      ++(*faces_evaluated);
    }
  }
  return derivative;
}

void MetalDiffusionModel::conservativeLimitedEuler(
    std::span<const MetalDiffusionCell> cells,
    std::span<const MetalDiffusionFace> faces,
    std::span<double> metal_mass_code,
    double dt_code,
    std::uint64_t* faces_evaluated,
    std::uint64_t* limited_faces) const {
  struct RequestedTransfer {
    std::uint32_t donor = 0;
    std::uint32_t receiver = 0;
    double amount = 0.0;
  };
  std::vector<RequestedTransfer> transfers;
  transfers.reserve(faces.size());
  std::vector<double> outgoing(cells.size(), 0.0);
  std::vector<double> incoming(cells.size(), 0.0);

  for (const auto& face : faces) {
    const double conductance = faceConductance(m_config, cells, face);
    if (!(conductance > 0.0)) {
      continue;
    }
    const std::uint32_t left = face.left_cell;
    const std::uint32_t right = face.right_cell;
    const double z_left = cells[left].gas_mass_code > 0.0
        ? metal_mass_code[left] / cells[left].gas_mass_code : 0.0;
    const double z_right = cells[right].gas_mass_code > 0.0
        ? metal_mass_code[right] / cells[right].gas_mass_code : 0.0;
    const double signed_amount = dt_code * conductance * (z_left - z_right);
    if (signed_amount == 0.0) {
      continue;
    }
    RequestedTransfer transfer;
    if (signed_amount > 0.0) {
      transfer = {left, right, signed_amount};
    } else {
      transfer = {right, left, -signed_amount};
    }
    outgoing[transfer.donor] += transfer.amount;
    incoming[transfer.receiver] += transfer.amount;
    transfers.push_back(transfer);
    if (faces_evaluated != nullptr) {
      ++(*faces_evaluated);
    }
  }

  std::vector<double> donor_scale(cells.size(), 1.0);
  std::vector<double> receiver_scale(cells.size(), 1.0);
  for (std::size_t i = 0; i < cells.size(); ++i) {
    if (outgoing[i] > metal_mass_code[i] && outgoing[i] > 0.0) {
      donor_scale[i] = metal_mass_code[i] / outgoing[i];
    }
    const double capacity = std::max(0.0, cells[i].gas_mass_code - metal_mass_code[i]);
    if (incoming[i] > capacity && incoming[i] > 0.0) {
      receiver_scale[i] = capacity / incoming[i];
    }
  }

  for (const auto& transfer : transfers) {
    const double scale = std::min(
        donor_scale[transfer.donor], receiver_scale[transfer.receiver]);
    const double amount = transfer.amount * scale;
    if (scale < 1.0 && limited_faces != nullptr) {
      ++(*limited_faces);
    }
    metal_mass_code[transfer.donor] -= amount;
    metal_mass_code[transfer.receiver] += amount;
  }
}

void MetalDiffusionModel::explicitSspRk2(
    std::span<const MetalDiffusionCell> cells,
    std::span<const MetalDiffusionFace> faces,
    std::span<double> metal_mass_code,
    double dt_code,
    std::uint64_t* faces_evaluated,
    std::uint64_t* limited_faces) const {
  std::vector<double> initial(metal_mass_code.begin(), metal_mass_code.end());
  std::vector<double> stage(initial);
  conservativeLimitedEuler(
      cells, faces, stage, dt_code, faces_evaluated, limited_faces);
  conservativeLimitedEuler(
      cells, faces, stage, dt_code, faces_evaluated, limited_faces);
  for (std::size_t i = 0; i < stage.size(); ++i) {
    metal_mass_code[i] = 0.5 * initial[i] + 0.5 * stage[i];
  }
}

void MetalDiffusionModel::rkl2(
    std::span<const MetalDiffusionCell> cells,
    std::span<const MetalDiffusionFace> faces,
    std::span<double> metal_mass_code,
    double dt_code,
    std::uint32_t stages,
    std::uint64_t* faces_evaluated) const {
  if (stages <= 1U) {
    const auto derivative = evaluateDerivative(cells, faces, metal_mass_code, faces_evaluated);
    for (std::size_t i = 0; i < metal_mass_code.size(); ++i) {
      metal_mass_code[i] += dt_code * derivative[i];
    }
    return;
  }

  const double s = static_cast<double>(stages);
  const double w1 = 4.0 / (s * s + s - 2.0);
  std::vector<double> y0(metal_mass_code.begin(), metal_mass_code.end());
  const std::vector<double> f0 = evaluateDerivative(cells, faces, y0, faces_evaluated);
  std::vector<double> yjm2(y0);
  std::vector<double> yjm1(y0.size());
  constexpr double b1 = 1.0 / 3.0;
  for (std::size_t i = 0; i < y0.size(); ++i) {
    yjm1[i] = y0[i] + b1 * w1 * dt_code * f0[i];
  }

  auto b = [](std::uint32_t j) {
    if (j <= 1U) {
      return 1.0 / 3.0;
    }
    const double jd = static_cast<double>(j);
    return (jd * jd + jd - 2.0) / (2.0 * jd * (jd + 1.0));
  };
  auto a = [&b](std::uint32_t j) { return 1.0 - b(j); };

  for (std::uint32_t j = 2U; j <= stages; ++j) {
    const double jd = static_cast<double>(j);
    const double mu = ((2.0 * jd - 1.0) / jd) * b(j) / b(j - 1U);
    const double nu = -((jd - 1.0) / jd) * b(j) / b(j - 2U);
    const double mu_tilde = mu * w1;
    const double gamma_tilde = -a(j - 1U) * mu_tilde;
    const std::vector<double> fjm1 =
        evaluateDerivative(cells, faces, yjm1, faces_evaluated);
    std::vector<double> yj(y0.size(), 0.0);
    for (std::size_t i = 0; i < y0.size(); ++i) {
      yj[i] = mu * yjm1[i] + nu * yjm2[i] + (1.0 - mu - nu) * y0[i] +
          mu_tilde * dt_code * fjm1[i] + gamma_tilde * dt_code * f0[i];
    }
    yjm2.swap(yjm1);
    yjm1.swap(yj);
  }
  std::copy(yjm1.begin(), yjm1.end(), metal_mass_code.begin());
}


std::uint64_t MetalDiffusionModel::requiredWorkspaceBytes(
    std::size_t cell_count) const {
  // Five scalar lanes are sufficient for both integrators. Explicit SSPRK2
  // uses the fifth lane as an immutable face-transfer source snapshot so the
  // two-pass limiter is algebraically identical to the legacy transfer list.
  return checkedWorkspaceBytes(cell_count, 5U);
}

MetalDiffusionStepReport MetalDiffusionModel::advanceFromView(
    MetalDiffusionFieldView view,
    std::span<const MetalDiffusionFace> faces,
    double dt_code,
    MetalDiffusionWorkspace& workspace) const {
  validateCompactInputs(m_config, view, faces, dt_code);
  MetalDiffusionStepReport report;
  const std::size_t count = view.gas_mass_code.size();
  const std::size_t required_buffer_count = 5U;
  for (std::size_t i = 0U; i < required_buffer_count; ++i) {
    workspace.m_buffers[i].resize(count);
  }
  workspace.m_high_water_bytes = std::max(
      workspace.m_high_water_bytes, workspace.ownedCapacityBytes());

  auto& working = workspace.m_buffers[0];
  std::copy(view.metal_mass_code.begin(), view.metal_mass_code.end(), working.begin());
  report.metal_mass_before_code = totalMetalMass(working);

  auto& scratch_a = workspace.m_buffers[2];
  std::fill(scratch_a.begin(), scratch_a.end(), 0.0);
  for (const auto& face : faces) {
    const double conductance = faceConductanceFromRhoKappa(view.rho_kappa_code, face);
    if (!(conductance > 0.0)) {
      continue;
    }
    scratch_a[face.left_cell] += conductance;
    scratch_a[face.right_cell] += conductance;
  }
  report.stable_dt_code = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0U; i < count; ++i) {
    if (scratch_a[i] > 0.0 && view.gas_mass_code[i] > 0.0) {
      report.stable_dt_code = std::min(
          report.stable_dt_code,
          m_config.parabolic_cfl * view.gas_mass_code[i] / scratch_a[i]);
    }
  }

  const auto limitedEuler = [&](std::span<double> metal_mass, double dt) {
    auto& outgoing = workspace.m_buffers[2];
    auto& incoming = workspace.m_buffers[3];
    auto& transfer_source = workspace.m_buffers[4];
    std::copy(metal_mass.begin(), metal_mass.end(), transfer_source.begin());
    std::fill(outgoing.begin(), outgoing.end(), 0.0);
    std::fill(incoming.begin(), incoming.end(), 0.0);
    for (const auto& face : faces) {
      const double conductance = faceConductanceFromRhoKappa(
          view.rho_kappa_code, face);
      if (!(conductance > 0.0)) {
        continue;
      }
      const std::uint32_t left = face.left_cell;
      const std::uint32_t right = face.right_cell;
      const double z_left = view.gas_mass_code[left] > 0.0
          ? transfer_source[left] / view.gas_mass_code[left] : 0.0;
      const double z_right = view.gas_mass_code[right] > 0.0
          ? transfer_source[right] / view.gas_mass_code[right] : 0.0;
      const double signed_amount = dt * conductance * (z_left - z_right);
      if (signed_amount > 0.0) {
        outgoing[left] += signed_amount;
        incoming[right] += signed_amount;
        ++report.faces_evaluated;
      } else if (signed_amount < 0.0) {
        outgoing[right] -= signed_amount;
        incoming[left] -= signed_amount;
        ++report.faces_evaluated;
      }
    }
    for (std::size_t i = 0U; i < count; ++i) {
      const double donor_scale =
          outgoing[i] > metal_mass[i] && outgoing[i] > 0.0
          ? metal_mass[i] / outgoing[i] : 1.0;
      const double capacity = std::max(0.0, view.gas_mass_code[i] - metal_mass[i]);
      const double receiver_scale =
          incoming[i] > capacity && incoming[i] > 0.0
          ? capacity / incoming[i] : 1.0;
      outgoing[i] = donor_scale;
      incoming[i] = receiver_scale;
    }
    for (const auto& face : faces) {
      const double conductance = faceConductanceFromRhoKappa(
          view.rho_kappa_code, face);
      if (!(conductance > 0.0)) {
        continue;
      }
      const std::uint32_t left = face.left_cell;
      const std::uint32_t right = face.right_cell;
      const double z_left = view.gas_mass_code[left] > 0.0
          ? transfer_source[left] / view.gas_mass_code[left] : 0.0;
      const double z_right = view.gas_mass_code[right] > 0.0
          ? transfer_source[right] / view.gas_mass_code[right] : 0.0;
      const double signed_amount = dt * conductance * (z_left - z_right);
      if (signed_amount == 0.0) {
        continue;
      }
      const std::uint32_t donor = signed_amount > 0.0 ? left : right;
      const std::uint32_t receiver = signed_amount > 0.0 ? right : left;
      const double requested = std::abs(signed_amount);
      const double scale = std::min(outgoing[donor], incoming[receiver]);
      if (scale < 1.0) {
        ++report.limited_faces;
      }
      const double amount = requested * scale;
      metal_mass[donor] -= amount;
      metal_mass[receiver] += amount;
    }
  };

  const auto derivativeInto = [&](std::span<const double> metal_mass,
                                  std::span<double> derivative) {
    std::fill(derivative.begin(), derivative.end(), 0.0);
    for (const auto& face : faces) {
      const double conductance = faceConductanceFromRhoKappa(
          view.rho_kappa_code, face);
      if (!(conductance > 0.0)) {
        continue;
      }
      const std::uint32_t left = face.left_cell;
      const std::uint32_t right = face.right_cell;
      const double z_left = view.gas_mass_code[left] > 0.0
          ? metal_mass[left] / view.gas_mass_code[left] : 0.0;
      const double z_right = view.gas_mass_code[right] > 0.0
          ? metal_mass[right] / view.gas_mass_code[right] : 0.0;
      const double rate = conductance * (z_left - z_right);
      derivative[left] -= rate;
      derivative[right] += rate;
      ++report.faces_evaluated;
    }
  };

  if (!m_config.enabled || dt_code == 0.0 ||
      !std::isfinite(report.stable_dt_code)) {
    report.metal_mass_after_code = report.metal_mass_before_code;
  } else if (m_config.time_integrator ==
             core::MetalDiffusionTimeIntegrator::kExplicitSubcycling) {
    const double ratio = dt_code / report.stable_dt_code;
    const std::uint32_t subcycles = std::max(
        1U, static_cast<std::uint32_t>(std::ceil(std::max(1.0, ratio))));
    if (subcycles > m_config.max_subcycles) {
      throw std::runtime_error(
          "MetalDiffusionModel: explicit subcycling exceeds configured maximum");
    }
    report.subcycles = subcycles;
    const double sub_dt = dt_code / static_cast<double>(subcycles);
    auto& initial = workspace.m_buffers[1];
    for (std::uint32_t subcycle = 0U; subcycle < subcycles; ++subcycle) {
      std::copy(working.begin(), working.end(), initial.begin());
      limitedEuler(working, sub_dt);
      limitedEuler(working, sub_dt);
      for (std::size_t i = 0U; i < count; ++i) {
        working[i] = 0.5 * initial[i] + 0.5 * working[i];
      }
    }
  } else {
    const std::uint32_t stages = rklStageCount(
        dt_code, report.stable_dt_code, m_config.max_rkl_stages);
    const double capacity = report.stable_dt_code *
        (static_cast<double>(stages) * stages + stages - 2.0) / 4.0;
    if (capacity + 64.0 * std::numeric_limits<double>::epsilon() * dt_code <
        dt_code) {
      throw std::runtime_error("MetalDiffusionModel: RKL2 stage limit is insufficient");
    }
    report.rkl_stages = stages;
    if (stages <= 1U) {
      auto& derivative = workspace.m_buffers[4];
      derivativeInto(working, derivative);
      for (std::size_t i = 0U; i < count; ++i) {
        working[i] += dt_code * derivative[i];
      }
    } else {
      auto& f0 = workspace.m_buffers[1];
      auto& yjm2 = workspace.m_buffers[2];
      auto& yjm1 = workspace.m_buffers[3];
      auto& derivative = workspace.m_buffers[4];
      derivativeInto(working, f0);
      std::copy(working.begin(), working.end(), yjm2.begin());
      constexpr double b1 = 1.0 / 3.0;
      const double sd = static_cast<double>(stages);
      const double w1 = 4.0 / (sd * sd + sd - 2.0);
      for (std::size_t i = 0U; i < count; ++i) {
        yjm1[i] = working[i] + b1 * w1 * dt_code * f0[i];
      }
      const auto b = [](std::uint32_t j) {
        if (j <= 1U) {
          return 1.0 / 3.0;
        }
        const double jd = static_cast<double>(j);
        return (jd * jd + jd - 2.0) / (2.0 * jd * (jd + 1.0));
      };
      const auto a = [&b](std::uint32_t j) { return 1.0 - b(j); };
      for (std::uint32_t j = 2U; j <= stages; ++j) {
        const double jd = static_cast<double>(j);
        const double mu = ((2.0 * jd - 1.0) / jd) * b(j) / b(j - 1U);
        const double nu = -((jd - 1.0) / jd) * b(j) / b(j - 2U);
        const double mu_tilde = mu * w1;
        const double gamma_tilde = -a(j - 1U) * mu_tilde;
        derivativeInto(yjm1, derivative);
        for (std::size_t i = 0U; i < count; ++i) {
          yjm2[i] = mu * yjm1[i] + nu * yjm2[i] +
              (1.0 - mu - nu) * working[i] +
              mu_tilde * dt_code * derivative[i] +
              gamma_tilde * dt_code * f0[i];
        }
        yjm2.swap(yjm1);
      }
      std::copy(yjm1.begin(), yjm1.end(), working.begin());
    }
    for (std::size_t i = 0U; i < count; ++i) {
      const double tolerance = 256.0 * std::numeric_limits<double>::epsilon() *
          std::max(1.0, view.gas_mass_code[i]);
      if (working[i] < -tolerance ||
          working[i] > view.gas_mass_code[i] + tolerance) {
        throw std::runtime_error(
            "MetalDiffusionModel: RKL2 produced an unbounded state; use explicit_subcycling");
      }
      working[i] = std::clamp(working[i], 0.0, view.gas_mass_code[i]);
    }
  }

  report.metal_mass_after_code = totalMetalMass(working);
  report.conservation_residual_code = report.metal_mass_after_code -
      report.metal_mass_before_code + report.open_boundary_loss_code;
  report.minimum_metallicity = std::numeric_limits<double>::infinity();
  report.maximum_metallicity = 0.0;
  for (std::size_t i = 0U; i < count; ++i) {
    const double metallicity = view.gas_mass_code[i] > 0.0
        ? working[i] / view.gas_mass_code[i] : 0.0;
    report.minimum_metallicity = std::min(report.minimum_metallicity, metallicity);
    report.maximum_metallicity = std::max(report.maximum_metallicity, metallicity);
  }
  if (count == 0U) {
    report.minimum_metallicity = 0.0;
  }
  std::copy(working.begin(), working.end(), view.metal_mass_code.begin());
  return report;
}

MetalDiffusionStepReport MetalDiffusionModel::advance(
    std::span<MetalDiffusionCell> cells,
    std::span<const MetalDiffusionFace> faces,
    double dt_code) const {
  validateInputs(m_config, cells, faces, dt_code);
  MetalDiffusionStepReport report;
  std::vector<double> metal_mass(cells.size(), 0.0);
  for (std::size_t i = 0; i < cells.size(); ++i) {
    metal_mass[i] = cells[i].metal_mass_code;
  }
  report.metal_mass_before_code = totalMetalMass(metal_mass);
  report.stable_dt_code = stableExplicitTimestepCode(cells, faces);

  if (!m_config.enabled || dt_code == 0.0 || !std::isfinite(report.stable_dt_code)) {
    report.metal_mass_after_code = report.metal_mass_before_code;
  } else if (m_config.time_integrator ==
             core::MetalDiffusionTimeIntegrator::kExplicitSubcycling) {
    const double ratio = dt_code / report.stable_dt_code;
    const std::uint32_t subcycles = std::max(
        1U, static_cast<std::uint32_t>(std::ceil(std::max(1.0, ratio))));
    if (subcycles > m_config.max_subcycles) {
      throw std::runtime_error(
          "MetalDiffusionModel: explicit subcycling exceeds configured maximum");
    }
    report.subcycles = subcycles;
    const double sub_dt = dt_code / static_cast<double>(subcycles);
    for (std::uint32_t subcycle = 0; subcycle < subcycles; ++subcycle) {
      explicitSspRk2(
          cells, faces, metal_mass, sub_dt,
          &report.faces_evaluated, &report.limited_faces);
    }
  } else {
    const std::uint32_t stages = rklStageCount(
        dt_code, report.stable_dt_code, m_config.max_rkl_stages);
    const double capacity = report.stable_dt_code *
        (static_cast<double>(stages) * stages + stages - 2.0) / 4.0;
    if (capacity + 64.0 * std::numeric_limits<double>::epsilon() * dt_code < dt_code) {
      throw std::runtime_error("MetalDiffusionModel: RKL2 stage limit is insufficient");
    }
    report.rkl_stages = stages;
    rkl2(cells, faces, metal_mass, dt_code, stages, &report.faces_evaluated);
    // RKL2 is used only while its linear stage result remains physically bounded.
    // The code fails closed rather than silently clipping a non-monotone state.
    for (std::size_t i = 0; i < cells.size(); ++i) {
      const double tolerance = 256.0 * std::numeric_limits<double>::epsilon() *
          std::max(1.0, cells[i].gas_mass_code);
      if (metal_mass[i] < -tolerance || metal_mass[i] > cells[i].gas_mass_code + tolerance) {
        throw std::runtime_error(
            "MetalDiffusionModel: RKL2 produced an unbounded state; use explicit_subcycling");
      }
      metal_mass[i] = std::clamp(metal_mass[i], 0.0, cells[i].gas_mass_code);
    }
  }

  report.metal_mass_after_code = totalMetalMass(metal_mass);
  report.conservation_residual_code = report.metal_mass_after_code -
      report.metal_mass_before_code + report.open_boundary_loss_code;
  report.minimum_metallicity = std::numeric_limits<double>::infinity();
  report.maximum_metallicity = 0.0;
  for (std::size_t i = 0; i < cells.size(); ++i) {
    cells[i].metal_mass_code = metal_mass[i];
    const double metallicity = cells[i].gas_mass_code > 0.0
        ? metal_mass[i] / cells[i].gas_mass_code : 0.0;
    report.minimum_metallicity = std::min(report.minimum_metallicity, metallicity);
    report.maximum_metallicity = std::max(report.maximum_metallicity, metallicity);
  }
  if (cells.empty()) {
    report.minimum_metallicity = 0.0;
  }
  return report;
}

MetalDiffusionConfig makeMetalDiffusionConfig(
    const core::PhysicsConfig& physics_config) {
  return MetalDiffusionConfig{
      .enabled = physics_config.enable_metal_diffusion,
      .model = physics_config.metal_diffusion_model,
      .time_integrator = physics_config.metal_diffusion_time_integrator,
      .smagorinsky_coefficient = physics_config.metal_diffusion_coefficient,
      .parabolic_cfl = physics_config.metal_diffusion_cfl,
      .max_subcycles = physics_config.metal_diffusion_max_subcycles,
      .max_rkl_stages = physics_config.metal_diffusion_max_rkl_stages,
      .diffusivity_floor_code = physics_config.metal_diffusion_coefficient_floor_code,
      .diffusivity_ceiling_code = physics_config.metal_diffusion_coefficient_ceiling_code,
  };
}

}  // namespace cosmosim::physics
