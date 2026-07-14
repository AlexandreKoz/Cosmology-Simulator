#include "cosmosim/gravity/pm_solver.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <numeric>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>
#include <limits>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/device_buffer.hpp"
#include "cosmosim/gravity/tree_pm_split_kernel.hpp"
#if COSMOSIM_ENABLE_CUDA
#include <cuda_runtime.h>
#include "cosmosim/gravity/pm_cuda_kernels.hpp"
#endif

#if COSMOSIM_ENABLE_FFTW
#include <fftw3.h>
#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#include <fftw3-mpi.h>
#endif
#endif

namespace cosmosim::gravity {
namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;

#if COSMOSIM_ENABLE_MPI
bool queryActiveMpiWorld(int& world_size, int& world_rank) noexcept {
  world_size = 1;
  world_rank = 0;
  int initialized = 0;
  MPI_Initialized(&initialized);
  if (initialized == 0) {
    return false;
  }
  int finalized = 0;
  MPI_Finalized(&finalized);
  if (finalized != 0) {
    return false;
  }
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  return true;
}
#endif

[[nodiscard]] std::size_t wrapIndex(std::ptrdiff_t i, std::size_t n) {
  const std::ptrdiff_t n_signed = static_cast<std::ptrdiff_t>(n);
  std::ptrdiff_t wrapped = i % n_signed;
  if (wrapped < 0) {
    wrapped += n_signed;
  }
  return static_cast<std::size_t>(wrapped);
}

[[nodiscard]] double wrapPosition(double x, double box_size) {
  const double wrapped = std::fmod(x, box_size);
  return wrapped < 0.0 ? wrapped + box_size : wrapped;
}

[[nodiscard]] bool positionInsideOpenDomain(double x, double box_size) {
  return x >= 0.0 && x < box_size;
}

[[nodiscard]] double sinc(double x) {
  if (std::abs(x) < 1.0e-12) {
    return 1.0;
  }
  return std::sin(x) / x;
}

struct PmAxisStencil1d {
  std::array<std::ptrdiff_t, 3> offsets{};
  std::array<double, 3> weights{};
  std::size_t count = 0;
};

struct PmDensityContributionRecord {
  std::uint32_t origin_rank = 0;
  std::uint32_t destination_rank = 0;
  std::uint32_t record_sequence = 0;
  std::uint32_t global_ix = 0;
  std::uint32_t global_iy = 0;
  std::uint32_t global_iz = 0;
  std::uint64_t exchange_epoch = 0;
  double mass_contribution = 0.0;
};

struct PmInterpolationRequestRecord {
  std::uint32_t origin_rank = 0;
  std::uint32_t destination_rank = 0;
  std::uint32_t request_sequence = 0;
  std::uint32_t origin_particle_index = 0;
  std::uint32_t global_ix = 0;
  std::uint32_t global_iy = 0;
  std::uint32_t global_iz = 0;
  std::uint64_t exchange_epoch = 0;
  double weight = 0.0;
};

struct PmForceContributionRecord {
  std::uint32_t source_rank = 0;
  std::uint32_t origin_rank = 0;
  std::uint32_t request_sequence = 0;
  std::uint32_t origin_particle_index = 0;
  std::uint64_t exchange_epoch = 0;
  double accel_x = 0.0;
  double accel_y = 0.0;
  double accel_z = 0.0;
};

struct PmPotentialContributionRecord {
  std::uint32_t source_rank = 0;
  std::uint32_t origin_rank = 0;
  std::uint32_t request_sequence = 0;
  std::uint32_t origin_particle_index = 0;
  std::uint64_t exchange_epoch = 0;
  double potential = 0.0;
};

struct PmInterpolationRequestRegistryEntry {
  std::uint32_t origin_particle_index = 0;
  std::uint64_t exchange_epoch = 0;
  int expected_sender_rank = -1;
};

[[nodiscard]] std::string pmRoutingDiagnostic(
    std::string_view stage,
    int receiver_rank,
    int sender_rank,
    std::uint64_t exchange_epoch,
    std::uint32_t request_sequence,
    std::uint32_t origin_particle_index,
    std::string_view detail) {
  return std::string(stage) + " routing failure local_receiver_rank=" + std::to_string(receiver_rank) +
      " sender_rank=" + std::to_string(sender_rank) + " exchange_epoch=" + std::to_string(exchange_epoch) +
      " request_sequence=" + std::to_string(request_sequence) +
      " origin_particle_index=" + std::to_string(origin_particle_index) + ": " + std::string(detail);
}

constexpr std::uint32_t k_pm_wire_magic = 0x31574d50U;  // "PMW1" in little-endian byte order.
constexpr std::uint32_t k_pm_wire_version = 1U;

enum class PmWireRecordKind : std::uint32_t {
  kDensityContribution = 1U,
  kForceRequest = 2U,
  kForceResponse = 3U,
  kPotentialRequest = 4U,
  kPotentialResponse = 5U,
};

constexpr std::size_t k_pm_density_wire_bytes = 56U;
constexpr std::size_t k_pm_interpolation_request_wire_bytes = 56U;
constexpr std::size_t k_pm_force_response_wire_bytes = 64U;
constexpr std::size_t k_pm_potential_response_wire_bytes = 48U;
static_assert(sizeof(double) == sizeof(std::uint64_t), "PM wire protocol requires binary64-sized doubles");
static_assert(
    std::numeric_limits<double>::is_iec559 && std::numeric_limits<double>::digits == 53,
    "PM wire protocol requires IEEE-754 binary64 doubles");

void writePmWireU32(std::span<std::uint8_t> bytes, std::size_t& offset, std::uint32_t value) {
  if (offset > bytes.size() || bytes.size() - offset < sizeof(value)) {
    throw std::overflow_error("PM wire encoder exceeded its fixed record size");
  }
  for (unsigned shift = 0U; shift < 32U; shift += 8U) {
    bytes[offset++] = static_cast<std::uint8_t>((value >> shift) & 0xffU);
  }
}

void writePmWireU64(std::span<std::uint8_t> bytes, std::size_t& offset, std::uint64_t value) {
  if (offset > bytes.size() || bytes.size() - offset < sizeof(value)) {
    throw std::overflow_error("PM wire encoder exceeded its fixed record size");
  }
  for (unsigned shift = 0U; shift < 64U; shift += 8U) {
    bytes[offset++] = static_cast<std::uint8_t>((value >> shift) & 0xffU);
  }
}

void writePmWireDouble(std::span<std::uint8_t> bytes, std::size_t& offset, double value) {
  writePmWireU64(bytes, offset, std::bit_cast<std::uint64_t>(value));
}

[[nodiscard]] std::uint32_t readPmWireU32(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset,
    std::string_view context) {
  if (offset > bytes.size() || bytes.size() - offset < sizeof(std::uint32_t)) {
    throw std::runtime_error(std::string(context) + " PM wire record is truncated");
  }
  std::uint32_t value = 0U;
  for (unsigned shift = 0U; shift < 32U; shift += 8U) {
    value |= static_cast<std::uint32_t>(bytes[offset++]) << shift;
  }
  return value;
}

[[nodiscard]] std::uint64_t readPmWireU64(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset,
    std::string_view context) {
  if (offset > bytes.size() || bytes.size() - offset < sizeof(std::uint64_t)) {
    throw std::runtime_error(std::string(context) + " PM wire record is truncated");
  }
  std::uint64_t value = 0U;
  for (unsigned shift = 0U; shift < 64U; shift += 8U) {
    value |= static_cast<std::uint64_t>(bytes[offset++]) << shift;
  }
  return value;
}

[[nodiscard]] double readPmWireDouble(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset,
    std::string_view context) {
  return std::bit_cast<double>(readPmWireU64(bytes, offset, context));
}

void writePmWireHeader(
    std::span<std::uint8_t> bytes,
    std::size_t& offset,
    PmWireRecordKind kind) {
  writePmWireU32(bytes, offset, k_pm_wire_magic);
  writePmWireU32(bytes, offset, k_pm_wire_version);
  writePmWireU32(bytes, offset, static_cast<std::uint32_t>(kind));
}

void readAndValidatePmWireHeader(
    std::span<const std::uint8_t> bytes,
    std::size_t& offset,
    PmWireRecordKind expected_kind,
    std::string_view context) {
  const std::uint32_t magic = readPmWireU32(bytes, offset, context);
  const std::uint32_t version = readPmWireU32(bytes, offset, context);
  const std::uint32_t kind = readPmWireU32(bytes, offset, context);
  if (magic != k_pm_wire_magic) {
    throw std::runtime_error(std::string(context) + " PM wire magic is invalid");
  }
  if (version != k_pm_wire_version) {
    throw std::runtime_error(
        std::string(context) + " PM wire version is unsupported: received=" + std::to_string(version) +
        " expected=" + std::to_string(k_pm_wire_version));
  }
  if (kind != static_cast<std::uint32_t>(expected_kind)) {
    throw std::runtime_error(
        std::string(context) + " PM wire record kind is invalid: received=" + std::to_string(kind) +
        " expected=" + std::to_string(static_cast<std::uint32_t>(expected_kind)));
  }
}

void encodePmDensityRecord(
    const PmDensityContributionRecord& record,
    std::span<std::uint8_t> bytes) {
  if (bytes.size() != k_pm_density_wire_bytes) {
    throw std::invalid_argument("PM density wire encoder requires an exact 56-byte record");
  }
  std::size_t offset = 0U;
  writePmWireHeader(bytes, offset, PmWireRecordKind::kDensityContribution);
  writePmWireU32(bytes, offset, record.origin_rank);
  writePmWireU32(bytes, offset, record.destination_rank);
  writePmWireU32(bytes, offset, record.record_sequence);
  writePmWireU32(bytes, offset, record.global_ix);
  writePmWireU32(bytes, offset, record.global_iy);
  writePmWireU32(bytes, offset, record.global_iz);
  writePmWireU32(bytes, offset, 0U);
  writePmWireU64(bytes, offset, record.exchange_epoch);
  writePmWireDouble(bytes, offset, record.mass_contribution);
  if (offset != bytes.size()) {
    throw std::logic_error("PM density wire encoder size contract drifted");
  }
}

[[nodiscard]] PmDensityContributionRecord decodePmDensityRecord(
    std::span<const std::uint8_t> bytes) {
  constexpr std::string_view context = "PM density contribution";
  if (bytes.size() != k_pm_density_wire_bytes) {
    throw std::runtime_error("PM density wire decoder requires an exact 56-byte record");
  }
  std::size_t offset = 0U;
  readAndValidatePmWireHeader(bytes, offset, PmWireRecordKind::kDensityContribution, context);
  PmDensityContributionRecord record;
  record.origin_rank = readPmWireU32(bytes, offset, context);
  record.destination_rank = readPmWireU32(bytes, offset, context);
  record.record_sequence = readPmWireU32(bytes, offset, context);
  record.global_ix = readPmWireU32(bytes, offset, context);
  record.global_iy = readPmWireU32(bytes, offset, context);
  record.global_iz = readPmWireU32(bytes, offset, context);
  const std::uint32_t reserved = readPmWireU32(bytes, offset, context);
  record.exchange_epoch = readPmWireU64(bytes, offset, context);
  record.mass_contribution = readPmWireDouble(bytes, offset, context);
  if (reserved != 0U || offset != bytes.size()) {
    throw std::runtime_error("PM density wire record has nonzero reserved flags or trailing data");
  }
  return record;
}

void encodePmInterpolationRequest(
    const PmInterpolationRequestRecord& record,
    PmWireRecordKind kind,
    std::span<std::uint8_t> bytes) {
  if ((kind != PmWireRecordKind::kForceRequest && kind != PmWireRecordKind::kPotentialRequest) ||
      bytes.size() != k_pm_interpolation_request_wire_bytes) {
    throw std::invalid_argument("PM interpolation request wire encoder received an invalid kind or record size");
  }
  std::size_t offset = 0U;
  writePmWireHeader(bytes, offset, kind);
  writePmWireU32(bytes, offset, record.origin_rank);
  writePmWireU32(bytes, offset, record.destination_rank);
  writePmWireU32(bytes, offset, record.request_sequence);
  writePmWireU32(bytes, offset, record.origin_particle_index);
  writePmWireU32(bytes, offset, record.global_ix);
  writePmWireU32(bytes, offset, record.global_iy);
  writePmWireU32(bytes, offset, record.global_iz);
  writePmWireU64(bytes, offset, record.exchange_epoch);
  writePmWireDouble(bytes, offset, record.weight);
  if (offset != bytes.size()) {
    throw std::logic_error("PM interpolation request wire encoder size contract drifted");
  }
}

[[nodiscard]] PmInterpolationRequestRecord decodePmInterpolationRequest(
    std::span<const std::uint8_t> bytes,
    PmWireRecordKind expected_kind,
    std::string_view context) {
  if ((expected_kind != PmWireRecordKind::kForceRequest &&
       expected_kind != PmWireRecordKind::kPotentialRequest) ||
      bytes.size() != k_pm_interpolation_request_wire_bytes) {
    throw std::runtime_error(std::string(context) + " PM request wire record has an invalid kind or size");
  }
  std::size_t offset = 0U;
  readAndValidatePmWireHeader(bytes, offset, expected_kind, context);
  PmInterpolationRequestRecord record;
  record.origin_rank = readPmWireU32(bytes, offset, context);
  record.destination_rank = readPmWireU32(bytes, offset, context);
  record.request_sequence = readPmWireU32(bytes, offset, context);
  record.origin_particle_index = readPmWireU32(bytes, offset, context);
  record.global_ix = readPmWireU32(bytes, offset, context);
  record.global_iy = readPmWireU32(bytes, offset, context);
  record.global_iz = readPmWireU32(bytes, offset, context);
  record.exchange_epoch = readPmWireU64(bytes, offset, context);
  record.weight = readPmWireDouble(bytes, offset, context);
  if (offset != bytes.size()) {
    throw std::runtime_error(std::string(context) + " PM request wire record has trailing data");
  }
  return record;
}

void encodePmForceResponse(
    const PmForceContributionRecord& record,
    std::span<std::uint8_t> bytes) {
  if (bytes.size() != k_pm_force_response_wire_bytes) {
    throw std::invalid_argument("PM force response wire encoder requires an exact 64-byte record");
  }
  std::size_t offset = 0U;
  writePmWireHeader(bytes, offset, PmWireRecordKind::kForceResponse);
  writePmWireU32(bytes, offset, record.source_rank);
  writePmWireU32(bytes, offset, record.origin_rank);
  writePmWireU32(bytes, offset, record.request_sequence);
  writePmWireU32(bytes, offset, record.origin_particle_index);
  writePmWireU32(bytes, offset, 0U);
  writePmWireU64(bytes, offset, record.exchange_epoch);
  writePmWireDouble(bytes, offset, record.accel_x);
  writePmWireDouble(bytes, offset, record.accel_y);
  writePmWireDouble(bytes, offset, record.accel_z);
  if (offset != bytes.size()) {
    throw std::logic_error("PM force response wire encoder size contract drifted");
  }
}

[[nodiscard]] PmForceContributionRecord decodePmForceResponse(
    std::span<const std::uint8_t> bytes) {
  constexpr std::string_view context = "PM force response";
  if (bytes.size() != k_pm_force_response_wire_bytes) {
    throw std::runtime_error("PM force response wire decoder requires an exact 64-byte record");
  }
  std::size_t offset = 0U;
  readAndValidatePmWireHeader(bytes, offset, PmWireRecordKind::kForceResponse, context);
  PmForceContributionRecord record;
  record.source_rank = readPmWireU32(bytes, offset, context);
  record.origin_rank = readPmWireU32(bytes, offset, context);
  record.request_sequence = readPmWireU32(bytes, offset, context);
  record.origin_particle_index = readPmWireU32(bytes, offset, context);
  const std::uint32_t reserved = readPmWireU32(bytes, offset, context);
  record.exchange_epoch = readPmWireU64(bytes, offset, context);
  record.accel_x = readPmWireDouble(bytes, offset, context);
  record.accel_y = readPmWireDouble(bytes, offset, context);
  record.accel_z = readPmWireDouble(bytes, offset, context);
  if (reserved != 0U || offset != bytes.size()) {
    throw std::runtime_error("PM force response wire record has nonzero reserved flags or trailing data");
  }
  return record;
}

void encodePmPotentialResponse(
    const PmPotentialContributionRecord& record,
    std::span<std::uint8_t> bytes) {
  if (bytes.size() != k_pm_potential_response_wire_bytes) {
    throw std::invalid_argument("PM potential response wire encoder requires an exact 48-byte record");
  }
  std::size_t offset = 0U;
  writePmWireHeader(bytes, offset, PmWireRecordKind::kPotentialResponse);
  writePmWireU32(bytes, offset, record.source_rank);
  writePmWireU32(bytes, offset, record.origin_rank);
  writePmWireU32(bytes, offset, record.request_sequence);
  writePmWireU32(bytes, offset, record.origin_particle_index);
  writePmWireU32(bytes, offset, 0U);
  writePmWireU64(bytes, offset, record.exchange_epoch);
  writePmWireDouble(bytes, offset, record.potential);
  if (offset != bytes.size()) {
    throw std::logic_error("PM potential response wire encoder size contract drifted");
  }
}

[[nodiscard]] PmPotentialContributionRecord decodePmPotentialResponse(
    std::span<const std::uint8_t> bytes) {
  constexpr std::string_view context = "PM potential response";
  if (bytes.size() != k_pm_potential_response_wire_bytes) {
    throw std::runtime_error("PM potential response wire decoder requires an exact 48-byte record");
  }
  std::size_t offset = 0U;
  readAndValidatePmWireHeader(bytes, offset, PmWireRecordKind::kPotentialResponse, context);
  PmPotentialContributionRecord record;
  record.source_rank = readPmWireU32(bytes, offset, context);
  record.origin_rank = readPmWireU32(bytes, offset, context);
  record.request_sequence = readPmWireU32(bytes, offset, context);
  record.origin_particle_index = readPmWireU32(bytes, offset, context);
  const std::uint32_t reserved = readPmWireU32(bytes, offset, context);
  record.exchange_epoch = readPmWireU64(bytes, offset, context);
  record.potential = readPmWireDouble(bytes, offset, context);
  if (reserved != 0U || offset != bytes.size()) {
    throw std::runtime_error("PM potential response wire record has nonzero reserved flags or trailing data");
  }
  return record;
}

[[nodiscard]] std::size_t checkedPmWireByteCount(
    std::size_t record_count,
    std::size_t record_bytes,
    std::string_view context) {
  if (record_bytes == 0U || record_count > std::numeric_limits<std::size_t>::max() / record_bytes) {
    throw std::overflow_error(std::string(context) + " PM wire byte count overflows size_t");
  }
  return record_count * record_bytes;
}

template <typename Record, typename Encoder>
void encodePmWireRecords(
    std::span<const Record> records,
    std::size_t record_bytes,
    std::vector<std::uint8_t>& wire,
    Encoder&& encoder,
    std::string_view context) {
  wire.resize(checkedPmWireByteCount(records.size(), record_bytes, context));
  for (std::size_t i = 0; i < records.size(); ++i) {
    encoder(records[i], std::span<std::uint8_t>(wire).subspan(i * record_bytes, record_bytes));
  }
}

template <typename Record, typename Decoder>
void decodePmWireRecords(
    std::span<const std::uint8_t> wire,
    std::size_t record_bytes,
    std::vector<Record>& records,
    Decoder&& decoder,
    std::string_view context) {
  if (record_bytes == 0U || wire.size() % record_bytes != 0U) {
    throw std::runtime_error(std::string(context) + " PM wire payload is not record-aligned");
  }
  const std::size_t record_count = wire.size() / record_bytes;
  records.resize(record_count);
  for (std::size_t i = 0; i < record_count; ++i) {
    records[i] = decoder(wire.subspan(i * record_bytes, record_bytes));
  }
}

#if COSMOSIM_ENABLE_MPI
[[nodiscard]] const std::uint8_t* nonNullPmWireData(const std::vector<std::uint8_t>& bytes) {
  static constexpr std::uint8_t empty_payload = 0U;
  return bytes.empty() ? &empty_payload : bytes.data();
}

[[nodiscard]] std::uint8_t* nonNullPmWireData(std::vector<std::uint8_t>& bytes) {
  static std::uint8_t empty_payload = 0U;
  return bytes.empty() ? &empty_payload : bytes.data();
}

void validatePmExchangeEpochConsensus(std::uint64_t exchange_epoch, std::string_view context) {
  std::uint64_t minimum_epoch = 0U;
  std::uint64_t maximum_epoch = 0U;
  MPI_Allreduce(&exchange_epoch, &minimum_epoch, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&exchange_epoch, &maximum_epoch, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  if (minimum_epoch != maximum_epoch) {
    throw std::runtime_error(
        std::string(context) + " ranks disagree on PM exchange epoch: minimum=" +
        std::to_string(minimum_epoch) + " maximum=" + std::to_string(maximum_epoch));
  }
}

template <typename Callable>
void measurePmMpiWait(double& accumulated_ms, Callable&& call) {
  const auto start = std::chrono::steady_clock::now();
  call();
  const auto stop = std::chrono::steady_clock::now();
  accumulated_ms += std::chrono::duration<double, std::milli>(stop - start).count();
}

enum class PmCollectiveEntryKind : int {
  kAssignDensity = 1,
  kSolvePoissonPeriodic = 2,
  kInterpolateForces = 3,
  kInterpolatePotential = 4,
  kSolvePoissonIsolatedOpen = 5,
  kInterpolateForceTargetView = 6,
  kSolveForParticles = 7,
};

void validatePmCollectiveEntryConsensus(
    PmCollectiveEntryKind entry_kind,
    bool rank_local_serial_layout,
    const PmGridShape& shape,
    const PmSolveOptions& options,
    double& accumulated_mpi_wait_ms,
    int control_lane_0 = 0,
    int control_lane_1 = 0) {
  // Public PM entry is collective over MPI_COMM_WORLD. Bind every mesh and
  // solver control that can select a different exchange, FFT plan, or kernel
  // before any rank enters a layout-specific phase.
  const std::array<std::uint64_t, 21> local_fingerprint{
      static_cast<std::uint64_t>(entry_kind),
      rank_local_serial_layout ? 1U : 0U,
      static_cast<std::uint64_t>(shape.nx),
      static_cast<std::uint64_t>(shape.ny),
      static_cast<std::uint64_t>(shape.nz),
      std::bit_cast<std::uint64_t>(options.box_size_mpc_comoving),
      std::bit_cast<std::uint64_t>(options.box_size_x_mpc_comoving),
      std::bit_cast<std::uint64_t>(options.box_size_y_mpc_comoving),
      std::bit_cast<std::uint64_t>(options.box_size_z_mpc_comoving),
      std::bit_cast<std::uint64_t>(options.scale_factor),
      std::bit_cast<std::uint64_t>(options.gravitational_constant_code),
      static_cast<std::uint64_t>(options.assignment_scheme),
      options.enable_window_deconvolution ? 1U : 0U,
      static_cast<std::uint64_t>(options.execution_policy),
      static_cast<std::uint64_t>(options.data_residency),
      static_cast<std::uint64_t>(options.decomposition_mode),
      static_cast<std::uint64_t>(options.boundary_condition),
      std::bit_cast<std::uint64_t>(options.tree_pm_split_scale_comoving),
      options.isolated_open_root_workspace_limit_bytes,
      static_cast<std::uint64_t>(control_lane_0),
      static_cast<std::uint64_t>(control_lane_1),
  };
  std::array<std::uint64_t, 21> minimum_fingerprint{};
  std::array<std::uint64_t, 21> maximum_fingerprint{};
  measurePmMpiWait(accumulated_mpi_wait_ms, [&]() {
    MPI_Allreduce(
        local_fingerprint.data(),
        minimum_fingerprint.data(),
        static_cast<int>(local_fingerprint.size()),
        MPI_UINT64_T,
        MPI_MIN,
        MPI_COMM_WORLD);
  });
  measurePmMpiWait(accumulated_mpi_wait_ms, [&]() {
    MPI_Allreduce(
        local_fingerprint.data(),
        maximum_fingerprint.data(),
        static_cast<int>(local_fingerprint.size()),
        MPI_UINT64_T,
        MPI_MAX,
        MPI_COMM_WORLD);
  });
  if (minimum_fingerprint != maximum_fingerprint) {
    throw std::runtime_error(
        "PM collective entry disagreement: api_kind_min=" + std::to_string(minimum_fingerprint[0]) +
        " api_kind_max=" + std::to_string(maximum_fingerprint[0]) +
        " rank_local_serial_min=" + std::to_string(minimum_fingerprint[1]) +
        " rank_local_serial_max=" + std::to_string(maximum_fingerprint[1]) +
        " mesh_min=(" + std::to_string(minimum_fingerprint[2]) + "," +
        std::to_string(minimum_fingerprint[3]) + "," +
        std::to_string(minimum_fingerprint[4]) + ") mesh_max=(" +
        std::to_string(maximum_fingerprint[2]) + "," +
        std::to_string(maximum_fingerprint[3]) + "," +
        std::to_string(maximum_fingerprint[4]) + ") control_0_min=" +
        std::to_string(minimum_fingerprint[19]) + " control_0_max=" +
        std::to_string(maximum_fingerprint[19]) + " control_1_min=" +
        std::to_string(minimum_fingerprint[20]) + " control_1_max=" +
        std::to_string(maximum_fingerprint[20]));
  }
}

void throwIfPmPayloadValidationFailed(
    std::string_view stage,
    int world_rank,
    int world_size,
    const std::string& local_error,
    double& accumulated_mpi_wait_ms) {
  const int local_failure = local_error.empty() ? 0 : 1;
  int global_failure_count = 0;
  measurePmMpiWait(accumulated_mpi_wait_ms, [&]() {
    MPI_Allreduce(&local_failure, &global_failure_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  });
  if (global_failure_count == 0) {
    return;
  }

  const int local_failure_rank = local_failure != 0 ? world_rank : world_size;
  int first_failure_rank = world_size;
  measurePmMpiWait(accumulated_mpi_wait_ms, [&]() {
    MPI_Allreduce(&local_failure_rank, &first_failure_rank, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  });
  const std::string detail = local_failure != 0
      ? " local_detail=" + local_error
      : " local_detail=peer_rank_rejected_payload";
  throw std::runtime_error(
      std::string(stage) + " coordinated PM payload validation failed: failure_count=" +
      std::to_string(global_failure_count) + " first_failure_rank=" + std::to_string(first_failure_rank) +
      detail);
}

template <typename Callable>
void runPmCoordinatedPhase(
    std::string_view stage,
    int world_rank,
    int world_size,
    double& accumulated_mpi_wait_ms,
    Callable&& phase) {
  std::string local_error;
  try {
    phase();
  } catch (const std::exception& error) {
    local_error = error.what();
  } catch (...) {
    local_error = "non-standard exception";
  }
  throwIfPmPayloadValidationFailed(
      stage,
      world_rank,
      world_size,
      local_error,
      accumulated_mpi_wait_ms);
}

[[nodiscard]] int checkedMpiRecordLimitForByteTransport(
    std::size_t record_size,
    std::string_view context) {
  const std::uint64_t max_mpi_int = static_cast<std::uint64_t>(std::numeric_limits<int>::max());
  if (record_size == 0U || record_size > max_mpi_int) {
    throw std::invalid_argument(
        std::string(context) + " record size=" + std::to_string(record_size) +
        " cannot be represented by MPI_BYTE int arguments");
  }
  return static_cast<int>(max_mpi_int / static_cast<std::uint64_t>(record_size));
}

[[nodiscard]] int checkedMpiRecordCountOrDisplacement(
    std::size_t record_value,
    std::size_t record_size,
    std::string_view context,
    std::string_view quantity,
    int rank) {
  const int record_limit = checkedMpiRecordLimitForByteTransport(record_size, context);
  if (record_value > static_cast<std::size_t>(record_limit)) {
    throw std::invalid_argument(
        std::string(context) + " " + std::string(quantity) + " exceeds MPI_BYTE record limit on rank " +
        std::to_string(rank) + ": records=" + std::to_string(record_value) +
        " record_size=" + std::to_string(record_size) +
        " max_records=" + std::to_string(record_limit));
  }
  return static_cast<int>(record_value);
}

[[nodiscard]] std::size_t checkedMpiReceivedRecordCount(
    int record_count,
    std::size_t record_size,
    std::string_view context,
    int rank) {
  if (record_count < 0) {
    throw std::invalid_argument(
        std::string(context) + " received negative record count from rank " + std::to_string(rank) +
        ": count=" + std::to_string(record_count));
  }
  const std::size_t count = static_cast<std::size_t>(record_count);
  static_cast<void>(checkedMpiRecordCountOrDisplacement(
      count, record_size, context, "received record count", rank));
  return count;
}

[[nodiscard]] std::size_t checkedMpiRecordTotal(
    std::size_t current_total,
    std::size_t record_count,
    std::size_t record_size,
    std::string_view context,
    int rank) {
  if (record_count > std::numeric_limits<std::size_t>::max() - current_total) {
    throw std::invalid_argument(
        std::string(context) + " record total overflows size_t while processing rank " + std::to_string(rank));
  }
  const std::size_t next_total = current_total + record_count;
  static_cast<void>(checkedMpiRecordCountOrDisplacement(
      next_total, record_size, context, "cumulative record total", rank));
  return next_total;
}

void checkedMpiRecordLayoutToByteLayout(
    std::span<const int> record_counts,
    std::span<const int> record_displacements,
    std::size_t record_size,
    std::span<int> byte_counts,
    std::span<int> byte_displacements,
    std::string_view context) {
  if (record_counts.size() != record_displacements.size() ||
      record_counts.size() != byte_counts.size() ||
      record_counts.size() != byte_displacements.size()) {
    throw std::invalid_argument(
        std::string(context) + " MPI record/byte count-displacement span sizes do not match");
  }

  const int record_limit = checkedMpiRecordLimitForByteTransport(record_size, context);
  const std::uint64_t max_mpi_int = static_cast<std::uint64_t>(std::numeric_limits<int>::max());
  const std::uint64_t record_size_u64 = static_cast<std::uint64_t>(record_size);
  for (std::size_t rank = 0; rank < record_counts.size(); ++rank) {
    const int record_count = record_counts[rank];
    const int record_displacement = record_displacements[rank];
    if (record_count < 0 || record_displacement < 0) {
      throw std::invalid_argument(
          std::string(context) + " contains negative record count/displacement for rank " +
          std::to_string(rank) + ": count=" + std::to_string(record_count) +
          " displacement=" + std::to_string(record_displacement));
    }

    const std::uint64_t count_u64 = static_cast<std::uint64_t>(record_count);
    const std::uint64_t displacement_u64 = static_cast<std::uint64_t>(record_displacement);
    if (count_u64 > static_cast<std::uint64_t>(record_limit) ||
        displacement_u64 > static_cast<std::uint64_t>(record_limit)) {
      throw std::invalid_argument(
          std::string(context) + " record count/displacement exceeds MPI_BYTE record limit for rank " +
          std::to_string(rank) + ": count=" + std::to_string(record_count) +
          " displacement=" + std::to_string(record_displacement) +
          " record_size=" + std::to_string(record_size) +
          " max_records=" + std::to_string(record_limit));
    }

    const std::uint64_t byte_count = count_u64 * record_size_u64;
    const std::uint64_t byte_displacement = displacement_u64 * record_size_u64;
    if (byte_count > max_mpi_int || byte_displacement > max_mpi_int ||
        byte_count > max_mpi_int - byte_displacement) {
      throw std::invalid_argument(
          std::string(context) + " byte count/displacement exceeds MPI int range for rank " +
          std::to_string(rank) + ": byte_count=" + std::to_string(byte_count) +
          " byte_displacement=" + std::to_string(byte_displacement));
    }
    byte_counts[rank] = static_cast<int>(byte_count);
    byte_displacements[rank] = static_cast<int>(byte_displacement);
  }
}
#endif

[[nodiscard]] PmAxisStencil1d makeAxisStencil(double grid_position, PmAssignmentScheme scheme) {
  PmAxisStencil1d stencil{};
  if (scheme == PmAssignmentScheme::kCic) {
    const double base = std::floor(grid_position);
    const std::ptrdiff_t i0 = static_cast<std::ptrdiff_t>(base);
    const double t = grid_position - base;
    stencil.offsets = {i0, i0 + 1, 0};
    stencil.weights = {1.0 - t, t, 0.0};
    stencil.count = 2;
    return stencil;
  }

  const double center = std::floor(grid_position + 0.5);
  const std::ptrdiff_t ic = static_cast<std::ptrdiff_t>(center);
  const double delta = grid_position - center;
  const double w_m1 = 0.5 * std::pow(0.5 - delta, 2.0);
  const double w_0 = 0.75 - delta * delta;
  const double w_p1 = 0.5 * std::pow(0.5 + delta, 2.0);
  stencil.offsets = {ic - 1, ic, ic + 1};
  stencil.weights = {w_m1, w_0, w_p1};
  stencil.count = 3;
  return stencil;
}

[[nodiscard]] int assignmentWindowExponent(PmAssignmentScheme scheme) {
  switch (scheme) {
    case PmAssignmentScheme::kCic:
      return 2;
    case PmAssignmentScheme::kTsc:
      return 3;
  }
  throw std::invalid_argument("Unknown PM assignment scheme in assignmentWindowExponent");
}

[[nodiscard]] std::uint64_t bytesForGridSweep(std::size_t cell_count) {
  return static_cast<std::uint64_t>(cell_count * sizeof(double));
}

[[nodiscard]] std::uint64_t bytesForParticles(std::size_t particle_count) {
  return static_cast<std::uint64_t>(particle_count * sizeof(double) * 4U);
}

[[nodiscard]] std::size_t checkedProduct(std::size_t a, std::size_t b, std::string_view context) {
  if (a != 0U && b > std::numeric_limits<std::size_t>::max() / a) {
    throw std::overflow_error(std::string(context) + " size product overflows size_t");
  }
  return a * b;
}

[[nodiscard]] std::uint64_t checkedBytesForCount(
    std::size_t count,
    std::size_t element_size,
    std::string_view context) {
  if (element_size != 0U && count > std::numeric_limits<std::uint64_t>::max() / element_size) {
    throw std::overflow_error(std::string(context) + " byte estimate overflows uint64_t");
  }
  return static_cast<std::uint64_t>(count) * static_cast<std::uint64_t>(element_size);
}

[[nodiscard]] std::uint64_t checkedAddBytes(
    std::uint64_t a,
    std::uint64_t b,
    std::string_view context) {
  if (b > std::numeric_limits<std::uint64_t>::max() - a) {
    throw std::overflow_error(std::string(context) + " byte estimate overflows uint64_t");
  }
  return a + b;
}

[[nodiscard]] std::uint64_t estimateIsolatedRootWorkspaceBytes(
    const PmGridShape& shape,
    bool distributed_slabs) {
  const std::size_t physical_cells =
      checkedProduct(checkedProduct(shape.nx, shape.ny, "isolated PM physical grid"),
                     shape.nz,
                     "isolated PM physical grid");
  const std::size_t pad_nx = checkedProduct(2U, shape.nx, "isolated PM padded nx");
  const std::size_t pad_ny = checkedProduct(2U, shape.ny, "isolated PM padded ny");
  const std::size_t pad_nz = checkedProduct(2U, shape.nz, "isolated PM padded nz");
  const std::size_t padded_cells =
      checkedProduct(checkedProduct(pad_nx, pad_ny, "isolated PM padded grid"),
                     pad_nz,
                     "isolated PM padded grid");
  std::uint64_t bytes = 0;
  if (distributed_slabs) {
    bytes = checkedAddBytes(
        bytes,
        checkedBytesForCount(physical_cells, sizeof(double), "isolated PM root density gather"),
        "isolated PM root workspace");
  }
  bytes = checkedAddBytes(
      bytes,
      checkedBytesForCount(padded_cells, sizeof(std::complex<double>) * 3U, "isolated PM padded complex workspace"),
      "isolated PM root workspace");
  bytes = checkedAddBytes(
      bytes,
      checkedBytesForCount(padded_cells, sizeof(double), "isolated PM padded real workspace"),
      "isolated PM root workspace");
  bytes = checkedAddBytes(
      bytes,
      checkedBytesForCount(physical_cells, sizeof(double) * 4U, "isolated PM extracted fields"),
      "isolated PM root workspace");
  return bytes;
}

[[nodiscard]] std::string isolatedPmGuardMessage(
    const PmGridShape& shape,
    int world_size,
    int world_rank,
    std::uint64_t estimated_bytes,
    std::uint64_t limit_bytes) {
  return "isolated/open PM root-gather guard exceeded on rank " + std::to_string(world_rank) +
      ": grid_shape=(" + std::to_string(shape.nx) + "," + std::to_string(shape.ny) + "," +
      std::to_string(shape.nz) + ") ranks=" + std::to_string(world_size) +
      " estimated_root_bytes=" + std::to_string(estimated_bytes) +
      " configured_limit_bytes=" + std::to_string(limit_bytes) +
      " route=root_gather_scatter policy=bounded_small_grid_only";
}

struct BoxLengths {
  double lx = 0.0;
  double ly = 0.0;
  double lz = 0.0;
};

[[nodiscard]] BoxLengths effectiveBoxLengths(const PmSolveOptions& options) {
  const double scalar = options.box_size_mpc_comoving;
  const double lx = options.box_size_x_mpc_comoving > 0.0 ? options.box_size_x_mpc_comoving : scalar;
  const double ly = options.box_size_y_mpc_comoving > 0.0 ? options.box_size_y_mpc_comoving : scalar;
  const double lz = options.box_size_z_mpc_comoving > 0.0 ? options.box_size_z_mpc_comoving : scalar;
  return BoxLengths{.lx = lx, .ly = ly, .lz = lz};
}

void validateOptions(const PmGridShape& shape, const PmSolveOptions& options) {
  if (!shape.isValid()) {
    throw std::invalid_argument("PM grid shape must be non-zero in all dimensions");
  }
  if (!std::isfinite(options.box_size_mpc_comoving) ||
      !std::isfinite(options.box_size_x_mpc_comoving) ||
      !std::isfinite(options.box_size_y_mpc_comoving) ||
      !std::isfinite(options.box_size_z_mpc_comoving)) {
    throw std::invalid_argument("PM solve requires finite scalar and axis-aware box lengths");
  }
  const BoxLengths lengths = effectiveBoxLengths(options);
  if (!std::isfinite(lengths.lx) || !std::isfinite(lengths.ly) || !std::isfinite(lengths.lz) ||
      lengths.lx <= 0.0 || lengths.ly <= 0.0 || lengths.lz <= 0.0) {
    throw std::invalid_argument("PM solve requires finite positive effective box lengths on all axes");
  }
  if (!std::isfinite(options.scale_factor) || options.scale_factor <= 0.0) {
    throw std::invalid_argument("PM solve requires finite scale_factor > 0");
  }
  if (!std::isfinite(options.gravitational_constant_code) || options.gravitational_constant_code <= 0.0) {
    throw std::invalid_argument("PM solve requires finite gravitational_constant_code > 0");
  }
  if (!std::isfinite(options.tree_pm_split_scale_comoving)) {
    throw std::invalid_argument("PM solve requires a finite TreePM split scale");
  }
  if (options.execution_policy == core::ExecutionPolicy::kCuda && options.data_residency == PmDataResidencyPolicy::kHostOnly) {
    throw std::invalid_argument(
        "execution_policy=cuda requires data_residency=kPreferDevice for explicit host/device ownership");
  }
  if (options.execution_policy == core::ExecutionPolicy::kCuda && options.assignment_scheme != PmAssignmentScheme::kCic) {
    throw std::invalid_argument(
        "execution_policy=cuda currently supports only assignment_scheme=cic in this build");
  }
  if (options.execution_policy == core::ExecutionPolicy::kCuda &&
      options.boundary_condition != PmBoundaryCondition::kPeriodic) {
    throw std::invalid_argument(
        "execution_policy=cuda currently supports only boundary_condition=periodic in this build");
  }
  if (options.boundary_condition == PmBoundaryCondition::kIsolatedOpen &&
      options.enable_window_deconvolution) {
    throw std::invalid_argument("isolated PM currently does not support window deconvolution");
  }
}

void validateSingleRankFullDomainGridContract(const PmGridStorage& grid, std::string_view callsite) {
  if (!grid.ownsFullDomain()) {
    throw std::invalid_argument(
        std::string(callsite) +
        " requires full-domain PM slab ownership in single-rank mode; use a valid distributed slab layout for multi-rank PM execution");
  }
}

}  // namespace

class PmSolver::Impl {
 public:
  struct PlanKey {
    int world_size = 1;
    int world_rank = 0;
    std::size_t owned_begin_x = 0;
    std::size_t owned_end_x = 0;
    core::PmDecompositionMode decomposition_mode = core::PmDecompositionMode::kSlab;

    [[nodiscard]] bool operator==(const PlanKey& other) const noexcept {
      return world_size == other.world_size && world_rank == other.world_rank && owned_begin_x == other.owned_begin_x &&
          owned_end_x == other.owned_end_x && decomposition_mode == other.decomposition_mode;
    }
  };

  struct PlanKeyHasher {
    [[nodiscard]] std::size_t operator()(const PlanKey& key) const noexcept {
      std::size_t seed = static_cast<std::size_t>(key.world_size * 1315423911U + key.world_rank);
      seed ^= key.owned_begin_x + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
      seed ^= key.owned_end_x + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
      seed ^= static_cast<std::size_t>(key.decomposition_mode) + 0x9e3779b97f4a7c15ULL + (seed << 6U) + (seed >> 2U);
      return seed;
    }
  };

  struct PlanResources {
    parallel::PmSlabLayout layout{};
    std::vector<double> real;
    std::vector<std::complex<double>> fourier;
    std::vector<std::complex<double>> potential_k;
    std::vector<std::complex<double>> working_k;
    std::vector<double> poisson_kernel;
    std::vector<double> grad_kx;
    std::vector<double> grad_ky;
    std::vector<double> grad_kz;
#if COSMOSIM_ENABLE_FFTW
    fftw_plan forward_plan = nullptr;
    fftw_plan inverse_plan = nullptr;
#endif
    bool is_distributed = false;
    bool spectral_transposed = false;
    std::size_t real_z_stride = 0;
    std::size_t transposed_local_ny = 0;
    std::size_t transposed_begin_y = 0;
    bool spectral_operators_ready = false;
    double cached_lx = 0.0;
    double cached_ly = 0.0;
    double cached_lz = 0.0;
    double cached_split_scale = -1.0;
    double cached_scale_factor = -1.0;
    double cached_gravitational_constant_code = -1.0;
    bool cached_window_deconvolution = false;
    PmAssignmentScheme cached_assignment_scheme = PmAssignmentScheme::kCic;
  };

  struct DensityExchangeBuffers {
    std::vector<std::vector<PmDensityContributionRecord>> send_records_by_rank;
    std::vector<PmDensityContributionRecord> send_flat_records;
    std::vector<PmDensityContributionRecord> recv_flat_records;
    std::vector<std::uint8_t> send_wire;
    std::vector<std::uint8_t> recv_wire;
    std::vector<int> send_counts;
    std::vector<int> send_displs;
    std::vector<int> recv_counts;
    std::vector<int> recv_displs;
    std::vector<int> send_counts_bytes;
    std::vector<int> send_displs_bytes;
    std::vector<int> recv_counts_bytes;
    std::vector<int> recv_displs_bytes;
  };

  template <typename ContributionRecord>
  struct InterpolationExchangeBuffers {
    std::vector<std::vector<PmInterpolationRequestRecord>> send_requests_by_rank;
    std::vector<PmInterpolationRequestRecord> send_requests_flat;
    std::vector<PmInterpolationRequestRecord> recv_requests_flat;
    std::vector<std::uint8_t> send_request_wire;
    std::vector<std::uint8_t> recv_request_wire;
    std::vector<std::vector<ContributionRecord>> send_contribs_by_rank;
    std::vector<ContributionRecord> send_contribs_flat;
    std::vector<ContributionRecord> recv_contribs_flat;
    std::vector<std::uint8_t> send_contrib_wire;
    std::vector<std::uint8_t> recv_contrib_wire;
    std::vector<int> send_counts;
    std::vector<int> send_displs;
    std::vector<int> recv_counts;
    std::vector<int> recv_displs;
    std::vector<int> send_counts_bytes;
    std::vector<int> send_displs_bytes;
    std::vector<int> recv_counts_bytes;
    std::vector<int> recv_displs_bytes;
    std::vector<int> send_contrib_counts;
    std::vector<int> send_contrib_displs;
    std::vector<int> recv_contrib_counts;
    std::vector<int> recv_contrib_displs;
    std::vector<int> send_contrib_counts_bytes;
    std::vector<int> send_contrib_displs_bytes;
    std::vector<int> recv_contrib_counts_bytes;
    std::vector<int> recv_contrib_displs_bytes;
  };

  explicit Impl(PmGridShape shape) : m_shape(shape) {}

  ~Impl() {
#if COSMOSIM_ENABLE_FFTW
    for (auto& [_, plan] : m_plan_cache) {
      if (plan.forward_plan != nullptr) {
        fftw_destroy_plan(plan.forward_plan);
      }
      if (plan.inverse_plan != nullptr) {
        fftw_destroy_plan(plan.inverse_plan);
      }
    }
    if (m_isolated_workspace.forward_plan != nullptr) {
      fftw_destroy_plan(m_isolated_workspace.forward_plan);
    }
    if (m_isolated_workspace.inverse_plan != nullptr) {
      fftw_destroy_plan(m_isolated_workspace.inverse_plan);
    }
#endif
  }

  struct IsolatedOpenWorkspace {
    std::size_t nx = 0;
    std::size_t ny = 0;
    std::size_t nz = 0;
    std::vector<std::complex<double>> rho_k;
    std::vector<std::complex<double>> kernel_k;
    std::vector<std::complex<double>> scratch;
    std::vector<double> potential_real;
    bool kernel_ready = false;
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double split_scale = -1.0;
#if COSMOSIM_ENABLE_FFTW
    fftw_plan forward_plan = nullptr;
    fftw_plan inverse_plan = nullptr;
#endif
  };

  void ensureIsolatedWorkspace(std::size_t nx, std::size_t ny, std::size_t nz) {
    const bool shape_changed = nx != m_isolated_workspace.nx || ny != m_isolated_workspace.ny || nz != m_isolated_workspace.nz;
    if (!shape_changed) {
      return;
    }
#if COSMOSIM_ENABLE_FFTW
    if (m_isolated_workspace.forward_plan != nullptr) {
      fftw_destroy_plan(m_isolated_workspace.forward_plan);
      m_isolated_workspace.forward_plan = nullptr;
    }
    if (m_isolated_workspace.inverse_plan != nullptr) {
      fftw_destroy_plan(m_isolated_workspace.inverse_plan);
      m_isolated_workspace.inverse_plan = nullptr;
    }
#endif
    m_isolated_workspace.nx = nx;
    m_isolated_workspace.ny = ny;
    m_isolated_workspace.nz = nz;
    const std::size_t n_tot = nx * ny * nz;
    m_isolated_workspace.rho_k.assign(n_tot, {0.0, 0.0});
    m_isolated_workspace.kernel_k.assign(n_tot, {0.0, 0.0});
    m_isolated_workspace.scratch.assign(n_tot, {0.0, 0.0});
    m_isolated_workspace.potential_real.assign(n_tot, 0.0);
    m_isolated_workspace.kernel_ready = false;
    m_isolated_workspace.dx = 0.0;
    m_isolated_workspace.dy = 0.0;
    m_isolated_workspace.dz = 0.0;
    m_isolated_workspace.split_scale = -1.0;
#if COSMOSIM_ENABLE_FFTW
    m_isolated_workspace.forward_plan = fftw_plan_dft_3d(
        static_cast<int>(nx),
        static_cast<int>(ny),
        static_cast<int>(nz),
        reinterpret_cast<fftw_complex*>(m_isolated_workspace.scratch.data()),
        reinterpret_cast<fftw_complex*>(m_isolated_workspace.scratch.data()),
        FFTW_FORWARD,
        FFTW_ESTIMATE);
    m_isolated_workspace.inverse_plan = fftw_plan_dft_3d(
        static_cast<int>(nx),
        static_cast<int>(ny),
        static_cast<int>(nz),
        reinterpret_cast<fftw_complex*>(m_isolated_workspace.scratch.data()),
        reinterpret_cast<fftw_complex*>(m_isolated_workspace.scratch.data()),
        FFTW_BACKWARD,
        FFTW_ESTIMATE);
    if (m_isolated_workspace.forward_plan == nullptr || m_isolated_workspace.inverse_plan == nullptr) {
      throw std::runtime_error("failed to create isolated PM FFTW plans");
    }
#endif
  }

  void isolatedForward(std::span<std::complex<double>> field) {
    std::copy(field.begin(), field.end(), m_isolated_workspace.scratch.begin());
#if COSMOSIM_ENABLE_FFTW
    fftw_execute(m_isolated_workspace.forward_plan);
    std::copy(m_isolated_workspace.scratch.begin(), m_isolated_workspace.scratch.end(), field.begin());
#else
    naiveComplexDft(field, /*forward=*/true);
#endif
  }

  void isolatedInverse(std::span<std::complex<double>> field) {
    std::copy(field.begin(), field.end(), m_isolated_workspace.scratch.begin());
#if COSMOSIM_ENABLE_FFTW
    fftw_execute(m_isolated_workspace.inverse_plan);
    std::copy(m_isolated_workspace.scratch.begin(), m_isolated_workspace.scratch.end(), field.begin());
#else
    naiveComplexDft(field, /*forward=*/false);
#endif
  }

  [[nodiscard]] IsolatedOpenWorkspace& isolatedWorkspace() { return m_isolated_workspace; }

  [[nodiscard]] PlanResources& planForLayout(const parallel::PmSlabLayout& layout, core::PmDecompositionMode decomposition_mode) {
    const PlanKey key{
        .world_size = layout.world_size,
        .world_rank = layout.world_rank,
        .owned_begin_x = layout.owned_x.begin_x,
        .owned_end_x = layout.owned_x.end_x,
        .decomposition_mode = decomposition_mode,
    };
    auto it = m_plan_cache.find(key);
#if COSMOSIM_ENABLE_FFTW && COSMOSIM_ENABLE_MPI
    if (layout.world_size > 1) {
      const int local_cache_hit = it != m_plan_cache.end() ? 1 : 0;
      int global_cache_hit_count = 0;
      MPI_Allreduce(
          &local_cache_hit,
          &global_cache_hit_count,
          1,
          MPI_INT,
          MPI_SUM,
          MPI_COMM_WORLD);
      if (global_cache_hit_count != 0 && global_cache_hit_count != layout.world_size) {
        throw std::runtime_error(
            "Distributed PM FFT plan cache state diverged across ranks; all ranks must build or reuse together");
      }
      if (global_cache_hit_count == layout.world_size) {
        m_active_key = key;
        return it->second;
      }
    } else
#endif
    if (it != m_plan_cache.end()) {
      m_active_key = key;
      return it->second;
    }

    PlanResources plan{};
    plan.layout = layout;
    const std::size_t nz_complex = m_shape.nz / 2U + 1U;
    plan.real_z_stride = m_shape.nz;
    std::size_t expected_local_complex_size = 0U;
#if COSMOSIM_ENABLE_FFTW && COSMOSIM_ENABLE_MPI
    double plan_mpi_wait_ms = 0.0;
    if (layout.world_size > 1) {
      runPmCoordinatedPhase(
          "PmSolver distributed FFT logical-size preparation",
          layout.world_rank,
          layout.world_size,
          plan_mpi_wait_ms,
          [&]() {
            expected_local_complex_size = checkedProduct(
                checkedProduct(layout.local_nx(), m_shape.ny, "PM local Fourier extent"),
                nz_complex,
                "PM local Fourier extent");
          });
    } else
#endif
    {
      expected_local_complex_size = checkedProduct(
          checkedProduct(layout.local_nx(), m_shape.ny, "PM local Fourier extent"),
          nz_complex,
          "PM local Fourier extent");
    }
    std::size_t allocated_local_complex_size = expected_local_complex_size;

#if COSMOSIM_ENABLE_FFTW
#if COSMOSIM_ENABLE_MPI
    if (layout.world_size > 1) {
      int mpi_world_size = 1;
      int mpi_world_rank = 0;
      queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
      runPmCoordinatedPhase(
          "PmSolver distributed FFT plan metadata preflight",
          mpi_world_rank,
          mpi_world_size,
          plan_mpi_wait_ms,
          [&]() {
            if (mpi_world_size != layout.world_size || mpi_world_rank != layout.world_rank) {
              throw std::invalid_argument(
                  "PM slab layout world metadata must match MPI_COMM_WORLD for distributed FFT path");
            }
          });
      fftw_mpi_init();
      ptrdiff_t backend_local_nx = 0;
      ptrdiff_t backend_begin_x = 0;
      ptrdiff_t backend_alloc_local = 0;
      ptrdiff_t backend_local_ny = 0;
      ptrdiff_t backend_begin_y = 0;
      if (decomposition_mode == core::PmDecompositionMode::kPencil) {
        backend_alloc_local = fftw_mpi_local_size_3d_transposed(
            static_cast<ptrdiff_t>(m_shape.nx),
            static_cast<ptrdiff_t>(m_shape.ny),
            static_cast<ptrdiff_t>(nz_complex),
            MPI_COMM_WORLD,
            &backend_local_nx,
            &backend_begin_x,
            &backend_local_ny,
            &backend_begin_y);
      } else {
        backend_alloc_local = fftw_mpi_local_size_3d(
            static_cast<ptrdiff_t>(m_shape.nx),
            static_cast<ptrdiff_t>(m_shape.ny),
            static_cast<ptrdiff_t>(nz_complex),
            MPI_COMM_WORLD,
            &backend_local_nx,
            &backend_begin_x);
      }
      runPmCoordinatedPhase(
          "PmSolver distributed FFT backend/storage preparation",
          mpi_world_rank,
          mpi_world_size,
          plan_mpi_wait_ms,
          [&]() {
      const bool local_extent_mismatch =
          backend_local_nx != static_cast<ptrdiff_t>(layout.local_nx());
      // FFTW leaves local_0_start unspecified for a zero-width input slab and
      // commonly reports zero, whereas CHUI's canonical empty range is
      // [global_nx, global_nx).  The origin has no ownership meaning when the
      // extent is zero, so compare it only for non-empty slabs.
      const bool nonempty_origin_mismatch =
          backend_local_nx > 0 &&
          backend_begin_x != static_cast<ptrdiff_t>(layout.owned_x.begin_x);
      if (local_extent_mismatch || nonempty_origin_mismatch) {
        std::ostringstream message;
        message << "PM slab layout is incompatible with FFTW MPI ownership for this communicator: rank="
                << layout.world_rank << ", configured=[" << layout.owned_x.begin_x << ',' << layout.owned_x.end_x
                << "), backend=[" << backend_begin_x << ',' << (backend_begin_x + backend_local_nx) << ')';
        throw std::invalid_argument(message.str());
      }
      if (backend_alloc_local < 0) {
        throw std::runtime_error("FFTW MPI reported a negative local allocation size for distributed PM plan");
      }
      // FFTW permits a participating rank to own no real-space slab when
      // world_size > nx and may then report alloc_local == 0.  Keep the
      // authoritative logical extents at zero while providing non-null dummy
      // storage to planners/backends that still require valid pointer values.
      allocated_local_complex_size = std::max<std::size_t>(
          1U,
          static_cast<std::size_t>(backend_alloc_local));
      plan.is_distributed = true;
      plan.real_z_stride = 2U * nz_complex;
      plan.real.assign(
          checkedProduct(2U, allocated_local_complex_size, "FFTW MPI local real allocation"),
          0.0);
      if (decomposition_mode == core::PmDecompositionMode::kPencil) {
        plan.spectral_transposed = true;
        plan.transposed_local_ny = static_cast<std::size_t>(backend_local_ny);
        plan.transposed_begin_y = static_cast<std::size_t>(backend_begin_y);
      }
      if (allocated_local_complex_size < expected_local_complex_size) {
        throw std::runtime_error("FFTW MPI local allocation is smaller than the expected slab-local Fourier extent");
      }
      plan.fourier.assign(allocated_local_complex_size, std::complex<double>(0.0, 0.0));
      plan.potential_k.assign(allocated_local_complex_size, std::complex<double>(0.0, 0.0));
      plan.working_k.assign(allocated_local_complex_size, std::complex<double>(0.0, 0.0));
      plan.poisson_kernel.assign(allocated_local_complex_size, 0.0);
      plan.grad_kx.assign(allocated_local_complex_size, 0.0);
      plan.grad_ky.assign(allocated_local_complex_size, 0.0);
      plan.grad_kz.assign(allocated_local_complex_size, 0.0);
          });
    }
#endif
#if COSMOSIM_ENABLE_MPI
    if (layout.world_size > 1) {
      plan.forward_plan = fftw_mpi_plan_dft_r2c_3d(
          static_cast<ptrdiff_t>(m_shape.nx),
          static_cast<ptrdiff_t>(m_shape.ny),
          static_cast<ptrdiff_t>(m_shape.nz),
          plan.real.data(),
          reinterpret_cast<fftw_complex*>(plan.fourier.data()),
          MPI_COMM_WORLD,
          decomposition_mode == core::PmDecompositionMode::kPencil
              ? (COSMOSIM_FFTW_PLANNER_FLAGS | FFTW_MPI_TRANSPOSED_OUT)
              : COSMOSIM_FFTW_PLANNER_FLAGS);
      plan.inverse_plan = fftw_mpi_plan_dft_c2r_3d(
          static_cast<ptrdiff_t>(m_shape.nx),
          static_cast<ptrdiff_t>(m_shape.ny),
          static_cast<ptrdiff_t>(m_shape.nz),
          reinterpret_cast<fftw_complex*>(plan.fourier.data()),
          plan.real.data(),
          MPI_COMM_WORLD,
          decomposition_mode == core::PmDecompositionMode::kPencil
              ? (COSMOSIM_FFTW_PLANNER_FLAGS | FFTW_MPI_TRANSPOSED_IN)
              : COSMOSIM_FFTW_PLANNER_FLAGS);
      runPmCoordinatedPhase(
          "PmSolver distributed FFT plan creation",
          layout.world_rank,
          layout.world_size,
          plan_mpi_wait_ms,
          [&]() {
            if (plan.forward_plan == nullptr || plan.inverse_plan == nullptr) {
              throw std::runtime_error("Failed to create distributed FFTW plans for PM solver");
            }
          });
    } else
#endif
    {
      plan.real_z_stride = m_shape.nz;
      plan.real.assign(
          checkedProduct(
              checkedProduct(layout.local_nx(), m_shape.ny, "PM local real extent"),
              plan.real_z_stride,
              "PM local real extent"),
          0.0);
      plan.fourier.assign(expected_local_complex_size, std::complex<double>(0.0, 0.0));
      plan.potential_k.assign(expected_local_complex_size, std::complex<double>(0.0, 0.0));
      plan.working_k.assign(expected_local_complex_size, std::complex<double>(0.0, 0.0));
      plan.poisson_kernel.assign(expected_local_complex_size, 0.0);
      plan.grad_kx.assign(expected_local_complex_size, 0.0);
      plan.grad_ky.assign(expected_local_complex_size, 0.0);
      plan.grad_kz.assign(expected_local_complex_size, 0.0);
      plan.forward_plan = fftw_plan_dft_r2c_3d(
          static_cast<int>(m_shape.nx),
          static_cast<int>(m_shape.ny),
          static_cast<int>(m_shape.nz),
          plan.real.data(),
          reinterpret_cast<fftw_complex*>(plan.fourier.data()),
          COSMOSIM_FFTW_PLANNER_FLAGS);
      plan.inverse_plan = fftw_plan_dft_c2r_3d(
          static_cast<int>(m_shape.nx),
          static_cast<int>(m_shape.ny),
          static_cast<int>(m_shape.nz),
          reinterpret_cast<fftw_complex*>(plan.fourier.data()),
          plan.real.data(),
          COSMOSIM_FFTW_PLANNER_FLAGS);
    }
    if (plan.forward_plan == nullptr || plan.inverse_plan == nullptr) {
      throw std::runtime_error("Failed to create FFTW plans for PM solver");
    }
#else
    if (!layout.ownsFullDomain()) {
      throw std::invalid_argument(
          "PM solver naive DFT fallback requires full-domain slab ownership; distributed PM requires COSMOSIM_ENABLE_FFTW=ON and COSMOSIM_ENABLE_MPI=ON");
    }
    plan.real.assign(
        checkedProduct(
            checkedProduct(layout.local_nx(), m_shape.ny, "PM local real extent"),
            plan.real_z_stride,
            "PM local real extent"),
        0.0);
    plan.fourier.assign(expected_local_complex_size, std::complex<double>(0.0, 0.0));
    plan.potential_k.assign(expected_local_complex_size, std::complex<double>(0.0, 0.0));
    plan.working_k.assign(expected_local_complex_size, std::complex<double>(0.0, 0.0));
    plan.poisson_kernel.assign(expected_local_complex_size, 0.0);
    plan.grad_kx.assign(expected_local_complex_size, 0.0);
    plan.grad_ky.assign(expected_local_complex_size, 0.0);
    plan.grad_kz.assign(expected_local_complex_size, 0.0);
#endif

#if COSMOSIM_ENABLE_FFTW && COSMOSIM_ENABLE_MPI
    if (layout.world_size > 1) {
      PlanResources* committed_plan = nullptr;
      runPmCoordinatedPhase(
          "PmSolver distributed FFT plan cache commit",
          layout.world_rank,
          layout.world_size,
          plan_mpi_wait_ms,
          [&]() {
            auto [insert_it, inserted] = m_plan_cache.emplace(key, std::move(plan));
            if (!inserted) {
              throw std::logic_error("Distributed PM FFT plan cache changed during coordinated build");
            }
            ++m_plan_build_count;
            m_active_key = key;
            committed_plan = &insert_it->second;
          });
      return *committed_plan;
    }
#endif
    auto [insert_it, inserted] = m_plan_cache.emplace(key, std::move(plan));
    (void)inserted;
    ++m_plan_build_count;
    m_active_key = key;
    return insert_it->second;
  }

  [[nodiscard]] std::span<double> realGrid() { return activePlan().real; }
  [[nodiscard]] std::span<std::complex<double>> fourierGrid() { return activePlan().fourier; }
  [[nodiscard]] std::span<std::complex<double>> potentialScratch() { return activePlan().potential_k; }
  [[nodiscard]] std::span<std::complex<double>> workingScratch() { return activePlan().working_k; }

  void ensureSpectralOperators(
      PlanResources& plan,
      const BoxLengths& lengths,
      const PmSolveOptions& options,
      const PmGridShape& shape) {
    if (plan.spectral_operators_ready &&
        plan.cached_lx == lengths.lx &&
        plan.cached_ly == lengths.ly &&
        plan.cached_lz == lengths.lz &&
        plan.cached_split_scale == options.tree_pm_split_scale_comoving &&
        plan.cached_scale_factor == options.scale_factor &&
        plan.cached_gravitational_constant_code == options.gravitational_constant_code &&
        plan.cached_window_deconvolution == options.enable_window_deconvolution &&
        plan.cached_assignment_scheme == options.assignment_scheme) {
      return;
    }

    std::fill(plan.poisson_kernel.begin(), plan.poisson_kernel.end(), 0.0);
    std::fill(plan.grad_kx.begin(), plan.grad_kx.end(), 0.0);
    std::fill(plan.grad_ky.begin(), plan.grad_ky.end(), 0.0);
    std::fill(plan.grad_kz.begin(), plan.grad_kz.end(), 0.0);

    const std::size_t nz_complex = shape.nz / 2U + 1U;
    const double prefactor = -4.0 * k_pi * options.gravitational_constant_code;
    const double dkx = 2.0 * k_pi / lengths.lx;
    const double dky = 2.0 * k_pi / lengths.ly;
    const double dkz = 2.0 * k_pi / lengths.lz;

    auto set_entry = [&](std::size_t index, double kx, double ky, double kz) {
      const double k2 = kx * kx + ky * ky + kz * kz;
      if (k2 == 0.0) {
        return;
      }
      double window_correction = 1.0;
      if (options.enable_window_deconvolution) {
        const int window_exponent = assignmentWindowExponent(options.assignment_scheme);
        const double wx = std::pow(sinc(0.5 * kx * lengths.lx / static_cast<double>(shape.nx)), static_cast<double>(window_exponent));
        const double wy = std::pow(sinc(0.5 * ky * lengths.ly / static_cast<double>(shape.ny)), static_cast<double>(window_exponent));
        const double wz = std::pow(sinc(0.5 * kz * lengths.lz / static_cast<double>(shape.nz)), static_cast<double>(window_exponent));
        const double transfer_window = wx * wy * wz;
        window_correction = 1.0 / std::max(transfer_window * transfer_window, 1.0e-12);
      }
      double split_filter = 1.0;
      if (options.tree_pm_split_scale_comoving > 0.0) {
        split_filter = treePmGaussianFourierLongRangeFilter(std::sqrt(k2), options.tree_pm_split_scale_comoving);
      }
      plan.poisson_kernel[index] = prefactor * window_correction * split_filter / k2;
      plan.grad_kx[index] = kx;
      plan.grad_ky[index] = ky;
      plan.grad_kz[index] = kz;
    };

    if (plan.spectral_transposed) {
      for (std::size_t local_iy = 0; local_iy < plan.transposed_local_ny; ++local_iy) {
        const std::size_t iy = plan.transposed_begin_y + local_iy;
        const std::ptrdiff_t ny_mode = iy <= shape.ny / 2U ? static_cast<std::ptrdiff_t>(iy)
                                                           : static_cast<std::ptrdiff_t>(iy) - static_cast<std::ptrdiff_t>(shape.ny);
        const double ky = dky * static_cast<double>(ny_mode);
        for (std::size_t ix = 0; ix < shape.nx; ++ix) {
          const std::ptrdiff_t nx_mode = ix <= shape.nx / 2U ? static_cast<std::ptrdiff_t>(ix)
                                                             : static_cast<std::ptrdiff_t>(ix) - static_cast<std::ptrdiff_t>(shape.nx);
          const double kx = dkx * static_cast<double>(nx_mode);
          for (std::size_t iz = 0; iz < nz_complex; ++iz) {
            const double kz = dkz * static_cast<double>(iz);
            const std::size_t index = (local_iy * shape.nx + ix) * nz_complex + iz;
            set_entry(index, kx, ky, kz);
          }
        }
      }
    } else {
      const std::size_t global_x_begin = plan.layout.owned_x.begin_x;
      for (std::size_t local_ix = 0; local_ix < plan.layout.local_nx(); ++local_ix) {
        const std::size_t ix = global_x_begin + local_ix;
        const std::ptrdiff_t nx_mode = ix <= shape.nx / 2U ? static_cast<std::ptrdiff_t>(ix)
                                                           : static_cast<std::ptrdiff_t>(ix) - static_cast<std::ptrdiff_t>(shape.nx);
        const double kx = dkx * static_cast<double>(nx_mode);
        for (std::size_t iy = 0; iy < shape.ny; ++iy) {
          const std::ptrdiff_t ny_mode = iy <= shape.ny / 2U ? static_cast<std::ptrdiff_t>(iy)
                                                             : static_cast<std::ptrdiff_t>(iy) - static_cast<std::ptrdiff_t>(shape.ny);
          const double ky = dky * static_cast<double>(ny_mode);
          for (std::size_t iz = 0; iz < nz_complex; ++iz) {
            const double kz = dkz * static_cast<double>(iz);
            const std::size_t index = (local_ix * shape.ny + iy) * nz_complex + iz;
            set_entry(index, kx, ky, kz);
          }
        }
      }
    }

    plan.spectral_operators_ready = true;
    plan.cached_lx = lengths.lx;
    plan.cached_ly = lengths.ly;
    plan.cached_lz = lengths.lz;
    plan.cached_split_scale = options.tree_pm_split_scale_comoving;
    plan.cached_scale_factor = options.scale_factor;
    plan.cached_gravitational_constant_code = options.gravitational_constant_code;
    plan.cached_window_deconvolution = options.enable_window_deconvolution;
    plan.cached_assignment_scheme = options.assignment_scheme;
  }

  double forwardFft() {
    const auto start = std::chrono::steady_clock::now();
#if COSMOSIM_ENABLE_FFTW
    fftw_execute(activePlan().forward_plan);
#else
    naiveForwardDft();
#endif
    const auto stop = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(stop - start).count();
  }

  double inverseFft() {
    const auto start = std::chrono::steady_clock::now();
#if COSMOSIM_ENABLE_FFTW
    fftw_execute(activePlan().inverse_plan);
#else
    naiveInverseDft();
#endif
    const auto stop = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(stop - start).count();
  }

  [[nodiscard]] std::size_t planCount() const noexcept { return m_plan_cache.size(); }
  [[nodiscard]] std::size_t planBuildCount() const noexcept { return m_plan_build_count; }

  [[nodiscard]] std::uint64_t nextDistributedExchangeEpoch() {
    if (m_next_distributed_exchange_epoch == std::numeric_limits<std::uint64_t>::max()) {
      throw std::overflow_error(
          "PmSolver distributed exchange epoch would overflow; recreate the solver before another exchange");
    }
    const std::uint64_t epoch = m_next_distributed_exchange_epoch;
    ++m_next_distributed_exchange_epoch;
    return epoch;
  }

  [[nodiscard]] DensityExchangeBuffers& densityExchangeBuffersForLayout(const parallel::PmSlabLayout& layout) {
    if (m_density_exchange.world_size != layout.world_size || m_density_exchange.world_rank != layout.world_rank) {
      m_density_exchange.world_size = layout.world_size;
      m_density_exchange.world_rank = layout.world_rank;
      m_density_exchange.buffers.send_records_by_rank.assign(static_cast<std::size_t>(layout.world_size), {});
      m_density_exchange.buffers.send_counts.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.send_displs.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.recv_counts.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.recv_displs.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.send_counts_bytes.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.send_displs_bytes.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.recv_counts_bytes.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.recv_displs_bytes.assign(static_cast<std::size_t>(layout.world_size), 0);
    }
    return m_density_exchange.buffers;
  }

  [[nodiscard]] InterpolationExchangeBuffers<PmForceContributionRecord>&
  forceInterpolationExchangeBuffersForLayout(const parallel::PmSlabLayout& layout) {
    if (m_force_exchange.world_size != layout.world_size || m_force_exchange.world_rank != layout.world_rank) {
      m_force_exchange.world_size = layout.world_size;
      m_force_exchange.world_rank = layout.world_rank;
      resetInterpolationExchangeForWorld(layout.world_size, m_force_exchange.buffers);
    }
    return m_force_exchange.buffers;
  }

  [[nodiscard]] InterpolationExchangeBuffers<PmPotentialContributionRecord>&
  potentialInterpolationExchangeBuffersForLayout(const parallel::PmSlabLayout& layout) {
    if (m_potential_exchange.world_size != layout.world_size ||
        m_potential_exchange.world_rank != layout.world_rank) {
      m_potential_exchange.world_size = layout.world_size;
      m_potential_exchange.world_rank = layout.world_rank;
      resetInterpolationExchangeForWorld(layout.world_size, m_potential_exchange.buffers);
    }
    return m_potential_exchange.buffers;
  }

 private:
  template <typename ContributionRecord>
  static void resetInterpolationExchangeForWorld(
      int world_size,
      InterpolationExchangeBuffers<ContributionRecord>& buffers) {
    buffers.send_requests_by_rank.assign(static_cast<std::size_t>(world_size), {});
    buffers.send_contribs_by_rank.assign(static_cast<std::size_t>(world_size), {});
    buffers.send_counts.assign(static_cast<std::size_t>(world_size), 0);
    buffers.send_displs.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_counts.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_displs.assign(static_cast<std::size_t>(world_size), 0);
    buffers.send_counts_bytes.assign(static_cast<std::size_t>(world_size), 0);
    buffers.send_displs_bytes.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_counts_bytes.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_displs_bytes.assign(static_cast<std::size_t>(world_size), 0);
    buffers.send_contrib_counts.assign(static_cast<std::size_t>(world_size), 0);
    buffers.send_contrib_displs.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_contrib_counts.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_contrib_displs.assign(static_cast<std::size_t>(world_size), 0);
    buffers.send_contrib_counts_bytes.assign(static_cast<std::size_t>(world_size), 0);
    buffers.send_contrib_displs_bytes.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_contrib_counts_bytes.assign(static_cast<std::size_t>(world_size), 0);
    buffers.recv_contrib_displs_bytes.assign(static_cast<std::size_t>(world_size), 0);
  }

  [[nodiscard]] PlanResources& activePlan() {
    if (!m_active_key.has_value()) {
      throw std::logic_error("PM solver plan has not been initialized for the active slab layout");
    }
    return m_plan_cache.at(*m_active_key);
  }

#if !COSMOSIM_ENABLE_FFTW
  void naiveComplexDft(std::span<std::complex<double>> field, bool forward) {
    const std::size_t nx = m_isolated_workspace.nx;
    const std::size_t ny = m_isolated_workspace.ny;
    const std::size_t nz = m_isolated_workspace.nz;
    const std::size_t total = nx * ny * nz;
    std::vector<std::complex<double>> out(total, {0.0, 0.0});
    const double sign = forward ? -1.0 : 1.0;
    for (std::size_t kx = 0; kx < nx; ++kx) {
      for (std::size_t ky = 0; ky < ny; ++ky) {
        for (std::size_t kz = 0; kz < nz; ++kz) {
          std::complex<double> acc(0.0, 0.0);
          for (std::size_t x = 0; x < nx; ++x) {
            for (std::size_t y = 0; y < ny; ++y) {
              for (std::size_t z = 0; z < nz; ++z) {
                const double phase = sign * 2.0 * k_pi *
                    (static_cast<double>(kx * x) / static_cast<double>(nx) +
                     static_cast<double>(ky * y) / static_cast<double>(ny) +
                     static_cast<double>(kz * z) / static_cast<double>(nz));
                const std::complex<double> euler(std::cos(phase), std::sin(phase));
                acc += field[(x * ny + y) * nz + z] * euler;
              }
            }
          }
          out[(kx * ny + ky) * nz + kz] = acc;
        }
      }
    }
    if (!forward) {
      const double inv_total = 1.0 / static_cast<double>(total);
      for (auto& v : out) {
        v *= inv_total;
      }
    }
    std::copy(out.begin(), out.end(), field.begin());
  }

  void naiveForwardDft() {
    if (!activePlan().layout.ownsFullDomain()) {
      throw std::logic_error(
          "PM solver naiveForwardDft requires full-domain ownership; distributed PM is unavailable without FFTW/MPI");
    }
    auto& real = activePlan().real;
    auto& fourier = activePlan().fourier;
    const std::size_t nz_complex = m_shape.nz / 2U + 1U;
    for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
      for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
        for (std::size_t iz = 0; iz < nz_complex; ++iz) {
          std::complex<double> acc(0.0, 0.0);
          for (std::size_t x = 0; x < m_shape.nx; ++x) {
            for (std::size_t y = 0; y < m_shape.ny; ++y) {
              for (std::size_t z = 0; z < m_shape.nz; ++z) {
                const double phase = -2.0 * k_pi *
                    (static_cast<double>(ix * x) / static_cast<double>(m_shape.nx) +
                     static_cast<double>(iy * y) / static_cast<double>(m_shape.ny) +
                     static_cast<double>(iz * z) / static_cast<double>(m_shape.nz));
                const std::complex<double> euler(std::cos(phase), std::sin(phase));
                const std::size_t rindex = (x * m_shape.ny + y) * m_shape.nz + z;
                acc += real[rindex] * euler;
              }
            }
          }
          const std::size_t cindex = (ix * m_shape.ny + iy) * nz_complex + iz;
          fourier[cindex] = acc;
        }
      }
    }
  }

  void naiveInverseDft() {
    if (!activePlan().layout.ownsFullDomain()) {
      throw std::logic_error(
          "PM solver naiveInverseDft requires full-domain ownership; distributed PM is unavailable without FFTW/MPI");
    }
    auto& real = activePlan().real;
    auto& fourier = activePlan().fourier;
    const std::size_t total = m_shape.cellCount();
    std::fill(real.begin(), real.end(), 0.0);
    const std::size_t nz_complex = m_shape.nz / 2U + 1U;

    for (std::size_t x = 0; x < m_shape.nx; ++x) {
      for (std::size_t y = 0; y < m_shape.ny; ++y) {
        for (std::size_t z = 0; z < m_shape.nz; ++z) {
          std::complex<double> acc(0.0, 0.0);
          for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
            for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
              for (std::size_t iz = 0; iz < nz_complex; ++iz) {
                const double phase = 2.0 * k_pi *
                    (static_cast<double>(ix * x) / static_cast<double>(m_shape.nx) +
                     static_cast<double>(iy * y) / static_cast<double>(m_shape.ny) +
                     static_cast<double>(iz * z) / static_cast<double>(m_shape.nz));
                const std::complex<double> euler(std::cos(phase), std::sin(phase));
                const std::size_t cindex = (ix * m_shape.ny + iy) * nz_complex + iz;
                if (iz == 0 || (m_shape.nz % 2U == 0U && iz == m_shape.nz / 2U)) {
                  acc += fourier[cindex] * euler;
                } else {
                  acc += fourier[cindex] * euler + std::conj(fourier[cindex]) * std::conj(euler);
                }
              }
            }
          }
          real[(x * m_shape.ny + y) * m_shape.nz + z] = acc.real() / static_cast<double>(total);
        }
      }
    }
  }
#endif

  PmGridShape m_shape;
  std::unordered_map<PlanKey, PlanResources, PlanKeyHasher> m_plan_cache;
  std::optional<PlanKey> m_active_key;
  std::size_t m_plan_build_count = 0;
  std::uint64_t m_next_distributed_exchange_epoch = 1;
  struct {
    int world_size = 1;
    int world_rank = 0;
    DensityExchangeBuffers buffers;
  } m_density_exchange{};
  struct {
    int world_size = 1;
    int world_rank = 0;
    InterpolationExchangeBuffers<PmForceContributionRecord> buffers;
  } m_force_exchange{};
  struct {
    int world_size = 1;
    int world_rank = 0;
    InterpolationExchangeBuffers<PmPotentialContributionRecord> buffers;
  } m_potential_exchange{};
  IsolatedOpenWorkspace m_isolated_workspace{};
};

std::size_t PmGridShape::cellCount() const {
  return nx * ny * nz;
}

bool PmGridShape::isValid() const {
  return nx > 0 && ny > 0 && nz > 0;
}

void PmProfiler::reset() {
  m_totals = {};
}

void PmProfiler::append(const PmProfileEvent& event) {
  m_totals.bytes_moved += event.bytes_moved;
  m_totals.routed_density_records += event.routed_density_records;
  m_totals.routed_force_requests += event.routed_force_requests;
  m_totals.routed_potential_requests += event.routed_potential_requests;
  m_totals.routed_density_peer_count += event.routed_density_peer_count;
  m_totals.routed_force_peer_count += event.routed_force_peer_count;
  m_totals.routed_potential_peer_count += event.routed_potential_peer_count;
  m_totals.routed_mpi_bytes_sent += event.routed_mpi_bytes_sent;
  m_totals.routed_mpi_bytes_received += event.routed_mpi_bytes_received;
  m_totals.force_halo_cache_hits += event.force_halo_cache_hits;
  m_totals.isolated_open_root_workspace_estimate_bytes =
      std::max(m_totals.isolated_open_root_workspace_estimate_bytes,
               event.isolated_open_root_workspace_estimate_bytes);
  m_totals.isolated_open_root_workspace_limit_bytes =
      std::max(m_totals.isolated_open_root_workspace_limit_bytes,
               event.isolated_open_root_workspace_limit_bytes);
  m_totals.isolated_open_gather_bytes += event.isolated_open_gather_bytes;
  m_totals.assign_ms += event.assign_ms;
  m_totals.fft_forward_ms += event.fft_forward_ms;
  m_totals.poisson_ms += event.poisson_ms;
  m_totals.gradient_ms += event.gradient_ms;
  m_totals.fft_inverse_ms += event.fft_inverse_ms;
  m_totals.fft_transpose_ms += event.fft_transpose_ms;
  m_totals.interpolate_ms += event.interpolate_ms;
  m_totals.routed_mpi_wait_ms += event.routed_mpi_wait_ms;
  m_totals.transfer_h2d_ms += event.transfer_h2d_ms;
  m_totals.transfer_d2h_ms += event.transfer_d2h_ms;
  m_totals.device_kernel_ms += event.device_kernel_ms;
  m_totals.fft_transpose_bytes += event.fft_transpose_bytes;
}

const PmProfileEvent& PmProfiler::totals() const {
  return m_totals;
}

PmGridStorage::PmGridStorage(PmGridShape shape)
    : PmGridStorage(
          shape,
          parallel::makePmSlabLayout(shape.nx, shape.ny, shape.nz, /*world_size=*/1, /*world_rank=*/0)) {}

PmGridStorage::PmGridStorage(PmGridShape shape, parallel::PmSlabLayout layout)
    : m_shape(shape),
      m_layout(std::move(layout)),
      m_density(m_layout.localCellCount(), 0.0),
      m_potential(m_layout.localCellCount(), 0.0),
      m_force_x(m_layout.localCellCount(), 0.0),
      m_force_y(m_layout.localCellCount(), 0.0),
      m_force_z(m_layout.localCellCount(), 0.0) {
  if (!m_shape.isValid()) {
    throw std::invalid_argument("PM grid shape must be valid");
  }
  if (!m_layout.isValid()) {
    throw std::invalid_argument("PM slab layout must be valid");
  }
  if (m_layout.global_nx != m_shape.nx || m_layout.global_ny != m_shape.ny || m_layout.global_nz != m_shape.nz) {
    throw std::invalid_argument("PM slab layout global shape must match PM grid shape");
  }
}

const PmGridShape& PmGridStorage::shape() const {
  return m_shape;
}

const parallel::PmSlabLayout& PmGridStorage::slabLayout() const {
  return m_layout;
}

bool PmGridStorage::ownsFullDomain() const noexcept {
  return m_layout.ownsFullDomain();
}

std::size_t PmGridStorage::localCellCount() const noexcept {
  return m_layout.localCellCount();
}

std::span<double> PmGridStorage::density() {
  return m_density;
}

std::span<const double> PmGridStorage::density() const {
  return m_density;
}

std::span<double> PmGridStorage::potential() {
  return m_potential;
}

std::span<const double> PmGridStorage::potential() const {
  return m_potential;
}

std::span<double> PmGridStorage::force_x() {
  return m_force_x;
}

std::span<const double> PmGridStorage::force_x() const {
  return m_force_x;
}

std::span<double> PmGridStorage::force_y() {
  return m_force_y;
}

std::span<const double> PmGridStorage::force_y() const {
  return m_force_y;
}

std::span<double> PmGridStorage::force_z() {
  return m_force_z;
}

std::span<const double> PmGridStorage::force_z() const {
  return m_force_z;
}

void PmGridStorage::clearForceHaloCache() {
  m_force_halo_cache = ForceHaloCache{};
}

void PmGridStorage::setForceHaloCache(
    const parallel::PmSlabHaloExchangeResult& force_x_halo,
    const parallel::PmSlabHaloExchangeResult& force_y_halo,
    const parallel::PmSlabHaloExchangeResult& force_z_halo,
    std::uint64_t exchange_sequence) {
  const auto require_same_shape = [&](const parallel::PmSlabHaloExchangeResult& component, std::string_view label) {
    if (component.halo_depth_x != force_x_halo.halo_depth_x ||
        component.left_peer_rank != force_x_halo.left_peer_rank ||
        component.right_peer_rank != force_x_halo.right_peer_rank ||
        component.left_halo.size() != force_x_halo.left_halo.size() ||
        component.right_halo.size() != force_x_halo.right_halo.size()) {
      throw std::invalid_argument(std::string("PM force halo cache component shape mismatch for ") + std::string(label));
    }
  };
  require_same_shape(force_y_halo, "force_y");
  require_same_shape(force_z_halo, "force_z");

  const std::size_t plane_size = m_layout.global_ny * m_layout.global_nz;
  if (force_x_halo.halo_depth_x == 0 || plane_size == 0) {
    clearForceHaloCache();
    return;
  }
  const std::size_t expected_values = force_x_halo.halo_depth_x * plane_size;
  if (force_x_halo.left_halo.size() != expected_values || force_x_halo.right_halo.size() != expected_values) {
    throw std::invalid_argument("PM force halo cache received a halo payload with inconsistent plane count");
  }

  m_force_halo_cache.left_force_x = force_x_halo.left_halo;
  m_force_halo_cache.left_force_y = force_y_halo.left_halo;
  m_force_halo_cache.left_force_z = force_z_halo.left_halo;
  m_force_halo_cache.right_force_x = force_x_halo.right_halo;
  m_force_halo_cache.right_force_y = force_y_halo.right_halo;
  m_force_halo_cache.right_force_z = force_z_halo.right_halo;
  m_force_halo_cache.halo_depth_x = force_x_halo.halo_depth_x;
  m_force_halo_cache.left_peer_rank = force_x_halo.left_peer_rank;
  m_force_halo_cache.right_peer_rank = force_x_halo.right_peer_rank;
  m_force_halo_cache.exchange_sequence = exchange_sequence;
  m_force_halo_cache.valid = true;
}

bool PmGridStorage::hasForceHaloCache() const noexcept {
  return m_force_halo_cache.valid;
}

bool PmGridStorage::tryLoadForceFromHalo(
    std::size_t global_x,
    std::size_t global_y,
    std::size_t global_z,
    double& force_x_value,
    double& force_y_value,
    double& force_z_value) const {
  if (!m_force_halo_cache.valid || global_y >= m_layout.global_ny || global_z >= m_layout.global_nz ||
      m_force_halo_cache.halo_depth_x == 0 || m_layout.global_nx == 0) {
    return false;
  }
  const std::size_t depth = m_force_halo_cache.halo_depth_x;
  const std::size_t plane_size = m_layout.global_ny * m_layout.global_nz;
  const auto plane_offset = [&](std::size_t halo_x) {
    return halo_x * plane_size + global_y * m_layout.global_nz + global_z;
  };

  if (m_force_halo_cache.left_peer_rank >= 0) {
    for (std::size_t halo_x = 0; halo_x < depth; ++halo_x) {
      const std::size_t halo_global_x = (m_layout.owned_x.begin_x + m_layout.global_nx - depth + halo_x) %
          m_layout.global_nx;
      if (halo_global_x == global_x) {
        const std::size_t index = plane_offset(halo_x);
        force_x_value = m_force_halo_cache.left_force_x[index];
        force_y_value = m_force_halo_cache.left_force_y[index];
        force_z_value = m_force_halo_cache.left_force_z[index];
        return true;
      }
    }
  }

  if (m_force_halo_cache.right_peer_rank >= 0) {
    for (std::size_t halo_x = 0; halo_x < depth; ++halo_x) {
      const std::size_t halo_global_x = (m_layout.owned_x.end_x + halo_x) % m_layout.global_nx;
      if (halo_global_x == global_x) {
        const std::size_t index = plane_offset(halo_x);
        force_x_value = m_force_halo_cache.right_force_x[index];
        force_y_value = m_force_halo_cache.right_force_y[index];
        force_z_value = m_force_halo_cache.right_force_z[index];
        return true;
      }
    }
  }
  return false;
}

std::size_t PmGridStorage::linearIndex(std::size_t ix, std::size_t iy, std::size_t iz) const {
  return m_layout.localLinearIndex(ix, iy, iz);
}

void PmGridStorage::clear() {
  std::fill(m_density.begin(), m_density.end(), 0.0);
  std::fill(m_potential.begin(), m_potential.end(), 0.0);
  std::fill(m_force_x.begin(), m_force_x.end(), 0.0);
  std::fill(m_force_y.begin(), m_force_y.end(), 0.0);
  std::fill(m_force_z.begin(), m_force_z.end(), 0.0);
  clearForceHaloCache();
}

void PmGridStorage::appendMemoryReport(core::MemoryReportBuilder& builder) const {
  const auto add = [&builder](std::string label, const auto& container) {
    const std::uint64_t bytes = core::ownedCapacityBytesForContainer(container);
    builder.addEntry(core::MemoryEntry{.subsystem = core::MemorySubsystem::kPmMesh,
                                       .lifetime = core::MemoryLifetime::kTransient,
                                       .label = std::move(label),
                                       .owned_capacity_bytes = bytes,
                                       .high_water_bytes = bytes});
  };
  add("pm_mesh.density", m_density);
  add("pm_mesh.potential", m_potential);
  add("pm_mesh.force_x", m_force_x);
  add("pm_mesh.force_y", m_force_y);
  add("pm_mesh.force_z", m_force_z);
  add("pm_mesh.force_halo_left_x", m_force_halo_cache.left_force_x);
  add("pm_mesh.force_halo_left_y", m_force_halo_cache.left_force_y);
  add("pm_mesh.force_halo_left_z", m_force_halo_cache.left_force_z);
  add("pm_mesh.force_halo_right_x", m_force_halo_cache.right_force_x);
  add("pm_mesh.force_halo_right_y", m_force_halo_cache.right_force_y);
  add("pm_mesh.force_halo_right_z", m_force_halo_cache.right_force_z);
}

PmSolver::PmSolver(PmGridShape shape) : m_shape(shape), m_impl(std::make_unique<Impl>(shape)) {
  if (!shape.isValid()) {
    throw std::invalid_argument("PM solver requires valid shape");
  }
}

PmSolver::~PmSolver() = default;
PmSolver::PmSolver(PmSolver&&) noexcept = default;
PmSolver& PmSolver::operator=(PmSolver&&) noexcept = default;

const PmGridShape& PmSolver::shape() const {
  return m_shape;
}

void PmSolver::assignDensity(
    PmGridStorage& grid,
    const PmMassSourceView& source_view,
    const PmSolveOptions& options,
    PmProfileEvent* profile) const {
  assignDensity(
      grid,
      source_view.pos_x_comoving,
      source_view.pos_y_comoving,
      source_view.pos_z_comoving,
      source_view.mass_code,
      options,
      profile);
}

void PmSolver::assignDensity(
    PmGridStorage& grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    const PmSolveOptions& options,
    PmProfileEvent* profile) const {
#if COSMOSIM_ENABLE_MPI
  double distributed_preflight_mpi_wait_ms = 0.0;
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
  // Public PM entry is collective whenever the hosting MPI world has peers.
  // The vote precedes any layout-based branch, so divergent rank-local layout
  // metadata cannot let one rank throw while another enters an exchange.
  if (mpi_world_size > 1) {
    const bool rank_local_serial_layout =
        grid.slabLayout().world_size == 1 && grid.ownsFullDomain();
    validatePmCollectiveEntryConsensus(
        PmCollectiveEntryKind::kAssignDensity,
        rank_local_serial_layout,
        m_shape,
        options,
        distributed_preflight_mpi_wait_ms);
    runPmCoordinatedPhase(
        "PmSolver::assignDensity distributed API preflight",
        mpi_world_rank,
        mpi_world_size,
        distributed_preflight_mpi_wait_ms,
        [&]() {
          validateOptions(m_shape, options);
          if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
            throw std::invalid_argument("PM solver/grid shape mismatch in assignDensity");
          }
          if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() || pos_x.size() != mass.size()) {
            throw std::invalid_argument("Particle coordinate/mass spans must match in assignDensity");
          }
          if (!grid.slabLayout().isValid()) {
            throw std::invalid_argument("PmSolver::assignDensity requires a valid PM slab layout");
          }
          if (!rank_local_serial_layout &&
              (mpi_world_size != grid.slabLayout().world_size ||
               mpi_world_rank != grid.slabLayout().world_rank)) {
            throw std::invalid_argument(
                "PmSolver::assignDensity slab layout world metadata must match MPI_COMM_WORLD");
          }
          if (m_shape.nx > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
              m_shape.ny > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
              m_shape.nz > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
            throw std::invalid_argument("PmSolver::assignDensity mesh dimensions exceed fixed-width routing indices");
          }
        });
  }
#endif
  validateOptions(m_shape, options);
  if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
    throw std::invalid_argument("PM solver/grid shape mismatch in assignDensity");
  }
  if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() || pos_x.size() != mass.size()) {
    throw std::invalid_argument("Particle coordinate/mass spans must match in assignDensity");
  }
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::assignDensity requires a valid PM slab layout");
  }
  const bool distributed_slabs = grid.slabLayout().world_size > 1;
#if !COSMOSIM_ENABLE_MPI
  if (distributed_slabs) {
    throw std::invalid_argument(
        "PmSolver::assignDensity distributed slabs require COSMOSIM_ENABLE_MPI=ON");
  }
#else
  if (distributed_slabs) {
    int mpi_world_size = 1;
    int mpi_world_rank = 0;
    queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
    if (mpi_world_size != grid.slabLayout().world_size || mpi_world_rank != grid.slabLayout().world_rank) {
      throw std::invalid_argument(
          "PmSolver::assignDensity slab layout world metadata must match MPI_COMM_WORLD");
    }
  } else {
    validateSingleRankFullDomainGridContract(grid, "PmSolver::assignDensity");
  }
#endif

  const auto start = std::chrono::steady_clock::now();
  std::fill(grid.density().begin(), grid.density().end(), 0.0);

  const BoxLengths lengths = effectiveBoxLengths(options);
  const double inv_dx = static_cast<double>(m_shape.nx) / lengths.lx;
  const double inv_dy = static_cast<double>(m_shape.ny) / lengths.ly;
  const double inv_dz = static_cast<double>(m_shape.nz) / lengths.lz;
  if (m_shape.nx > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
      m_shape.ny > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
      m_shape.nz > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
    throw std::invalid_argument("PmSolver::assignDensity mesh dimensions exceed fixed-width routing indices");
  }
  std::optional<std::size_t> first_non_finite_source_index;
  for (std::size_t p = 0; p < mass.size(); ++p) {
    if (!std::isfinite(pos_x[p]) || !std::isfinite(pos_y[p]) || !std::isfinite(pos_z[p]) ||
        !std::isfinite(mass[p])) {
      first_non_finite_source_index = p;
      break;
    }
  }
  if (!distributed_slabs && first_non_finite_source_index.has_value()) {
    throw std::invalid_argument(
        "PmSolver::assignDensity requires finite particle coordinates and masses; particle_index=" +
        std::to_string(*first_non_finite_source_index));
  }

  const auto accumulate_owned = [&](const PmDensityContributionRecord& record) {
    if (record.global_ix >= m_shape.nx || record.global_iy >= m_shape.ny || record.global_iz >= m_shape.nz) {
      throw std::invalid_argument("PmSolver::assignDensity received out-of-range PM contribution record");
    }
    if (!grid.slabLayout().ownsGlobalX(record.global_ix)) {
      throw std::invalid_argument("PmSolver::assignDensity received contribution for non-owned PM slab x-index");
    }
    if (!std::isfinite(record.mass_contribution)) {
      throw std::invalid_argument("PmSolver::assignDensity received a non-finite mass contribution");
    }
    grid.density()[grid.linearIndex(record.global_ix, record.global_iy, record.global_iz)] += record.mass_contribution;
  };

  if (!distributed_slabs) {
    for (std::size_t p = 0; p < mass.size(); ++p) {
      const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
      if (!periodic &&
          (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
           !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
           !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
        continue;
      }
      const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
      const double y = (periodic ? wrapPosition(pos_y[p], lengths.ly) : pos_y[p]) * inv_dy;
      const double z = (periodic ? wrapPosition(pos_z[p], lengths.lz) : pos_z[p]) * inv_dz;

      const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
      const PmAxisStencil1d sy = makeAxisStencil(y, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(z, options.assignment_scheme);

      for (std::size_t dx = 0; dx < sx.count; ++dx) {
        if (!periodic && (sx.offsets[dx] < 0 || sx.offsets[dx] >= static_cast<int>(m_shape.nx))) {
          continue;
        }
        const std::size_t ix = periodic ? wrapIndex(sx.offsets[dx], m_shape.nx) : static_cast<std::size_t>(sx.offsets[dx]);
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          if (!periodic && (sy.offsets[dy] < 0 || sy.offsets[dy] >= static_cast<int>(m_shape.ny))) {
            continue;
          }
          const std::size_t iy = periodic ? wrapIndex(sy.offsets[dy], m_shape.ny) : static_cast<std::size_t>(sy.offsets[dy]);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            if (!periodic && (sz.offsets[dz] < 0 || sz.offsets[dz] >= static_cast<int>(m_shape.nz))) {
              continue;
            }
            const std::size_t iz = periodic ? wrapIndex(sz.offsets[dz], m_shape.nz) : static_cast<std::size_t>(sz.offsets[dz]);
            const double weight = sx.weights[dx] * sy.weights[dy] * sz.weights[dz];
            accumulate_owned(PmDensityContributionRecord{
                .global_ix = static_cast<std::uint32_t>(ix),
                .global_iy = static_cast<std::uint32_t>(iy),
                .global_iz = static_cast<std::uint32_t>(iz),
                .mass_contribution = mass[p] * weight,
            });
          }
        }
      }
    }
  } else {
#if COSMOSIM_ENABLE_MPI
    double routed_mpi_wait_ms = distributed_preflight_mpi_wait_ms;
    const int local_sources_valid = first_non_finite_source_index.has_value() ? 0 : 1;
    int all_sources_valid = 0;
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      MPI_Allreduce(&local_sources_valid, &all_sources_valid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    });
    if (all_sources_valid == 0) {
      const std::string local_detail = first_non_finite_source_index.has_value()
          ? " local_first_invalid_particle_index=" + std::to_string(*first_non_finite_source_index)
          : " local_sources_valid=true";
      throw std::invalid_argument(
          "PmSolver::assignDensity rejected non-finite coordinates or masses on at least one MPI rank;" +
          local_detail);
    }
    std::uint64_t exchange_epoch = 0U;
    runPmCoordinatedPhase(
        "PmSolver::assignDensity exchange epoch allocation",
        grid.slabLayout().world_rank,
        grid.slabLayout().world_size,
        routed_mpi_wait_ms,
        [&]() { exchange_epoch = m_impl->nextDistributedExchangeEpoch(); });
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      validatePmExchangeEpochConsensus(exchange_epoch, "PmSolver::assignDensity");
    });
    constexpr std::string_view density_exchange_context = "PmSolver::assignDensity density contribution exchange";
    Impl::DensityExchangeBuffers* exchange_ptr = nullptr;
    runPmCoordinatedPhase(
        "PmSolver::assignDensity density send preparation",
        grid.slabLayout().world_rank,
        grid.slabLayout().world_size,
        routed_mpi_wait_ms,
        [&]() {
    exchange_ptr = &m_impl->densityExchangeBuffersForLayout(grid.slabLayout());
    auto& exchange = *exchange_ptr;
    for (auto& per_rank : exchange.send_records_by_rank) {
      per_rank.clear();
    }
    exchange.send_flat_records.clear();
    exchange.recv_flat_records.clear();
    exchange.send_wire.clear();
    exchange.recv_wire.clear();
    std::fill(exchange.send_counts.begin(), exchange.send_counts.end(), 0);
    std::fill(exchange.send_displs.begin(), exchange.send_displs.end(), 0);
    std::fill(exchange.recv_counts.begin(), exchange.recv_counts.end(), 0);
    std::fill(exchange.recv_displs.begin(), exchange.recv_displs.end(), 0);
    std::fill(exchange.send_counts_bytes.begin(), exchange.send_counts_bytes.end(), 0);
    std::fill(exchange.send_displs_bytes.begin(), exchange.send_displs_bytes.end(), 0);
    std::fill(exchange.recv_counts_bytes.begin(), exchange.recv_counts_bytes.end(), 0);
    std::fill(exchange.recv_displs_bytes.begin(), exchange.recv_displs_bytes.end(), 0);

    std::uint32_t next_record_sequence = 0U;
    for (std::size_t p = 0; p < mass.size(); ++p) {
      const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
      if (!periodic &&
          (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
           !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
           !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
        continue;
      }
      const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
      const double y = (periodic ? wrapPosition(pos_y[p], lengths.ly) : pos_y[p]) * inv_dy;
      const double z = (periodic ? wrapPosition(pos_z[p], lengths.lz) : pos_z[p]) * inv_dz;
      const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
      const PmAxisStencil1d sy = makeAxisStencil(y, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(z, options.assignment_scheme);

      for (std::size_t dx = 0; dx < sx.count; ++dx) {
        if (!periodic && (sx.offsets[dx] < 0 || sx.offsets[dx] >= static_cast<std::ptrdiff_t>(m_shape.nx))) {
          continue;
        }
        const std::size_t ix = periodic
            ? wrapIndex(sx.offsets[dx], m_shape.nx)
            : static_cast<std::size_t>(sx.offsets[dx]);
        const int destination_rank = parallel::pmOwnerRankForGlobalX(m_shape.nx, grid.slabLayout().world_size, ix);
        auto& batch = exchange.send_records_by_rank[static_cast<std::size_t>(destination_rank)];
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          if (!periodic && (sy.offsets[dy] < 0 || sy.offsets[dy] >= static_cast<std::ptrdiff_t>(m_shape.ny))) {
            continue;
          }
          const std::size_t iy = periodic
              ? wrapIndex(sy.offsets[dy], m_shape.ny)
              : static_cast<std::size_t>(sy.offsets[dy]);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            if (!periodic && (sz.offsets[dz] < 0 || sz.offsets[dz] >= static_cast<std::ptrdiff_t>(m_shape.nz))) {
              continue;
            }
            const std::size_t iz = periodic
                ? wrapIndex(sz.offsets[dz], m_shape.nz)
                : static_cast<std::size_t>(sz.offsets[dz]);
            const double weight = sx.weights[dx] * sy.weights[dy] * sz.weights[dz];
            if (next_record_sequence == std::numeric_limits<std::uint32_t>::max()) {
              throw std::invalid_argument("PmSolver::assignDensity record sequence exceeds routing token limit");
            }
            batch.push_back(PmDensityContributionRecord{
                .origin_rank = static_cast<std::uint32_t>(grid.slabLayout().world_rank),
                .destination_rank = static_cast<std::uint32_t>(destination_rank),
                .record_sequence = next_record_sequence,
                .global_ix = static_cast<std::uint32_t>(ix),
                .global_iy = static_cast<std::uint32_t>(iy),
                .global_iz = static_cast<std::uint32_t>(iz),
                .exchange_epoch = exchange_epoch,
                .mass_contribution = mass[p] * weight,
            });
            ++next_record_sequence;
          }
        }
      }
    }

    std::size_t total_send_records = 0;
    for (int rank = 0; rank < grid.slabLayout().world_size; ++rank) {
      const std::size_t count = exchange.send_records_by_rank[static_cast<std::size_t>(rank)].size();
      const int mpi_count = checkedMpiRecordCountOrDisplacement(
          count, k_pm_density_wire_bytes, density_exchange_context, "send record count", rank);
      const int mpi_displacement = checkedMpiRecordCountOrDisplacement(
          total_send_records, k_pm_density_wire_bytes, density_exchange_context, "send record displacement", rank);
      const std::size_t next_total = checkedMpiRecordTotal(
          total_send_records, count, k_pm_density_wire_bytes, density_exchange_context, rank);
      exchange.send_counts[static_cast<std::size_t>(rank)] = mpi_count;
      exchange.send_displs[static_cast<std::size_t>(rank)] = mpi_displacement;
      total_send_records = next_total;
    }
    exchange.send_flat_records.reserve(total_send_records);
    for (int rank = 0; rank < grid.slabLayout().world_size; ++rank) {
      const auto& records = exchange.send_records_by_rank[static_cast<std::size_t>(rank)];
      exchange.send_flat_records.insert(exchange.send_flat_records.end(), records.begin(), records.end());
    }
    encodePmWireRecords<PmDensityContributionRecord>(
        exchange.send_flat_records,
        k_pm_density_wire_bytes,
        exchange.send_wire,
        [](const PmDensityContributionRecord& record, std::span<std::uint8_t> bytes) {
          encodePmDensityRecord(record, bytes);
        },
        density_exchange_context);
        });
    auto& exchange = *exchange_ptr;

    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      MPI_Alltoall(
          exchange.send_counts.data(),
          1,
          MPI_INT,
          exchange.recv_counts.data(),
          1,
          MPI_INT,
          MPI_COMM_WORLD);
    });

    runPmCoordinatedPhase(
        "PmSolver::assignDensity density receive-layout preparation",
        grid.slabLayout().world_rank,
        grid.slabLayout().world_size,
        routed_mpi_wait_ms,
        [&]() {
    std::size_t total_recv_records = 0;
    for (int rank = 0; rank < grid.slabLayout().world_size; ++rank) {
      const std::size_t count = checkedMpiReceivedRecordCount(
          exchange.recv_counts[static_cast<std::size_t>(rank)],
          k_pm_density_wire_bytes,
          density_exchange_context,
          rank);
      const int mpi_displacement = checkedMpiRecordCountOrDisplacement(
          total_recv_records, k_pm_density_wire_bytes, density_exchange_context, "received record displacement", rank);
      const std::size_t next_total = checkedMpiRecordTotal(
          total_recv_records, count, k_pm_density_wire_bytes, density_exchange_context, rank);
      exchange.recv_displs[static_cast<std::size_t>(rank)] = mpi_displacement;
      total_recv_records = next_total;
    }
    exchange.recv_wire.resize(checkedPmWireByteCount(
        total_recv_records, k_pm_density_wire_bytes, density_exchange_context));

    checkedMpiRecordLayoutToByteLayout(
        exchange.send_counts,
        exchange.send_displs,
        k_pm_density_wire_bytes,
        exchange.send_counts_bytes,
        exchange.send_displs_bytes,
        "PmSolver::assignDensity density contribution send MPI_Alltoallv");
    checkedMpiRecordLayoutToByteLayout(
        exchange.recv_counts,
        exchange.recv_displs,
        k_pm_density_wire_bytes,
        exchange.recv_counts_bytes,
        exchange.recv_displs_bytes,
        "PmSolver::assignDensity density contribution receive MPI_Alltoallv");
        });

    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      MPI_Alltoallv(
          nonNullPmWireData(exchange.send_wire),
          exchange.send_counts_bytes.data(),
          exchange.send_displs_bytes.data(),
          MPI_BYTE,
          nonNullPmWireData(exchange.recv_wire),
          exchange.recv_counts_bytes.data(),
          exchange.recv_displs_bytes.data(),
          MPI_BYTE,
          MPI_COMM_WORLD);
    });

    std::string local_density_validation_error;
    try {
      decodePmWireRecords<PmDensityContributionRecord>(
          exchange.recv_wire,
          k_pm_density_wire_bytes,
          exchange.recv_flat_records,
          [](std::span<const std::uint8_t> bytes) { return decodePmDensityRecord(bytes); },
          density_exchange_context);

      for (int source_rank = 0; source_rank < grid.slabLayout().world_size; ++source_rank) {
        const int begin = exchange.recv_displs[static_cast<std::size_t>(source_rank)];
        const int count = exchange.recv_counts[static_cast<std::size_t>(source_rank)];
        std::optional<std::uint32_t> previous_sequence;
        for (int i = 0; i < count; ++i) {
          const PmDensityContributionRecord& record =
              exchange.recv_flat_records[static_cast<std::size_t>(begin + i)];
          if (record.origin_rank != static_cast<std::uint32_t>(source_rank)) {
            throw std::invalid_argument(
                "density origin rank does not match MPI sender segment");
          }
          if (record.destination_rank != static_cast<std::uint32_t>(grid.slabLayout().world_rank)) {
            throw std::invalid_argument(
                "density destination rank does not match receiving slab rank");
          }
          if (record.exchange_epoch != exchange_epoch) {
            throw std::invalid_argument(
                "density exchange epoch is stale/mismatched; received=" +
                std::to_string(record.exchange_epoch) + " expected=" + std::to_string(exchange_epoch));
          }
          if (previous_sequence.has_value() && record.record_sequence <= *previous_sequence) {
            throw std::invalid_argument(
                "density record sequence is duplicate or non-monotonic for sender rank " +
                std::to_string(source_rank));
          }
          previous_sequence = record.record_sequence;
          if (record.global_ix >= m_shape.nx || record.global_iy >= m_shape.ny || record.global_iz >= m_shape.nz) {
            throw std::invalid_argument("density contribution mesh index is out of range");
          }
          const int expected_owner =
              parallel::pmOwnerRankForGlobalX(m_shape.nx, grid.slabLayout().world_size, record.global_ix);
          if (expected_owner != grid.slabLayout().world_rank ||
              expected_owner != static_cast<int>(record.destination_rank)) {
            throw std::invalid_argument(
                "density contribution expected owner does not match routed destination");
          }
          if (!std::isfinite(record.mass_contribution)) {
            throw std::invalid_argument("density contribution mass is not finite");
          }
        }
      }
    } catch (const std::exception& error) {
      local_density_validation_error = error.what();
    } catch (...) {
      local_density_validation_error = "non-standard exception while validating density payload";
    }
    throwIfPmPayloadValidationFailed(
        "PmSolver::assignDensity density receive",
        grid.slabLayout().world_rank,
        grid.slabLayout().world_size,
        local_density_validation_error,
        routed_mpi_wait_ms);

    for (const PmDensityContributionRecord& record : exchange.recv_flat_records) {
      grid.density()[grid.linearIndex(record.global_ix, record.global_iy, record.global_iz)] +=
          record.mass_contribution;
    }
    if (profile != nullptr) {
      profile->routed_density_records += static_cast<std::uint64_t>(exchange.send_flat_records.size());
      profile->routed_density_peer_count += static_cast<std::uint64_t>(std::count_if(
          exchange.send_counts.begin(),
          exchange.send_counts.end(),
          [](int count) { return count > 0; }));
      profile->routed_mpi_bytes_sent += static_cast<std::uint64_t>(exchange.send_wire.size());
      profile->routed_mpi_bytes_received += static_cast<std::uint64_t>(exchange.recv_wire.size());
      profile->routed_mpi_wait_ms += routed_mpi_wait_ms;
      profile->bytes_moved += checkedBytesForCount(
          exchange.send_flat_records.size() + exchange.recv_flat_records.size(),
          k_pm_density_wire_bytes,
          "PM density routed exchange profile");
    }
#endif
  }

  const double cell_volume =
      (lengths.lx * lengths.ly * lengths.lz) / static_cast<double>(m_shape.cellCount());
  for (double& density_cell : grid.density()) {
    density_cell /= cell_volume;
  }

  const auto stop = std::chrono::steady_clock::now();
  if (profile != nullptr) {
    profile->assign_ms += std::chrono::duration<double, std::milli>(stop - start).count();
    profile->bytes_moved += bytesForGridSweep(m_shape.cellCount());
    profile->bytes_moved += bytesForParticles(pos_x.size());
  }
}

void PmSolver::solvePoissonPeriodic(PmGridStorage& grid, const PmSolveOptions& options, PmProfileEvent* profile) {
#if COSMOSIM_ENABLE_MPI
  double distributed_plan_preflight_mpi_wait_ms = 0.0;
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
  if (mpi_world_size > 1) {
    const bool rank_local_serial_layout =
        grid.slabLayout().world_size == 1 && grid.ownsFullDomain();
    validatePmCollectiveEntryConsensus(
        PmCollectiveEntryKind::kSolvePoissonPeriodic,
        rank_local_serial_layout,
        m_shape,
        options,
        distributed_plan_preflight_mpi_wait_ms);
    runPmCoordinatedPhase(
        "PmSolver::solvePoissonPeriodic distributed API preflight",
        mpi_world_rank,
        mpi_world_size,
        distributed_plan_preflight_mpi_wait_ms,
        [&]() {
          validateOptions(m_shape, options);
          if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
            throw std::invalid_argument("PM solver/grid shape mismatch in solvePoissonPeriodic");
          }
          if (!grid.slabLayout().isValid()) {
            throw std::invalid_argument("PmSolver::solvePoissonPeriodic requires valid PM slab layout");
          }
          if (!rank_local_serial_layout &&
              (mpi_world_size != grid.slabLayout().world_size ||
               mpi_world_rank != grid.slabLayout().world_rank)) {
            throw std::invalid_argument(
                "PmSolver::solvePoissonPeriodic slab layout world metadata must match MPI_COMM_WORLD");
          }
        });
    if (profile != nullptr) {
      profile->routed_mpi_wait_ms += distributed_plan_preflight_mpi_wait_ms;
    }
  }
#endif
  validateOptions(m_shape, options);
  if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
    throw std::invalid_argument("PM solver/grid shape mismatch in solvePoissonPeriodic");
  }
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::solvePoissonPeriodic requires valid PM slab layout");
  }
  if (grid.slabLayout().world_size > 1) {
#if !(COSMOSIM_ENABLE_FFTW && COSMOSIM_ENABLE_MPI)
    throw std::invalid_argument(
        "PmSolver::solvePoissonPeriodic distributed slabs require COSMOSIM_ENABLE_FFTW=ON and COSMOSIM_ENABLE_MPI=ON");
#endif
  } else if (!grid.ownsFullDomain()) {
    throw std::invalid_argument("PmSolver::solvePoissonPeriodic single-rank path requires full-domain slab ownership");
  }

  auto& plan = m_impl->planForLayout(grid.slabLayout(), options.decomposition_mode);
  auto real = m_impl->realGrid();
  auto fourier = m_impl->fourierGrid();
  auto potential_k = m_impl->potentialScratch();
  auto working_k = m_impl->workingScratch();

  std::fill(real.begin(), real.end(), 0.0);
  if (plan.is_distributed) {
    for (std::size_t local_ix = 0; local_ix < grid.slabLayout().local_nx(); ++local_ix) {
      const std::size_t global_ix = grid.slabLayout().globalXFromLocal(local_ix);
      for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
        const std::size_t compact_base = grid.linearIndex(global_ix, iy, 0);
        const std::size_t fftw_base = (local_ix * m_shape.ny + iy) * plan.real_z_stride;
        std::copy_n(grid.density().begin() + static_cast<std::ptrdiff_t>(compact_base), m_shape.nz, real.begin() + static_cast<std::ptrdiff_t>(fftw_base));
      }
    }
  } else {
    std::copy(grid.density().begin(), grid.density().end(), real.begin());
  }

  double local_density_sum = std::accumulate(real.begin(), real.end(), 0.0);
  double global_density_sum = local_density_sum;
#if COSMOSIM_ENABLE_MPI
  if (grid.slabLayout().world_size > 1) {
    MPI_Allreduce(&local_density_sum, &global_density_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
#endif
  const double mean_density = global_density_sum / static_cast<double>(m_shape.cellCount());
  for (double& value : real) {
    value -= mean_density;
  }

  const double forward_fft_ms = m_impl->forwardFft();
  if (profile != nullptr) {
    profile->fft_forward_ms += forward_fft_ms;
    if (plan.spectral_transposed) {
      profile->fft_transpose_ms += forward_fft_ms;
    }
  }

  const auto poisson_start = std::chrono::steady_clock::now();
  const BoxLengths lengths = effectiveBoxLengths(options);
  m_impl->ensureSpectralOperators(plan, lengths, options, m_shape);
  for (std::size_t i = 0; i < fourier.size(); ++i) {
    fourier[i] *= plan.poisson_kernel[i];
  }
  const auto poisson_stop = std::chrono::steady_clock::now();

  if (profile != nullptr) {
    profile->poisson_ms += std::chrono::duration<double, std::milli>(poisson_stop - poisson_start).count();
  }

  std::copy(fourier.begin(), fourier.end(), potential_k.begin());

  auto inverse_into = [this, profile, &plan, &grid](std::span<const std::complex<double>> src, std::span<double> dst) {
    auto fourier_dst = m_impl->fourierGrid();
    std::copy(src.begin(), src.end(), fourier_dst.begin());
    const double fft_time = m_impl->inverseFft();
    auto real_values = m_impl->realGrid();
#if COSMOSIM_ENABLE_FFTW
    const double normalization = 1.0 / static_cast<double>(m_shape.cellCount());
    if (plan.is_distributed) {
      for (std::size_t local_ix = 0; local_ix < plan.layout.local_nx(); ++local_ix) {
        const std::size_t global_ix = plan.layout.globalXFromLocal(local_ix);
        for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
          const std::size_t compact_base = grid.linearIndex(global_ix, iy, 0);
          const std::size_t fftw_base = (local_ix * m_shape.ny + iy) * plan.real_z_stride;
          for (std::size_t iz = 0; iz < m_shape.nz; ++iz) {
            dst[compact_base + iz] = real_values[fftw_base + iz] * normalization;
          }
        }
      }
    } else {
      for (std::size_t i = 0; i < dst.size(); ++i) {
        dst[i] = real_values[i] * normalization;
      }
    }
#else
    std::copy(real_values.begin(), real_values.end(), dst.begin());
#endif
    if (profile != nullptr) {
      profile->fft_inverse_ms += fft_time;
      if (plan.spectral_transposed) {
        profile->fft_transpose_ms += fft_time;
      }
    }
  };

  const auto grad_start = std::chrono::steady_clock::now();
  inverse_into(potential_k, grid.potential());

  auto fill_gradient_k = [&](std::span<std::complex<double>> dst, int axis) {
    const std::span<const double> grad = axis == 0
        ? std::span<const double>(plan.grad_kx.data(), plan.grad_kx.size())
        : (axis == 1
            ? std::span<const double>(plan.grad_ky.data(), plan.grad_ky.size())
            : std::span<const double>(plan.grad_kz.data(), plan.grad_kz.size()));
    for (std::size_t index = 0; index < dst.size(); ++index) {
      dst[index] = std::complex<double>(0.0, -grad[index]) * potential_k[index];
    }
  };

  fill_gradient_k(working_k, 0);
  inverse_into(working_k, grid.force_x());
  fill_gradient_k(working_k, 1);
  inverse_into(working_k, grid.force_y());
  fill_gradient_k(working_k, 2);
  inverse_into(working_k, grid.force_z());
  const auto grad_stop = std::chrono::steady_clock::now();
  if (profile != nullptr) {
    profile->gradient_ms += std::chrono::duration<double, std::milli>(grad_stop - grad_start).count();
  }

  if (profile != nullptr) {
    profile->bytes_moved += bytesForGridSweep(m_shape.cellCount()) * 6U;
    if (plan.spectral_transposed) {
      profile->fft_transpose_bytes +=
          static_cast<std::uint64_t>(sizeof(std::complex<double>)) * static_cast<std::uint64_t>(fourier.size()) * 8ULL;
    }
  }
}

void PmSolver::solvePoissonIsolatedOpen(PmGridStorage& grid, const PmSolveOptions& options, PmProfileEvent* profile) {
#if COSMOSIM_ENABLE_MPI
  double isolated_preflight_mpi_wait_ms = 0.0;
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
  if (mpi_world_size > 1) {
    const bool rank_local_serial_layout =
        grid.slabLayout().world_size == 1 && grid.ownsFullDomain();
    validatePmCollectiveEntryConsensus(
        PmCollectiveEntryKind::kSolvePoissonIsolatedOpen,
        rank_local_serial_layout,
        m_shape,
        options,
        isolated_preflight_mpi_wait_ms);
    runPmCoordinatedPhase(
        "PmSolver::solvePoissonIsolatedOpen API preflight",
        mpi_world_rank,
        mpi_world_size,
        isolated_preflight_mpi_wait_ms,
        [&]() {
          validateOptions(m_shape, options);
          if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny ||
              grid.shape().nz != m_shape.nz) {
            throw std::invalid_argument("PM solver/grid shape mismatch in solvePoissonIsolatedOpen");
          }
          if (!grid.slabLayout().isValid()) {
            throw std::invalid_argument(
                "PmSolver::solvePoissonIsolatedOpen requires valid PM slab layout");
          }
          if (!rank_local_serial_layout &&
              (mpi_world_size != grid.slabLayout().world_size ||
               mpi_world_rank != grid.slabLayout().world_rank)) {
            throw std::invalid_argument(
                "PmSolver::solvePoissonIsolatedOpen slab layout world metadata must match MPI_COMM_WORLD");
          }
          if (m_shape.nx < 3 || m_shape.ny < 3 || m_shape.nz < 3) {
            throw std::invalid_argument(
                "PmSolver::solvePoissonIsolatedOpen requires nx,ny,nz >= 3 for one-sided boundary gradients");
          }
          const bool distributed_layout = !rank_local_serial_layout;
          const std::uint64_t workspace_estimate =
              estimateIsolatedRootWorkspaceBytes(m_shape, distributed_layout);
          if (options.isolated_open_root_workspace_limit_bytes == 0U ||
              workspace_estimate > options.isolated_open_root_workspace_limit_bytes) {
            throw std::runtime_error(isolatedPmGuardMessage(
                m_shape,
                distributed_layout ? grid.slabLayout().world_size : 1,
                distributed_layout ? grid.slabLayout().world_rank : 0,
                workspace_estimate,
                options.isolated_open_root_workspace_limit_bytes));
          }
          if (grid.localCellCount() >
              static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::overflow_error("isolated PM local slab count exceeds MPI int limit");
          }
          if (distributed_layout) {
            std::size_t displacement = 0U;
            for (int rank = 0; rank < grid.slabLayout().world_size; ++rank) {
              const auto owned = parallel::pmOwnedXRangeForRank(
                  grid.slabLayout().global_nx, grid.slabLayout().world_size, rank);
              const std::size_t count = checkedProduct(
                  checkedProduct(
                      owned.extentX(),
                      grid.slabLayout().global_ny,
                      "isolated PM preflight slab count"),
                  grid.slabLayout().global_nz,
                  "isolated PM preflight slab count");
              if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
                  displacement > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                throw std::overflow_error(
                    "isolated PM gather/scatter count or displacement exceeds MPI int limit");
              }
              if (count > std::numeric_limits<std::size_t>::max() - displacement) {
                throw std::overflow_error(
                    "isolated PM preflight gather/scatter displacement overflows size_t");
              }
              displacement += count;
            }
          }
        });
    if (profile != nullptr) {
      profile->routed_mpi_wait_ms += isolated_preflight_mpi_wait_ms;
    }
  }
#endif
  validateOptions(m_shape, options);
  if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
    throw std::invalid_argument("PM solver/grid shape mismatch in solvePoissonIsolatedOpen");
  }
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::solvePoissonIsolatedOpen requires valid PM slab layout");
  }
#if !COSMOSIM_ENABLE_MPI
  if (grid.slabLayout().world_size > 1) {
    throw std::invalid_argument("PmSolver::solvePoissonIsolatedOpen distributed slabs require COSMOSIM_ENABLE_MPI=ON");
  }
#endif

  const bool distributed_slabs = grid.slabLayout().world_size > 1;
#if COSMOSIM_ENABLE_MPI
  if (distributed_slabs) {
    int mpi_world_size = 1;
    int mpi_world_rank = 0;
    queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
    if (mpi_world_size != grid.slabLayout().world_size || mpi_world_rank != grid.slabLayout().world_rank) {
      throw std::invalid_argument("PmSolver::solvePoissonIsolatedOpen slab layout world metadata must match MPI_COMM_WORLD");
    }
  } else
#endif
  if (!grid.ownsFullDomain()) {
    throw std::invalid_argument("PmSolver::solvePoissonIsolatedOpen single-rank path requires full-domain ownership");
  }

  const BoxLengths lengths = effectiveBoxLengths(options);
  const double dx = lengths.lx / static_cast<double>(m_shape.nx);
  const double dy = lengths.ly / static_cast<double>(m_shape.ny);
  const double dz = lengths.lz / static_cast<double>(m_shape.nz);
  if (m_shape.nx < 3 || m_shape.ny < 3 || m_shape.nz < 3) {
    throw std::invalid_argument("PmSolver::solvePoissonIsolatedOpen requires nx,ny,nz >= 3 for one-sided boundary gradients");
  }

  constexpr int k_root = 0;
  const bool is_root = !distributed_slabs || grid.slabLayout().world_rank == k_root;
  const auto poisson_start = std::chrono::steady_clock::now();

  const int guard_world_size = distributed_slabs ? grid.slabLayout().world_size : 1;
  const int guard_world_rank = distributed_slabs ? grid.slabLayout().world_rank : 0;
  const std::uint64_t isolated_root_workspace_estimate_bytes =
      estimateIsolatedRootWorkspaceBytes(m_shape, distributed_slabs);
  if (options.isolated_open_root_workspace_limit_bytes == 0U ||
      isolated_root_workspace_estimate_bytes > options.isolated_open_root_workspace_limit_bytes) {
    throw std::runtime_error(isolatedPmGuardMessage(
        m_shape,
        guard_world_size,
        guard_world_rank,
        isolated_root_workspace_estimate_bytes,
        options.isolated_open_root_workspace_limit_bytes));
  }

  std::vector<double> global_density;
#if COSMOSIM_ENABLE_MPI
  if (distributed_slabs) {
    const std::size_t global_size = checkedProduct(
        checkedProduct(grid.slabLayout().global_nx, grid.slabLayout().global_ny, "isolated PM distributed gather"),
        grid.slabLayout().global_nz,
        "isolated PM distributed gather");
    std::vector<int> recv_counts(static_cast<std::size_t>(grid.slabLayout().world_size), 0);
    std::vector<int> recv_displs(static_cast<std::size_t>(grid.slabLayout().world_size), 0);
    std::size_t disp = 0;
    for (int rank = 0; rank < grid.slabLayout().world_size; ++rank) {
      const auto owned = parallel::pmOwnedXRangeForRank(grid.slabLayout().global_nx, grid.slabLayout().world_size, rank);
      const std::size_t count = owned.extentX() * grid.slabLayout().global_ny * grid.slabLayout().global_nz;
      if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("isolated PM distributed slab count exceeds MPI int limit");
      }
      if (disp > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("isolated PM distributed slab displacement exceeds MPI int limit");
      }
      recv_counts[static_cast<std::size_t>(rank)] = static_cast<int>(count);
      recv_displs[static_cast<std::size_t>(rank)] = static_cast<int>(disp);
      disp += count;
    }
    if (is_root) {
      global_density.assign(global_size, 0.0);
    }
    MPI_Gatherv(
        grid.density().data(),
        static_cast<int>(grid.density().size()),
        MPI_DOUBLE,
        is_root ? global_density.data() : nullptr,
        recv_counts.data(),
        recv_displs.data(),
        MPI_DOUBLE,
        k_root,
        MPI_COMM_WORLD);
  }
#endif

  const std::size_t pad_nx = 2U * m_shape.nx;
  const std::size_t pad_ny = 2U * m_shape.ny;
  const std::size_t pad_nz = 2U * m_shape.nz;
  const std::size_t pad_total =
      checkedProduct(checkedProduct(pad_nx, pad_ny, "isolated PM padded workspace"),
                     pad_nz,
                     "isolated PM padded workspace");
  std::vector<double> full_potential;
  std::vector<double> full_force_x;
  std::vector<double> full_force_y;
  std::vector<double> full_force_z;

  if (is_root) {
    m_impl->ensureIsolatedWorkspace(pad_nx, pad_ny, pad_nz);
    auto& ws = m_impl->isolatedWorkspace();
    std::fill(ws.rho_k.begin(), ws.rho_k.end(), std::complex<double>(0.0, 0.0));
    for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
      for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
        for (std::size_t iz = 0; iz < m_shape.nz; ++iz) {
          const std::size_t src = distributed_slabs ? (ix * m_shape.ny + iy) * m_shape.nz + iz : grid.linearIndex(ix, iy, iz);
          const std::size_t dst = (ix * pad_ny + iy) * pad_nz + iz;
          ws.rho_k[dst] = {(distributed_slabs ? global_density[src] : grid.density()[src]), 0.0};
        }
      }
    }
    m_impl->isolatedForward(ws.rho_k);

    if (!ws.kernel_ready || ws.dx != dx || ws.dy != dy || ws.dz != dz || ws.split_scale != options.tree_pm_split_scale_comoving) {
      ws.dx = dx;
      ws.dy = dy;
      ws.dz = dz;
      ws.split_scale = options.tree_pm_split_scale_comoving;
      std::fill(ws.kernel_k.begin(), ws.kernel_k.end(), std::complex<double>(0.0, 0.0));
      for (std::size_t ix = 0; ix < pad_nx; ++ix) {
        const std::ptrdiff_t sx = ix <= m_shape.nx ? static_cast<std::ptrdiff_t>(ix)
                                                   : static_cast<std::ptrdiff_t>(ix) - static_cast<std::ptrdiff_t>(pad_nx);
        const double rx = static_cast<double>(sx) * dx;
        for (std::size_t iy = 0; iy < pad_ny; ++iy) {
          const std::ptrdiff_t sy = iy <= m_shape.ny ? static_cast<std::ptrdiff_t>(iy)
                                                     : static_cast<std::ptrdiff_t>(iy) - static_cast<std::ptrdiff_t>(pad_ny);
          const double ry = static_cast<double>(sy) * dy;
          for (std::size_t iz = 0; iz < pad_nz; ++iz) {
            const std::ptrdiff_t sz = iz <= m_shape.nz ? static_cast<std::ptrdiff_t>(iz)
                                                       : static_cast<std::ptrdiff_t>(iz) - static_cast<std::ptrdiff_t>(pad_nz);
            const double rz = static_cast<double>(sz) * dz;
            const double r = std::sqrt(rx * rx + ry * ry + rz * rz);
            double kernel = 0.0;
            if (r > 0.0) {
              kernel = -1.0 / r;
              if (options.tree_pm_split_scale_comoving > 0.0) {
                kernel = -std::erf(0.5 * r / options.tree_pm_split_scale_comoving) / r;
              }
            }
            ws.kernel_k[(ix * pad_ny + iy) * pad_nz + iz] = {kernel, 0.0};
          }
        }
      }
      m_impl->isolatedForward(ws.kernel_k);
      ws.kernel_ready = true;
    }

    for (std::size_t i = 0; i < pad_total; ++i) {
      ws.rho_k[i] *= ws.kernel_k[i];
    }
    m_impl->isolatedInverse(ws.rho_k);

    const double cell_volume = dx * dy * dz;
    const double prefactor = options.gravitational_constant_code * cell_volume;
    full_potential.assign(m_shape.cellCount(), 0.0);
    full_force_x.assign(m_shape.cellCount(), 0.0);
    full_force_y.assign(m_shape.cellCount(), 0.0);
    full_force_z.assign(m_shape.cellCount(), 0.0);
    for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
      for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
        for (std::size_t iz = 0; iz < m_shape.nz; ++iz) {
          const std::size_t physical_index = (ix * m_shape.ny + iy) * m_shape.nz + iz;
          const std::size_t padded_index = (ix * pad_ny + iy) * pad_nz + iz;
          double value = ws.rho_k[padded_index].real();
#if COSMOSIM_ENABLE_FFTW
          value /= static_cast<double>(pad_total);
#endif
          full_potential[physical_index] = prefactor * value;
        }
      }
    }

    const auto fieldPotential = [&](std::size_t ix, std::size_t iy, std::size_t iz) -> double {
      return full_potential[(ix * m_shape.ny + iy) * m_shape.nz + iz];
    };
    for (std::size_t ix = 0; ix < m_shape.nx; ++ix) {
      for (std::size_t iy = 0; iy < m_shape.ny; ++iy) {
        for (std::size_t iz = 0; iz < m_shape.nz; ++iz) {
          const std::size_t center = (ix * m_shape.ny + iy) * m_shape.nz + iz;
          const auto grad_x = [&]() {
            if (ix == 0) return (-3.0 * fieldPotential(0, iy, iz) + 4.0 * fieldPotential(1, iy, iz) - fieldPotential(2, iy, iz)) / (2.0 * dx);
            if (ix + 1U == m_shape.nx) return (3.0 * fieldPotential(m_shape.nx - 1U, iy, iz) - 4.0 * fieldPotential(m_shape.nx - 2U, iy, iz) + fieldPotential(m_shape.nx - 3U, iy, iz)) / (2.0 * dx);
            return (fieldPotential(ix + 1U, iy, iz) - fieldPotential(ix - 1U, iy, iz)) / (2.0 * dx);
          }();
          const auto grad_y = [&]() {
            if (iy == 0) return (-3.0 * fieldPotential(ix, 0, iz) + 4.0 * fieldPotential(ix, 1, iz) - fieldPotential(ix, 2, iz)) / (2.0 * dy);
            if (iy + 1U == m_shape.ny) return (3.0 * fieldPotential(ix, m_shape.ny - 1U, iz) - 4.0 * fieldPotential(ix, m_shape.ny - 2U, iz) + fieldPotential(ix, m_shape.ny - 3U, iz)) / (2.0 * dy);
            return (fieldPotential(ix, iy + 1U, iz) - fieldPotential(ix, iy - 1U, iz)) / (2.0 * dy);
          }();
          const auto grad_z = [&]() {
            if (iz == 0) return (-3.0 * fieldPotential(ix, iy, 0) + 4.0 * fieldPotential(ix, iy, 1) - fieldPotential(ix, iy, 2)) / (2.0 * dz);
            if (iz + 1U == m_shape.nz) return (3.0 * fieldPotential(ix, iy, m_shape.nz - 1U) - 4.0 * fieldPotential(ix, iy, m_shape.nz - 2U) + fieldPotential(ix, iy, m_shape.nz - 3U)) / (2.0 * dz);
            return (fieldPotential(ix, iy, iz + 1U) - fieldPotential(ix, iy, iz - 1U)) / (2.0 * dz);
          }();
          full_force_x[center] = -grad_x;
          full_force_y[center] = -grad_y;
          full_force_z[center] = -grad_z;
        }
      }
    }
  }

#if COSMOSIM_ENABLE_MPI
  if (distributed_slabs) {
    std::vector<int> recv_counts(static_cast<std::size_t>(grid.slabLayout().world_size), 0);
    std::vector<int> recv_displs(static_cast<std::size_t>(grid.slabLayout().world_size), 0);
    std::size_t disp = 0;
    for (int rank = 0; rank < grid.slabLayout().world_size; ++rank) {
      const auto owned = parallel::pmOwnedXRangeForRank(grid.slabLayout().global_nx, grid.slabLayout().world_size, rank);
      const std::size_t count = owned.extentX() * grid.slabLayout().global_ny * grid.slabLayout().global_nz;
      if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("isolated PM scatter slab count exceeds MPI int limit");
      }
      if (disp > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument("isolated PM scatter slab displacement exceeds MPI int limit");
      }
      recv_counts[static_cast<std::size_t>(rank)] = static_cast<int>(count);
      recv_displs[static_cast<std::size_t>(rank)] = static_cast<int>(disp);
      disp += count;
    }
    MPI_Scatterv(is_root ? full_potential.data() : nullptr, recv_counts.data(), recv_displs.data(), MPI_DOUBLE, grid.potential().data(), static_cast<int>(grid.potential().size()), MPI_DOUBLE, k_root, MPI_COMM_WORLD);
    MPI_Scatterv(is_root ? full_force_x.data() : nullptr, recv_counts.data(), recv_displs.data(), MPI_DOUBLE, grid.force_x().data(), static_cast<int>(grid.force_x().size()), MPI_DOUBLE, k_root, MPI_COMM_WORLD);
    MPI_Scatterv(is_root ? full_force_y.data() : nullptr, recv_counts.data(), recv_displs.data(), MPI_DOUBLE, grid.force_y().data(), static_cast<int>(grid.force_y().size()), MPI_DOUBLE, k_root, MPI_COMM_WORLD);
    MPI_Scatterv(is_root ? full_force_z.data() : nullptr, recv_counts.data(), recv_displs.data(), MPI_DOUBLE, grid.force_z().data(), static_cast<int>(grid.force_z().size()), MPI_DOUBLE, k_root, MPI_COMM_WORLD);
  } else
#endif
  {
    std::copy(full_potential.begin(), full_potential.end(), grid.potential().begin());
    std::copy(full_force_x.begin(), full_force_x.end(), grid.force_x().begin());
    std::copy(full_force_y.begin(), full_force_y.end(), grid.force_y().begin());
    std::copy(full_force_z.begin(), full_force_z.end(), grid.force_z().begin());
  }

  if (profile != nullptr) {
    const auto stop = std::chrono::steady_clock::now();
    profile->poisson_ms += std::chrono::duration<double, std::milli>(stop - poisson_start).count();
    profile->isolated_open_root_workspace_estimate_bytes =
        std::max(profile->isolated_open_root_workspace_estimate_bytes, isolated_root_workspace_estimate_bytes);
    profile->isolated_open_root_workspace_limit_bytes = options.isolated_open_root_workspace_limit_bytes;
    if (distributed_slabs) {
      const std::uint64_t local_gather_scatter_bytes =
          checkedBytesForCount(
              checkedProduct(5U, grid.localCellCount(), "isolated PM gather/scatter profile"),
              sizeof(double),
              "isolated PM gather/scatter profile");
      profile->isolated_open_gather_bytes += local_gather_scatter_bytes;
      profile->bytes_moved += local_gather_scatter_bytes;
    }
  }
}

void PmSolver::interpolateForces(
    const PmGridStorage& grid,
    const PmForceTargetView& target_view,
    const PmSolveOptions& options,
    PmProfileEvent* profile) const {
#if COSMOSIM_ENABLE_MPI
  double target_view_preflight_mpi_wait_ms = 0.0;
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
  if (mpi_world_size > 1) {
    const bool rank_local_serial_layout =
        grid.slabLayout().world_size == 1 && grid.ownsFullDomain();
    validatePmCollectiveEntryConsensus(
        PmCollectiveEntryKind::kInterpolateForceTargetView,
        rank_local_serial_layout,
        m_shape,
        options,
        target_view_preflight_mpi_wait_ms,
        static_cast<int>(target_view.output_layout));
    runPmCoordinatedPhase(
        "PmSolver::interpolateForces target-view API preflight",
        mpi_world_rank,
        mpi_world_size,
        target_view_preflight_mpi_wait_ms,
        [&]() {
          validateOptions(m_shape, options);
          if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny ||
              grid.shape().nz != m_shape.nz) {
            throw std::invalid_argument("PM solver/grid shape mismatch in interpolateForces target view");
          }
          if (!grid.slabLayout().isValid()) {
            throw std::invalid_argument(
                "PmSolver::interpolateForces target view requires a valid PM slab layout");
          }
          if (!rank_local_serial_layout &&
              (mpi_world_size != grid.slabLayout().world_size ||
               mpi_world_rank != grid.slabLayout().world_rank)) {
            throw std::invalid_argument(
                "PmSolver::interpolateForces target-view slab metadata must match MPI_COMM_WORLD");
          }
          const std::size_t local_active_count = target_view.active_particle_index.size();
          if (local_active_count != target_view.pos_x_comoving.size() ||
              local_active_count != target_view.pos_y_comoving.size() ||
              local_active_count != target_view.pos_z_comoving.size()) {
            throw std::invalid_argument(
                "PmSolver::interpolateForces active target coordinate view extents must match");
          }
          switch (target_view.output_layout) {
            case PmForceOutputLayout::kCompactActive:
              if (local_active_count != target_view.accel_x_comoving.size() ||
                  local_active_count != target_view.accel_y_comoving.size() ||
                  local_active_count != target_view.accel_z_comoving.size()) {
                throw std::invalid_argument(
                    "PmSolver::interpolateForces compact active output extents must match active count");
              }
              break;
            case PmForceOutputLayout::kIndexedGlobal:
              for (const std::uint32_t particle_index : target_view.active_particle_index) {
                if (particle_index >= target_view.accel_x_comoving.size() ||
                    particle_index >= target_view.accel_y_comoving.size() ||
                    particle_index >= target_view.accel_z_comoving.size()) {
                  throw std::out_of_range(
                      "PmSolver::interpolateForces indexed active particle output index out of range");
                }
              }
              break;
            default:
              throw std::invalid_argument("PmSolver::interpolateForces unknown PM force output layout");
          }
        });
    if (profile != nullptr) {
      profile->routed_mpi_wait_ms += target_view_preflight_mpi_wait_ms;
    }
  }
#endif
  const std::size_t active_count = target_view.active_particle_index.size();
  if (active_count != target_view.pos_x_comoving.size() || active_count != target_view.pos_y_comoving.size() ||
      active_count != target_view.pos_z_comoving.size()) {
    throw std::invalid_argument("PmSolver::interpolateForces active target coordinate view extents must match");
  }

  switch (target_view.output_layout) {
    case PmForceOutputLayout::kCompactActive: {
      if (active_count != target_view.accel_x_comoving.size() ||
          active_count != target_view.accel_y_comoving.size() ||
          active_count != target_view.accel_z_comoving.size()) {
        throw std::invalid_argument("PmSolver::interpolateForces compact active output extents must match active count");
      }
      interpolateForces(
          grid,
          target_view.pos_x_comoving,
          target_view.pos_y_comoving,
          target_view.pos_z_comoving,
          target_view.accel_x_comoving,
          target_view.accel_y_comoving,
          target_view.accel_z_comoving,
          options,
          profile);
      return;
    }
    case PmForceOutputLayout::kIndexedGlobal: {
      for (const std::uint32_t particle_index : target_view.active_particle_index) {
        if (particle_index >= target_view.accel_x_comoving.size() ||
            particle_index >= target_view.accel_y_comoving.size() ||
            particle_index >= target_view.accel_z_comoving.size()) {
          throw std::out_of_range("PmSolver::interpolateForces indexed active particle output index out of range");
        }
      }
      std::vector<double> compact_ax(active_count, 0.0);
      std::vector<double> compact_ay(active_count, 0.0);
      std::vector<double> compact_az(active_count, 0.0);
      interpolateForces(
          grid,
          target_view.pos_x_comoving,
          target_view.pos_y_comoving,
          target_view.pos_z_comoving,
          compact_ax,
          compact_ay,
          compact_az,
          options,
          profile);
      for (std::size_t active_i = 0; active_i < active_count; ++active_i) {
        const std::uint32_t particle_index = target_view.active_particle_index[active_i];
        target_view.accel_x_comoving[particle_index] = compact_ax[active_i];
        target_view.accel_y_comoving[particle_index] = compact_ay[active_i];
        target_view.accel_z_comoving[particle_index] = compact_az[active_i];
      }
      return;
    }
  }
  throw std::invalid_argument("PmSolver::interpolateForces unknown PM force output layout");
}

void PmSolver::interpolateForces(
    const PmGridStorage& grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<double> accel_x,
    std::span<double> accel_y,
    std::span<double> accel_z,
    const PmSolveOptions& options,
    PmProfileEvent* profile) const {
#if COSMOSIM_ENABLE_MPI
  double distributed_preflight_mpi_wait_ms = 0.0;
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
  if (mpi_world_size > 1) {
    const bool rank_local_serial_layout =
        grid.slabLayout().world_size == 1 && grid.ownsFullDomain();
    validatePmCollectiveEntryConsensus(
        PmCollectiveEntryKind::kInterpolateForces,
        rank_local_serial_layout,
        m_shape,
        options,
        distributed_preflight_mpi_wait_ms);
    runPmCoordinatedPhase(
        "PmSolver::interpolateForces distributed API preflight",
        mpi_world_rank,
        mpi_world_size,
        distributed_preflight_mpi_wait_ms,
        [&]() {
          validateOptions(m_shape, options);
          if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
            throw std::invalid_argument("PM solver/grid shape mismatch in interpolateForces");
          }
          if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() || pos_x.size() != accel_x.size() ||
              pos_x.size() != accel_y.size() || pos_x.size() != accel_z.size()) {
            throw std::invalid_argument("Particle coordinate/acceleration spans must match in interpolateForces");
          }
          if (!grid.slabLayout().isValid()) {
            throw std::invalid_argument("PmSolver::interpolateForces requires a valid PM slab layout");
          }
          if (!rank_local_serial_layout &&
              (mpi_world_size != grid.slabLayout().world_size ||
               mpi_world_rank != grid.slabLayout().world_rank)) {
            throw std::invalid_argument(
                "PmSolver::interpolateForces slab layout world metadata must match MPI_COMM_WORLD");
          }
          if (m_shape.nx > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
              m_shape.ny > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
              m_shape.nz > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
            throw std::invalid_argument(
                "PmSolver::interpolateForces mesh dimensions exceed fixed-width routing indices");
          }
        });
  }
#endif
  validateOptions(m_shape, options);
  if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
    throw std::invalid_argument("PM solver/grid shape mismatch in interpolateForces");
  }
  if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() || pos_x.size() != accel_x.size() ||
      pos_x.size() != accel_y.size() || pos_x.size() != accel_z.size()) {
    throw std::invalid_argument("Particle coordinate/acceleration spans must match in interpolateForces");
  }
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::interpolateForces requires a valid PM slab layout");
  }

  const bool distributed_slabs = grid.slabLayout().world_size > 1;
#if !COSMOSIM_ENABLE_MPI
  if (distributed_slabs) {
    throw std::invalid_argument("PmSolver::interpolateForces distributed slabs require COSMOSIM_ENABLE_MPI=ON");
  }
  validateSingleRankFullDomainGridContract(grid, "PmSolver::interpolateForces");
#else
  if (distributed_slabs) {
    int mpi_world_size = 1;
    int mpi_world_rank = 0;
    queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
    if (mpi_world_size != grid.slabLayout().world_size || mpi_world_rank != grid.slabLayout().world_rank) {
      throw std::invalid_argument("PmSolver::interpolateForces slab layout world metadata must match MPI_COMM_WORLD");
    }
  } else {
    validateSingleRankFullDomainGridContract(grid, "PmSolver::interpolateForces");
  }
#endif

  const auto start = std::chrono::steady_clock::now();

  const BoxLengths lengths = effectiveBoxLengths(options);
  const double inv_dx = static_cast<double>(m_shape.nx) / lengths.lx;
  const double inv_dy = static_cast<double>(m_shape.ny) / lengths.ly;
  const double inv_dz = static_cast<double>(m_shape.nz) / lengths.lz;

  if (distributed_slabs &&
      (m_shape.nx > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
       m_shape.ny > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
       m_shape.nz > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))) {
    throw std::invalid_argument("PmSolver::interpolateForces mesh dimensions exceed fixed-width routing indices");
  }
  std::optional<std::size_t> first_non_finite_target_index;
  for (std::size_t p = 0; p < pos_x.size(); ++p) {
    if (!std::isfinite(pos_x[p]) || !std::isfinite(pos_y[p]) || !std::isfinite(pos_z[p])) {
      first_non_finite_target_index = p;
      break;
    }
  }
  if (!distributed_slabs && first_non_finite_target_index.has_value()) {
    throw std::invalid_argument(
        "PmSolver::interpolateForces requires finite particle coordinates; particle_index=" +
        std::to_string(*first_non_finite_target_index));
  }

  if (!distributed_slabs) {
    for (std::size_t p = 0; p < pos_x.size(); ++p) {
      const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
      if (!periodic &&
          (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
           !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
           !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
        accel_x[p] = 0.0;
        accel_y[p] = 0.0;
        accel_z[p] = 0.0;
        continue;
      }
      const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
      const double y = (periodic ? wrapPosition(pos_y[p], lengths.ly) : pos_y[p]) * inv_dy;
      const double z = (periodic ? wrapPosition(pos_z[p], lengths.lz) : pos_z[p]) * inv_dz;

      const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
      const PmAxisStencil1d sy = makeAxisStencil(y, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(z, options.assignment_scheme);

      double gx = 0.0;
      double gy = 0.0;
      double gz = 0.0;

      for (std::size_t dx = 0; dx < sx.count; ++dx) {
        if (!periodic && (sx.offsets[dx] < 0 || sx.offsets[dx] >= static_cast<int>(m_shape.nx))) {
          continue;
        }
        const std::size_t ix = periodic ? wrapIndex(sx.offsets[dx], m_shape.nx) : static_cast<std::size_t>(sx.offsets[dx]);
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          if (!periodic && (sy.offsets[dy] < 0 || sy.offsets[dy] >= static_cast<int>(m_shape.ny))) {
            continue;
          }
          const std::size_t iy = periodic ? wrapIndex(sy.offsets[dy], m_shape.ny) : static_cast<std::size_t>(sy.offsets[dy]);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            if (!periodic && (sz.offsets[dz] < 0 || sz.offsets[dz] >= static_cast<int>(m_shape.nz))) {
              continue;
            }
            const std::size_t iz = periodic ? wrapIndex(sz.offsets[dz], m_shape.nz) : static_cast<std::size_t>(sz.offsets[dz]);
            const double weight = sx.weights[dx] * sy.weights[dy] * sz.weights[dz];
            const std::size_t index = grid.linearIndex(ix, iy, iz);
            gx += weight * grid.force_x()[index];
            gy += weight * grid.force_y()[index];
            gz += weight * grid.force_z()[index];
          }
        }
      }

      accel_x[p] = gx;
      accel_y[p] = gy;
      accel_z[p] = gz;
    }
  } else {
#if COSMOSIM_ENABLE_MPI
    const int world_size = grid.slabLayout().world_size;
    double routed_mpi_wait_ms = distributed_preflight_mpi_wait_ms;
    const int local_targets_valid = first_non_finite_target_index.has_value() ? 0 : 1;
    int all_targets_valid = 0;
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      MPI_Allreduce(&local_targets_valid, &all_targets_valid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    });
    if (all_targets_valid == 0) {
      const std::string local_detail = first_non_finite_target_index.has_value()
          ? " local_first_invalid_particle_index=" + std::to_string(*first_non_finite_target_index)
          : " local_targets_valid=true";
      throw std::invalid_argument(
          "PmSolver::interpolateForces rejected non-finite target coordinates on at least one MPI rank;" +
          local_detail);
    }
    std::uint64_t exchange_epoch = 0U;
    runPmCoordinatedPhase(
        "PmSolver::interpolateForces exchange epoch allocation",
        grid.slabLayout().world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() { exchange_epoch = m_impl->nextDistributedExchangeEpoch(); });
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      validatePmExchangeEpochConsensus(exchange_epoch, "PmSolver::interpolateForces");
    });
    constexpr std::string_view request_exchange_context = "PmSolver::interpolateForces request exchange";
    Impl::InterpolationExchangeBuffers<PmForceContributionRecord>* exchange_ptr = nullptr;
    std::vector<PmInterpolationRequestRegistryEntry> request_registry;
    std::vector<std::uint8_t> response_received;
    std::uint64_t force_halo_cache_hits = 0;
    runPmCoordinatedPhase(
        "PmSolver::interpolateForces request send preparation",
        grid.slabLayout().world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
    exchange_ptr = &m_impl->forceInterpolationExchangeBuffersForLayout(grid.slabLayout());
    auto& exchange = *exchange_ptr;
    for (auto& per_rank : exchange.send_requests_by_rank) {
      per_rank.clear();
    }
    for (auto& per_rank : exchange.send_contribs_by_rank) {
      per_rank.clear();
    }
    exchange.send_requests_flat.clear();
    exchange.recv_requests_flat.clear();
    exchange.send_request_wire.clear();
    exchange.recv_request_wire.clear();
    exchange.send_contribs_flat.clear();
    exchange.recv_contribs_flat.clear();
    exchange.send_contrib_wire.clear();
    exchange.recv_contrib_wire.clear();
    std::fill(exchange.send_counts.begin(), exchange.send_counts.end(), 0);
    std::fill(exchange.send_displs.begin(), exchange.send_displs.end(), 0);
    std::fill(exchange.recv_counts.begin(), exchange.recv_counts.end(), 0);
    std::fill(exchange.recv_displs.begin(), exchange.recv_displs.end(), 0);
    std::fill(exchange.send_counts_bytes.begin(), exchange.send_counts_bytes.end(), 0);
    std::fill(exchange.send_displs_bytes.begin(), exchange.send_displs_bytes.end(), 0);
    std::fill(exchange.recv_counts_bytes.begin(), exchange.recv_counts_bytes.end(), 0);
    std::fill(exchange.recv_displs_bytes.begin(), exchange.recv_displs_bytes.end(), 0);
    std::fill(exchange.send_contrib_counts.begin(), exchange.send_contrib_counts.end(), 0);
    std::fill(exchange.send_contrib_displs.begin(), exchange.send_contrib_displs.end(), 0);
    std::fill(exchange.recv_contrib_counts.begin(), exchange.recv_contrib_counts.end(), 0);
    std::fill(exchange.recv_contrib_displs.begin(), exchange.recv_contrib_displs.end(), 0);
    std::fill(exchange.send_contrib_counts_bytes.begin(), exchange.send_contrib_counts_bytes.end(), 0);
    std::fill(exchange.send_contrib_displs_bytes.begin(), exchange.send_contrib_displs_bytes.end(), 0);
    std::fill(exchange.recv_contrib_counts_bytes.begin(), exchange.recv_contrib_counts_bytes.end(), 0);
    std::fill(exchange.recv_contrib_displs_bytes.begin(), exchange.recv_contrib_displs_bytes.end(), 0);

    std::fill(accel_x.begin(), accel_x.end(), 0.0);
    std::fill(accel_y.begin(), accel_y.end(), 0.0);
    std::fill(accel_z.begin(), accel_z.end(), 0.0);

    std::uint32_t next_request_sequence = 0;
    const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
    for (std::size_t p = 0; p < pos_x.size(); ++p) {
      if (p > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
        throw std::invalid_argument("PmSolver::interpolateForces origin particle index exceeds routing token limit");
      }
      if (!periodic &&
          (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
           !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
           !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
        continue;
      }
      const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
      const double y = (periodic ? wrapPosition(pos_y[p], lengths.ly) : pos_y[p]) * inv_dy;
      const double z = (periodic ? wrapPosition(pos_z[p], lengths.lz) : pos_z[p]) * inv_dz;
      const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
      const PmAxisStencil1d sy = makeAxisStencil(y, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(z, options.assignment_scheme);
      for (std::size_t dx = 0; dx < sx.count; ++dx) {
        if (!periodic && (sx.offsets[dx] < 0 || sx.offsets[dx] >= static_cast<std::ptrdiff_t>(m_shape.nx))) {
          continue;
        }
        const std::size_t ix = periodic
            ? wrapIndex(sx.offsets[dx], m_shape.nx)
            : static_cast<std::size_t>(sx.offsets[dx]);
        const int destination_rank = parallel::pmOwnerRankForGlobalX(m_shape.nx, world_size, ix);
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          if (!periodic && (sy.offsets[dy] < 0 || sy.offsets[dy] >= static_cast<std::ptrdiff_t>(m_shape.ny))) {
            continue;
          }
          const std::size_t iy = periodic
              ? wrapIndex(sy.offsets[dy], m_shape.ny)
              : static_cast<std::size_t>(sy.offsets[dy]);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            if (!periodic && (sz.offsets[dz] < 0 || sz.offsets[dz] >= static_cast<std::ptrdiff_t>(m_shape.nz))) {
              continue;
            }
            const std::size_t iz = periodic
                ? wrapIndex(sz.offsets[dz], m_shape.nz)
                : static_cast<std::size_t>(sz.offsets[dz]);
            const double weight = sx.weights[dx] * sy.weights[dy] * sz.weights[dz];
            if (destination_rank == grid.slabLayout().world_rank) {
              const std::size_t local_index = grid.linearIndex(ix, iy, iz);
              if (!std::isfinite(grid.force_x()[local_index]) || !std::isfinite(grid.force_y()[local_index]) ||
                  !std::isfinite(grid.force_z()[local_index])) {
                throw std::invalid_argument(
                    "PmSolver::interpolateForces encountered non-finite owner-local mesh force data");
              }
              accel_x[p] += weight * grid.force_x()[local_index];
              accel_y[p] += weight * grid.force_y()[local_index];
              accel_z[p] += weight * grid.force_z()[local_index];
              continue;
            }

            double halo_fx = 0.0;
            double halo_fy = 0.0;
            double halo_fz = 0.0;
            if (grid.tryLoadForceFromHalo(ix, iy, iz, halo_fx, halo_fy, halo_fz)) {
              if (!std::isfinite(halo_fx) || !std::isfinite(halo_fy) || !std::isfinite(halo_fz)) {
                throw std::invalid_argument(
                    "PmSolver::interpolateForces encountered non-finite PM force halo data");
              }
              accel_x[p] += weight * halo_fx;
              accel_y[p] += weight * halo_fy;
              accel_z[p] += weight * halo_fz;
              ++force_halo_cache_hits;
              continue;
            }

            auto& batch = exchange.send_requests_by_rank[static_cast<std::size_t>(destination_rank)];
            if (next_request_sequence == std::numeric_limits<std::uint32_t>::max()) {
              throw std::invalid_argument("PmSolver::interpolateForces request sequence exceeds routing token limit");
            }
            request_registry.push_back(PmInterpolationRequestRegistryEntry{
                .origin_particle_index = static_cast<std::uint32_t>(p),
                .exchange_epoch = exchange_epoch,
                .expected_sender_rank = destination_rank,
            });
            batch.push_back(PmInterpolationRequestRecord{
                .origin_rank = static_cast<std::uint32_t>(grid.slabLayout().world_rank),
                .destination_rank = static_cast<std::uint32_t>(destination_rank),
                .request_sequence = next_request_sequence,
                .origin_particle_index = static_cast<std::uint32_t>(p),
                .global_ix = static_cast<std::uint32_t>(ix),
                .global_iy = static_cast<std::uint32_t>(iy),
                .global_iz = static_cast<std::uint32_t>(iz),
                .exchange_epoch = exchange_epoch,
                .weight = weight,
            });
            ++next_request_sequence;
          }
        }
      }
    }
    response_received.assign(request_registry.size(), 0U);

    std::size_t total_send = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      const std::size_t count = exchange.send_requests_by_rank[static_cast<std::size_t>(rank)].size();
      const int mpi_count = checkedMpiRecordCountOrDisplacement(
          count, k_pm_interpolation_request_wire_bytes, request_exchange_context, "send record count", rank);
      const int mpi_displacement = checkedMpiRecordCountOrDisplacement(
          total_send, k_pm_interpolation_request_wire_bytes, request_exchange_context, "send record displacement", rank);
      const std::size_t next_total = checkedMpiRecordTotal(
          total_send, count, k_pm_interpolation_request_wire_bytes, request_exchange_context, rank);
      exchange.send_counts[static_cast<std::size_t>(rank)] = mpi_count;
      exchange.send_displs[static_cast<std::size_t>(rank)] = mpi_displacement;
      total_send = next_total;
    }

    exchange.send_requests_flat.reserve(total_send);
    for (int rank = 0; rank < world_size; ++rank) {
      const auto& source = exchange.send_requests_by_rank[static_cast<std::size_t>(rank)];
      exchange.send_requests_flat.insert(exchange.send_requests_flat.end(), source.begin(), source.end());
    }
    encodePmWireRecords<PmInterpolationRequestRecord>(
        exchange.send_requests_flat,
        k_pm_interpolation_request_wire_bytes,
        exchange.send_request_wire,
        [](const PmInterpolationRequestRecord& record, std::span<std::uint8_t> bytes) {
          encodePmInterpolationRequest(record, PmWireRecordKind::kForceRequest, bytes);
        },
        request_exchange_context);
        });
    auto& exchange = *exchange_ptr;

    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      MPI_Alltoall(
          exchange.send_counts.data(), 1, MPI_INT, exchange.recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    });

    runPmCoordinatedPhase(
        "PmSolver::interpolateForces request receive-layout preparation",
        grid.slabLayout().world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
    std::size_t total_recv = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      const std::size_t count = checkedMpiReceivedRecordCount(
          exchange.recv_counts[static_cast<std::size_t>(rank)],
          k_pm_interpolation_request_wire_bytes,
          request_exchange_context,
          rank);
      const int mpi_displacement = checkedMpiRecordCountOrDisplacement(
          total_recv, k_pm_interpolation_request_wire_bytes, request_exchange_context, "received record displacement", rank);
      const std::size_t next_total = checkedMpiRecordTotal(
          total_recv, count, k_pm_interpolation_request_wire_bytes, request_exchange_context, rank);
      exchange.recv_displs[static_cast<std::size_t>(rank)] = mpi_displacement;
      total_recv = next_total;
    }
    exchange.recv_request_wire.resize(checkedPmWireByteCount(
        total_recv, k_pm_interpolation_request_wire_bytes, request_exchange_context));

    checkedMpiRecordLayoutToByteLayout(
        exchange.send_counts,
        exchange.send_displs,
        k_pm_interpolation_request_wire_bytes,
        exchange.send_counts_bytes,
        exchange.send_displs_bytes,
        "PmSolver::interpolateForces request send MPI_Alltoallv");
    checkedMpiRecordLayoutToByteLayout(
        exchange.recv_counts,
        exchange.recv_displs,
        k_pm_interpolation_request_wire_bytes,
        exchange.recv_counts_bytes,
        exchange.recv_displs_bytes,
        "PmSolver::interpolateForces request receive MPI_Alltoallv");
        });
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      MPI_Alltoallv(
          nonNullPmWireData(exchange.send_request_wire),
          exchange.send_counts_bytes.data(),
          exchange.send_displs_bytes.data(),
          MPI_BYTE,
          nonNullPmWireData(exchange.recv_request_wire),
          exchange.recv_counts_bytes.data(),
          exchange.recv_displs_bytes.data(),
          MPI_BYTE,
          MPI_COMM_WORLD);
    });

    std::string local_force_request_validation_error;
    try {
      decodePmWireRecords<PmInterpolationRequestRecord>(
          exchange.recv_request_wire,
          k_pm_interpolation_request_wire_bytes,
          exchange.recv_requests_flat,
          [](std::span<const std::uint8_t> bytes) {
            return decodePmInterpolationRequest(
                bytes, PmWireRecordKind::kForceRequest, "PmSolver::interpolateForces request decode");
          },
          request_exchange_context);

      for (int source_rank = 0; source_rank < world_size; ++source_rank) {
      auto& batch = exchange.send_contribs_by_rank[static_cast<std::size_t>(source_rank)];
      const int begin = exchange.recv_displs[static_cast<std::size_t>(source_rank)];
      const int count = exchange.recv_counts[static_cast<std::size_t>(source_rank)];
      std::optional<std::uint32_t> previous_sequence;
      for (int i = 0; i < count; ++i) {
        const auto& request = exchange.recv_requests_flat[static_cast<std::size_t>(begin + i)];
        if (request.origin_rank != static_cast<std::uint32_t>(source_rank)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request origin rank does not match MPI sender segment"));
        }
        if (request.destination_rank != static_cast<std::uint32_t>(grid.slabLayout().world_rank)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request destination rank does not match receiving slab rank"));
        }
        if (request.exchange_epoch != exchange_epoch) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request exchange epoch is stale/mismatched; expected=" + std::to_string(exchange_epoch)));
        }
        if (previous_sequence.has_value() && request.request_sequence <= *previous_sequence) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request sequence is duplicate or non-monotonic for sender segment"));
        }
        previous_sequence = request.request_sequence;
        if (!std::isfinite(request.weight)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request interpolation weight is not finite"));
        }
        if (request.global_ix >= m_shape.nx || request.global_iy >= m_shape.ny || request.global_iz >= m_shape.nz) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request PM index out of range"));
        }
        if (!grid.slabLayout().ownsGlobalX(request.global_ix)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request x-index is not owned by receiving slab"));
        }
        const int expected_owner = parallel::pmOwnerRankForGlobalX(m_shape.nx, world_size, request.global_ix);
        if (expected_owner != grid.slabLayout().world_rank ||
            expected_owner != static_cast<int>(request.destination_rank)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request expected mesh owner does not match routed destination"));
        }
        const std::size_t index = grid.linearIndex(request.global_ix, request.global_iy, request.global_iz);
        if (!std::isfinite(grid.force_x()[index]) || !std::isfinite(grid.force_y()[index]) ||
            !std::isfinite(grid.force_z()[index])) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "mesh force contribution is not finite"));
        }
        batch.push_back(PmForceContributionRecord{
            .source_rank = static_cast<std::uint32_t>(grid.slabLayout().world_rank),
            .origin_rank = request.origin_rank,
            .request_sequence = request.request_sequence,
            .origin_particle_index = request.origin_particle_index,
            .exchange_epoch = request.exchange_epoch,
            .accel_x = request.weight * grid.force_x()[index],
            .accel_y = request.weight * grid.force_y()[index],
            .accel_z = request.weight * grid.force_z()[index],
        });
      }
      }
    } catch (const std::exception& error) {
      local_force_request_validation_error = error.what();
    } catch (...) {
      local_force_request_validation_error = "non-standard exception while validating force requests";
    }
    throwIfPmPayloadValidationFailed(
        "PmSolver::interpolateForces request receive",
        grid.slabLayout().world_rank,
        world_size,
        local_force_request_validation_error,
        routed_mpi_wait_ms);

    constexpr std::string_view response_exchange_context = "PmSolver::interpolateForces response exchange";
    runPmCoordinatedPhase(
        "PmSolver::interpolateForces response send preparation",
        grid.slabLayout().world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
    std::size_t total_send_contribs = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      const std::size_t count = exchange.send_contribs_by_rank[static_cast<std::size_t>(rank)].size();
      const int mpi_count = checkedMpiRecordCountOrDisplacement(
          count, k_pm_force_response_wire_bytes, response_exchange_context, "send record count", rank);
      const int mpi_displacement = checkedMpiRecordCountOrDisplacement(
          total_send_contribs, k_pm_force_response_wire_bytes, response_exchange_context, "send record displacement", rank);
      const std::size_t next_total = checkedMpiRecordTotal(
          total_send_contribs, count, k_pm_force_response_wire_bytes, response_exchange_context, rank);
      exchange.send_contrib_counts[static_cast<std::size_t>(rank)] = mpi_count;
      exchange.send_contrib_displs[static_cast<std::size_t>(rank)] = mpi_displacement;
      total_send_contribs = next_total;
    }

    exchange.send_contribs_flat.reserve(total_send_contribs);
    for (int rank = 0; rank < world_size; ++rank) {
      const auto& source = exchange.send_contribs_by_rank[static_cast<std::size_t>(rank)];
      exchange.send_contribs_flat.insert(exchange.send_contribs_flat.end(), source.begin(), source.end());
    }
    encodePmWireRecords<PmForceContributionRecord>(
        exchange.send_contribs_flat,
        k_pm_force_response_wire_bytes,
        exchange.send_contrib_wire,
        [](const PmForceContributionRecord& record, std::span<std::uint8_t> bytes) {
          encodePmForceResponse(record, bytes);
        },
        response_exchange_context);
        });

    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      MPI_Alltoall(
          exchange.send_contrib_counts.data(),
          1,
          MPI_INT,
          exchange.recv_contrib_counts.data(),
          1,
          MPI_INT,
          MPI_COMM_WORLD);
    });

    runPmCoordinatedPhase(
        "PmSolver::interpolateForces response receive-layout preparation",
        grid.slabLayout().world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
    std::size_t total_recv_contribs = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      const std::size_t count = checkedMpiReceivedRecordCount(
          exchange.recv_contrib_counts[static_cast<std::size_t>(rank)],
          k_pm_force_response_wire_bytes,
          response_exchange_context,
          rank);
      const int mpi_displacement = checkedMpiRecordCountOrDisplacement(
          total_recv_contribs, k_pm_force_response_wire_bytes, response_exchange_context, "received record displacement", rank);
      const std::size_t next_total = checkedMpiRecordTotal(
          total_recv_contribs, count, k_pm_force_response_wire_bytes, response_exchange_context, rank);
      exchange.recv_contrib_displs[static_cast<std::size_t>(rank)] = mpi_displacement;
      total_recv_contribs = next_total;
    }
    exchange.recv_contrib_wire.resize(checkedPmWireByteCount(
        total_recv_contribs, k_pm_force_response_wire_bytes, response_exchange_context));

    checkedMpiRecordLayoutToByteLayout(
        exchange.send_contrib_counts,
        exchange.send_contrib_displs,
        k_pm_force_response_wire_bytes,
        exchange.send_contrib_counts_bytes,
        exchange.send_contrib_displs_bytes,
        "PmSolver::interpolateForces response send MPI_Alltoallv");
    checkedMpiRecordLayoutToByteLayout(
        exchange.recv_contrib_counts,
        exchange.recv_contrib_displs,
        k_pm_force_response_wire_bytes,
        exchange.recv_contrib_counts_bytes,
        exchange.recv_contrib_displs_bytes,
        "PmSolver::interpolateForces response receive MPI_Alltoallv");
        });
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      MPI_Alltoallv(
          nonNullPmWireData(exchange.send_contrib_wire),
          exchange.send_contrib_counts_bytes.data(),
          exchange.send_contrib_displs_bytes.data(),
          MPI_BYTE,
          nonNullPmWireData(exchange.recv_contrib_wire),
          exchange.recv_contrib_counts_bytes.data(),
          exchange.recv_contrib_displs_bytes.data(),
          MPI_BYTE,
          MPI_COMM_WORLD);
    });

    std::string local_force_response_validation_error;
    try {
      decodePmWireRecords<PmForceContributionRecord>(
          exchange.recv_contrib_wire,
          k_pm_force_response_wire_bytes,
          exchange.recv_contribs_flat,
          [](std::span<const std::uint8_t> bytes) { return decodePmForceResponse(bytes); },
          response_exchange_context);

      for (int sender_rank = 0; sender_rank < world_size; ++sender_rank) {
      const int begin = exchange.recv_contrib_displs[static_cast<std::size_t>(sender_rank)];
      const int count = exchange.recv_contrib_counts[static_cast<std::size_t>(sender_rank)];
      for (int i = 0; i < count; ++i) {
        const auto& contribution = exchange.recv_contribs_flat[static_cast<std::size_t>(begin + i)];
        if (contribution.source_rank != static_cast<std::uint32_t>(sender_rank)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response source rank does not match MPI sender segment"));
        }
        if (contribution.origin_rank != static_cast<std::uint32_t>(grid.slabLayout().world_rank)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response origin rank does not match receiver rank"));
        }
        if (contribution.exchange_epoch != exchange_epoch) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response exchange epoch is stale/mismatched; expected=" + std::to_string(exchange_epoch)));
        }
        if (contribution.request_sequence >= request_registry.size()) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response request sequence is out of range"));
        }
        const auto& expected_request = request_registry[contribution.request_sequence];
        if (expected_request.exchange_epoch != exchange_epoch) {
          throw std::logic_error(pmRoutingDiagnostic(
              "PmSolver::interpolateForces origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              expected_request.exchange_epoch,
              contribution.request_sequence,
              expected_request.origin_particle_index,
              "origin request registry exchange epoch mismatches current exchange"));
        }
        if (sender_rank != expected_request.expected_sender_rank) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response sender rank mismatches request destination slab; expected_sender_rank=" +
                  std::to_string(expected_request.expected_sender_rank)));
        }
        if (contribution.origin_particle_index != expected_request.origin_particle_index) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response origin slot mismatches request sequence"));
        }
        if (contribution.origin_particle_index >= pos_x.size()) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response origin particle index is out of range"));
        }
        if (!std::isfinite(contribution.accel_x) || !std::isfinite(contribution.accel_y) ||
            !std::isfinite(contribution.accel_z)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response acceleration contribution is not finite"));
        }
        if (response_received[contribution.request_sequence] != 0U) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "duplicate response contribution for issued request"));
        }
        response_received[contribution.request_sequence] = 1U;
        const std::size_t p = static_cast<std::size_t>(contribution.origin_particle_index);
        accel_x[p] += contribution.accel_x;
        accel_y[p] += contribution.accel_y;
        accel_z[p] += contribution.accel_z;
      }
      }

      for (std::size_t request_index = 0; request_index < request_registry.size(); ++request_index) {
        if (response_received[request_index] == 0U) {
          const auto& expected_request = request_registry[request_index];
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolateForces origin accumulation",
              grid.slabLayout().world_rank,
              expected_request.expected_sender_rank,
              expected_request.exchange_epoch,
              static_cast<std::uint32_t>(request_index),
              expected_request.origin_particle_index,
              "missing response contribution for issued request"));
        }
      }
    } catch (const std::exception& error) {
      local_force_response_validation_error = error.what();
    } catch (...) {
      local_force_response_validation_error = "non-standard exception while validating force responses";
    }
    throwIfPmPayloadValidationFailed(
        "PmSolver::interpolateForces response receive",
        grid.slabLayout().world_rank,
        world_size,
        local_force_response_validation_error,
        routed_mpi_wait_ms);
    if (profile != nullptr) {
      profile->routed_force_requests += static_cast<std::uint64_t>(request_registry.size());
      profile->routed_force_peer_count += static_cast<std::uint64_t>(std::count_if(
          exchange.send_counts.begin(),
          exchange.send_counts.end(),
          [](int count) { return count > 0; }));
      profile->routed_mpi_bytes_sent += static_cast<std::uint64_t>(
          exchange.send_request_wire.size() + exchange.send_contrib_wire.size());
      profile->routed_mpi_bytes_received += static_cast<std::uint64_t>(
          exchange.recv_request_wire.size() + exchange.recv_contrib_wire.size());
      profile->routed_mpi_wait_ms += routed_mpi_wait_ms;
      profile->force_halo_cache_hits += force_halo_cache_hits;
      profile->bytes_moved += checkedAddBytes(
          checkedBytesForCount(
              exchange.send_requests_flat.size() + exchange.recv_requests_flat.size(),
              k_pm_interpolation_request_wire_bytes,
              "PM force routed request profile"),
          checkedBytesForCount(
              exchange.send_contribs_flat.size() + exchange.recv_contribs_flat.size(),
              k_pm_force_response_wire_bytes,
              "PM force routed response profile"),
          "PM force routed exchange profile");
    }
#endif
  }

  const auto stop = std::chrono::steady_clock::now();
  if (profile != nullptr) {
    profile->interpolate_ms += std::chrono::duration<double, std::milli>(stop - start).count();
    profile->bytes_moved += bytesForParticles(pos_x.size());
  }
}

void PmSolver::interpolatePotential(
    const PmGridStorage& grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<double> potential,
    const PmSolveOptions& options,
    PmProfileEvent* profile) const {
#if COSMOSIM_ENABLE_MPI
  double distributed_preflight_mpi_wait_ms = 0.0;
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
  if (mpi_world_size > 1) {
    const bool rank_local_serial_layout =
        grid.slabLayout().world_size == 1 && grid.ownsFullDomain();
    validatePmCollectiveEntryConsensus(
        PmCollectiveEntryKind::kInterpolatePotential,
        rank_local_serial_layout,
        m_shape,
        options,
        distributed_preflight_mpi_wait_ms);
    runPmCoordinatedPhase(
        "PmSolver::interpolatePotential distributed API preflight",
        mpi_world_rank,
        mpi_world_size,
        distributed_preflight_mpi_wait_ms,
        [&]() {
          validateOptions(m_shape, options);
          if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
            throw std::invalid_argument("PM solver/grid shape mismatch in interpolatePotential");
          }
          if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() || pos_x.size() != potential.size()) {
            throw std::invalid_argument("Particle coordinate/potential spans must match in interpolatePotential");
          }
          if (!grid.slabLayout().isValid()) {
            throw std::invalid_argument("PmSolver::interpolatePotential requires a valid PM slab layout");
          }
          if (!rank_local_serial_layout &&
              (mpi_world_size != grid.slabLayout().world_size ||
               mpi_world_rank != grid.slabLayout().world_rank)) {
            throw std::invalid_argument(
                "PmSolver::interpolatePotential slab layout world metadata must match MPI_COMM_WORLD");
          }
          if (m_shape.nx > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
              m_shape.ny > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
              m_shape.nz > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
            throw std::invalid_argument(
                "PmSolver::interpolatePotential mesh dimensions exceed fixed-width routing indices");
          }
        });
  }
#endif
  validateOptions(m_shape, options);
  if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny || grid.shape().nz != m_shape.nz) {
    throw std::invalid_argument("PM solver/grid shape mismatch in interpolatePotential");
  }
  if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() || pos_x.size() != potential.size()) {
    throw std::invalid_argument("Particle coordinate/potential spans must match in interpolatePotential");
  }
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::interpolatePotential requires a valid PM slab layout");
  }

  const bool distributed_slabs = grid.slabLayout().world_size > 1;
#if !COSMOSIM_ENABLE_MPI
  if (distributed_slabs) {
    throw std::invalid_argument("PmSolver::interpolatePotential distributed slabs require COSMOSIM_ENABLE_MPI=ON");
  }
  validateSingleRankFullDomainGridContract(grid, "PmSolver::interpolatePotential");
#else
  if (distributed_slabs) {
    int mpi_world_size = 1;
    int mpi_world_rank = 0;
    queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
    if (mpi_world_size != grid.slabLayout().world_size || mpi_world_rank != grid.slabLayout().world_rank) {
      throw std::invalid_argument(
          "PmSolver::interpolatePotential slab layout world metadata must match MPI_COMM_WORLD");
    }
  } else {
    validateSingleRankFullDomainGridContract(grid, "PmSolver::interpolatePotential");
  }
#endif

  const auto start = std::chrono::steady_clock::now();
  const BoxLengths lengths = effectiveBoxLengths(options);
  const double inv_dx = static_cast<double>(m_shape.nx) / lengths.lx;
  const double inv_dy = static_cast<double>(m_shape.ny) / lengths.ly;
  const double inv_dz = static_cast<double>(m_shape.nz) / lengths.lz;
  if (distributed_slabs &&
      (m_shape.nx > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
       m_shape.ny > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()) ||
       m_shape.nz > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max()))) {
    throw std::invalid_argument("PmSolver::interpolatePotential mesh dimensions exceed fixed-width routing indices");
  }
  std::optional<std::size_t> first_non_finite_target_index;
  for (std::size_t p = 0; p < pos_x.size(); ++p) {
    if (!std::isfinite(pos_x[p]) || !std::isfinite(pos_y[p]) || !std::isfinite(pos_z[p])) {
      first_non_finite_target_index = p;
      break;
    }
  }
  if (!distributed_slabs && first_non_finite_target_index.has_value()) {
    throw std::invalid_argument(
        "PmSolver::interpolatePotential requires finite particle coordinates; particle_index=" +
        std::to_string(*first_non_finite_target_index));
  }

  if (!distributed_slabs) {
    for (std::size_t p = 0; p < pos_x.size(); ++p) {
      const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
      if (!periodic &&
          (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
           !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
           !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
        potential[p] = 0.0;
        continue;
      }
      const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
      const double y = (periodic ? wrapPosition(pos_y[p], lengths.ly) : pos_y[p]) * inv_dy;
      const double z = (periodic ? wrapPosition(pos_z[p], lengths.lz) : pos_z[p]) * inv_dz;

      const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
      const PmAxisStencil1d sy = makeAxisStencil(y, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(z, options.assignment_scheme);

      double phi = 0.0;
      for (std::size_t dx = 0; dx < sx.count; ++dx) {
        if (!periodic && (sx.offsets[dx] < 0 || sx.offsets[dx] >= static_cast<std::ptrdiff_t>(m_shape.nx))) {
          continue;
        }
        const std::size_t ix = periodic
            ? wrapIndex(sx.offsets[dx], m_shape.nx)
            : static_cast<std::size_t>(sx.offsets[dx]);
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          if (!periodic && (sy.offsets[dy] < 0 || sy.offsets[dy] >= static_cast<std::ptrdiff_t>(m_shape.ny))) {
            continue;
          }
          const std::size_t iy = periodic
              ? wrapIndex(sy.offsets[dy], m_shape.ny)
              : static_cast<std::size_t>(sy.offsets[dy]);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            if (!periodic && (sz.offsets[dz] < 0 || sz.offsets[dz] >= static_cast<std::ptrdiff_t>(m_shape.nz))) {
              continue;
            }
            const std::size_t iz = periodic
                ? wrapIndex(sz.offsets[dz], m_shape.nz)
                : static_cast<std::size_t>(sz.offsets[dz]);
            const double weight = sx.weights[dx] * sy.weights[dy] * sz.weights[dz];
            phi += weight * grid.potential()[grid.linearIndex(ix, iy, iz)];
          }
        }
      }
      potential[p] = phi;
    }
  } else {
#if COSMOSIM_ENABLE_MPI
    const int world_size = grid.slabLayout().world_size;
    double routed_mpi_wait_ms = distributed_preflight_mpi_wait_ms;
    const int local_targets_valid = first_non_finite_target_index.has_value() ? 0 : 1;
    int all_targets_valid = 0;
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      MPI_Allreduce(&local_targets_valid, &all_targets_valid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    });
    if (all_targets_valid == 0) {
      const std::string local_detail = first_non_finite_target_index.has_value()
          ? " local_first_invalid_particle_index=" + std::to_string(*first_non_finite_target_index)
          : " local_targets_valid=true";
      throw std::invalid_argument(
          "PmSolver::interpolatePotential rejected non-finite target coordinates on at least one MPI rank;" +
          local_detail);
    }
    std::uint64_t exchange_epoch = 0U;
    runPmCoordinatedPhase(
        "PmSolver::interpolatePotential exchange epoch allocation",
        grid.slabLayout().world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() { exchange_epoch = m_impl->nextDistributedExchangeEpoch(); });
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      validatePmExchangeEpochConsensus(exchange_epoch, "PmSolver::interpolatePotential");
    });
    constexpr std::string_view request_exchange_context = "PmSolver::interpolatePotential request exchange";
    Impl::InterpolationExchangeBuffers<PmPotentialContributionRecord>* exchange_ptr = nullptr;
    std::vector<PmInterpolationRequestRegistryEntry> request_registry;
    std::vector<std::uint8_t> response_received;
    runPmCoordinatedPhase(
        "PmSolver::interpolatePotential request send preparation",
        grid.slabLayout().world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
    exchange_ptr = &m_impl->potentialInterpolationExchangeBuffersForLayout(grid.slabLayout());
    auto& exchange = *exchange_ptr;
    for (auto& per_rank : exchange.send_requests_by_rank) {
      per_rank.clear();
    }
    for (auto& per_rank : exchange.send_contribs_by_rank) {
      per_rank.clear();
    }
    exchange.send_requests_flat.clear();
    exchange.recv_requests_flat.clear();
    exchange.send_request_wire.clear();
    exchange.recv_request_wire.clear();
    exchange.send_contribs_flat.clear();
    exchange.recv_contribs_flat.clear();
    exchange.send_contrib_wire.clear();
    exchange.recv_contrib_wire.clear();
    std::fill(exchange.send_counts.begin(), exchange.send_counts.end(), 0);
    std::fill(exchange.send_displs.begin(), exchange.send_displs.end(), 0);
    std::fill(exchange.recv_counts.begin(), exchange.recv_counts.end(), 0);
    std::fill(exchange.recv_displs.begin(), exchange.recv_displs.end(), 0);
    std::fill(exchange.send_counts_bytes.begin(), exchange.send_counts_bytes.end(), 0);
    std::fill(exchange.send_displs_bytes.begin(), exchange.send_displs_bytes.end(), 0);
    std::fill(exchange.recv_counts_bytes.begin(), exchange.recv_counts_bytes.end(), 0);
    std::fill(exchange.recv_displs_bytes.begin(), exchange.recv_displs_bytes.end(), 0);
    std::fill(exchange.send_contrib_counts.begin(), exchange.send_contrib_counts.end(), 0);
    std::fill(exchange.send_contrib_displs.begin(), exchange.send_contrib_displs.end(), 0);
    std::fill(exchange.recv_contrib_counts.begin(), exchange.recv_contrib_counts.end(), 0);
    std::fill(exchange.recv_contrib_displs.begin(), exchange.recv_contrib_displs.end(), 0);
    std::fill(exchange.send_contrib_counts_bytes.begin(), exchange.send_contrib_counts_bytes.end(), 0);
    std::fill(exchange.send_contrib_displs_bytes.begin(), exchange.send_contrib_displs_bytes.end(), 0);
    std::fill(exchange.recv_contrib_counts_bytes.begin(), exchange.recv_contrib_counts_bytes.end(), 0);
    std::fill(exchange.recv_contrib_displs_bytes.begin(), exchange.recv_contrib_displs_bytes.end(), 0);

    std::uint32_t next_request_sequence = 0;
    const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
    for (std::size_t p = 0; p < pos_x.size(); ++p) {
      if (p > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
        throw std::invalid_argument("PmSolver::interpolatePotential origin particle index exceeds routing token limit");
      }
      if (!periodic &&
          (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
           !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
           !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
        continue;
      }
      const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
      const double y = (periodic ? wrapPosition(pos_y[p], lengths.ly) : pos_y[p]) * inv_dy;
      const double z = (periodic ? wrapPosition(pos_z[p], lengths.lz) : pos_z[p]) * inv_dz;
      const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
      const PmAxisStencil1d sy = makeAxisStencil(y, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(z, options.assignment_scheme);
      for (std::size_t dx = 0; dx < sx.count; ++dx) {
        if (!periodic && (sx.offsets[dx] < 0 || sx.offsets[dx] >= static_cast<std::ptrdiff_t>(m_shape.nx))) {
          continue;
        }
        const std::size_t ix = periodic
            ? wrapIndex(sx.offsets[dx], m_shape.nx)
            : static_cast<std::size_t>(sx.offsets[dx]);
        const int destination_rank = parallel::pmOwnerRankForGlobalX(m_shape.nx, world_size, ix);
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          if (!periodic && (sy.offsets[dy] < 0 || sy.offsets[dy] >= static_cast<std::ptrdiff_t>(m_shape.ny))) {
            continue;
          }
          const std::size_t iy = periodic
              ? wrapIndex(sy.offsets[dy], m_shape.ny)
              : static_cast<std::size_t>(sy.offsets[dy]);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            if (!periodic && (sz.offsets[dz] < 0 || sz.offsets[dz] >= static_cast<std::ptrdiff_t>(m_shape.nz))) {
              continue;
            }
            const std::size_t iz = periodic
                ? wrapIndex(sz.offsets[dz], m_shape.nz)
                : static_cast<std::size_t>(sz.offsets[dz]);
            const double weight = sx.weights[dx] * sy.weights[dy] * sz.weights[dz];
            if (next_request_sequence == std::numeric_limits<std::uint32_t>::max()) {
              throw std::invalid_argument("PmSolver::interpolatePotential request sequence exceeds routing token limit");
            }
            auto& batch = exchange.send_requests_by_rank[static_cast<std::size_t>(destination_rank)];
            request_registry.push_back(PmInterpolationRequestRegistryEntry{
                .origin_particle_index = static_cast<std::uint32_t>(p),
                .exchange_epoch = exchange_epoch,
                .expected_sender_rank = destination_rank,
            });
            batch.push_back(PmInterpolationRequestRecord{
                .origin_rank = static_cast<std::uint32_t>(grid.slabLayout().world_rank),
                .destination_rank = static_cast<std::uint32_t>(destination_rank),
                .request_sequence = next_request_sequence,
                .origin_particle_index = static_cast<std::uint32_t>(p),
                .global_ix = static_cast<std::uint32_t>(ix),
                .global_iy = static_cast<std::uint32_t>(iy),
                .global_iz = static_cast<std::uint32_t>(iz),
                .exchange_epoch = exchange_epoch,
                .weight = weight,
            });
            ++next_request_sequence;
          }
        }
      }
    }
    response_received.assign(request_registry.size(), 0U);

    std::size_t total_send = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      const std::size_t count = exchange.send_requests_by_rank[static_cast<std::size_t>(rank)].size();
      const int mpi_count = checkedMpiRecordCountOrDisplacement(
          count, k_pm_interpolation_request_wire_bytes, request_exchange_context, "send record count", rank);
      const int mpi_displacement = checkedMpiRecordCountOrDisplacement(
          total_send, k_pm_interpolation_request_wire_bytes, request_exchange_context, "send record displacement", rank);
      const std::size_t next_total = checkedMpiRecordTotal(
          total_send, count, k_pm_interpolation_request_wire_bytes, request_exchange_context, rank);
      exchange.send_counts[static_cast<std::size_t>(rank)] = mpi_count;
      exchange.send_displs[static_cast<std::size_t>(rank)] = mpi_displacement;
      total_send = next_total;
    }

    exchange.send_requests_flat.reserve(total_send);
    for (int rank = 0; rank < world_size; ++rank) {
      const auto& source = exchange.send_requests_by_rank[static_cast<std::size_t>(rank)];
      exchange.send_requests_flat.insert(exchange.send_requests_flat.end(), source.begin(), source.end());
    }
    encodePmWireRecords<PmInterpolationRequestRecord>(
        exchange.send_requests_flat,
        k_pm_interpolation_request_wire_bytes,
        exchange.send_request_wire,
        [](const PmInterpolationRequestRecord& record, std::span<std::uint8_t> bytes) {
          encodePmInterpolationRequest(record, PmWireRecordKind::kPotentialRequest, bytes);
        },
        request_exchange_context);
        });
    auto& exchange = *exchange_ptr;

    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      MPI_Alltoall(
          exchange.send_counts.data(), 1, MPI_INT, exchange.recv_counts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    });

    runPmCoordinatedPhase(
        "PmSolver::interpolatePotential request receive-layout preparation",
        grid.slabLayout().world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
    std::size_t total_recv = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      const std::size_t count = checkedMpiReceivedRecordCount(
          exchange.recv_counts[static_cast<std::size_t>(rank)],
          k_pm_interpolation_request_wire_bytes,
          request_exchange_context,
          rank);
      const int mpi_displacement = checkedMpiRecordCountOrDisplacement(
          total_recv, k_pm_interpolation_request_wire_bytes, request_exchange_context, "received record displacement", rank);
      const std::size_t next_total = checkedMpiRecordTotal(
          total_recv, count, k_pm_interpolation_request_wire_bytes, request_exchange_context, rank);
      exchange.recv_displs[static_cast<std::size_t>(rank)] = mpi_displacement;
      total_recv = next_total;
    }
    exchange.recv_request_wire.resize(checkedPmWireByteCount(
        total_recv, k_pm_interpolation_request_wire_bytes, request_exchange_context));

    checkedMpiRecordLayoutToByteLayout(
        exchange.send_counts,
        exchange.send_displs,
        k_pm_interpolation_request_wire_bytes,
        exchange.send_counts_bytes,
        exchange.send_displs_bytes,
        "PmSolver::interpolatePotential request send MPI_Alltoallv");
    checkedMpiRecordLayoutToByteLayout(
        exchange.recv_counts,
        exchange.recv_displs,
        k_pm_interpolation_request_wire_bytes,
        exchange.recv_counts_bytes,
        exchange.recv_displs_bytes,
        "PmSolver::interpolatePotential request receive MPI_Alltoallv");
        });
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      MPI_Alltoallv(
          nonNullPmWireData(exchange.send_request_wire),
          exchange.send_counts_bytes.data(),
          exchange.send_displs_bytes.data(),
          MPI_BYTE,
          nonNullPmWireData(exchange.recv_request_wire),
          exchange.recv_counts_bytes.data(),
          exchange.recv_displs_bytes.data(),
          MPI_BYTE,
          MPI_COMM_WORLD);
    });

    std::string local_potential_request_validation_error;
    try {
      decodePmWireRecords<PmInterpolationRequestRecord>(
          exchange.recv_request_wire,
          k_pm_interpolation_request_wire_bytes,
          exchange.recv_requests_flat,
          [](std::span<const std::uint8_t> bytes) {
            return decodePmInterpolationRequest(
                bytes, PmWireRecordKind::kPotentialRequest, "PmSolver::interpolatePotential request decode");
          },
          request_exchange_context);

      for (int source_rank = 0; source_rank < world_size; ++source_rank) {
      auto& batch = exchange.send_contribs_by_rank[static_cast<std::size_t>(source_rank)];
      const int begin = exchange.recv_displs[static_cast<std::size_t>(source_rank)];
      const int count = exchange.recv_counts[static_cast<std::size_t>(source_rank)];
      std::optional<std::uint32_t> previous_sequence;
      for (int i = 0; i < count; ++i) {
        const auto& request = exchange.recv_requests_flat[static_cast<std::size_t>(begin + i)];
        if (request.origin_rank != static_cast<std::uint32_t>(source_rank)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request origin rank does not match MPI sender segment"));
        }
        if (request.destination_rank != static_cast<std::uint32_t>(grid.slabLayout().world_rank)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request destination rank does not match receiving slab rank"));
        }
        if (request.exchange_epoch != exchange_epoch) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request exchange epoch is stale/mismatched; expected=" + std::to_string(exchange_epoch)));
        }
        if (previous_sequence.has_value() && request.request_sequence <= *previous_sequence) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request sequence is duplicate or non-monotonic for sender segment"));
        }
        previous_sequence = request.request_sequence;
        if (!std::isfinite(request.weight)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request interpolation weight is not finite"));
        }
        if (request.global_ix >= m_shape.nx || request.global_iy >= m_shape.ny || request.global_iz >= m_shape.nz) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request PM index out of range"));
        }
        if (!grid.slabLayout().ownsGlobalX(request.global_ix)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request x-index is not owned by receiving slab"));
        }
        const int expected_owner = parallel::pmOwnerRankForGlobalX(m_shape.nx, world_size, request.global_ix);
        if (expected_owner != grid.slabLayout().world_rank ||
            expected_owner != static_cast<int>(request.destination_rank)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "request expected mesh owner does not match routed destination"));
        }
        const std::size_t index = grid.linearIndex(request.global_ix, request.global_iy, request.global_iz);
        if (!std::isfinite(grid.potential()[index])) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential slab evaluation",
              grid.slabLayout().world_rank,
              source_rank,
              request.exchange_epoch,
              request.request_sequence,
              request.origin_particle_index,
              "mesh potential contribution is not finite"));
        }
        batch.push_back(PmPotentialContributionRecord{
            .source_rank = static_cast<std::uint32_t>(grid.slabLayout().world_rank),
            .origin_rank = request.origin_rank,
            .request_sequence = request.request_sequence,
            .origin_particle_index = request.origin_particle_index,
            .exchange_epoch = request.exchange_epoch,
            .potential = request.weight * grid.potential()[index],
        });
      }
      }
    } catch (const std::exception& error) {
      local_potential_request_validation_error = error.what();
    } catch (...) {
      local_potential_request_validation_error = "non-standard exception while validating potential requests";
    }
    throwIfPmPayloadValidationFailed(
        "PmSolver::interpolatePotential request receive",
        grid.slabLayout().world_rank,
        world_size,
        local_potential_request_validation_error,
        routed_mpi_wait_ms);

    constexpr std::string_view response_exchange_context = "PmSolver::interpolatePotential response exchange";
    runPmCoordinatedPhase(
        "PmSolver::interpolatePotential response send preparation",
        grid.slabLayout().world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
    std::size_t total_send_contribs = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      const std::size_t count = exchange.send_contribs_by_rank[static_cast<std::size_t>(rank)].size();
      const int mpi_count = checkedMpiRecordCountOrDisplacement(
          count, k_pm_potential_response_wire_bytes, response_exchange_context, "send record count", rank);
      const int mpi_displacement = checkedMpiRecordCountOrDisplacement(
          total_send_contribs, k_pm_potential_response_wire_bytes, response_exchange_context, "send record displacement", rank);
      const std::size_t next_total = checkedMpiRecordTotal(
          total_send_contribs, count, k_pm_potential_response_wire_bytes, response_exchange_context, rank);
      exchange.send_contrib_counts[static_cast<std::size_t>(rank)] = mpi_count;
      exchange.send_contrib_displs[static_cast<std::size_t>(rank)] = mpi_displacement;
      total_send_contribs = next_total;
    }

    exchange.send_contribs_flat.reserve(total_send_contribs);
    for (int rank = 0; rank < world_size; ++rank) {
      const auto& source = exchange.send_contribs_by_rank[static_cast<std::size_t>(rank)];
      exchange.send_contribs_flat.insert(exchange.send_contribs_flat.end(), source.begin(), source.end());
    }
    encodePmWireRecords<PmPotentialContributionRecord>(
        exchange.send_contribs_flat,
        k_pm_potential_response_wire_bytes,
        exchange.send_contrib_wire,
        [](const PmPotentialContributionRecord& record, std::span<std::uint8_t> bytes) {
          encodePmPotentialResponse(record, bytes);
        },
        response_exchange_context);
        });

    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      MPI_Alltoall(
          exchange.send_contrib_counts.data(),
          1,
          MPI_INT,
          exchange.recv_contrib_counts.data(),
          1,
          MPI_INT,
          MPI_COMM_WORLD);
    });

    runPmCoordinatedPhase(
        "PmSolver::interpolatePotential response receive-layout preparation",
        grid.slabLayout().world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
    std::size_t total_recv_contribs = 0;
    for (int rank = 0; rank < world_size; ++rank) {
      const std::size_t count = checkedMpiReceivedRecordCount(
          exchange.recv_contrib_counts[static_cast<std::size_t>(rank)],
          k_pm_potential_response_wire_bytes,
          response_exchange_context,
          rank);
      const int mpi_displacement = checkedMpiRecordCountOrDisplacement(
          total_recv_contribs, k_pm_potential_response_wire_bytes, response_exchange_context, "received record displacement", rank);
      const std::size_t next_total = checkedMpiRecordTotal(
          total_recv_contribs, count, k_pm_potential_response_wire_bytes, response_exchange_context, rank);
      exchange.recv_contrib_displs[static_cast<std::size_t>(rank)] = mpi_displacement;
      total_recv_contribs = next_total;
    }
    exchange.recv_contrib_wire.resize(checkedPmWireByteCount(
        total_recv_contribs, k_pm_potential_response_wire_bytes, response_exchange_context));

    checkedMpiRecordLayoutToByteLayout(
        exchange.send_contrib_counts,
        exchange.send_contrib_displs,
        k_pm_potential_response_wire_bytes,
        exchange.send_contrib_counts_bytes,
        exchange.send_contrib_displs_bytes,
        "PmSolver::interpolatePotential response send MPI_Alltoallv");
    checkedMpiRecordLayoutToByteLayout(
        exchange.recv_contrib_counts,
        exchange.recv_contrib_displs,
        k_pm_potential_response_wire_bytes,
        exchange.recv_contrib_counts_bytes,
        exchange.recv_contrib_displs_bytes,
        "PmSolver::interpolatePotential response receive MPI_Alltoallv");
        });
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      MPI_Alltoallv(
          nonNullPmWireData(exchange.send_contrib_wire),
          exchange.send_contrib_counts_bytes.data(),
          exchange.send_contrib_displs_bytes.data(),
          MPI_BYTE,
          nonNullPmWireData(exchange.recv_contrib_wire),
          exchange.recv_contrib_counts_bytes.data(),
          exchange.recv_contrib_displs_bytes.data(),
          MPI_BYTE,
          MPI_COMM_WORLD);
    });

    std::string local_potential_response_validation_error;
    try {
      decodePmWireRecords<PmPotentialContributionRecord>(
          exchange.recv_contrib_wire,
          k_pm_potential_response_wire_bytes,
          exchange.recv_contribs_flat,
          [](std::span<const std::uint8_t> bytes) { return decodePmPotentialResponse(bytes); },
          response_exchange_context);

      std::fill(potential.begin(), potential.end(), 0.0);
      for (int sender_rank = 0; sender_rank < world_size; ++sender_rank) {
      const int begin = exchange.recv_contrib_displs[static_cast<std::size_t>(sender_rank)];
      const int count = exchange.recv_contrib_counts[static_cast<std::size_t>(sender_rank)];
      for (int i = 0; i < count; ++i) {
        const auto& contribution = exchange.recv_contribs_flat[static_cast<std::size_t>(begin + i)];
        if (contribution.source_rank != static_cast<std::uint32_t>(sender_rank)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response source rank does not match MPI sender segment"));
        }
        if (contribution.origin_rank != static_cast<std::uint32_t>(grid.slabLayout().world_rank)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response origin rank does not match receiver rank"));
        }
        if (contribution.exchange_epoch != exchange_epoch) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response exchange epoch is stale/mismatched; expected=" + std::to_string(exchange_epoch)));
        }
        if (contribution.request_sequence >= request_registry.size()) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response request sequence is out of range"));
        }
        const auto& expected_request = request_registry[contribution.request_sequence];
        if (expected_request.exchange_epoch != exchange_epoch) {
          throw std::logic_error(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              expected_request.exchange_epoch,
              contribution.request_sequence,
              expected_request.origin_particle_index,
              "origin request registry exchange epoch mismatches current exchange"));
        }
        if (sender_rank != expected_request.expected_sender_rank) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response sender rank mismatches request destination slab; expected_sender_rank=" +
                  std::to_string(expected_request.expected_sender_rank)));
        }
        if (contribution.origin_particle_index != expected_request.origin_particle_index) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response origin slot mismatches request sequence"));
        }
        if (contribution.origin_particle_index >= pos_x.size()) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response origin particle index is out of range"));
        }
        if (!std::isfinite(contribution.potential)) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "response potential contribution is not finite"));
        }
        if (response_received[contribution.request_sequence] != 0U) {
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential origin accumulation",
              grid.slabLayout().world_rank,
              sender_rank,
              contribution.exchange_epoch,
              contribution.request_sequence,
              contribution.origin_particle_index,
              "duplicate response contribution for issued request"));
        }
        response_received[contribution.request_sequence] = 1U;
        potential[static_cast<std::size_t>(contribution.origin_particle_index)] += contribution.potential;
      }
      }

      for (std::size_t request_index = 0; request_index < request_registry.size(); ++request_index) {
        if (response_received[request_index] == 0U) {
          const auto& expected_request = request_registry[request_index];
          throw std::invalid_argument(pmRoutingDiagnostic(
              "PmSolver::interpolatePotential origin accumulation",
              grid.slabLayout().world_rank,
              expected_request.expected_sender_rank,
              expected_request.exchange_epoch,
              static_cast<std::uint32_t>(request_index),
              expected_request.origin_particle_index,
              "missing response contribution for issued request"));
        }
      }
    } catch (const std::exception& error) {
      local_potential_response_validation_error = error.what();
    } catch (...) {
      local_potential_response_validation_error = "non-standard exception while validating potential responses";
    }
    throwIfPmPayloadValidationFailed(
        "PmSolver::interpolatePotential response receive",
        grid.slabLayout().world_rank,
        world_size,
        local_potential_response_validation_error,
        routed_mpi_wait_ms);
    if (profile != nullptr) {
      profile->routed_potential_requests += static_cast<std::uint64_t>(request_registry.size());
      profile->routed_potential_peer_count += static_cast<std::uint64_t>(std::count_if(
          exchange.send_counts.begin(),
          exchange.send_counts.end(),
          [](int count) { return count > 0; }));
      profile->routed_mpi_bytes_sent += static_cast<std::uint64_t>(
          exchange.send_request_wire.size() + exchange.send_contrib_wire.size());
      profile->routed_mpi_bytes_received += static_cast<std::uint64_t>(
          exchange.recv_request_wire.size() + exchange.recv_contrib_wire.size());
      profile->routed_mpi_wait_ms += routed_mpi_wait_ms;
      profile->bytes_moved += checkedAddBytes(
          checkedBytesForCount(
              exchange.send_requests_flat.size() + exchange.recv_requests_flat.size(),
              k_pm_interpolation_request_wire_bytes,
              "PM potential routed request profile"),
          checkedBytesForCount(
              exchange.send_contribs_flat.size() + exchange.recv_contribs_flat.size(),
              k_pm_potential_response_wire_bytes,
              "PM potential routed response profile"),
          "PM potential routed exchange profile");
    }
#endif
  }

  const auto stop = std::chrono::steady_clock::now();
  if (profile != nullptr) {
    profile->interpolate_ms += std::chrono::duration<double, std::milli>(stop - start).count();
    profile->bytes_moved += bytesForParticles(pos_x.size());
  }
}

void PmSolver::solveForParticles(
    PmGridStorage& grid,
    std::span<const double> pos_x,
    std::span<const double> pos_y,
    std::span<const double> pos_z,
    std::span<const double> mass,
    std::span<double> accel_x,
    std::span<double> accel_y,
    std::span<double> accel_z,
    const PmSolveOptions& options,
    PmProfileEvent* profile) {
#if COSMOSIM_ENABLE_MPI
  double solve_preflight_mpi_wait_ms = 0.0;
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  queryActiveMpiWorld(mpi_world_size, mpi_world_rank);
  if (mpi_world_size > 1) {
    const bool rank_local_serial_layout =
        grid.slabLayout().world_size == 1 && grid.ownsFullDomain();
    validatePmCollectiveEntryConsensus(
        PmCollectiveEntryKind::kSolveForParticles,
        rank_local_serial_layout,
        m_shape,
        options,
        solve_preflight_mpi_wait_ms,
        static_cast<int>(options.execution_policy),
        static_cast<int>(options.boundary_condition));
    runPmCoordinatedPhase(
        "PmSolver::solveForParticles API preflight",
        mpi_world_rank,
        mpi_world_size,
        solve_preflight_mpi_wait_ms,
        [&]() {
          validateOptions(m_shape, options);
          if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny ||
              grid.shape().nz != m_shape.nz) {
            throw std::invalid_argument("PM solver/grid shape mismatch in solveForParticles");
          }
          if (!grid.slabLayout().isValid()) {
            throw std::invalid_argument("PmSolver::solveForParticles requires a valid PM slab layout");
          }
          if (!rank_local_serial_layout &&
              (mpi_world_size != grid.slabLayout().world_size ||
               mpi_world_rank != grid.slabLayout().world_rank)) {
            throw std::invalid_argument(
                "PmSolver::solveForParticles slab layout world metadata must match MPI_COMM_WORLD");
          }
          if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() ||
              pos_x.size() != mass.size() || pos_x.size() != accel_x.size() ||
              pos_x.size() != accel_y.size() || pos_x.size() != accel_z.size()) {
            throw std::invalid_argument(
                "Particle coordinate/mass/acceleration spans must match in solveForParticles");
          }
#if !COSMOSIM_ENABLE_CUDA
          if (options.execution_policy == core::ExecutionPolicy::kCuda) {
            throw std::runtime_error(
                "PM solve requested execution_policy=cuda, but this build has COSMOSIM_ENABLE_CUDA=OFF");
          }
#else
          if (options.execution_policy == core::ExecutionPolicy::kCuda &&
              !rank_local_serial_layout) {
            throw std::invalid_argument(
                "Distributed PM solve does not support execution_policy=cuda");
          }
#endif
        });
    if (profile != nullptr) {
      profile->routed_mpi_wait_ms += solve_preflight_mpi_wait_ms;
    }
  }
#endif
  validateOptions(m_shape, options);

  if (options.execution_policy == core::ExecutionPolicy::kCuda) {
#if COSMOSIM_ENABLE_CUDA
    cudaStream_t stream = nullptr;
    if (cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking) != cudaSuccess) {
      throw std::runtime_error("Failed to create CUDA stream for PM solve");
    }

    try {
      const auto copy_h2d_start = std::chrono::steady_clock::now();
      core::DeviceBufferDouble pos_x_device(pos_x.size());
      core::DeviceBufferDouble pos_y_device(pos_y.size());
      core::DeviceBufferDouble pos_z_device(pos_z.size());
      core::DeviceBufferDouble mass_device(mass.size());
      core::DeviceBufferDouble density_device(m_shape.cellCount());

      core::DeviceBufferDouble force_x_device(m_shape.cellCount());
      core::DeviceBufferDouble force_y_device(m_shape.cellCount());
      core::DeviceBufferDouble force_z_device(m_shape.cellCount());
      core::DeviceBufferDouble accel_x_device(accel_x.size());
      core::DeviceBufferDouble accel_y_device(accel_y.size());
      core::DeviceBufferDouble accel_z_device(accel_z.size());

      pos_x_device.copyFromHost(pos_x, stream);
      pos_y_device.copyFromHost(pos_y, stream);
      pos_z_device.copyFromHost(pos_z, stream);
      mass_device.copyFromHost(mass, stream);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing H2D particle copy");
      }
      const auto copy_h2d_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->transfer_h2d_ms += std::chrono::duration<double, std::milli>(copy_h2d_stop - copy_h2d_start).count();
      }

      const BoxLengths lengths = effectiveBoxLengths(options);

      const auto kernel_assign_start = std::chrono::steady_clock::now();
      pmCudaAssignDensityCic(
          PmCudaAssignLaunch{
              pos_x.size(),
              m_shape.nx,
              m_shape.ny,
              m_shape.nz,
              lengths.lx,
              lengths.ly,
              lengths.lz},
          pos_x_device.data(),
          pos_y_device.data(),
          pos_z_device.data(),
          mass_device.data(),
          density_device.data(),
          stream);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing PM assignment kernel");
      }
      const auto kernel_assign_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->device_kernel_ms +=
            std::chrono::duration<double, std::milli>(kernel_assign_stop - kernel_assign_start).count();
      }

      const auto copy_density_start = std::chrono::steady_clock::now();
      density_device.copyToHost(grid.density(), stream);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing D2H density copy");
      }
      const auto copy_density_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->transfer_d2h_ms += std::chrono::duration<double, std::milli>(copy_density_stop - copy_density_start).count();
      }
      const double cell_volume =
          (lengths.lx * lengths.ly * lengths.lz) / static_cast<double>(m_shape.cellCount());
      for (double& density_cell : grid.density()) {
        density_cell /= cell_volume;
      }

      solvePoissonPeriodic(grid, options, profile);

      const auto copy_forces_start = std::chrono::steady_clock::now();
      force_x_device.copyFromHost(grid.force_x(), stream);
      force_y_device.copyFromHost(grid.force_y(), stream);
      force_z_device.copyFromHost(grid.force_z(), stream);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing H2D force copy");
      }
      const auto copy_forces_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->transfer_h2d_ms += std::chrono::duration<double, std::milli>(copy_forces_stop - copy_forces_start).count();
      }

      const auto kernel_interp_start = std::chrono::steady_clock::now();
      pmCudaInterpolateForcesCic(
          PmCudaInterpLaunch{
              pos_x.size(),
              m_shape.nx,
              m_shape.ny,
              m_shape.nz,
              lengths.lx,
              lengths.ly,
              lengths.lz},
          pos_x_device.data(),
          pos_y_device.data(),
          pos_z_device.data(),
          force_x_device.data(),
          force_y_device.data(),
          force_z_device.data(),
          accel_x_device.data(),
          accel_y_device.data(),
          accel_z_device.data(),
          stream);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing PM interpolation kernel");
      }
      const auto kernel_interp_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->device_kernel_ms +=
            std::chrono::duration<double, std::milli>(kernel_interp_stop - kernel_interp_start).count();
      }

      const auto copy_accel_start = std::chrono::steady_clock::now();
      accel_x_device.copyToHost(accel_x, stream);
      accel_y_device.copyToHost(accel_y, stream);
      accel_z_device.copyToHost(accel_z, stream);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing D2H acceleration copy");
      }
      const auto copy_accel_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->transfer_d2h_ms += std::chrono::duration<double, std::milli>(copy_accel_stop - copy_accel_start).count();
      }
    } catch (...) {
      (void)cudaStreamDestroy(stream);
      throw;
    }

    (void)cudaStreamDestroy(stream);
    if (profile != nullptr) {
      profile->bytes_moved += bytesForParticles(pos_x.size()) * 2U;
      profile->bytes_moved += bytesForGridSweep(m_shape.cellCount()) * 4U;
    }
    return;
#else
    throw std::runtime_error("PM solve requested execution_policy=cuda, but this build has COSMOSIM_ENABLE_CUDA=OFF");
#endif
  }

  assignDensity(grid, pos_x, pos_y, pos_z, mass, options, profile);
  if (options.boundary_condition == PmBoundaryCondition::kPeriodic) {
    solvePoissonPeriodic(grid, options, profile);
  } else {
    solvePoissonIsolatedOpen(grid, options, profile);
  }
  interpolateForces(grid, pos_x, pos_y, pos_z, accel_x, accel_y, accel_z, options, profile);
}

bool PmSolver::fftBackendAvailable() {
#if COSMOSIM_ENABLE_FFTW
  return true;
#else
  return false;
#endif
}

bool PmSolver::cudaBackendAvailable() {
#if COSMOSIM_ENABLE_CUDA
  int device_count = 0;
  return cudaGetDeviceCount(&device_count) == cudaSuccess && device_count > 0;
#else
  return false;
#endif
}

std::string PmSolver::fftBackendName() {
#if COSMOSIM_ENABLE_FFTW
  return "fftw";
#else
  return "naive_dft";
#endif
}

std::size_t PmSolver::cachedPlanCount() const {
  return m_impl->planCount();
}

std::size_t PmSolver::planBuildCount() const {
  return m_impl->planBuildCount();
}

bool treePmSupportedByBuild() {
  return true;
}

void requireTreePmSupportOrThrow(core::GravitySolver gravity_solver) {
  static_cast<void>(gravity_solver);
}

}  // namespace cosmosim::gravity
