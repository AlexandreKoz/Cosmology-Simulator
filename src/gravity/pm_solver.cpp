#include "cosmosim/gravity/pm_solver.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <numeric>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>
#include <limits>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"
#include "cosmosim/core/device_buffer.hpp"
#include "cosmosim/gravity/tree_pm_split_kernel.hpp"
#if COSMOSIM_ENABLE_CUDA
#include <cuda_runtime.h>
#include "cosmosim/gravity/pm_cuda_kernels.hpp"
#endif

#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif
#if COSMOSIM_ENABLE_FFTW
#include <fftw3.h>
#if COSMOSIM_ENABLE_MPI
#include <fftw3-mpi.h>
#endif
#endif

namespace cosmosim::gravity {
namespace {

constexpr double k_pi = 3.141592653589793238462643383279502884;

#if COSMOSIM_ENABLE_MPI
[[nodiscard]] std::string pmMpiErrorText(int error_code) {
  char message[MPI_MAX_ERROR_STRING]{};
  int length = 0;
  if (MPI_Error_string(error_code, message, &length) != MPI_SUCCESS || length <= 0) {
    return "MPI error code " + std::to_string(error_code);
  }
  return std::string(message, static_cast<std::size_t>(length));
}

void requirePmMpiSuccess(int error_code, std::string_view context) {
  if (error_code != MPI_SUCCESS) {
    throw std::runtime_error(
        std::string(context) + " failed: " + pmMpiErrorText(error_code));
  }
}

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

struct PmDensityPlaneRecord {
  std::uint32_t origin_rank = 0;
  std::uint32_t destination_rank = 0;
  std::uint32_t record_sequence = 0;
  std::uint32_t plane_count = 0;
  std::array<std::uint32_t, 3> global_ix{};
  std::uint64_t exchange_epoch = 0;
  double y_grid = 0.0;
  double z_grid = 0.0;
  double mass_code = 0.0;
  std::array<double, 3> x_weight{};
};

struct PmPlaneInterpolationRequestRecord {
  std::uint32_t origin_rank = 0;
  std::uint32_t destination_rank = 0;
  std::uint32_t request_sequence = 0;
  std::uint32_t origin_particle_index = 0;
  std::uint32_t plane_count = 0;
  std::array<std::uint32_t, 3> global_ix{};
  std::uint64_t exchange_epoch = 0;
  double y_grid = 0.0;
  double z_grid = 0.0;
  std::array<double, 3> x_weight{};
};

[[nodiscard, maybe_unused]] std::string pmRoutingDiagnostic(
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
  kForceResponse = 3U,
  kPotentialResponse = 5U,
  kDensityPlane = 6U,
  kForcePlaneRequest = 7U,
  kPotentialPlaneRequest = 8U,
};

constexpr std::size_t k_pm_density_wire_bytes = 56U;
constexpr std::size_t k_pm_density_plane_wire_bytes = 96U;
constexpr std::size_t k_pm_plane_interpolation_request_wire_bytes = 96U;
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

[[maybe_unused]] void encodePmDensityRecord(
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

[[nodiscard, maybe_unused]] PmDensityContributionRecord decodePmDensityRecord(
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

void encodePmDensityPlaneRecord(
    const PmDensityPlaneRecord& record,
    std::span<std::uint8_t> bytes) {
  if (bytes.size() != k_pm_density_plane_wire_bytes) {
    throw std::invalid_argument("PM density-plane wire encoder requires an exact 96-byte record");
  }
  std::size_t offset = 0U;
  writePmWireHeader(bytes, offset, PmWireRecordKind::kDensityPlane);
  writePmWireU32(bytes, offset, record.origin_rank);
  writePmWireU32(bytes, offset, record.destination_rank);
  writePmWireU32(bytes, offset, record.record_sequence);
  writePmWireU32(bytes, offset, record.plane_count);
  for (const std::uint32_t ix : record.global_ix) {
    writePmWireU32(bytes, offset, ix);
  }
  writePmWireU64(bytes, offset, record.exchange_epoch);
  writePmWireDouble(bytes, offset, record.y_grid);
  writePmWireDouble(bytes, offset, record.z_grid);
  writePmWireDouble(bytes, offset, record.mass_code);
  for (const double weight : record.x_weight) {
    writePmWireDouble(bytes, offset, weight);
  }
  if (offset != bytes.size()) {
    throw std::logic_error("PM density-plane wire encoder size contract drifted");
  }
}

[[nodiscard]] PmDensityPlaneRecord decodePmDensityPlaneRecord(
    std::span<const std::uint8_t> bytes) {
  constexpr std::string_view context = "PM density-plane request";
  if (bytes.size() != k_pm_density_plane_wire_bytes) {
    throw std::runtime_error("PM density-plane wire decoder requires an exact 96-byte record");
  }
  std::size_t offset = 0U;
  readAndValidatePmWireHeader(bytes, offset, PmWireRecordKind::kDensityPlane, context);
  PmDensityPlaneRecord record;
  record.origin_rank = readPmWireU32(bytes, offset, context);
  record.destination_rank = readPmWireU32(bytes, offset, context);
  record.record_sequence = readPmWireU32(bytes, offset, context);
  record.plane_count = readPmWireU32(bytes, offset, context);
  for (std::uint32_t& ix : record.global_ix) {
    ix = readPmWireU32(bytes, offset, context);
  }
  record.exchange_epoch = readPmWireU64(bytes, offset, context);
  record.y_grid = readPmWireDouble(bytes, offset, context);
  record.z_grid = readPmWireDouble(bytes, offset, context);
  record.mass_code = readPmWireDouble(bytes, offset, context);
  for (double& weight : record.x_weight) {
    weight = readPmWireDouble(bytes, offset, context);
  }
  if (offset != bytes.size()) {
    throw std::runtime_error("PM density-plane wire record has trailing data");
  }
  return record;
}

void encodePmPlaneInterpolationRequest(
    const PmPlaneInterpolationRequestRecord& record,
    PmWireRecordKind kind,
    std::span<std::uint8_t> bytes) {
  if ((kind != PmWireRecordKind::kForcePlaneRequest &&
       kind != PmWireRecordKind::kPotentialPlaneRequest) ||
      bytes.size() != k_pm_plane_interpolation_request_wire_bytes) {
    throw std::invalid_argument(
        "PM plane-interpolation wire encoder received an invalid kind or record size");
  }
  std::size_t offset = 0U;
  writePmWireHeader(bytes, offset, kind);
  writePmWireU32(bytes, offset, record.origin_rank);
  writePmWireU32(bytes, offset, record.destination_rank);
  writePmWireU32(bytes, offset, record.request_sequence);
  writePmWireU32(bytes, offset, record.origin_particle_index);
  writePmWireU32(bytes, offset, record.plane_count);
  for (const std::uint32_t ix : record.global_ix) {
    writePmWireU32(bytes, offset, ix);
  }
  writePmWireU32(bytes, offset, 0U);
  writePmWireU64(bytes, offset, record.exchange_epoch);
  writePmWireDouble(bytes, offset, record.y_grid);
  writePmWireDouble(bytes, offset, record.z_grid);
  for (const double weight : record.x_weight) {
    writePmWireDouble(bytes, offset, weight);
  }
  if (offset != bytes.size()) {
    throw std::logic_error("PM plane-interpolation wire encoder size contract drifted");
  }
}

[[nodiscard]] PmPlaneInterpolationRequestRecord decodePmPlaneInterpolationRequest(
    std::span<const std::uint8_t> bytes,
    PmWireRecordKind expected_kind,
    std::string_view context) {
  if (expected_kind != PmWireRecordKind::kForcePlaneRequest &&
      expected_kind != PmWireRecordKind::kPotentialPlaneRequest) {
    throw std::invalid_argument("PM plane-interpolation decoder received an invalid record kind");
  }
  if (bytes.size() != k_pm_plane_interpolation_request_wire_bytes) {
    throw std::runtime_error("PM plane-interpolation wire decoder requires an exact 96-byte record");
  }
  std::size_t offset = 0U;
  readAndValidatePmWireHeader(bytes, offset, expected_kind, context);
  PmPlaneInterpolationRequestRecord record;
  record.origin_rank = readPmWireU32(bytes, offset, context);
  record.destination_rank = readPmWireU32(bytes, offset, context);
  record.request_sequence = readPmWireU32(bytes, offset, context);
  record.origin_particle_index = readPmWireU32(bytes, offset, context);
  record.plane_count = readPmWireU32(bytes, offset, context);
  for (std::uint32_t& ix : record.global_ix) {
    ix = readPmWireU32(bytes, offset, context);
  }
  const std::uint32_t reserved = readPmWireU32(bytes, offset, context);
  record.exchange_epoch = readPmWireU64(bytes, offset, context);
  record.y_grid = readPmWireDouble(bytes, offset, context);
  record.z_grid = readPmWireDouble(bytes, offset, context);
  for (double& weight : record.x_weight) {
    weight = readPmWireDouble(bytes, offset, context);
  }
  if (reserved != 0U || offset != bytes.size()) {
    throw std::runtime_error(
        "PM plane-interpolation wire record has nonzero reserved flags or trailing data");
  }
  return record;
}

[[maybe_unused]] void encodePmForceResponse(
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

[[nodiscard, maybe_unused]] PmForceContributionRecord decodePmForceResponse(
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

[[maybe_unused]] void encodePmPotentialResponse(
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

[[nodiscard, maybe_unused]] PmPotentialContributionRecord decodePmPotentialResponse(
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

[[nodiscard, maybe_unused]] std::size_t checkedPmWireByteCount(
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
  requirePmMpiSuccess(
      MPI_Allreduce(&exchange_epoch, &minimum_epoch, 1, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD),
      "PM exchange-epoch minimum MPI_Allreduce");
  requirePmMpiSuccess(
      MPI_Allreduce(&exchange_epoch, &maximum_epoch, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD),
      "PM exchange-epoch maximum MPI_Allreduce");
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
    double& accumulated_mpi_wait_ms) {
  // Public PM entry is collective over MPI_COMM_WORLD. Bind every
  // communicator-global mesh/solver control that can select a different
  // exchange, FFT plan, or kernel before any rank enters a layout-specific
  // phase. Rank-local target coordinate/output layouts are deliberately not
  // fingerprinted: empty/uneven ranks may legitimately use indexed-source or
  // compact-active target views while participating in the same PM protocol.
  // Those representation contracts are validated in the coordinated local
  // preflight below rather than promoted to false communicator invariants.
  const std::array<std::uint64_t, 20> local_fingerprint{
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
      options.routing_exchange_batch_bytes,
      options.isolated_open_root_workspace_limit_bytes,
  };
  std::array<std::uint64_t, 20> minimum_fingerprint{};
  std::array<std::uint64_t, 20> maximum_fingerprint{};
  measurePmMpiWait(accumulated_mpi_wait_ms, [&]() {
    requirePmMpiSuccess(
        MPI_Allreduce(
            local_fingerprint.data(),
            minimum_fingerprint.data(),
            static_cast<int>(local_fingerprint.size()),
            MPI_UINT64_T,
            MPI_MIN,
            MPI_COMM_WORLD),
        "PM collective fingerprint minimum MPI_Allreduce");
  });
  measurePmMpiWait(accumulated_mpi_wait_ms, [&]() {
    requirePmMpiSuccess(
        MPI_Allreduce(
            local_fingerprint.data(),
            maximum_fingerprint.data(),
            static_cast<int>(local_fingerprint.size()),
            MPI_UINT64_T,
            MPI_MAX,
            MPI_COMM_WORLD),
        "PM collective fingerprint maximum MPI_Allreduce");
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
        std::to_string(maximum_fingerprint[4]) + ")");
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
    requirePmMpiSuccess(
        MPI_Allreduce(&local_failure, &global_failure_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD),
        "PM payload-validation failure-count MPI_Allreduce");
  });
  if (global_failure_count == 0) {
    return;
  }

  const int local_failure_rank = local_failure != 0 ? world_rank : world_size;
  int first_failure_rank = world_size;
  measurePmMpiWait(accumulated_mpi_wait_ms, [&]() {
    requirePmMpiSuccess(
        MPI_Allreduce(&local_failure_rank, &first_failure_rank, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD),
        "PM payload-validation first-rank MPI_Allreduce");
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

[[nodiscard]] int checkedMpiByteCount(
    std::uint64_t byte_value,
    std::string_view context) {
  const std::uint64_t max_mpi_int = static_cast<std::uint64_t>(std::numeric_limits<int>::max());
  if (byte_value > max_mpi_int) {
    throw std::invalid_argument(
        std::string(context) + " exceeds MPI int byte-count/displacement capacity: bytes=" +
        std::to_string(byte_value));
  }
  return static_cast<int>(byte_value);
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
  switch (scheme) {
    case PmAssignmentScheme::kCic: {
      const double base = std::floor(grid_position);
      const std::ptrdiff_t i0 = static_cast<std::ptrdiff_t>(base);
      const double t = grid_position - base;
      stencil.offsets = {i0, i0 + 1, 0};
      stencil.weights = {1.0 - t, t, 0.0};
      stencil.count = 2;
      return stencil;
    }
    case PmAssignmentScheme::kTsc: {
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
  }
  throw std::invalid_argument("Unknown PM assignment scheme in makeAxisStencil");
}

struct PmXPlaneGroup {
  int destination_rank = -1;
  std::uint32_t plane_count = 0U;
  std::array<std::uint32_t, 3> global_ix{};
  std::array<double, 3> x_weight{};
};

[[nodiscard]] std::size_t makePmXPlaneGroups(
    const PmAxisStencil1d& sx,
    bool periodic,
    const PmGridShape& shape,
    const PmDecompositionView& decomposition_view,
    std::array<PmXPlaneGroup, 3>& groups) {
  groups = {};
  std::size_t group_count = 0U;
  for (std::size_t dx = 0; dx < sx.count; ++dx) {
    if (!periodic &&
        (sx.offsets[dx] < 0 || sx.offsets[dx] >= static_cast<std::ptrdiff_t>(shape.nx))) {
      continue;
    }
    const std::size_t ix = periodic
        ? wrapIndex(sx.offsets[dx], shape.nx)
        : static_cast<std::size_t>(sx.offsets[dx]);
    const int destination_rank = decomposition_view.realOwnerRank(PmGlobalCell{ix, 0U, 0U});
    std::size_t group_index = group_count;
    for (std::size_t i = 0; i < group_count; ++i) {
      if (groups[i].destination_rank == destination_rank) {
        group_index = i;
        break;
      }
    }
    if (group_index == group_count) {
      if (group_count >= groups.size()) {
        throw std::logic_error("PM x-stencil generated more destination groups than x planes");
      }
      groups[group_index].destination_rank = destination_rank;
      ++group_count;
    }
    auto& group = groups[group_index];
    if (group.plane_count >= group.global_ix.size()) {
      throw std::logic_error("PM x-plane destination group exceeded fixed stencil capacity");
    }
    const std::size_t slot = static_cast<std::size_t>(group.plane_count);
    group.global_ix[slot] = static_cast<std::uint32_t>(ix);
    group.x_weight[slot] = sx.weights[dx];
    ++group.plane_count;
  }
  return group_count;
}

void resizePmWireBufferBounded(
    std::vector<std::uint8_t>& buffer,
    std::size_t required_bytes,
    std::size_t policy_limit_bytes,
    std::string_view context) {
  if (required_bytes > policy_limit_bytes) {
    throw std::runtime_error(
        std::string(context) + " requires " + std::to_string(required_bytes) +
        " bytes, exceeding aggregate PM routing buffer policy " +
        std::to_string(policy_limit_bytes));
  }
  if (required_bytes > buffer.capacity() || buffer.capacity() > policy_limit_bytes) {
    // Grow (or shrink after a policy reduction) via a fresh exact-size vector
    // rather than resize() growth. This avoids implementation growth heuristics
    // retaining up to ~2x the logical payload and makes owned-capacity high water
    // track the configured routing policy.
    std::vector<std::uint8_t> replacement(required_bytes);
    buffer.swap(replacement);
  } else {
    buffer.resize(required_bytes);
  }
  if (buffer.capacity() > policy_limit_bytes) {
    throw std::runtime_error(
        std::string(context) + " retained vector capacity exceeds the PM routing policy");
  }
}

[[nodiscard]] std::size_t pmRoutingParticlesPerRound(
    std::uint64_t per_peer_batch_bytes,
    std::size_t wire_record_bytes) {
  if (wire_record_bytes == 0U || per_peer_batch_bytes < wire_record_bytes) {
    throw std::invalid_argument(
        "PM routing exchange batch must hold at least one fixed wire record");
  }
  const std::uint64_t records = per_peer_batch_bytes / wire_record_bytes;
  return static_cast<std::size_t>(std::min<std::uint64_t>(
      records,
      static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())));
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
  const std::size_t bytes = core::checkedSizeMultiply(
      cell_count, sizeof(double), "PM grid sweep bytes");
  return core::checkedIntegralNarrow<std::uint64_t>(bytes, "PM grid sweep bytes");
}

[[nodiscard]] std::uint64_t bytesForParticles(std::size_t particle_count) {
  const std::size_t bytes_per_particle = core::checkedSizeMultiply(
      sizeof(double), 4U, "PM particle profile bytes per particle");
  const std::size_t bytes = core::checkedSizeMultiply(
      particle_count, bytes_per_particle, "PM particle profile bytes");
  return core::checkedIntegralNarrow<std::uint64_t>(bytes, "PM particle profile bytes");
}

[[nodiscard]] std::size_t checkedProduct(std::size_t a, std::size_t b, std::string_view context) {
  return core::checkedSizeMultiply(a, b, context);
}

struct PmPlanLocalStorageLayout {
  std::size_t logical_local_complex_cells = 0U;
  std::size_t allocated_local_complex_cells = 0U;
  std::size_t real_element_count = 0U;
  std::size_t real_z_stride = 0U;
  std::size_t transposed_local_ny = 0U;
  std::size_t transposed_begin_y = 0U;
  bool is_distributed = false;
  bool spectral_transposed = false;
  bool used_backend_allocation_query = false;
};

[[nodiscard]] PmPlanLocalStorageLayout determinePmPlanLocalStorageLayout(
    PmGridShape shape,
    const parallel::PmSlabLayout& layout,
    core::PmDecompositionMode decomposition_mode,
    bool require_active_distributed_backend) {
  if (!shape.isValid()) {
    throw std::invalid_argument("PM plan storage layout requires a valid grid shape");
  }
  if (!layout.isValid() ||
      layout.global_nx != shape.nx ||
      layout.global_ny != shape.ny ||
      layout.global_nz != shape.nz) {
    throw std::invalid_argument(
        "PM plan storage layout requires a slab layout matching the PM grid shape");
  }

  PmPlanLocalStorageLayout storage;
  const std::size_t nz_complex = shape.nz / 2U + 1U;
  storage.logical_local_complex_cells = checkedProduct(
      checkedProduct(layout.local_nx(), shape.ny, "PM plan local Fourier extent"),
      nz_complex,
      "PM plan local Fourier extent");
  storage.allocated_local_complex_cells = storage.logical_local_complex_cells;
  storage.real_z_stride = shape.nz;

  if (layout.world_size <= 1) {
    storage.real_element_count = checkedProduct(
        checkedProduct(layout.local_nx(), shape.ny, "PM plan serial real extent"),
        shape.nz,
        "PM plan serial real extent");
    return storage;
  }

  storage.is_distributed = true;
  storage.real_z_stride = checkedProduct(2U, nz_complex, "PM plan distributed real stride");

#if COSMOSIM_ENABLE_FFTW && COSMOSIM_ENABLE_MPI
  int mpi_world_size = 1;
  int mpi_world_rank = 0;
  if (queryActiveMpiWorld(mpi_world_size, mpi_world_rank)) {
    if (mpi_world_size != layout.world_size || mpi_world_rank != layout.world_rank) {
      throw std::invalid_argument(
          "PM slab layout world metadata must match active MPI_COMM_WORLD for distributed FFT storage");
    }

    fftw_mpi_init();
    ptrdiff_t backend_local_nx = 0;
    ptrdiff_t backend_begin_x = 0;
    ptrdiff_t backend_local_ny = 0;
    ptrdiff_t backend_begin_y = 0;
    ptrdiff_t backend_alloc_local = 0;
    if (decomposition_mode == core::PmDecompositionMode::kPencil) {
      backend_alloc_local = fftw_mpi_local_size_3d_transposed(
          static_cast<ptrdiff_t>(shape.nx),
          static_cast<ptrdiff_t>(shape.ny),
          static_cast<ptrdiff_t>(nz_complex),
          MPI_COMM_WORLD,
          &backend_local_nx,
          &backend_begin_x,
          &backend_local_ny,
          &backend_begin_y);
    } else {
      backend_alloc_local = fftw_mpi_local_size_3d(
          static_cast<ptrdiff_t>(shape.nx),
          static_cast<ptrdiff_t>(shape.ny),
          static_cast<ptrdiff_t>(nz_complex),
          MPI_COMM_WORLD,
          &backend_local_nx,
          &backend_begin_x);
    }
    if (backend_alloc_local < 0) {
      throw std::runtime_error(
          "FFTW MPI reported a negative local allocation size for distributed PM storage");
    }
    const bool local_extent_mismatch =
        backend_local_nx != static_cast<ptrdiff_t>(layout.local_nx());
    // FFTW leaves local_0_start unspecified for a zero-width input slab. The
    // origin carries no ownership meaning in that legal case.
    const bool nonempty_origin_mismatch =
        backend_local_nx > 0 &&
        backend_begin_x != static_cast<ptrdiff_t>(layout.owned_x.begin_x);
    if (local_extent_mismatch || nonempty_origin_mismatch) {
      std::ostringstream message;
      message << "PM slab layout is incompatible with FFTW MPI ownership for this communicator: rank="
              << layout.world_rank << ", configured=[" << layout.owned_x.begin_x << ','
              << layout.owned_x.end_x << "), backend=[" << backend_begin_x << ','
              << (backend_begin_x + backend_local_nx) << ')';
      throw std::invalid_argument(message.str());
    }

    storage.allocated_local_complex_cells = std::max<std::size_t>(
        1U, static_cast<std::size_t>(backend_alloc_local));
    if (storage.allocated_local_complex_cells < storage.logical_local_complex_cells) {
      throw std::runtime_error(
          "FFTW MPI local allocation is smaller than CHUI's logical PM Fourier extent");
    }
    if (decomposition_mode == core::PmDecompositionMode::kPencil) {
      if (backend_local_ny < 0 || backend_begin_y < 0) {
        throw std::runtime_error(
            "FFTW MPI reported negative transposed PM ownership metadata");
      }
      storage.spectral_transposed = true;
      storage.transposed_local_ny = static_cast<std::size_t>(backend_local_ny);
      storage.transposed_begin_y = static_cast<std::size_t>(backend_begin_y);
    }
    storage.used_backend_allocation_query = true;
  }
#endif

  if (!storage.used_backend_allocation_query) {
    if (require_active_distributed_backend) {
      throw std::runtime_error(
          "distributed PM plan storage requires an active FFTW-MPI session");
    }
    // Conservative backend-independent fallback for preflight/compile-only
    // environments. A transposed FFT can require the larger of input x-slab
    // or output y-slab spectral storage.
    std::size_t conservative_complex = storage.logical_local_complex_cells;
    if (decomposition_mode == core::PmDecompositionMode::kPencil) {
      const parallel::PmSlabRange y_range = parallel::pmOwnedXRangeForRank(
          shape.ny, layout.world_size, layout.world_rank);
      const std::size_t transposed_complex = checkedProduct(
          checkedProduct(shape.nx, y_range.extentX(),
                         "PM plan transposed Fourier extent"),
          nz_complex,
          "PM plan transposed Fourier extent");
      conservative_complex = std::max(conservative_complex, transposed_complex);
    }
    storage.allocated_local_complex_cells =
        std::max<std::size_t>(1U, conservative_complex);
  }

  storage.real_element_count = checkedProduct(
      2U,
      storage.allocated_local_complex_cells,
      "PM plan distributed padded real extent");
  return storage;
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
  switch (options.assignment_scheme) {
    case PmAssignmentScheme::kCic:
    case PmAssignmentScheme::kTsc:
      break;
    default:
      throw std::invalid_argument("PM solve has invalid assignment_scheme");
  }
  switch (options.boundary_condition) {
    case PmBoundaryCondition::kPeriodic:
    case PmBoundaryCondition::kIsolatedOpen:
      break;
    default:
      throw std::invalid_argument("PM solve has invalid boundary_condition");
  }
  switch (options.data_residency) {
    case PmDataResidencyPolicy::kHostOnly:
    case PmDataResidencyPolicy::kPreferDevice:
      break;
    default:
      throw std::invalid_argument("PM solve has invalid data_residency");
  }
  switch (options.execution_policy) {
    case core::ExecutionPolicy::kHostSerial:
    case core::ExecutionPolicy::kCuda:
      break;
    default:
      throw std::invalid_argument("PM solve has invalid execution_policy");
  }
  switch (options.decomposition_mode) {
    case core::PmDecompositionMode::kSlab:
    case core::PmDecompositionMode::kPencil:
      break;
    default:
      throw std::invalid_argument("PM solve has invalid decomposition_mode");
  }
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
  if (options.routing_exchange_batch_bytes <
      std::max(k_pm_density_plane_wire_bytes, k_pm_plane_interpolation_request_wire_bytes)) {
    throw std::invalid_argument(
        "PM solve routing_exchange_batch_bytes is too small for one routing record");
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

PmRoutingCapacityModel modelPmRoutingCapacity(
    int world_size,
    std::uint64_t configured_per_peer_max_bytes,
    std::uint64_t aggregate_workspace_limit_bytes,
    std::uint64_t fixed_workspace_bytes,
    std::size_t wire_record_bytes) {
  if (world_size < 1) {
    throw std::invalid_argument("PM routing capacity model requires world_size >= 1");
  }
  if (wire_record_bytes == 0U) {
    throw std::invalid_argument("PM routing capacity model requires a nonzero wire record size");
  }
  const std::uint64_t record_bytes = static_cast<std::uint64_t>(wire_record_bytes);
  if (fixed_workspace_bytes > aggregate_workspace_limit_bytes) {
    throw std::invalid_argument(
        "PM routing fixed workspace exceeds aggregate workspace limit");
  }

  PmRoutingCapacityModel model{
      .configured_per_peer_max_bytes = configured_per_peer_max_bytes,
      .fixed_workspace_bytes = fixed_workspace_bytes,
  };
  if (world_size == 1) {
    model.max_simultaneous_workspace_bytes = fixed_workspace_bytes;
    return model;
  }
  if (configured_per_peer_max_bytes < record_bytes) {
    throw std::invalid_argument(
        "PM routing configured per-peer maximum cannot hold one wire record");
  }

  const std::uint64_t remote_peers = static_cast<std::uint64_t>(world_size - 1);
  if (remote_peers > std::numeric_limits<std::uint64_t>::max() / 2U) {
    throw std::overflow_error("PM routing peer-slot count overflows uint64 capacity");
  }
  const std::uint64_t simultaneous_peer_slots = 2U * remote_peers;
  const std::uint64_t payload_budget =
      aggregate_workspace_limit_bytes - fixed_workspace_bytes;
  const std::uint64_t aggregate_limited_peer_bytes =
      payload_budget / simultaneous_peer_slots;
  std::uint64_t effective_peer_bytes =
      std::min(configured_per_peer_max_bytes, aggregate_limited_peer_bytes);
  effective_peer_bytes -= effective_peer_bytes % record_bytes;
  if (effective_peer_bytes < record_bytes) {
    throw std::invalid_argument(
        "PM routing aggregate workspace cannot hold one record per active peer direction");
  }
  if (remote_peers > std::numeric_limits<std::uint64_t>::max() / effective_peer_bytes) {
    throw std::overflow_error("PM routing directional payload capacity overflows uint64");
  }
  const std::uint64_t directional_payload = remote_peers * effective_peer_bytes;
  if (directional_payload > std::numeric_limits<std::uint64_t>::max() / 2U) {
    throw std::overflow_error("PM routing combined payload capacity overflows uint64");
  }
  const std::uint64_t combined_payload = 2U * directional_payload;
  if (combined_payload > aggregate_workspace_limit_bytes - fixed_workspace_bytes) {
    throw std::logic_error("PM routing capacity model violated aggregate workspace limit");
  }

  model.effective_per_peer_payload_bytes = effective_peer_bytes;
  model.max_send_payload_bytes = directional_payload;
  model.max_receive_payload_bytes = directional_payload;
  model.max_simultaneous_workspace_bytes = fixed_workspace_bytes + combined_payload;
  return model;
}

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
      constexpr std::size_t k_hash_magic = static_cast<std::size_t>(0x9e3779b97f4a7c15ULL);
      std::size_t seed = static_cast<std::size_t>(key.world_size) * static_cast<std::size_t>(1315423911U) +
          static_cast<std::size_t>(key.world_rank);
      seed ^= key.owned_begin_x + k_hash_magic + (seed << 6U) + (seed >> 2U);
      seed ^= key.owned_end_x + k_hash_magic + (seed << 6U) + (seed >> 2U);
      seed ^= static_cast<std::size_t>(key.decomposition_mode) + k_hash_magic + (seed << 6U) + (seed >> 2U);
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
    std::uint64_t workspace_high_water_bytes = 0U;
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
    // Reused for per-peer packing offsets so bounded rounds do not allocate
    // O(world_size) scratch repeatedly. Capacity is part of the routing budget.
    std::vector<std::size_t> cursor;
  };

  struct PlaneInterpolationExchangeBuffers {
    std::uint64_t workspace_high_water_bytes = 0U;
    // Two physical wire buffers are reused across request/response phases. Requests
    // are 96-byte plane records and responses are no larger than 64 bytes.
    // Responses compact in-place, avoiding four simultaneous payload copies.
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
    std::vector<int> send_response_counts_bytes;
    std::vector<int> send_response_displs_bytes;
    std::vector<int> recv_response_counts_bytes;
    std::vector<int> recv_response_displs_bytes;
    // Shared by request packing and later response-stream validation; the two
    // uses are phase-disjoint, so one retained rank-scale vector is sufficient.
    std::vector<std::size_t> cursor;
  };


  struct IndexedTargetWorkspace {
    std::vector<double> gathered_x;
    std::vector<double> gathered_y;
    std::vector<double> gathered_z;
    std::vector<double> compact_ax;
    std::vector<double> compact_ay;
    std::vector<double> compact_az;
  };

#if COSMOSIM_ENABLE_CUDA
  // Persistent owner for the CUDA PM staging backend. The currently validated
  // path still performs FFT/Poisson/gradient on the host, but stream and device
  // storage no longer have per-call ownership. Reserved FFT/filter workspaces
  // and an opaque plan slot make the resident backend boundary explicit without
  // claiming cuFFT execution in builds where it has not been implemented.
  struct CudaBackendWorkspace {
    cudaStream_t stream = nullptr;
    core::DeviceBufferDouble pos_x;
    core::DeviceBufferDouble pos_y;
    core::DeviceBufferDouble pos_z;
    core::DeviceBufferDouble mass;
    core::DeviceBufferDouble density;
    core::DeviceBufferDouble force_x;
    core::DeviceBufferDouble force_y;
    core::DeviceBufferDouble force_z;
    core::DeviceBufferDouble accel_x;
    core::DeviceBufferDouble accel_y;
    core::DeviceBufferDouble accel_z;
    core::DeviceBufferDouble fft_workspace;
    core::DeviceBufferDouble green_filter_workspace;
    std::uintptr_t fft_plan_handle = 0U;

    ~CudaBackendWorkspace() {
      if (stream != nullptr) {
        (void)cudaStreamDestroy(stream);
      }
    }

    [[nodiscard]] core::CudaStreamView prepare(
        std::size_t particle_count, std::size_t grid_cell_count) {
      if (stream == nullptr &&
          cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking) != cudaSuccess) {
        throw std::runtime_error("Failed to create persistent CUDA PM stream");
      }
      pos_x.resize(particle_count);
      pos_y.resize(particle_count);
      pos_z.resize(particle_count);
      mass.resize(particle_count);
      density.resize(grid_cell_count);
      force_x.resize(grid_cell_count);
      force_y.resize(grid_cell_count);
      force_z.resize(grid_cell_count);
      accel_x.resize(particle_count);
      accel_y.resize(particle_count);
      accel_z.resize(particle_count);
      return core::CudaStreamView{stream};
    }
  };
#endif

  explicit Impl(PmGridShape shape) : m_shape(shape) {}

  void shutdownBackendResources() {
#if COSMOSIM_ENABLE_FFTW
#if COSMOSIM_ENABLE_MPI
    bool has_distributed_plan = false;
    for (const auto& [_, plan] : m_plan_cache) {
      has_distributed_plan = has_distributed_plan ||
          (plan.is_distributed &&
           (plan.forward_plan != nullptr || plan.inverse_plan != nullptr));
    }
    if (has_distributed_plan) {
      int world_size = 1;
      int world_rank = 0;
      if (!queryActiveMpiWorld(world_size, world_rank)) {
        throw std::runtime_error(
            "PM FFTW-MPI plan shutdown requires an active MPI session");
      }
    }
#endif
    for (auto& [_, plan] : m_plan_cache) {
      if (plan.forward_plan != nullptr) {
        fftw_destroy_plan(plan.forward_plan);
        plan.forward_plan = nullptr;
      }
      if (plan.inverse_plan != nullptr) {
        fftw_destroy_plan(plan.inverse_plan);
        plan.inverse_plan = nullptr;
      }
    }
    if (m_isolated_workspace.forward_plan != nullptr) {
      fftw_destroy_plan(m_isolated_workspace.forward_plan);
      m_isolated_workspace.forward_plan = nullptr;
    }
    if (m_isolated_workspace.inverse_plan != nullptr) {
      fftw_destroy_plan(m_isolated_workspace.inverse_plan);
      m_isolated_workspace.inverse_plan = nullptr;
    }
#endif
    m_plan_cache.clear();
    m_active_key.reset();
  }

  void shutdownBackendResourcesNoexcept() noexcept {
#if COSMOSIM_ENABLE_FFTW && COSMOSIM_ENABLE_MPI
    bool has_distributed_plan = false;
    for (const auto& [_, plan] : m_plan_cache) {
      has_distributed_plan = has_distributed_plan ||
          (plan.is_distributed &&
           (plan.forward_plan != nullptr || plan.inverse_plan != nullptr));
    }
    if (has_distributed_plan) {
      int world_size = 1;
      int world_rank = 0;
      if (!queryActiveMpiWorld(world_size, world_rank)) {
        std::fprintf(stderr,
            "CHUI MPI LIFETIME ERROR: PM FFTW-MPI plans reached destructor outside an active MPI session; explicit backend shutdown ordering was violated\n");
        // The process is already outside legal MPI teardown. Do not call FFTW
        // destroy routines for MPI plans here; process teardown reclaims them.
        for (auto& [_, plan] : m_plan_cache) {
          if (plan.is_distributed) {
            plan.forward_plan = nullptr;
            plan.inverse_plan = nullptr;
          }
        }
      }
    }
#endif
    try {
      shutdownBackendResources();
    } catch (...) {
      // Destructors must stay noexcept. Any late distributed plan was surfaced
      // above; for other backend teardown failures process exit remains the
      // final containment boundary.
      m_plan_cache.clear();
      m_active_key.reset();
    }
  }

  ~Impl() { shutdownBackendResourcesNoexcept(); }

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
      requirePmMpiSuccess(
          MPI_Allreduce(
              &local_cache_hit,
              &global_cache_hit_count,
              1,
              MPI_INT,
              MPI_SUM,
              MPI_COMM_WORLD),
          "PmSolver distributed FFT plan-cache agreement");
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
    PmPlanLocalStorageLayout storage_layout;
#if COSMOSIM_ENABLE_FFTW && COSMOSIM_ENABLE_MPI
    double plan_mpi_wait_ms = 0.0;
    if (layout.world_size > 1) {
      runPmCoordinatedPhase(
          "PmSolver distributed FFT storage-layout preparation",
          layout.world_rank,
          layout.world_size,
          plan_mpi_wait_ms,
          [&]() {
            storage_layout = determinePmPlanLocalStorageLayout(
                m_shape, layout, decomposition_mode, true);
          });
    } else
#endif
    {
      storage_layout = determinePmPlanLocalStorageLayout(
          m_shape, layout, decomposition_mode, false);
    }
    const std::size_t expected_local_complex_size =
        storage_layout.logical_local_complex_cells;
    const std::size_t allocated_local_complex_size =
        storage_layout.allocated_local_complex_cells;
    plan.real_z_stride = storage_layout.real_z_stride;


#if COSMOSIM_ENABLE_FFTW
#if COSMOSIM_ENABLE_MPI
    if (layout.world_size > 1) {
      plan.is_distributed = storage_layout.is_distributed;
      plan.spectral_transposed = storage_layout.spectral_transposed;
      plan.transposed_local_ny = storage_layout.transposed_local_ny;
      plan.transposed_begin_y = storage_layout.transposed_begin_y;
      plan.real.assign(storage_layout.real_element_count, 0.0);
      plan.fourier.assign(allocated_local_complex_size, std::complex<double>(0.0, 0.0));
      plan.potential_k.assign(allocated_local_complex_size, std::complex<double>(0.0, 0.0));
      plan.working_k.assign(allocated_local_complex_size, std::complex<double>(0.0, 0.0));
      plan.poisson_kernel.assign(allocated_local_complex_size, 0.0);
      plan.grad_kx.assign(allocated_local_complex_size, 0.0);
      plan.grad_ky.assign(allocated_local_complex_size, 0.0);
      plan.grad_kz.assign(allocated_local_complex_size, 0.0);
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
      plan.real_z_stride = storage_layout.real_z_stride;
      plan.real.assign(storage_layout.real_element_count, 0.0);
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
    plan.real.assign(storage_layout.real_element_count, 0.0);
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

    auto set_entry = [&](std::size_t index,
                         double kx,
                         double ky,
                         double kz,
                         double grad_kx,
                         double grad_ky,
                         double grad_kz) {
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
        split_filter = treePmGaussianFourierLongRangeFilterUnchecked(std::sqrt(k2), options.tree_pm_split_scale_comoving);
      }
      plan.poisson_kernel[index] = prefactor * window_correction * split_filter / k2;
      plan.grad_kx[index] = grad_kx;
      plan.grad_ky[index] = grad_ky;
      plan.grad_kz[index] = grad_kz;
    };

    const auto first_derivative_mode = [](std::size_t mode_index,
                                          std::size_t mode_count,
                                          double wave_number) {
      // In an even-length real DFT the Nyquist coefficient is self-conjugate.
      // A first derivative is odd, so assigning +/-i*k to that self-conjugate
      // coefficient would violate the Hermitian contract required by c2r.
      // The standard pseudospectral convention is therefore a zero Nyquist
      // multiplier on the differentiated axis only.
      return mode_count % 2U == 0U && mode_index == mode_count / 2U
          ? 0.0
          : wave_number;
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
            set_entry(
                index, kx, ky, kz,
                first_derivative_mode(ix, shape.nx, kx),
                first_derivative_mode(iy, shape.ny, ky),
                first_derivative_mode(iz, shape.nz, kz));
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
            set_entry(
                index, kx, ky, kz,
                first_derivative_mode(ix, shape.nx, kx),
                first_derivative_mode(iy, shape.ny, ky),
                first_derivative_mode(iz, shape.nz, kz));
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
      m_density_exchange.buffers.send_counts.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.send_displs.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.recv_counts.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.recv_displs.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.send_counts_bytes.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.send_displs_bytes.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.recv_counts_bytes.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.recv_displs_bytes.assign(static_cast<std::size_t>(layout.world_size), 0);
      m_density_exchange.buffers.cursor.assign(static_cast<std::size_t>(layout.world_size), 0U);
    }
    return m_density_exchange.buffers;
  }

  [[nodiscard]] PlaneInterpolationExchangeBuffers&
  planeInterpolationExchangeBuffersForLayout(const parallel::PmSlabLayout& layout) {
    if (m_plane_interpolation_exchange.world_size != layout.world_size || m_plane_interpolation_exchange.world_rank != layout.world_rank) {
      m_plane_interpolation_exchange.world_size = layout.world_size;
      m_plane_interpolation_exchange.world_rank = layout.world_rank;
      auto& b = m_plane_interpolation_exchange.buffers;
      const std::size_t ranks = static_cast<std::size_t>(layout.world_size);
      b.send_counts.assign(ranks, 0);
      b.send_displs.assign(ranks, 0);
      b.recv_counts.assign(ranks, 0);
      b.recv_displs.assign(ranks, 0);
      b.send_counts_bytes.assign(ranks, 0);
      b.send_displs_bytes.assign(ranks, 0);
      b.recv_counts_bytes.assign(ranks, 0);
      b.recv_displs_bytes.assign(ranks, 0);
      b.send_response_counts_bytes.assign(ranks, 0);
      b.send_response_displs_bytes.assign(ranks, 0);
      b.recv_response_counts_bytes.assign(ranks, 0);
      b.recv_response_displs_bytes.assign(ranks, 0);
      b.cursor.assign(ranks, 0U);
    }
    return m_plane_interpolation_exchange.buffers;
  }

  [[nodiscard]] IndexedTargetWorkspace& indexedTargetWorkspace() noexcept {
    return m_indexed_target_workspace;
  }
#if COSMOSIM_ENABLE_CUDA
  [[nodiscard]] CudaBackendWorkspace& cudaWorkspace() noexcept {
    return m_cuda_workspace;
  }
#endif
  void appendOwnedMemoryReport(core::MemoryReportBuilder& builder) const {
    appendMemoryReport(builder);
  }

#if COSMOSIM_ENABLE_CUDA
  // Current CUDA staging strategy. All transfer/deposit/host-field-solve/
  // interpolation choreography lives behind Impl so a future resident strategy
  // can replace these stages without changing solveForParticles().
  void solveCudaStagingBackend(
      PmSolver& solver,
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
    auto& cuda_workspace = m_cuda_workspace;
    const core::CudaStreamView stream_view = cuda_workspace.prepare(
        pos_x.size(), m_shape.cellCount());
    const cudaStream_t stream = cuda_workspace.stream;
    auto& pos_x_device = cuda_workspace.pos_x;
    auto& pos_y_device = cuda_workspace.pos_y;
    auto& pos_z_device = cuda_workspace.pos_z;
    auto& mass_device = cuda_workspace.mass;
    auto& density_device = cuda_workspace.density;
    auto& force_x_device = cuda_workspace.force_x;
    auto& force_y_device = cuda_workspace.force_y;
    auto& force_z_device = cuda_workspace.force_z;
    auto& accel_x_device = cuda_workspace.accel_x;
    auto& accel_y_device = cuda_workspace.accel_y;
    auto& accel_z_device = cuda_workspace.accel_z;

    {
      const auto copy_h2d_start = std::chrono::steady_clock::now();
      pos_x_device.copyFromHostAsync(pos_x, stream_view);
      pos_y_device.copyFromHostAsync(pos_y, stream_view);
      pos_z_device.copyFromHostAsync(pos_z, stream_view);
      mass_device.copyFromHostAsync(mass, stream_view);
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
          stream_view);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing PM assignment kernel");
      }
      const auto kernel_assign_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->device_kernel_ms +=
            std::chrono::duration<double, std::milli>(kernel_assign_stop - kernel_assign_start).count();
      }

      const auto copy_density_start = std::chrono::steady_clock::now();
      density_device.copyToHostAsync(grid.density(), stream_view);
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

      solver.solvePoissonPeriodic(grid, options, profile);

      const auto copy_forces_start = std::chrono::steady_clock::now();
      force_x_device.copyFromHostAsync(grid.force_x(), stream_view);
      force_y_device.copyFromHostAsync(grid.force_y(), stream_view);
      force_z_device.copyFromHostAsync(grid.force_z(), stream_view);
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
          stream_view);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing PM interpolation kernel");
      }
      const auto kernel_interp_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->device_kernel_ms +=
            std::chrono::duration<double, std::milli>(kernel_interp_stop - kernel_interp_start).count();
      }

      const auto copy_accel_start = std::chrono::steady_clock::now();
      accel_x_device.copyToHostAsync(accel_x, stream_view);
      accel_y_device.copyToHostAsync(accel_y, stream_view);
      accel_z_device.copyToHostAsync(accel_z, stream_view);
      if (cudaStreamSynchronize(stream) != cudaSuccess) {
        throw std::runtime_error("Failed while synchronizing D2H acceleration copy");
      }
      const auto copy_accel_stop = std::chrono::steady_clock::now();
      if (profile != nullptr) {
        profile->transfer_d2h_ms += std::chrono::duration<double, std::milli>(copy_accel_stop - copy_accel_start).count();
      }
    }

    if (profile != nullptr) {
      profile->bytes_moved += bytesForParticles(pos_x.size()) * 2U;
      profile->bytes_moved += bytesForGridSweep(m_shape.cellCount()) * 4U;
    }
    return;
  }
#endif

 private:
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

  void appendMemoryReport(core::MemoryReportBuilder& builder) const {
    const auto checked_add = [](std::uint64_t lhs, std::uint64_t rhs) {
      if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
        throw std::overflow_error("PM memory report byte sum overflow");
      }
      return lhs + rhs;
    };
    const auto accumulate = [&checked_add](
                                std::uint64_t& current,
                                std::uint64_t& capacity,
                                const auto& container) {
      current = checked_add(current, core::currentSizeBytesForContainer(container));
      capacity = checked_add(capacity, core::ownedCapacityBytesForContainer(container));
    };
    const auto emit = [&builder](
                          core::MemorySubsystem subsystem,
                          core::MemoryLifetime lifetime,
                          std::string label,
                          std::uint64_t current,
                          std::uint64_t capacity,
                          std::uint64_t high_water,
                          std::string note) {
      builder.addEntry(core::MemoryEntry{
          .subsystem = subsystem,
          .lifetime = lifetime,
          .label = std::move(label),
          .current_size_bytes = current,
          .owned_capacity_bytes = capacity,
          .high_water_bytes = std::max(capacity, high_water),
          .estimated_next_step_bytes = 0U,
          .uncertainty_note = std::move(note),
      });
    };

    {
      std::uint64_t current = 0U;
      std::uint64_t capacity = 0U;
      accumulate(current, capacity, m_indexed_target_workspace.gathered_x);
      accumulate(current, capacity, m_indexed_target_workspace.gathered_y);
      accumulate(current, capacity, m_indexed_target_workspace.gathered_z);
      accumulate(current, capacity, m_indexed_target_workspace.compact_ax);
      accumulate(current, capacity, m_indexed_target_workspace.compact_ay);
      accumulate(current, capacity, m_indexed_target_workspace.compact_az);
      emit(core::MemorySubsystem::kScratch, core::MemoryLifetime::kTransient,
           "pm_solver.indexed_target_workspace", current, capacity, capacity,
           "retained indexed-target coordinate/output scratch; capacity is its allocation high-water");
    }

    {
      std::uint64_t current = 0U;
      std::uint64_t capacity = 0U;
      for (const auto& [_, resources] : m_plan_cache) {
        accumulate(current, capacity, resources.real);
        accumulate(current, capacity, resources.fourier);
        accumulate(current, capacity, resources.potential_k);
        accumulate(current, capacity, resources.working_k);
        accumulate(current, capacity, resources.poisson_kernel);
        accumulate(current, capacity, resources.grad_kx);
        accumulate(current, capacity, resources.grad_ky);
        accumulate(current, capacity, resources.grad_kz);
      }
      emit(core::MemorySubsystem::kPmMesh, core::MemoryLifetime::kPersistent,
           "pm_solver.plan_cache_owned_arrays", current, capacity, capacity,
           "owned FFT/Poisson arrays only; FFTW/cuFFT plan-internal allocations remain external/unknown");
    }

    {
      std::uint64_t current = 0U;
      std::uint64_t capacity = 0U;
      accumulate(current, capacity, m_isolated_workspace.rho_k);
      accumulate(current, capacity, m_isolated_workspace.kernel_k);
      accumulate(current, capacity, m_isolated_workspace.scratch);
      accumulate(current, capacity, m_isolated_workspace.potential_real);
      emit(core::MemorySubsystem::kPmMesh, core::MemoryLifetime::kTransient,
           "pm_solver.isolated_open_workspace", current, capacity, capacity,
           "retained owned arrays; distributed root/backend allocations are reported separately");
    }

    {
      const auto& b = m_density_exchange.buffers;
      std::uint64_t current = 0U;
      std::uint64_t capacity = 0U;
      accumulate(current, capacity, b.send_wire);
      accumulate(current, capacity, b.recv_wire);
      accumulate(current, capacity, b.send_counts);
      accumulate(current, capacity, b.send_displs);
      accumulate(current, capacity, b.recv_counts);
      accumulate(current, capacity, b.recv_displs);
      accumulate(current, capacity, b.send_counts_bytes);
      accumulate(current, capacity, b.send_displs_bytes);
      accumulate(current, capacity, b.recv_counts_bytes);
      accumulate(current, capacity, b.recv_displs_bytes);
      accumulate(current, capacity, b.cursor);
      emit(core::MemorySubsystem::kMpiBuffers, core::MemoryLifetime::kTransient,
           "pm_solver.density_exchange_workspace", current, capacity, b.workspace_high_water_bytes,
           "bounded PM density communication workspace; logical bytes, retained capacity, and historical high-water are distinct");
    }

    {
      const auto& b = m_plane_interpolation_exchange.buffers;
      std::uint64_t current = 0U;
      std::uint64_t capacity = 0U;
      accumulate(current, capacity, b.send_wire);
      accumulate(current, capacity, b.recv_wire);
      accumulate(current, capacity, b.send_counts);
      accumulate(current, capacity, b.send_displs);
      accumulate(current, capacity, b.recv_counts);
      accumulate(current, capacity, b.recv_displs);
      accumulate(current, capacity, b.send_counts_bytes);
      accumulate(current, capacity, b.send_displs_bytes);
      accumulate(current, capacity, b.recv_counts_bytes);
      accumulate(current, capacity, b.recv_displs_bytes);
      accumulate(current, capacity, b.send_response_counts_bytes);
      accumulate(current, capacity, b.send_response_displs_bytes);
      accumulate(current, capacity, b.recv_response_counts_bytes);
      accumulate(current, capacity, b.recv_response_displs_bytes);
      accumulate(current, capacity, b.cursor);
      emit(core::MemorySubsystem::kMpiBuffers, core::MemoryLifetime::kTransient,
           "pm_solver.plane_interpolation_exchange_workspace", current, capacity, b.workspace_high_water_bytes,
           "shared bounded plane-request PM force/potential workspace; no population-scale decoded copies");
    }

#if COSMOSIM_ENABLE_CUDA
    const auto emit_device = [&builder](std::string label, const core::DeviceBufferDouble& buffer) {
      builder.addEntry(core::MemoryEntry{
          .subsystem = core::MemorySubsystem::kPmMesh,
          .lifetime = core::MemoryLifetime::kPersistent,
          .label = std::move(label),
          .current_size_bytes = static_cast<std::uint64_t>(buffer.sizeBytes()),
          .owned_capacity_bytes = static_cast<std::uint64_t>(buffer.sizeBytes()),
          .high_water_bytes = static_cast<std::uint64_t>(buffer.highWaterBytes()),
          .estimated_next_step_bytes = 0U,
          .uncertainty_note = "persistent CUDA allocation owned by PmSolver; CUDA/cuFFT runtime internals excluded",
      });
    };
    emit_device("pm_solver.cuda.pos_x", m_cuda_workspace.pos_x);
    emit_device("pm_solver.cuda.pos_y", m_cuda_workspace.pos_y);
    emit_device("pm_solver.cuda.pos_z", m_cuda_workspace.pos_z);
    emit_device("pm_solver.cuda.mass", m_cuda_workspace.mass);
    emit_device("pm_solver.cuda.density", m_cuda_workspace.density);
    emit_device("pm_solver.cuda.force_x", m_cuda_workspace.force_x);
    emit_device("pm_solver.cuda.force_y", m_cuda_workspace.force_y);
    emit_device("pm_solver.cuda.force_z", m_cuda_workspace.force_z);
    emit_device("pm_solver.cuda.accel_x", m_cuda_workspace.accel_x);
    emit_device("pm_solver.cuda.accel_y", m_cuda_workspace.accel_y);
    emit_device("pm_solver.cuda.accel_z", m_cuda_workspace.accel_z);
    emit_device("pm_solver.cuda.fft_workspace", m_cuda_workspace.fft_workspace);
    emit_device("pm_solver.cuda.green_filter_workspace", m_cuda_workspace.green_filter_workspace);
    builder.addEntry(core::MemoryEntry{
        .subsystem = core::MemorySubsystem::kPmMesh,
        .lifetime = core::MemoryLifetime::kUnknown,
        .label = "pm_solver.cuda.plan_internal_workspace",
        .estimate_only = true,
        .uncertainty_note = "CUDA/cuFFT plan/runtime allocations are backend-owned and not measurable through the current API",
    });
#endif
  }

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
    PlaneInterpolationExchangeBuffers buffers;
  } m_plane_interpolation_exchange{};
  IndexedTargetWorkspace m_indexed_target_workspace{};
  IsolatedOpenWorkspace m_isolated_workspace{};
#if COSMOSIM_ENABLE_CUDA
  CudaBackendWorkspace m_cuda_workspace{};
#endif
};

std::size_t PmGridShape::cellCount() const {
  return checkedProduct(checkedProduct(nx, ny, "PmGridShape::cellCount"), nz, "PmGridShape::cellCount");
}

bool PmGridShape::isValid() const {
  return nx > 0 && ny > 0 && nz > 0;
}

PmPlanResourcesMemoryEstimate estimatePmPlanResourcesMemory(
    PmGridShape shape,
    const parallel::PmSlabLayout& layout,
    core::PmDecompositionMode decomposition_mode) {
  const PmPlanLocalStorageLayout storage_layout =
      determinePmPlanLocalStorageLayout(
          shape, layout, decomposition_mode, false);

  const auto checked_bytes = [](std::size_t count, std::size_t element_bytes,
                                std::string_view context) -> std::uint64_t {
    const std::size_t bytes = checkedProduct(count, element_bytes, context);
    return static_cast<std::uint64_t>(bytes);
  };

  PmPlanResourcesMemoryEstimate estimate;
  estimate.logical_local_complex_cells = static_cast<std::uint64_t>(storage_layout.logical_local_complex_cells);
  estimate.allocated_local_complex_cells = static_cast<std::uint64_t>(storage_layout.allocated_local_complex_cells);
  estimate.used_backend_allocation_query = storage_layout.used_backend_allocation_query;

  estimate.real_array_bytes = checked_bytes(
      storage_layout.real_element_count,
      sizeof(double),
      "PM plan-memory real bytes");

  estimate.complex_spectral_array_bytes = checked_bytes(
      checkedProduct(3U, storage_layout.allocated_local_complex_cells,
                     "PM plan-memory complex spectral arrays"),
      sizeof(std::complex<double>),
      "PM plan-memory complex spectral bytes");
  estimate.scalar_spectral_array_bytes = checked_bytes(
      checkedProduct(4U, storage_layout.allocated_local_complex_cells,
                     "PM plan-memory scalar spectral arrays"),
      sizeof(double),
      "PM plan-memory scalar spectral bytes");
  const auto checked_add_bytes = [](std::uint64_t lhs, std::uint64_t rhs,
                                    std::string_view context) {
    if (rhs > std::numeric_limits<std::uint64_t>::max() - lhs) {
      throw std::overflow_error(std::string(context) + " overflow");
    }
    return lhs + rhs;
  };
  const std::uint64_t spectral_bytes = checked_add_bytes(
      estimate.complex_spectral_array_bytes,
      estimate.scalar_spectral_array_bytes,
      "PM plan-memory spectral byte sum");
  estimate.total_owned_bytes = checked_add_bytes(
      estimate.real_array_bytes,
      spectral_bytes,
      "PM plan-memory total byte sum");
  return estimate;
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
  m_totals.routed_send_buffer_high_water_bytes = std::max(
      m_totals.routed_send_buffer_high_water_bytes, event.routed_send_buffer_high_water_bytes);
  m_totals.routed_receive_buffer_high_water_bytes = std::max(
      m_totals.routed_receive_buffer_high_water_bytes, event.routed_receive_buffer_high_water_bytes);
  m_totals.routed_combined_buffer_high_water_bytes = std::max(
      m_totals.routed_combined_buffer_high_water_bytes, event.routed_combined_buffer_high_water_bytes);
  m_totals.routed_workspace_high_water_bytes = std::max(
      m_totals.routed_workspace_high_water_bytes, event.routed_workspace_high_water_bytes);
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

std::size_t PmGridStorage::localCellCount() const {
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
                                       .current_size_bytes = core::currentSizeBytesForContainer(container),
                                       .owned_capacity_bytes = bytes,
                                       .high_water_bytes = bytes,
                                       .estimated_next_step_bytes = 0U,
                                       .uncertainty_note = "retained PM-grid capacity is an observed allocation high-water; next-step requirement is not predicted"});
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

void PmSolver::appendMemoryReport(core::MemoryReportBuilder& builder) const {
  m_impl->appendOwnedMemoryReport(builder);
}

void PmSolver::shutdownBackendResources() {
  m_impl->shutdownBackendResources();
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
  const PmDecompositionView decomposition_view(
      m_shape, grid.slabLayout(), options.decomposition_mode);
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
        !std::isfinite(mass[p]) || mass[p] < 0.0) {
      first_non_finite_source_index = p;
      break;
    }
  }
  if (!distributed_slabs && first_non_finite_source_index.has_value()) {
    throw std::invalid_argument(
        "PmSolver::assignDensity requires finite particle coordinates and finite non-negative masses; particle_index=" +
        std::to_string(*first_non_finite_source_index));
  }

  const auto accumulate_owned = [&](const PmDensityContributionRecord& record) {
    if (record.global_ix >= m_shape.nx || record.global_iy >= m_shape.ny || record.global_iz >= m_shape.nz) {
      throw std::invalid_argument("PmSolver::assignDensity received out-of-range PM contribution record");
    }
    if (!grid.slabLayout().ownsGlobalX(record.global_ix)) {
      throw std::invalid_argument("PmSolver::assignDensity received contribution for non-owned PM slab x-index");
    }
    if (!std::isfinite(record.mass_contribution) || record.mass_contribution < 0.0) {
      throw std::invalid_argument(
          "PmSolver::assignDensity received a non-finite or negative mass contribution");
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
    const int world_size = grid.slabLayout().world_size;
    const int world_rank = grid.slabLayout().world_rank;
    double routed_mpi_wait_ms = distributed_preflight_mpi_wait_ms;
    const int local_sources_valid = first_non_finite_source_index.has_value() ? 0 : 1;
    int all_sources_valid = 0;
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      requirePmMpiSuccess(
          MPI_Allreduce(&local_sources_valid, &all_sources_valid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD),
          "PmSolver::assignDensity source validation MPI_Allreduce");
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
        world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() { exchange_epoch = m_impl->nextDistributedExchangeEpoch(); });
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      validatePmExchangeEpochConsensus(exchange_epoch, "PmSolver::assignDensity");
    });

    Impl::DensityExchangeBuffers* exchange_ptr = nullptr;
    std::uint64_t routing_metadata_capacity = 0U;
    PmRoutingCapacityModel routing_capacity{};
    std::size_t routing_buffer_limit = 0U;
    std::size_t particles_per_round = 0U;
    runPmCoordinatedPhase(
        "PmSolver::assignDensity routing preflight",
        world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
          exchange_ptr = &m_impl->densityExchangeBuffersForLayout(grid.slabLayout());
          auto& prepared_exchange = *exchange_ptr;
          const auto add_routing_metadata = [&](const auto& values) {
            routing_metadata_capacity = checkedAddBytes(
                routing_metadata_capacity,
                core::ownedCapacityBytesForContainer(values),
                "PM density routing metadata capacity");
          };
          add_routing_metadata(prepared_exchange.send_counts);
          add_routing_metadata(prepared_exchange.send_displs);
          add_routing_metadata(prepared_exchange.recv_counts);
          add_routing_metadata(prepared_exchange.recv_displs);
          add_routing_metadata(prepared_exchange.send_counts_bytes);
          add_routing_metadata(prepared_exchange.send_displs_bytes);
          add_routing_metadata(prepared_exchange.recv_counts_bytes);
          add_routing_metadata(prepared_exchange.recv_displs_bytes);
          add_routing_metadata(prepared_exchange.cursor);
          routing_capacity = modelPmRoutingCapacity(
              world_size,
              options.routing_exchange_batch_bytes,
              k_pm_routing_modeled_workspace_limit_bytes,
              routing_metadata_capacity,
              k_pm_density_plane_wire_bytes);
          if (routing_capacity.max_send_payload_bytes >
              static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
            throw std::overflow_error("PM density routing buffer limit exceeds size_t capacity");
          }
          routing_buffer_limit =
              static_cast<std::size_t>(routing_capacity.max_send_payload_bytes);
          particles_per_round = pmRoutingParticlesPerRound(
              routing_capacity.effective_per_peer_payload_bytes,
              k_pm_density_plane_wire_bytes);
        });
    auto& exchange = *exchange_ptr;
    const std::uint64_t local_source_count = static_cast<std::uint64_t>(mass.size());
    std::uint64_t global_max_source_count = 0U;
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      requirePmMpiSuccess(
          MPI_Allreduce(
              &local_source_count,
              &global_max_source_count,
              1,
              MPI_UINT64_T,
              MPI_MAX,
              MPI_COMM_WORLD),
          "PmSolver::assignDensity max-source MPI_Allreduce");
    });
    const std::uint64_t round_capacity = static_cast<std::uint64_t>(particles_per_round);
    const std::uint64_t round_count = global_max_source_count == 0U
        ? 0U
        : 1U + (global_max_source_count - 1U) / round_capacity;
    const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
    std::uint32_t next_record_sequence = 0U;

    const auto accumulate_plane_record = [&](const PmDensityPlaneRecord& record) {
      if (record.plane_count == 0U || record.plane_count > 3U ||
          record.destination_rank != static_cast<std::uint32_t>(world_rank) ||
          record.exchange_epoch != exchange_epoch ||
          !std::isfinite(record.y_grid) || !std::isfinite(record.z_grid) ||
          !std::isfinite(record.mass_code) || record.mass_code < 0.0) {
        throw std::invalid_argument("PmSolver::assignDensity received invalid bounded density-plane record");
      }
      const PmAxisStencil1d sy = makeAxisStencil(record.y_grid, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(record.z_grid, options.assignment_scheme);
      for (std::size_t plane = 0; plane < static_cast<std::size_t>(record.plane_count); ++plane) {
        const std::size_t ix = static_cast<std::size_t>(record.global_ix[plane]);
        if (ix >= m_shape.nx || !grid.slabLayout().ownsGlobalX(ix) ||
            !std::isfinite(record.x_weight[plane]) || record.x_weight[plane] < 0.0) {
          throw std::invalid_argument("PmSolver::assignDensity density-plane x ownership/weight is invalid");
        }
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          if (!periodic &&
              (sy.offsets[dy] < 0 || sy.offsets[dy] >= static_cast<std::ptrdiff_t>(m_shape.ny))) {
            continue;
          }
          const std::size_t iy = periodic
              ? wrapIndex(sy.offsets[dy], m_shape.ny)
              : static_cast<std::size_t>(sy.offsets[dy]);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            if (!periodic &&
                (sz.offsets[dz] < 0 || sz.offsets[dz] >= static_cast<std::ptrdiff_t>(m_shape.nz))) {
              continue;
            }
            const std::size_t iz = periodic
                ? wrapIndex(sz.offsets[dz], m_shape.nz)
                : static_cast<std::size_t>(sz.offsets[dz]);
            const double weight = record.x_weight[plane] * sy.weights[dy] * sz.weights[dz];
            grid.density()[grid.linearIndex(ix, iy, iz)] += record.mass_code * weight;
          }
        }
      }
    };

    std::uint64_t logical_records_sent = 0U;
    std::uint64_t peers_touched = 0U;
    std::uint64_t total_wire_sent = 0U;
    std::uint64_t total_wire_received = 0U;
    std::uint64_t send_high_water = 0U;
    std::uint64_t receive_high_water = 0U;
    std::uint64_t combined_high_water = 0U;
    std::uint64_t workspace_high_water = routing_metadata_capacity;

    for (std::uint64_t round = 0U; round < round_count; ++round) {
      const std::uint64_t begin_u64 = round * round_capacity;
      const std::uint64_t end_u64 = std::min(
          local_source_count,
          begin_u64 > std::numeric_limits<std::uint64_t>::max() - round_capacity
              ? local_source_count
              : begin_u64 + round_capacity);
      const std::size_t begin = static_cast<std::size_t>(begin_u64);
      const std::size_t end = static_cast<std::size_t>(end_u64);

      std::fill(exchange.send_counts.begin(), exchange.send_counts.end(), 0);
      std::fill(exchange.send_displs.begin(), exchange.send_displs.end(), 0);
      std::fill(exchange.recv_counts.begin(), exchange.recv_counts.end(), 0);
      std::fill(exchange.recv_displs.begin(), exchange.recv_displs.end(), 0);
      std::fill(exchange.send_counts_bytes.begin(), exchange.send_counts_bytes.end(), 0);
      std::fill(exchange.send_displs_bytes.begin(), exchange.send_displs_bytes.end(), 0);
      std::fill(exchange.recv_counts_bytes.begin(), exchange.recv_counts_bytes.end(), 0);
      std::fill(exchange.recv_displs_bytes.begin(), exchange.recv_displs_bytes.end(), 0);
      exchange.send_wire.clear();
      exchange.recv_wire.clear();

      runPmCoordinatedPhase(
          "PmSolver::assignDensity bounded density send preparation",
          world_rank,
          world_size,
          routed_mpi_wait_ms,
          [&]() {
            // Pass 1 counts one compact record per particle/destination x-plane group.
            for (std::size_t p = begin; p < end; ++p) {
              if (!periodic &&
                  (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
                   !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
                   !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
                continue;
              }
              const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
              const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
              std::array<PmXPlaneGroup, 3> groups{};
              const std::size_t group_count = makePmXPlaneGroups(
                  sx, periodic, m_shape, decomposition_view, groups);
              for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
                const int destination_rank = groups[group_i].destination_rank;
                if (destination_rank == world_rank) {
                  continue;
                }
                int& count = exchange.send_counts[static_cast<std::size_t>(destination_rank)];
                if (count == std::numeric_limits<int>::max()) {
                  throw std::overflow_error("PM density per-peer record count exceeds MPI int capacity");
                }
                ++count;
              }
            }

            std::size_t total_send_records = 0U;
            for (int rank = 0; rank < world_size; ++rank) {
              const std::size_t count = static_cast<std::size_t>(
                  exchange.send_counts[static_cast<std::size_t>(rank)]);
              const std::uint64_t peer_bytes = checkedBytesForCount(
                  count, k_pm_density_plane_wire_bytes, "PM density per-peer routing batch");
              if (peer_bytes > routing_capacity.effective_per_peer_payload_bytes) {
                throw std::runtime_error("PM density routing exceeded effective per-peer batch high-water");
              }
              exchange.send_displs[static_cast<std::size_t>(rank)] =
                  checkedMpiRecordCountOrDisplacement(
                      total_send_records,
                      k_pm_density_plane_wire_bytes,
                      "PmSolver::assignDensity bounded density exchange",
                      "send record displacement",
                      rank);
              total_send_records = checkedMpiRecordTotal(
                  total_send_records,
                  count,
                  k_pm_density_plane_wire_bytes,
                  "PmSolver::assignDensity bounded density exchange",
                  rank);
            }
            resizePmWireBufferBounded(
                exchange.send_wire,
                checkedPmWireByteCount(
                    total_send_records,
                    k_pm_density_plane_wire_bytes,
                    "PmSolver::assignDensity bounded density exchange"),
                routing_buffer_limit,
                "PmSolver::assignDensity send buffer");

            for (int rank = 0; rank < world_size; ++rank) {
              exchange.cursor[static_cast<std::size_t>(rank)] = static_cast<std::size_t>(
                  exchange.send_displs[static_cast<std::size_t>(rank)]);
            }

            // Pass 2 deposits owner-local x planes immediately and encodes only
            // remote x-plane groups directly into the bounded wire buffer.
            for (std::size_t p = begin; p < end; ++p) {
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
              std::array<PmXPlaneGroup, 3> groups{};
              const std::size_t group_count = makePmXPlaneGroups(
                  sx, periodic, m_shape, decomposition_view, groups);
              for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
                const auto& group = groups[group_i];
                PmDensityPlaneRecord record{
                    .origin_rank = static_cast<std::uint32_t>(world_rank),
                    .destination_rank = static_cast<std::uint32_t>(group.destination_rank),
                    .record_sequence = next_record_sequence,
                    .plane_count = group.plane_count,
                    .global_ix = group.global_ix,
                    .exchange_epoch = exchange_epoch,
                    .y_grid = y,
                    .z_grid = z,
                    .mass_code = mass[p],
                    .x_weight = group.x_weight,
                };
                if (group.destination_rank == world_rank) {
                  accumulate_plane_record(record);
                  continue;
                }
                if (next_record_sequence == std::numeric_limits<std::uint32_t>::max()) {
                  throw std::overflow_error("PmSolver::assignDensity routing sequence exceeds uint32 capacity");
                }
                auto& cursor = exchange.cursor[static_cast<std::size_t>(group.destination_rank)];
                const std::size_t byte_offset = cursor * k_pm_density_plane_wire_bytes;
                encodePmDensityPlaneRecord(
                    record,
                    std::span<std::uint8_t>(exchange.send_wire).subspan(
                        byte_offset, k_pm_density_plane_wire_bytes));
                ++cursor;
                ++next_record_sequence;
              }
            }
          });

      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoall(
                exchange.send_counts.data(), 1, MPI_INT,
                exchange.recv_counts.data(), 1, MPI_INT,
                MPI_COMM_WORLD),
            "PmSolver::assignDensity count MPI_Alltoall");
      });

      runPmCoordinatedPhase(
          "PmSolver::assignDensity bounded density receive layout",
          world_rank,
          world_size,
          routed_mpi_wait_ms,
          [&]() {
            std::size_t total_recv_records = 0U;
            for (int rank = 0; rank < world_size; ++rank) {
              const std::size_t count = checkedMpiReceivedRecordCount(
                  exchange.recv_counts[static_cast<std::size_t>(rank)],
                  k_pm_density_plane_wire_bytes,
                  "PmSolver::assignDensity bounded density exchange",
                  rank);
              exchange.recv_displs[static_cast<std::size_t>(rank)] =
                  checkedMpiRecordCountOrDisplacement(
                      total_recv_records,
                      k_pm_density_plane_wire_bytes,
                      "PmSolver::assignDensity bounded density exchange",
                      "received record displacement",
                      rank);
              total_recv_records = checkedMpiRecordTotal(
                  total_recv_records,
                  count,
                  k_pm_density_plane_wire_bytes,
                  "PmSolver::assignDensity bounded density exchange",
                  rank);
            }
            resizePmWireBufferBounded(
                exchange.recv_wire,
                checkedPmWireByteCount(
                    total_recv_records,
                    k_pm_density_plane_wire_bytes,
                    "PmSolver::assignDensity bounded density exchange"),
                routing_buffer_limit,
                "PmSolver::assignDensity receive buffer");
            checkedMpiRecordLayoutToByteLayout(
                exchange.send_counts,
                exchange.send_displs,
                k_pm_density_plane_wire_bytes,
                exchange.send_counts_bytes,
                exchange.send_displs_bytes,
                "PmSolver::assignDensity bounded density send MPI_Alltoallv");
            checkedMpiRecordLayoutToByteLayout(
                exchange.recv_counts,
                exchange.recv_displs,
                k_pm_density_plane_wire_bytes,
                exchange.recv_counts_bytes,
                exchange.recv_displs_bytes,
                "PmSolver::assignDensity bounded density receive MPI_Alltoallv");
          });

      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoallv(
                nonNullPmWireData(exchange.send_wire),
                exchange.send_counts_bytes.data(),
                exchange.send_displs_bytes.data(),
                MPI_BYTE,
                nonNullPmWireData(exchange.recv_wire),
                exchange.recv_counts_bytes.data(),
                exchange.recv_displs_bytes.data(),
                MPI_BYTE,
                MPI_COMM_WORLD),
            "PmSolver::assignDensity payload MPI_Alltoallv");
      });

      std::string local_density_validation_error;
      try {
        for (int source_rank = 0; source_rank < world_size; ++source_rank) {
          const int begin_record = exchange.recv_displs[static_cast<std::size_t>(source_rank)];
          const int count = exchange.recv_counts[static_cast<std::size_t>(source_rank)];
          std::optional<std::uint32_t> previous_sequence;
          for (int i = 0; i < count; ++i) {
            const std::size_t record_index = static_cast<std::size_t>(begin_record + i);
            const PmDensityPlaneRecord record = decodePmDensityPlaneRecord(
                std::span<const std::uint8_t>(exchange.recv_wire).subspan(
                    record_index * k_pm_density_plane_wire_bytes,
                    k_pm_density_plane_wire_bytes));
            if (record.origin_rank != static_cast<std::uint32_t>(source_rank) ||
                record.destination_rank != static_cast<std::uint32_t>(world_rank) ||
                record.exchange_epoch != exchange_epoch) {
              throw std::invalid_argument("PM density-plane routing metadata does not match MPI segment");
            }
            if (previous_sequence.has_value() && record.record_sequence <= *previous_sequence) {
              throw std::invalid_argument("PM density-plane sequence is non-monotonic within sender segment");
            }
            previous_sequence = record.record_sequence;
            accumulate_plane_record(record);
          }
        }
      } catch (const std::exception& error) {
        local_density_validation_error = error.what();
      } catch (...) {
        local_density_validation_error = "non-standard exception while consuming bounded density payload";
      }
      throwIfPmPayloadValidationFailed(
          "PmSolver::assignDensity bounded density receive",
          world_rank,
          world_size,
          local_density_validation_error,
          routed_mpi_wait_ms);

      const std::uint64_t round_send_bytes = static_cast<std::uint64_t>(exchange.send_wire.size());
      const std::uint64_t round_recv_bytes = static_cast<std::uint64_t>(exchange.recv_wire.size());
      const std::uint64_t round_send_capacity = static_cast<std::uint64_t>(exchange.send_wire.capacity());
      const std::uint64_t round_recv_capacity = static_cast<std::uint64_t>(exchange.recv_wire.capacity());
      const std::uint64_t round_combined_capacity = checkedAddBytes(
          round_send_capacity,
          round_recv_capacity,
          "PM density combined routing buffer capacity");
      const std::uint64_t round_workspace_capacity = checkedAddBytes(
          routing_metadata_capacity,
          round_combined_capacity,
          "PM density routing workspace capacity");
      if (round_workspace_capacity > k_pm_routing_workspace_target_bytes) {
        throw std::logic_error("PM density routing workspace exceeded the M1A aggregate target");
      }
      logical_records_sent += static_cast<std::uint64_t>(
          exchange.send_wire.size() / k_pm_density_plane_wire_bytes);
      peers_touched += static_cast<std::uint64_t>(std::count_if(
          exchange.send_counts.begin(), exchange.send_counts.end(), [](int count) { return count > 0; }));
      total_wire_sent += round_send_bytes;
      total_wire_received += round_recv_bytes;
      send_high_water = std::max(send_high_water, round_send_capacity);
      receive_high_water = std::max(receive_high_water, round_recv_capacity);
      combined_high_water = std::max(combined_high_water, round_combined_capacity);
      workspace_high_water = std::max(workspace_high_water, round_workspace_capacity);
      exchange.workspace_high_water_bytes = std::max(
          exchange.workspace_high_water_bytes, round_workspace_capacity);
    }

    if (profile != nullptr) {
      profile->routed_density_records += logical_records_sent;
      profile->routed_density_peer_count += peers_touched;
      profile->routed_mpi_bytes_sent += total_wire_sent;
      profile->routed_mpi_bytes_received += total_wire_received;
      profile->routed_send_buffer_high_water_bytes = std::max(
          profile->routed_send_buffer_high_water_bytes, send_high_water);
      profile->routed_receive_buffer_high_water_bytes = std::max(
          profile->routed_receive_buffer_high_water_bytes, receive_high_water);
      profile->routed_combined_buffer_high_water_bytes = std::max(
          profile->routed_combined_buffer_high_water_bytes, combined_high_water);
      profile->routed_workspace_high_water_bytes = std::max(
          profile->routed_workspace_high_water_bytes, workspace_high_water);
      profile->routed_mpi_wait_ms += routed_mpi_wait_ms;
      profile->bytes_moved += total_wire_sent + total_wire_received;
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
        target_view_preflight_mpi_wait_ms);
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
          const bool indexed_coordinates =
              target_view.coordinate_layout == PmForceCoordinateLayout::kIndexedSource;
          if (target_view.pos_x_comoving.size() != target_view.pos_y_comoving.size() ||
              target_view.pos_x_comoving.size() != target_view.pos_z_comoving.size()) {
            throw std::invalid_argument(
                "PmSolver::interpolateForces coordinate source extents must match");
          }
          if (indexed_coordinates) {
            if (target_view.coordinate_source_index.size() != local_active_count) {
              throw std::invalid_argument(
                  "PmSolver::interpolateForces indexed coordinate map extent must match active count");
            }
            for (const TreeLocalIndex source_index : target_view.coordinate_source_index) {
              if (static_cast<std::size_t>(source_index) >= target_view.pos_x_comoving.size()) {
                throw std::out_of_range(
                    "PmSolver::interpolateForces coordinate source index out of range");
              }
            }
          } else {
            if (!target_view.coordinate_source_index.empty()) {
              throw std::invalid_argument(
                  "PmSolver::interpolateForces compact coordinates must not provide source indices");
            }
            if (local_active_count != target_view.pos_x_comoving.size()) {
              throw std::invalid_argument(
                  "PmSolver::interpolateForces compact coordinate extent must match active count");
            }
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
  const bool indexed_coordinates =
      target_view.coordinate_layout == PmForceCoordinateLayout::kIndexedSource;
  if (target_view.pos_x_comoving.size() != target_view.pos_y_comoving.size() ||
      target_view.pos_x_comoving.size() != target_view.pos_z_comoving.size()) {
    throw std::invalid_argument(
        "PmSolver::interpolateForces coordinate source extents must match");
  }
  if (indexed_coordinates) {
    if (target_view.coordinate_source_index.size() != active_count) {
      throw std::invalid_argument(
          "PmSolver::interpolateForces indexed coordinate map extent must match active count");
    }
  } else {
    if (!target_view.coordinate_source_index.empty()) {
      throw std::invalid_argument(
          "PmSolver::interpolateForces compact coordinates must not provide source indices");
    }
    if (active_count != target_view.pos_x_comoving.size()) {
      throw std::invalid_argument(
          "PmSolver::interpolateForces compact coordinate extent must match active count");
    }
  }

  auto& indexed_workspace = m_impl->indexedTargetWorkspace();
  std::span<const double> target_x = target_view.pos_x_comoving;
  std::span<const double> target_y = target_view.pos_y_comoving;
  std::span<const double> target_z = target_view.pos_z_comoving;
  if (indexed_coordinates) {
    indexed_workspace.gathered_x.resize(active_count);
    indexed_workspace.gathered_y.resize(active_count);
    indexed_workspace.gathered_z.resize(active_count);
    for (std::size_t active_i = 0; active_i < active_count; ++active_i) {
      const std::size_t source_i = static_cast<std::size_t>(
          target_view.coordinate_source_index[active_i]);
      if (source_i >= target_view.pos_x_comoving.size()) {
        throw std::out_of_range(
            "PmSolver::interpolateForces coordinate source index out of range");
      }
      indexed_workspace.gathered_x[active_i] = target_view.pos_x_comoving[source_i];
      indexed_workspace.gathered_y[active_i] = target_view.pos_y_comoving[source_i];
      indexed_workspace.gathered_z[active_i] = target_view.pos_z_comoving[source_i];
    }
    target_x = indexed_workspace.gathered_x;
    target_y = indexed_workspace.gathered_y;
    target_z = indexed_workspace.gathered_z;
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
          target_x,
          target_y,
          target_z,
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
      indexed_workspace.compact_ax.assign(active_count, 0.0);
      indexed_workspace.compact_ay.assign(active_count, 0.0);
      indexed_workspace.compact_az.assign(active_count, 0.0);
      interpolateForces(
          grid,
          target_x,
          target_y,
          target_z,
          indexed_workspace.compact_ax,
          indexed_workspace.compact_ay,
          indexed_workspace.compact_az,
          options,
          profile);
      for (std::size_t active_i = 0; active_i < active_count; ++active_i) {
        const std::uint32_t particle_index = target_view.active_particle_index[active_i];
        target_view.accel_x_comoving[particle_index] = indexed_workspace.compact_ax[active_i];
        target_view.accel_y_comoving[particle_index] = indexed_workspace.compact_ay[active_i];
        target_view.accel_z_comoving[particle_index] = indexed_workspace.compact_az[active_i];
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
  const PmDecompositionView decomposition_view(
      m_shape, grid.slabLayout(), options.decomposition_mode);
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
    const int world_rank = grid.slabLayout().world_rank;
    double routed_mpi_wait_ms = distributed_preflight_mpi_wait_ms;
    const int local_targets_valid = first_non_finite_target_index.has_value() ? 0 : 1;
    int all_targets_valid = 0;
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      requirePmMpiSuccess(
          MPI_Allreduce(&local_targets_valid, &all_targets_valid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD),
          "PmSolver::interpolateForces target validation MPI_Allreduce");
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
        world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() { exchange_epoch = m_impl->nextDistributedExchangeEpoch(); });
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      validatePmExchangeEpochConsensus(exchange_epoch, "PmSolver::interpolateForces");
    });

    Impl::PlaneInterpolationExchangeBuffers* exchange_ptr = nullptr;
    std::uint64_t routing_metadata_capacity = 0U;
    PmRoutingCapacityModel routing_capacity{};
    std::size_t routing_buffer_limit = 0U;
    std::size_t particles_per_round = 0U;
    runPmCoordinatedPhase(
        "PmSolver::interpolateForces routing preflight",
        world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
          exchange_ptr = &m_impl->planeInterpolationExchangeBuffersForLayout(grid.slabLayout());
          auto& prepared_exchange = *exchange_ptr;
          const auto add_routing_metadata = [&](const auto& values) {
            routing_metadata_capacity = checkedAddBytes(
                routing_metadata_capacity,
                core::ownedCapacityBytesForContainer(values),
                "PM force routing metadata capacity");
          };
          add_routing_metadata(prepared_exchange.send_counts);
          add_routing_metadata(prepared_exchange.send_displs);
          add_routing_metadata(prepared_exchange.recv_counts);
          add_routing_metadata(prepared_exchange.recv_displs);
          add_routing_metadata(prepared_exchange.send_counts_bytes);
          add_routing_metadata(prepared_exchange.send_displs_bytes);
          add_routing_metadata(prepared_exchange.recv_counts_bytes);
          add_routing_metadata(prepared_exchange.recv_displs_bytes);
          add_routing_metadata(prepared_exchange.send_response_counts_bytes);
          add_routing_metadata(prepared_exchange.send_response_displs_bytes);
          add_routing_metadata(prepared_exchange.recv_response_counts_bytes);
          add_routing_metadata(prepared_exchange.recv_response_displs_bytes);
          add_routing_metadata(prepared_exchange.cursor);
          routing_capacity = modelPmRoutingCapacity(
              world_size,
              options.routing_exchange_batch_bytes,
              k_pm_routing_modeled_workspace_limit_bytes,
              routing_metadata_capacity,
              k_pm_plane_interpolation_request_wire_bytes);
          if (routing_capacity.max_send_payload_bytes >
              static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
            throw std::overflow_error("PM force routing buffer limit exceeds size_t capacity");
          }
          routing_buffer_limit =
              static_cast<std::size_t>(routing_capacity.max_send_payload_bytes);
          particles_per_round = pmRoutingParticlesPerRound(
              routing_capacity.effective_per_peer_payload_bytes,
              k_pm_plane_interpolation_request_wire_bytes);
        });
    auto& exchange = *exchange_ptr;
    const std::uint64_t local_target_count = static_cast<std::uint64_t>(pos_x.size());
    std::uint64_t global_max_target_count = 0U;
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      requirePmMpiSuccess(
          MPI_Allreduce(
              &local_target_count,
              &global_max_target_count,
              1,
              MPI_UINT64_T,
              MPI_MAX,
              MPI_COMM_WORLD),
          "PmSolver::interpolateForces max-target MPI_Allreduce");
    });
    const std::uint64_t round_capacity = static_cast<std::uint64_t>(particles_per_round);
    const std::uint64_t round_count = global_max_target_count == 0U
        ? 0U
        : 1U + (global_max_target_count - 1U) / round_capacity;
    const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
    std::uint32_t next_request_sequence = 0U;

    std::fill(accel_x.begin(), accel_x.end(), 0.0);
    std::fill(accel_y.begin(), accel_y.end(), 0.0);
    std::fill(accel_z.begin(), accel_z.end(), 0.0);

    const auto evaluate_plane_request = [&](const PmPlaneInterpolationRequestRecord& request) {
      if (request.plane_count == 0U || request.plane_count > 3U ||
          request.destination_rank != static_cast<std::uint32_t>(world_rank) ||
          request.exchange_epoch != exchange_epoch ||
          !std::isfinite(request.y_grid) || !std::isfinite(request.z_grid)) {
        throw std::invalid_argument(
            "PmSolver::interpolateForces received invalid bounded force-plane request");
      }
      const PmAxisStencil1d sy = makeAxisStencil(request.y_grid, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(request.z_grid, options.assignment_scheme);
      PmForceContributionRecord response{
          .source_rank = static_cast<std::uint32_t>(world_rank),
          .origin_rank = request.origin_rank,
          .request_sequence = request.request_sequence,
          .origin_particle_index = request.origin_particle_index,
          .exchange_epoch = request.exchange_epoch,
          .accel_x = 0.0,
          .accel_y = 0.0,
          .accel_z = 0.0,
      };
      for (std::size_t plane = 0; plane < static_cast<std::size_t>(request.plane_count); ++plane) {
        const std::size_t ix = static_cast<std::size_t>(request.global_ix[plane]);
        if (ix >= m_shape.nx || !grid.slabLayout().ownsGlobalX(ix) ||
            !std::isfinite(request.x_weight[plane])) {
          throw std::invalid_argument(
              "PmSolver::interpolateForces force-plane x ownership/weight is invalid");
        }
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          if (!periodic &&
              (sy.offsets[dy] < 0 || sy.offsets[dy] >= static_cast<std::ptrdiff_t>(m_shape.ny))) {
            continue;
          }
          const std::size_t iy = periodic
              ? wrapIndex(sy.offsets[dy], m_shape.ny)
              : static_cast<std::size_t>(sy.offsets[dy]);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            if (!periodic &&
                (sz.offsets[dz] < 0 || sz.offsets[dz] >= static_cast<std::ptrdiff_t>(m_shape.nz))) {
              continue;
            }
            const std::size_t iz = periodic
                ? wrapIndex(sz.offsets[dz], m_shape.nz)
                : static_cast<std::size_t>(sz.offsets[dz]);
            const std::size_t local_index = decomposition_view.globalToLocalRealIndex(ix, iy, iz);
            if (!std::isfinite(grid.force_x()[local_index]) ||
                !std::isfinite(grid.force_y()[local_index]) ||
                !std::isfinite(grid.force_z()[local_index])) {
              throw std::invalid_argument(
                  "PmSolver::interpolateForces encountered non-finite owner-local mesh force data");
            }
            const double weight = request.x_weight[plane] * sy.weights[dy] * sz.weights[dz];
            response.accel_x += weight * grid.force_x()[local_index];
            response.accel_y += weight * grid.force_y()[local_index];
            response.accel_z += weight * grid.force_z()[local_index];
          }
        }
      }
      return response;
    };

    std::uint64_t logical_requests_sent = 0U;
    std::uint64_t peers_touched = 0U;
    std::uint64_t total_wire_sent = 0U;
    std::uint64_t total_wire_received = 0U;
    std::uint64_t send_high_water = 0U;
    std::uint64_t receive_high_water = 0U;
    std::uint64_t combined_high_water = 0U;
    std::uint64_t workspace_high_water = routing_metadata_capacity;

    for (std::uint64_t round = 0U; round < round_count; ++round) {
      const std::uint64_t begin_u64 = round * round_capacity;
      const std::uint64_t end_u64 = std::min(
          local_target_count,
          begin_u64 > std::numeric_limits<std::uint64_t>::max() - round_capacity
              ? local_target_count
              : begin_u64 + round_capacity);
      const std::size_t begin = static_cast<std::size_t>(begin_u64);
      const std::size_t end = static_cast<std::size_t>(end_u64);
      const std::uint32_t round_sequence_begin = next_request_sequence;

      std::fill(exchange.send_counts.begin(), exchange.send_counts.end(), 0);
      std::fill(exchange.send_displs.begin(), exchange.send_displs.end(), 0);
      std::fill(exchange.recv_counts.begin(), exchange.recv_counts.end(), 0);
      std::fill(exchange.recv_displs.begin(), exchange.recv_displs.end(), 0);
      std::fill(exchange.send_counts_bytes.begin(), exchange.send_counts_bytes.end(), 0);
      std::fill(exchange.send_displs_bytes.begin(), exchange.send_displs_bytes.end(), 0);
      std::fill(exchange.recv_counts_bytes.begin(), exchange.recv_counts_bytes.end(), 0);
      std::fill(exchange.recv_displs_bytes.begin(), exchange.recv_displs_bytes.end(), 0);
      std::fill(exchange.send_response_counts_bytes.begin(), exchange.send_response_counts_bytes.end(), 0);
      std::fill(exchange.send_response_displs_bytes.begin(), exchange.send_response_displs_bytes.end(), 0);
      std::fill(exchange.recv_response_counts_bytes.begin(), exchange.recv_response_counts_bytes.end(), 0);
      std::fill(exchange.recv_response_displs_bytes.begin(), exchange.recv_response_displs_bytes.end(), 0);
      exchange.send_wire.clear();
      exchange.recv_wire.clear();

      runPmCoordinatedPhase(
          "PmSolver::interpolateForces bounded request preparation",
          world_rank,
          world_size,
          routed_mpi_wait_ms,
          [&]() {
            // Pass 1 counts at most one request per particle/destination x-plane group.
            for (std::size_t p = begin; p < end; ++p) {
              if (p > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
                throw std::invalid_argument(
                    "PmSolver::interpolateForces origin particle index exceeds routing token limit");
              }
              if (!periodic &&
                  (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
                   !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
                   !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
                continue;
              }
              const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
              const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
              std::array<PmXPlaneGroup, 3> groups{};
              const std::size_t group_count = makePmXPlaneGroups(
                  sx, periodic, m_shape, decomposition_view, groups);
              for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
                const int destination_rank = groups[group_i].destination_rank;
                if (destination_rank == world_rank) {
                  continue;
                }
                int& count = exchange.send_counts[static_cast<std::size_t>(destination_rank)];
                if (count == std::numeric_limits<int>::max()) {
                  throw std::overflow_error("PM force per-peer request count exceeds MPI int capacity");
                }
                ++count;
              }
            }

            std::size_t total_send_records = 0U;
            for (int rank = 0; rank < world_size; ++rank) {
              const std::size_t count = static_cast<std::size_t>(
                  exchange.send_counts[static_cast<std::size_t>(rank)]);
              const std::uint64_t peer_bytes = checkedBytesForCount(
                  count, k_pm_plane_interpolation_request_wire_bytes, "PM force per-peer routing batch");
              if (peer_bytes > routing_capacity.effective_per_peer_payload_bytes) {
                throw std::runtime_error(
                    "PM force routing exceeded effective per-peer batch high-water");
              }
              exchange.send_displs[static_cast<std::size_t>(rank)] =
                  checkedMpiRecordCountOrDisplacement(
                      total_send_records,
                      k_pm_plane_interpolation_request_wire_bytes,
                      "PmSolver::interpolateForces bounded request exchange",
                      "send record displacement",
                      rank);
              total_send_records = checkedMpiRecordTotal(
                  total_send_records,
                  count,
                  k_pm_plane_interpolation_request_wire_bytes,
                  "PmSolver::interpolateForces bounded request exchange",
                  rank);
            }
            resizePmWireBufferBounded(
                exchange.send_wire,
                checkedPmWireByteCount(
                    total_send_records,
                    k_pm_plane_interpolation_request_wire_bytes,
                    "PmSolver::interpolateForces bounded request exchange"),
                routing_buffer_limit,
                "PmSolver::interpolateForces request send buffer");

            for (int rank = 0; rank < world_size; ++rank) {
              exchange.cursor[static_cast<std::size_t>(rank)] = static_cast<std::size_t>(
                  exchange.send_displs[static_cast<std::size_t>(rank)]);
            }

            // Pass 2 evaluates local x planes immediately and encodes only remote
            // x-plane groups. The receiver expands the y/z stencil and returns one
            // accumulated force response per request.
            for (std::size_t p = begin; p < end; ++p) {
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
              std::array<PmXPlaneGroup, 3> groups{};
              const std::size_t group_count = makePmXPlaneGroups(
                  sx, periodic, m_shape, decomposition_view, groups);
              for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
                const auto& group = groups[group_i];
                PmPlaneInterpolationRequestRecord request{
                    .origin_rank = static_cast<std::uint32_t>(world_rank),
                    .destination_rank = static_cast<std::uint32_t>(group.destination_rank),
                    .request_sequence = next_request_sequence,
                    .origin_particle_index = static_cast<std::uint32_t>(p),
                    .plane_count = group.plane_count,
                    .global_ix = group.global_ix,
                    .exchange_epoch = exchange_epoch,
                    .y_grid = y,
                    .z_grid = z,
                    .x_weight = group.x_weight,
                };
                if (group.destination_rank == world_rank) {
                  const PmForceContributionRecord local = evaluate_plane_request(request);
                  accel_x[p] += local.accel_x;
                  accel_y[p] += local.accel_y;
                  accel_z[p] += local.accel_z;
                  continue;
                }
                if (next_request_sequence == std::numeric_limits<std::uint32_t>::max()) {
                  throw std::overflow_error(
                      "PmSolver::interpolateForces routing sequence exceeds uint32 capacity");
                }
                auto& cursor = exchange.cursor[static_cast<std::size_t>(group.destination_rank)];
                encodePmPlaneInterpolationRequest(
                    request,
                    PmWireRecordKind::kForcePlaneRequest,
                    std::span<std::uint8_t>(exchange.send_wire).subspan(
                        cursor * k_pm_plane_interpolation_request_wire_bytes,
                        k_pm_plane_interpolation_request_wire_bytes));
                ++cursor;
                ++next_request_sequence;
              }
            }
          });

      const std::uint64_t request_send_bytes = static_cast<std::uint64_t>(exchange.send_wire.size());
      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoall(
                exchange.send_counts.data(), 1, MPI_INT,
                exchange.recv_counts.data(), 1, MPI_INT,
                MPI_COMM_WORLD),
            "PmSolver::interpolateForces request-count MPI_Alltoall");
      });

      runPmCoordinatedPhase(
          "PmSolver::interpolateForces bounded request receive layout",
          world_rank,
          world_size,
          routed_mpi_wait_ms,
          [&]() {
            std::size_t total_recv_records = 0U;
            for (int rank = 0; rank < world_size; ++rank) {
              const std::size_t count = checkedMpiReceivedRecordCount(
                  exchange.recv_counts[static_cast<std::size_t>(rank)],
                  k_pm_plane_interpolation_request_wire_bytes,
                  "PmSolver::interpolateForces bounded request exchange",
                  rank);
              exchange.recv_displs[static_cast<std::size_t>(rank)] =
                  checkedMpiRecordCountOrDisplacement(
                      total_recv_records,
                      k_pm_plane_interpolation_request_wire_bytes,
                      "PmSolver::interpolateForces bounded request exchange",
                      "received record displacement",
                      rank);
              total_recv_records = checkedMpiRecordTotal(
                  total_recv_records,
                  count,
                  k_pm_plane_interpolation_request_wire_bytes,
                  "PmSolver::interpolateForces bounded request exchange",
                  rank);
            }
            resizePmWireBufferBounded(
                exchange.recv_wire,
                checkedPmWireByteCount(
                    total_recv_records,
                    k_pm_plane_interpolation_request_wire_bytes,
                    "PmSolver::interpolateForces bounded request exchange"),
                routing_buffer_limit,
                "PmSolver::interpolateForces request receive buffer");
            checkedMpiRecordLayoutToByteLayout(
                exchange.send_counts,
                exchange.send_displs,
                k_pm_plane_interpolation_request_wire_bytes,
                exchange.send_counts_bytes,
                exchange.send_displs_bytes,
                "PmSolver::interpolateForces bounded request send MPI_Alltoallv");
            checkedMpiRecordLayoutToByteLayout(
                exchange.recv_counts,
                exchange.recv_displs,
                k_pm_plane_interpolation_request_wire_bytes,
                exchange.recv_counts_bytes,
                exchange.recv_displs_bytes,
                "PmSolver::interpolateForces bounded request receive MPI_Alltoallv");
          });
      const std::uint64_t request_recv_bytes = static_cast<std::uint64_t>(exchange.recv_wire.size());

      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoallv(
                nonNullPmWireData(exchange.send_wire),
                exchange.send_counts_bytes.data(),
                exchange.send_displs_bytes.data(),
                MPI_BYTE,
                nonNullPmWireData(exchange.recv_wire),
                exchange.recv_counts_bytes.data(),
                exchange.recv_displs_bytes.data(),
                MPI_BYTE,
                MPI_COMM_WORLD),
            "PmSolver::interpolateForces request MPI_Alltoallv");
      });

      std::string local_force_request_validation_error;
      try {
        // Compact responses in-place into the receive-request buffer. Since a
        // response (64 B) is smaller than a request (96 B), forward compaction
        // cannot overwrite a request that has not yet been decoded.
        for (int source_rank = 0; source_rank < world_size; ++source_rank) {
          const int request_begin = exchange.recv_displs[static_cast<std::size_t>(source_rank)];
          const int request_count = exchange.recv_counts[static_cast<std::size_t>(source_rank)];
          exchange.send_response_counts_bytes[static_cast<std::size_t>(source_rank)] =
              checkedMpiByteCount(
                  checkedBytesForCount(
                      static_cast<std::size_t>(request_count),
                      k_pm_force_response_wire_bytes,
                      "PM force response sender count"),
                  "PM force response sender count");
          exchange.send_response_displs_bytes[static_cast<std::size_t>(source_rank)] =
              checkedMpiByteCount(
                  checkedBytesForCount(
                      static_cast<std::size_t>(request_begin),
                      k_pm_force_response_wire_bytes,
                      "PM force response sender displacement"),
                  "PM force response sender displacement");
          std::optional<std::uint32_t> previous_sequence;
          for (int i = 0; i < request_count; ++i) {
            const std::size_t request_index = static_cast<std::size_t>(request_begin + i);
            const PmPlaneInterpolationRequestRecord request = decodePmPlaneInterpolationRequest(
                std::span<const std::uint8_t>(exchange.recv_wire).subspan(
                    request_index * k_pm_plane_interpolation_request_wire_bytes,
                    k_pm_plane_interpolation_request_wire_bytes),
                PmWireRecordKind::kForcePlaneRequest,
                "PmSolver::interpolateForces force-plane request");
            if (request.origin_rank != static_cast<std::uint32_t>(source_rank) ||
                request.destination_rank != static_cast<std::uint32_t>(world_rank) ||
                request.exchange_epoch != exchange_epoch) {
              throw std::invalid_argument(
                  "PM force-plane routing metadata does not match MPI sender/receiver segment");
            }
            if (previous_sequence.has_value() && request.request_sequence <= *previous_sequence) {
              throw std::invalid_argument(
                  "PM force-plane request sequence is non-monotonic within sender segment");
            }
            previous_sequence = request.request_sequence;
            const PmForceContributionRecord response = evaluate_plane_request(request);
            const std::size_t response_index =
                static_cast<std::size_t>(exchange.send_response_displs_bytes[static_cast<std::size_t>(source_rank)]) /
                    k_pm_force_response_wire_bytes +
                static_cast<std::size_t>(i);
            encodePmForceResponse(
                response,
                std::span<std::uint8_t>(exchange.recv_wire).subspan(
                    response_index * k_pm_force_response_wire_bytes,
                    k_pm_force_response_wire_bytes));
          }
        }
      } catch (const std::exception& error) {
        local_force_request_validation_error = error.what();
      } catch (...) {
        local_force_request_validation_error =
            "non-standard exception while consuming bounded force-plane requests";
      }
      throwIfPmPayloadValidationFailed(
          "PmSolver::interpolateForces bounded request receive",
          world_rank,
          world_size,
          local_force_request_validation_error,
          routed_mpi_wait_ms);

      const std::size_t total_received_requests = exchange.recv_wire.size() /
          k_pm_plane_interpolation_request_wire_bytes;
      resizePmWireBufferBounded(
          exchange.recv_wire,
          checkedPmWireByteCount(
              total_received_requests,
              k_pm_force_response_wire_bytes,
              "PmSolver::interpolateForces compact response send"),
          routing_buffer_limit,
          "PmSolver::interpolateForces response send buffer");

      const std::size_t total_sent_requests = exchange.send_wire.size() /
          k_pm_plane_interpolation_request_wire_bytes;
      for (int rank = 0; rank < world_size; ++rank) {
        exchange.recv_response_counts_bytes[static_cast<std::size_t>(rank)] =
            checkedMpiByteCount(
                checkedBytesForCount(
                    static_cast<std::size_t>(exchange.send_counts[static_cast<std::size_t>(rank)]),
                    k_pm_force_response_wire_bytes,
                    "PM force response receive count"),
                "PM force response receive count");
        exchange.recv_response_displs_bytes[static_cast<std::size_t>(rank)] =
            checkedMpiByteCount(
                checkedBytesForCount(
                    static_cast<std::size_t>(exchange.send_displs[static_cast<std::size_t>(rank)]),
                    k_pm_force_response_wire_bytes,
                    "PM force response receive displacement"),
                "PM force response receive displacement");
      }
      resizePmWireBufferBounded(
          exchange.send_wire,
          checkedPmWireByteCount(
              total_sent_requests,
              k_pm_force_response_wire_bytes,
              "PmSolver::interpolateForces compact response receive"),
          routing_buffer_limit,
          "PmSolver::interpolateForces response receive buffer");

      const std::uint64_t response_send_bytes = static_cast<std::uint64_t>(exchange.recv_wire.size());
      const std::uint64_t response_recv_bytes = static_cast<std::uint64_t>(exchange.send_wire.size());
      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoallv(
                nonNullPmWireData(exchange.recv_wire),
                exchange.send_response_counts_bytes.data(),
                exchange.send_response_displs_bytes.data(),
                MPI_BYTE,
                nonNullPmWireData(exchange.send_wire),
                exchange.recv_response_counts_bytes.data(),
                exchange.recv_response_displs_bytes.data(),
                MPI_BYTE,
                MPI_COMM_WORLD),
            "PmSolver::interpolateForces response MPI_Alltoallv");
      });

      std::string local_force_response_validation_error;
      try {
        // Re-generate the compact request stream instead of retaining an O(N)
        // request registry. This provides exact response identity validation with
        // only O(world_size) cursor state.
        for (int rank = 0; rank < world_size; ++rank) {
          exchange.cursor[static_cast<std::size_t>(rank)] = static_cast<std::size_t>(
              exchange.send_displs[static_cast<std::size_t>(rank)]);
        }
        std::uint32_t expected_sequence = round_sequence_begin;
        for (std::size_t p = begin; p < end; ++p) {
          if (!periodic &&
              (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
               !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
               !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
            continue;
          }
          const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
          const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
          std::array<PmXPlaneGroup, 3> groups{};
          const std::size_t group_count = makePmXPlaneGroups(
              sx, periodic, m_shape, decomposition_view, groups);
          for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
            const int sender_rank = groups[group_i].destination_rank;
            if (sender_rank == world_rank) {
              continue;
            }
            auto& cursor = exchange.cursor[static_cast<std::size_t>(sender_rank)];
            const std::size_t segment_end = static_cast<std::size_t>(
                exchange.send_displs[static_cast<std::size_t>(sender_rank)] +
                exchange.send_counts[static_cast<std::size_t>(sender_rank)]);
            if (cursor >= segment_end) {
              throw std::invalid_argument("PM force response segment ended before issued requests");
            }
            const PmForceContributionRecord response = decodePmForceResponse(
                std::span<const std::uint8_t>(exchange.send_wire).subspan(
                    cursor * k_pm_force_response_wire_bytes,
                    k_pm_force_response_wire_bytes));
            if (response.source_rank != static_cast<std::uint32_t>(sender_rank) ||
                response.origin_rank != static_cast<std::uint32_t>(world_rank) ||
                response.exchange_epoch != exchange_epoch ||
                response.request_sequence != expected_sequence ||
                response.origin_particle_index != static_cast<std::uint32_t>(p) ||
                !std::isfinite(response.accel_x) ||
                !std::isfinite(response.accel_y) ||
                !std::isfinite(response.accel_z)) {
              throw std::invalid_argument(
                  "PM force response does not match the regenerated bounded request stream");
            }
            accel_x[p] += response.accel_x;
            accel_y[p] += response.accel_y;
            accel_z[p] += response.accel_z;
            ++cursor;
            ++expected_sequence;
          }
        }
        if (expected_sequence != next_request_sequence) {
          throw std::logic_error("PM force response validation sequence accounting drifted");
        }
        for (int rank = 0; rank < world_size; ++rank) {
          const std::size_t expected_end = static_cast<std::size_t>(
              exchange.send_displs[static_cast<std::size_t>(rank)] +
              exchange.send_counts[static_cast<std::size_t>(rank)]);
          if (exchange.cursor[static_cast<std::size_t>(rank)] != expected_end) {
            throw std::invalid_argument("PM force response segment contains unconsumed records");
          }
        }
      } catch (const std::exception& error) {
        local_force_response_validation_error = error.what();
      } catch (...) {
        local_force_response_validation_error =
            "non-standard exception while validating bounded force responses";
      }
      throwIfPmPayloadValidationFailed(
          "PmSolver::interpolateForces bounded response receive",
          world_rank,
          world_size,
          local_force_response_validation_error,
          routed_mpi_wait_ms);

      logical_requests_sent += static_cast<std::uint64_t>(total_sent_requests);
      peers_touched += static_cast<std::uint64_t>(std::count_if(
          exchange.send_counts.begin(), exchange.send_counts.end(), [](int count) { return count > 0; }));
      total_wire_sent += request_send_bytes + response_send_bytes;
      total_wire_received += request_recv_bytes + response_recv_bytes;
      const std::uint64_t send_wire_capacity =
          static_cast<std::uint64_t>(exchange.send_wire.capacity());
      const std::uint64_t recv_wire_capacity =
          static_cast<std::uint64_t>(exchange.recv_wire.capacity());
      // The two physical buffers swap direction for the response phase, so
      // each directional high-water must consider both retained capacities.
      const std::uint64_t directional_capacity =
          std::max(send_wire_capacity, recv_wire_capacity);
      const std::uint64_t combined_capacity = checkedAddBytes(
          send_wire_capacity,
          recv_wire_capacity,
          "PM force combined routing buffer capacity");
      const std::uint64_t workspace_capacity = checkedAddBytes(
          routing_metadata_capacity,
          combined_capacity,
          "PM force routing workspace capacity");
      if (workspace_capacity > k_pm_routing_workspace_target_bytes) {
        throw std::logic_error("PM force routing workspace exceeded the M1A aggregate target");
      }
      send_high_water = std::max(send_high_water, directional_capacity);
      receive_high_water = std::max(receive_high_water, directional_capacity);
      combined_high_water = std::max(combined_high_water, combined_capacity);
      workspace_high_water = std::max(workspace_high_water, workspace_capacity);
      exchange.workspace_high_water_bytes = std::max(
          exchange.workspace_high_water_bytes, workspace_capacity);
    }

    if (profile != nullptr) {
      profile->routed_force_requests += logical_requests_sent;
      profile->routed_force_peer_count += peers_touched;
      profile->routed_mpi_bytes_sent += total_wire_sent;
      profile->routed_mpi_bytes_received += total_wire_received;
      profile->routed_send_buffer_high_water_bytes = std::max(
          profile->routed_send_buffer_high_water_bytes, send_high_water);
      profile->routed_receive_buffer_high_water_bytes = std::max(
          profile->routed_receive_buffer_high_water_bytes, receive_high_water);
      profile->routed_combined_buffer_high_water_bytes = std::max(
          profile->routed_combined_buffer_high_water_bytes, combined_high_water);
      profile->routed_workspace_high_water_bytes = std::max(
          profile->routed_workspace_high_water_bytes, workspace_high_water);
      profile->routed_mpi_wait_ms += routed_mpi_wait_ms;
      profile->bytes_moved += total_wire_sent + total_wire_received;
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
  const PmDecompositionView decomposition_view(
      m_shape, grid.slabLayout(), options.decomposition_mode);
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
    const int world_rank = grid.slabLayout().world_rank;
    double routed_mpi_wait_ms = distributed_preflight_mpi_wait_ms;
    const int local_targets_valid = first_non_finite_target_index.has_value() ? 0 : 1;
    int all_targets_valid = 0;
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      requirePmMpiSuccess(
          MPI_Allreduce(&local_targets_valid, &all_targets_valid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD),
          "PmSolver::interpolatePotential target validation MPI_Allreduce");
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
        world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() { exchange_epoch = m_impl->nextDistributedExchangeEpoch(); });
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      validatePmExchangeEpochConsensus(exchange_epoch, "PmSolver::interpolatePotential");
    });

    Impl::PlaneInterpolationExchangeBuffers* exchange_ptr = nullptr;
    std::uint64_t routing_metadata_capacity = 0U;
    PmRoutingCapacityModel routing_capacity{};
    std::size_t routing_buffer_limit = 0U;
    std::size_t particles_per_round = 0U;
    runPmCoordinatedPhase(
        "PmSolver::interpolatePotential routing preflight",
        world_rank,
        world_size,
        routed_mpi_wait_ms,
        [&]() {
          exchange_ptr = &m_impl->planeInterpolationExchangeBuffersForLayout(grid.slabLayout());
          auto& prepared_exchange = *exchange_ptr;
          const auto add_routing_metadata = [&](const auto& values) {
            routing_metadata_capacity = checkedAddBytes(
                routing_metadata_capacity,
                core::ownedCapacityBytesForContainer(values),
                "PM potential routing metadata capacity");
          };
          add_routing_metadata(prepared_exchange.send_counts);
          add_routing_metadata(prepared_exchange.send_displs);
          add_routing_metadata(prepared_exchange.recv_counts);
          add_routing_metadata(prepared_exchange.recv_displs);
          add_routing_metadata(prepared_exchange.send_counts_bytes);
          add_routing_metadata(prepared_exchange.send_displs_bytes);
          add_routing_metadata(prepared_exchange.recv_counts_bytes);
          add_routing_metadata(prepared_exchange.recv_displs_bytes);
          add_routing_metadata(prepared_exchange.send_response_counts_bytes);
          add_routing_metadata(prepared_exchange.send_response_displs_bytes);
          add_routing_metadata(prepared_exchange.recv_response_counts_bytes);
          add_routing_metadata(prepared_exchange.recv_response_displs_bytes);
          add_routing_metadata(prepared_exchange.cursor);
          routing_capacity = modelPmRoutingCapacity(
              world_size,
              options.routing_exchange_batch_bytes,
              k_pm_routing_modeled_workspace_limit_bytes,
              routing_metadata_capacity,
              k_pm_plane_interpolation_request_wire_bytes);
          if (routing_capacity.max_send_payload_bytes >
              static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
            throw std::overflow_error("PM potential routing buffer limit exceeds size_t capacity");
          }
          routing_buffer_limit =
              static_cast<std::size_t>(routing_capacity.max_send_payload_bytes);
          particles_per_round = pmRoutingParticlesPerRound(
              routing_capacity.effective_per_peer_payload_bytes,
              k_pm_plane_interpolation_request_wire_bytes);
        });
    auto& exchange = *exchange_ptr;
    const std::uint64_t local_target_count = static_cast<std::uint64_t>(pos_x.size());
    std::uint64_t global_max_target_count = 0U;
    measurePmMpiWait(routed_mpi_wait_ms, [&]() {
      requirePmMpiSuccess(
          MPI_Allreduce(
              &local_target_count,
              &global_max_target_count,
              1,
              MPI_UINT64_T,
              MPI_MAX,
              MPI_COMM_WORLD),
          "PmSolver::interpolatePotential max-target MPI_Allreduce");
    });
    const std::uint64_t round_capacity = static_cast<std::uint64_t>(particles_per_round);
    const std::uint64_t round_count = global_max_target_count == 0U
        ? 0U
        : 1U + (global_max_target_count - 1U) / round_capacity;
    const bool periodic = options.boundary_condition == PmBoundaryCondition::kPeriodic;
    std::uint32_t next_request_sequence = 0U;

    std::fill(potential.begin(), potential.end(), 0.0);

    const auto evaluate_plane_request = [&](const PmPlaneInterpolationRequestRecord& request) {
      if (request.plane_count == 0U || request.plane_count > 3U ||
          request.destination_rank != static_cast<std::uint32_t>(world_rank) ||
          request.exchange_epoch != exchange_epoch ||
          !std::isfinite(request.y_grid) || !std::isfinite(request.z_grid)) {
        throw std::invalid_argument(
            "PmSolver::interpolatePotential received invalid bounded potential-plane request");
      }
      const PmAxisStencil1d sy = makeAxisStencil(request.y_grid, options.assignment_scheme);
      const PmAxisStencil1d sz = makeAxisStencil(request.z_grid, options.assignment_scheme);
      PmPotentialContributionRecord response{
          .source_rank = static_cast<std::uint32_t>(world_rank),
          .origin_rank = request.origin_rank,
          .request_sequence = request.request_sequence,
          .origin_particle_index = request.origin_particle_index,
          .exchange_epoch = request.exchange_epoch,
          .potential = 0.0,
      };
      for (std::size_t plane = 0; plane < static_cast<std::size_t>(request.plane_count); ++plane) {
        const std::size_t ix = static_cast<std::size_t>(request.global_ix[plane]);
        if (ix >= m_shape.nx || !grid.slabLayout().ownsGlobalX(ix) ||
            !std::isfinite(request.x_weight[plane])) {
          throw std::invalid_argument(
              "PmSolver::interpolatePotential potential-plane x ownership/weight is invalid");
        }
        for (std::size_t dy = 0; dy < sy.count; ++dy) {
          if (!periodic &&
              (sy.offsets[dy] < 0 || sy.offsets[dy] >= static_cast<std::ptrdiff_t>(m_shape.ny))) {
            continue;
          }
          const std::size_t iy = periodic
              ? wrapIndex(sy.offsets[dy], m_shape.ny)
              : static_cast<std::size_t>(sy.offsets[dy]);
          for (std::size_t dz = 0; dz < sz.count; ++dz) {
            if (!periodic &&
                (sz.offsets[dz] < 0 || sz.offsets[dz] >= static_cast<std::ptrdiff_t>(m_shape.nz))) {
              continue;
            }
            const std::size_t iz = periodic
                ? wrapIndex(sz.offsets[dz], m_shape.nz)
                : static_cast<std::size_t>(sz.offsets[dz]);
            const std::size_t local_index = decomposition_view.globalToLocalRealIndex(ix, iy, iz);
            if (!std::isfinite(grid.potential()[local_index])) {
              throw std::invalid_argument(
                  "PmSolver::interpolatePotential encountered non-finite owner-local mesh potential data");
            }
            const double weight = request.x_weight[plane] * sy.weights[dy] * sz.weights[dz];
            response.potential += weight * grid.potential()[local_index];
          }
        }
      }
      return response;
    };

    std::uint64_t logical_requests_sent = 0U;
    std::uint64_t peers_touched = 0U;
    std::uint64_t total_wire_sent = 0U;
    std::uint64_t total_wire_received = 0U;
    std::uint64_t send_high_water = 0U;
    std::uint64_t receive_high_water = 0U;
    std::uint64_t combined_high_water = 0U;
    std::uint64_t workspace_high_water = routing_metadata_capacity;

    for (std::uint64_t round = 0U; round < round_count; ++round) {
      const std::uint64_t begin_u64 = round * round_capacity;
      const std::uint64_t end_u64 = std::min(
          local_target_count,
          begin_u64 > std::numeric_limits<std::uint64_t>::max() - round_capacity
              ? local_target_count
              : begin_u64 + round_capacity);
      const std::size_t begin = static_cast<std::size_t>(begin_u64);
      const std::size_t end = static_cast<std::size_t>(end_u64);
      const std::uint32_t round_sequence_begin = next_request_sequence;

      std::fill(exchange.send_counts.begin(), exchange.send_counts.end(), 0);
      std::fill(exchange.send_displs.begin(), exchange.send_displs.end(), 0);
      std::fill(exchange.recv_counts.begin(), exchange.recv_counts.end(), 0);
      std::fill(exchange.recv_displs.begin(), exchange.recv_displs.end(), 0);
      std::fill(exchange.send_counts_bytes.begin(), exchange.send_counts_bytes.end(), 0);
      std::fill(exchange.send_displs_bytes.begin(), exchange.send_displs_bytes.end(), 0);
      std::fill(exchange.recv_counts_bytes.begin(), exchange.recv_counts_bytes.end(), 0);
      std::fill(exchange.recv_displs_bytes.begin(), exchange.recv_displs_bytes.end(), 0);
      std::fill(exchange.send_response_counts_bytes.begin(), exchange.send_response_counts_bytes.end(), 0);
      std::fill(exchange.send_response_displs_bytes.begin(), exchange.send_response_displs_bytes.end(), 0);
      std::fill(exchange.recv_response_counts_bytes.begin(), exchange.recv_response_counts_bytes.end(), 0);
      std::fill(exchange.recv_response_displs_bytes.begin(), exchange.recv_response_displs_bytes.end(), 0);
      exchange.send_wire.clear();
      exchange.recv_wire.clear();

      runPmCoordinatedPhase(
          "PmSolver::interpolatePotential bounded request preparation",
          world_rank,
          world_size,
          routed_mpi_wait_ms,
          [&]() {
            // Pass 1 counts at most one request per particle/destination x-plane group.
            for (std::size_t p = begin; p < end; ++p) {
              if (p > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
                throw std::invalid_argument(
                    "PmSolver::interpolatePotential origin particle index exceeds routing token limit");
              }
              if (!periodic &&
                  (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
                   !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
                   !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
                continue;
              }
              const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
              const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
              std::array<PmXPlaneGroup, 3> groups{};
              const std::size_t group_count = makePmXPlaneGroups(
                  sx, periodic, m_shape, decomposition_view, groups);
              for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
                const int destination_rank = groups[group_i].destination_rank;
                if (destination_rank == world_rank) {
                  continue;
                }
                int& count = exchange.send_counts[static_cast<std::size_t>(destination_rank)];
                if (count == std::numeric_limits<int>::max()) {
                  throw std::overflow_error("PM potential per-peer request count exceeds MPI int capacity");
                }
                ++count;
              }
            }

            std::size_t total_send_records = 0U;
            for (int rank = 0; rank < world_size; ++rank) {
              const std::size_t count = static_cast<std::size_t>(
                  exchange.send_counts[static_cast<std::size_t>(rank)]);
              const std::uint64_t peer_bytes = checkedBytesForCount(
                  count, k_pm_plane_interpolation_request_wire_bytes, "PM potential per-peer routing batch");
              if (peer_bytes > routing_capacity.effective_per_peer_payload_bytes) {
                throw std::runtime_error(
                    "PM potential routing exceeded effective per-peer batch high-water");
              }
              exchange.send_displs[static_cast<std::size_t>(rank)] =
                  checkedMpiRecordCountOrDisplacement(
                      total_send_records,
                      k_pm_plane_interpolation_request_wire_bytes,
                      "PmSolver::interpolatePotential bounded request exchange",
                      "send record displacement",
                      rank);
              total_send_records = checkedMpiRecordTotal(
                  total_send_records,
                  count,
                  k_pm_plane_interpolation_request_wire_bytes,
                  "PmSolver::interpolatePotential bounded request exchange",
                  rank);
            }
            resizePmWireBufferBounded(
                exchange.send_wire,
                checkedPmWireByteCount(
                    total_send_records,
                    k_pm_plane_interpolation_request_wire_bytes,
                    "PmSolver::interpolatePotential bounded request exchange"),
                routing_buffer_limit,
                "PmSolver::interpolatePotential request send buffer");

            for (int rank = 0; rank < world_size; ++rank) {
              exchange.cursor[static_cast<std::size_t>(rank)] = static_cast<std::size_t>(
                  exchange.send_displs[static_cast<std::size_t>(rank)]);
            }

            // Pass 2 evaluates local x planes immediately and encodes only remote
            // x-plane groups. The receiver expands the y/z stencil and returns one
            // accumulated potential response per request.
            for (std::size_t p = begin; p < end; ++p) {
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
              std::array<PmXPlaneGroup, 3> groups{};
              const std::size_t group_count = makePmXPlaneGroups(
                  sx, periodic, m_shape, decomposition_view, groups);
              for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
                const auto& group = groups[group_i];
                PmPlaneInterpolationRequestRecord request{
                    .origin_rank = static_cast<std::uint32_t>(world_rank),
                    .destination_rank = static_cast<std::uint32_t>(group.destination_rank),
                    .request_sequence = next_request_sequence,
                    .origin_particle_index = static_cast<std::uint32_t>(p),
                    .plane_count = group.plane_count,
                    .global_ix = group.global_ix,
                    .exchange_epoch = exchange_epoch,
                    .y_grid = y,
                    .z_grid = z,
                    .x_weight = group.x_weight,
                };
                if (group.destination_rank == world_rank) {
                  const PmPotentialContributionRecord local = evaluate_plane_request(request);
                  potential[p] += local.potential;
                  continue;
                }
                if (next_request_sequence == std::numeric_limits<std::uint32_t>::max()) {
                  throw std::overflow_error(
                      "PmSolver::interpolatePotential routing sequence exceeds uint32 capacity");
                }
                auto& cursor = exchange.cursor[static_cast<std::size_t>(group.destination_rank)];
                encodePmPlaneInterpolationRequest(
                    request,
                    PmWireRecordKind::kPotentialPlaneRequest,
                    std::span<std::uint8_t>(exchange.send_wire).subspan(
                        cursor * k_pm_plane_interpolation_request_wire_bytes,
                        k_pm_plane_interpolation_request_wire_bytes));
                ++cursor;
                ++next_request_sequence;
              }
            }
          });

      const std::uint64_t request_send_bytes = static_cast<std::uint64_t>(exchange.send_wire.size());
      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoall(
                exchange.send_counts.data(), 1, MPI_INT,
                exchange.recv_counts.data(), 1, MPI_INT,
                MPI_COMM_WORLD),
            "PmSolver::interpolatePotential request-count MPI_Alltoall");
      });

      runPmCoordinatedPhase(
          "PmSolver::interpolatePotential bounded request receive layout",
          world_rank,
          world_size,
          routed_mpi_wait_ms,
          [&]() {
            std::size_t total_recv_records = 0U;
            for (int rank = 0; rank < world_size; ++rank) {
              const std::size_t count = checkedMpiReceivedRecordCount(
                  exchange.recv_counts[static_cast<std::size_t>(rank)],
                  k_pm_plane_interpolation_request_wire_bytes,
                  "PmSolver::interpolatePotential bounded request exchange",
                  rank);
              exchange.recv_displs[static_cast<std::size_t>(rank)] =
                  checkedMpiRecordCountOrDisplacement(
                      total_recv_records,
                      k_pm_plane_interpolation_request_wire_bytes,
                      "PmSolver::interpolatePotential bounded request exchange",
                      "received record displacement",
                      rank);
              total_recv_records = checkedMpiRecordTotal(
                  total_recv_records,
                  count,
                  k_pm_plane_interpolation_request_wire_bytes,
                  "PmSolver::interpolatePotential bounded request exchange",
                  rank);
            }
            resizePmWireBufferBounded(
                exchange.recv_wire,
                checkedPmWireByteCount(
                    total_recv_records,
                    k_pm_plane_interpolation_request_wire_bytes,
                    "PmSolver::interpolatePotential bounded request exchange"),
                routing_buffer_limit,
                "PmSolver::interpolatePotential request receive buffer");
            checkedMpiRecordLayoutToByteLayout(
                exchange.send_counts,
                exchange.send_displs,
                k_pm_plane_interpolation_request_wire_bytes,
                exchange.send_counts_bytes,
                exchange.send_displs_bytes,
                "PmSolver::interpolatePotential bounded request send MPI_Alltoallv");
            checkedMpiRecordLayoutToByteLayout(
                exchange.recv_counts,
                exchange.recv_displs,
                k_pm_plane_interpolation_request_wire_bytes,
                exchange.recv_counts_bytes,
                exchange.recv_displs_bytes,
                "PmSolver::interpolatePotential bounded request receive MPI_Alltoallv");
          });
      const std::uint64_t request_recv_bytes = static_cast<std::uint64_t>(exchange.recv_wire.size());

      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoallv(
                nonNullPmWireData(exchange.send_wire),
                exchange.send_counts_bytes.data(),
                exchange.send_displs_bytes.data(),
                MPI_BYTE,
                nonNullPmWireData(exchange.recv_wire),
                exchange.recv_counts_bytes.data(),
                exchange.recv_displs_bytes.data(),
                MPI_BYTE,
                MPI_COMM_WORLD),
            "PmSolver::interpolatePotential request MPI_Alltoallv");
      });

      std::string local_potential_request_validation_error;
      try {
        // Compact responses in-place into the receive-request buffer. Since a
        // response (48 B) is smaller than a request (96 B), forward compaction
        // cannot overwrite a request that has not yet been decoded.
        for (int source_rank = 0; source_rank < world_size; ++source_rank) {
          const int request_begin = exchange.recv_displs[static_cast<std::size_t>(source_rank)];
          const int request_count = exchange.recv_counts[static_cast<std::size_t>(source_rank)];
          exchange.send_response_counts_bytes[static_cast<std::size_t>(source_rank)] =
              checkedMpiByteCount(
                  checkedBytesForCount(
                      static_cast<std::size_t>(request_count),
                      k_pm_potential_response_wire_bytes,
                      "PM potential response sender count"),
                  "PM potential response sender count");
          exchange.send_response_displs_bytes[static_cast<std::size_t>(source_rank)] =
              checkedMpiByteCount(
                  checkedBytesForCount(
                      static_cast<std::size_t>(request_begin),
                      k_pm_potential_response_wire_bytes,
                      "PM potential response sender displacement"),
                  "PM potential response sender displacement");
          std::optional<std::uint32_t> previous_sequence;
          for (int i = 0; i < request_count; ++i) {
            const std::size_t request_index = static_cast<std::size_t>(request_begin + i);
            const PmPlaneInterpolationRequestRecord request = decodePmPlaneInterpolationRequest(
                std::span<const std::uint8_t>(exchange.recv_wire).subspan(
                    request_index * k_pm_plane_interpolation_request_wire_bytes,
                    k_pm_plane_interpolation_request_wire_bytes),
                PmWireRecordKind::kPotentialPlaneRequest,
                "PmSolver::interpolatePotential potential-plane request");
            if (request.origin_rank != static_cast<std::uint32_t>(source_rank) ||
                request.destination_rank != static_cast<std::uint32_t>(world_rank) ||
                request.exchange_epoch != exchange_epoch) {
              throw std::invalid_argument(
                  "PM potential-plane routing metadata does not match MPI sender/receiver segment");
            }
            if (previous_sequence.has_value() && request.request_sequence <= *previous_sequence) {
              throw std::invalid_argument(
                  "PM potential-plane request sequence is non-monotonic within sender segment");
            }
            previous_sequence = request.request_sequence;
            const PmPotentialContributionRecord response = evaluate_plane_request(request);
            const std::size_t response_index =
                static_cast<std::size_t>(exchange.send_response_displs_bytes[static_cast<std::size_t>(source_rank)]) /
                    k_pm_potential_response_wire_bytes +
                static_cast<std::size_t>(i);
            encodePmPotentialResponse(
                response,
                std::span<std::uint8_t>(exchange.recv_wire).subspan(
                    response_index * k_pm_potential_response_wire_bytes,
                    k_pm_potential_response_wire_bytes));
          }
        }
      } catch (const std::exception& error) {
        local_potential_request_validation_error = error.what();
      } catch (...) {
        local_potential_request_validation_error =
            "non-standard exception while consuming bounded potential-plane requests";
      }
      throwIfPmPayloadValidationFailed(
          "PmSolver::interpolatePotential bounded request receive",
          world_rank,
          world_size,
          local_potential_request_validation_error,
          routed_mpi_wait_ms);

      const std::size_t total_received_requests = exchange.recv_wire.size() /
          k_pm_plane_interpolation_request_wire_bytes;
      resizePmWireBufferBounded(
          exchange.recv_wire,
          checkedPmWireByteCount(
              total_received_requests,
              k_pm_potential_response_wire_bytes,
              "PmSolver::interpolatePotential compact response send"),
          routing_buffer_limit,
          "PmSolver::interpolatePotential response send buffer");

      const std::size_t total_sent_requests = exchange.send_wire.size() /
          k_pm_plane_interpolation_request_wire_bytes;
      for (int rank = 0; rank < world_size; ++rank) {
        exchange.recv_response_counts_bytes[static_cast<std::size_t>(rank)] =
            checkedMpiByteCount(
                checkedBytesForCount(
                    static_cast<std::size_t>(exchange.send_counts[static_cast<std::size_t>(rank)]),
                    k_pm_potential_response_wire_bytes,
                    "PM potential response receive count"),
                "PM potential response receive count");
        exchange.recv_response_displs_bytes[static_cast<std::size_t>(rank)] =
            checkedMpiByteCount(
                checkedBytesForCount(
                    static_cast<std::size_t>(exchange.send_displs[static_cast<std::size_t>(rank)]),
                    k_pm_potential_response_wire_bytes,
                    "PM potential response receive displacement"),
                "PM potential response receive displacement");
      }
      resizePmWireBufferBounded(
          exchange.send_wire,
          checkedPmWireByteCount(
              total_sent_requests,
              k_pm_potential_response_wire_bytes,
              "PmSolver::interpolatePotential compact response receive"),
          routing_buffer_limit,
          "PmSolver::interpolatePotential response receive buffer");

      const std::uint64_t response_send_bytes = static_cast<std::uint64_t>(exchange.recv_wire.size());
      const std::uint64_t response_recv_bytes = static_cast<std::uint64_t>(exchange.send_wire.size());
      measurePmMpiWait(routed_mpi_wait_ms, [&]() {
        requirePmMpiSuccess(
            MPI_Alltoallv(
                nonNullPmWireData(exchange.recv_wire),
                exchange.send_response_counts_bytes.data(),
                exchange.send_response_displs_bytes.data(),
                MPI_BYTE,
                nonNullPmWireData(exchange.send_wire),
                exchange.recv_response_counts_bytes.data(),
                exchange.recv_response_displs_bytes.data(),
                MPI_BYTE,
                MPI_COMM_WORLD),
            "PmSolver::interpolatePotential response MPI_Alltoallv");
      });

      std::string local_potential_response_validation_error;
      try {
        // Re-generate the compact request stream instead of retaining an O(N)
        // request registry. This provides exact response identity validation with
        // only O(world_size) cursor state.
        for (int rank = 0; rank < world_size; ++rank) {
          exchange.cursor[static_cast<std::size_t>(rank)] = static_cast<std::size_t>(
              exchange.send_displs[static_cast<std::size_t>(rank)]);
        }
        std::uint32_t expected_sequence = round_sequence_begin;
        for (std::size_t p = begin; p < end; ++p) {
          if (!periodic &&
              (!positionInsideOpenDomain(pos_x[p], lengths.lx) ||
               !positionInsideOpenDomain(pos_y[p], lengths.ly) ||
               !positionInsideOpenDomain(pos_z[p], lengths.lz))) {
            continue;
          }
          const double x = (periodic ? wrapPosition(pos_x[p], lengths.lx) : pos_x[p]) * inv_dx;
          const PmAxisStencil1d sx = makeAxisStencil(x, options.assignment_scheme);
          std::array<PmXPlaneGroup, 3> groups{};
          const std::size_t group_count = makePmXPlaneGroups(
              sx, periodic, m_shape, decomposition_view, groups);
          for (std::size_t group_i = 0; group_i < group_count; ++group_i) {
            const int sender_rank = groups[group_i].destination_rank;
            if (sender_rank == world_rank) {
              continue;
            }
            auto& cursor = exchange.cursor[static_cast<std::size_t>(sender_rank)];
            const std::size_t segment_end = static_cast<std::size_t>(
                exchange.send_displs[static_cast<std::size_t>(sender_rank)] +
                exchange.send_counts[static_cast<std::size_t>(sender_rank)]);
            if (cursor >= segment_end) {
              throw std::invalid_argument("PM potential response segment ended before issued requests");
            }
            const PmPotentialContributionRecord response = decodePmPotentialResponse(
                std::span<const std::uint8_t>(exchange.send_wire).subspan(
                    cursor * k_pm_potential_response_wire_bytes,
                    k_pm_potential_response_wire_bytes));
            if (response.source_rank != static_cast<std::uint32_t>(sender_rank) ||
                response.origin_rank != static_cast<std::uint32_t>(world_rank) ||
                response.exchange_epoch != exchange_epoch ||
                response.request_sequence != expected_sequence ||
                response.origin_particle_index != static_cast<std::uint32_t>(p) ||
                !std::isfinite(response.potential)) {
              throw std::invalid_argument(
                  "PM potential response does not match the regenerated bounded request stream");
            }
            potential[p] += response.potential;
            ++cursor;
            ++expected_sequence;
          }
        }
        if (expected_sequence != next_request_sequence) {
          throw std::logic_error("PM potential response validation sequence accounting drifted");
        }
        for (int rank = 0; rank < world_size; ++rank) {
          const std::size_t expected_end = static_cast<std::size_t>(
              exchange.send_displs[static_cast<std::size_t>(rank)] +
              exchange.send_counts[static_cast<std::size_t>(rank)]);
          if (exchange.cursor[static_cast<std::size_t>(rank)] != expected_end) {
            throw std::invalid_argument("PM potential response segment contains unconsumed records");
          }
        }
      } catch (const std::exception& error) {
        local_potential_response_validation_error = error.what();
      } catch (...) {
        local_potential_response_validation_error =
            "non-standard exception while validating bounded potential responses";
      }
      throwIfPmPayloadValidationFailed(
          "PmSolver::interpolatePotential bounded response receive",
          world_rank,
          world_size,
          local_potential_response_validation_error,
          routed_mpi_wait_ms);

      logical_requests_sent += static_cast<std::uint64_t>(total_sent_requests);
      peers_touched += static_cast<std::uint64_t>(std::count_if(
          exchange.send_counts.begin(), exchange.send_counts.end(), [](int count) { return count > 0; }));
      total_wire_sent += request_send_bytes + response_send_bytes;
      total_wire_received += request_recv_bytes + response_recv_bytes;
      const std::uint64_t send_wire_capacity =
          static_cast<std::uint64_t>(exchange.send_wire.capacity());
      const std::uint64_t recv_wire_capacity =
          static_cast<std::uint64_t>(exchange.recv_wire.capacity());
      // The two physical buffers swap direction for the response phase, so
      // each directional high-water must consider both retained capacities.
      const std::uint64_t directional_capacity =
          std::max(send_wire_capacity, recv_wire_capacity);
      const std::uint64_t combined_capacity = checkedAddBytes(
          send_wire_capacity,
          recv_wire_capacity,
          "PM potential combined routing buffer capacity");
      const std::uint64_t workspace_capacity = checkedAddBytes(
          routing_metadata_capacity,
          combined_capacity,
          "PM potential routing workspace capacity");
      if (workspace_capacity > k_pm_routing_workspace_target_bytes) {
        throw std::logic_error("PM potential routing workspace exceeded the M1A aggregate target");
      }
      send_high_water = std::max(send_high_water, directional_capacity);
      receive_high_water = std::max(receive_high_water, directional_capacity);
      combined_high_water = std::max(combined_high_water, combined_capacity);
      workspace_high_water = std::max(workspace_high_water, workspace_capacity);
      exchange.workspace_high_water_bytes = std::max(
          exchange.workspace_high_water_bytes, workspace_capacity);
    }

    if (profile != nullptr) {
      profile->routed_potential_requests += logical_requests_sent;
      profile->routed_potential_peer_count += peers_touched;
      profile->routed_mpi_bytes_sent += total_wire_sent;
      profile->routed_mpi_bytes_received += total_wire_received;
      profile->routed_send_buffer_high_water_bytes = std::max(
          profile->routed_send_buffer_high_water_bytes, send_high_water);
      profile->routed_receive_buffer_high_water_bytes = std::max(
          profile->routed_receive_buffer_high_water_bytes, receive_high_water);
      profile->routed_combined_buffer_high_water_bytes = std::max(
          profile->routed_combined_buffer_high_water_bytes, combined_high_water);
      profile->routed_workspace_high_water_bytes = std::max(
          profile->routed_workspace_high_water_bytes, workspace_high_water);
      profile->routed_mpi_wait_ms += routed_mpi_wait_ms;
      profile->bytes_moved += total_wire_sent + total_wire_received;
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
        solve_preflight_mpi_wait_ms);
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
  // Backend-independent API contract. Keep this outside the MPI branch so a
  // single-rank CUDA solve cannot bypass the source/target shape checks that
  // the host and distributed paths rely on. Multi-rank callers have already
  // passed the coordinated equivalent above, making these checks deterministic.
  validateOptions(m_shape, options);
  if (grid.shape().nx != m_shape.nx || grid.shape().ny != m_shape.ny ||
      grid.shape().nz != m_shape.nz) {
    throw std::invalid_argument("PM solver/grid shape mismatch in solveForParticles");
  }
  if (!grid.slabLayout().isValid()) {
    throw std::invalid_argument("PmSolver::solveForParticles requires a valid PM slab layout");
  }
  if (pos_x.size() != pos_y.size() || pos_x.size() != pos_z.size() ||
      pos_x.size() != mass.size() || pos_x.size() != accel_x.size() ||
      pos_x.size() != accel_y.size() || pos_x.size() != accel_z.size()) {
    throw std::invalid_argument(
        "Particle coordinate/mass/acceleration spans must match in solveForParticles");
  }
#if !COSMOSIM_ENABLE_MPI
  if (grid.slabLayout().world_size != 1) {
    throw std::invalid_argument(
        "PmSolver::solveForParticles distributed slabs require COSMOSIM_ENABLE_MPI=ON");
  }
#endif
  if (grid.slabLayout().world_size == 1) {
    validateSingleRankFullDomainGridContract(grid, "PmSolver::solveForParticles");
  }
#if !COSMOSIM_ENABLE_CUDA
  if (options.execution_policy == core::ExecutionPolicy::kCuda) {
    throw std::runtime_error(
        "PM solve requested execution_policy=cuda, but this build has COSMOSIM_ENABLE_CUDA=OFF");
  }
#endif

  if (options.execution_policy == core::ExecutionPolicy::kCuda) {
#if COSMOSIM_ENABLE_CUDA
    m_impl->solveCudaStagingBackend(
        *this, grid, pos_x, pos_y, pos_z, mass, accel_x, accel_y, accel_z,
        options, profile);
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
  return pmBackendProductionReady();
}

PmBackendCapability pmBackendCapability() noexcept {
#if COSMOSIM_ENABLE_FFTW
  return PmBackendCapability::kProductionFftw;
#else
  return PmBackendCapability::kDiagnosticNaiveDft;
#endif
}

std::string_view pmBackendCapabilityName(PmBackendCapability capability) noexcept {
  switch (capability) {
    case PmBackendCapability::kUnavailable:
      return "unavailable";
    case PmBackendCapability::kDiagnosticNaiveDft:
      return "diagnostic_naive_dft";
    case PmBackendCapability::kProductionFftw:
      return "production_fftw";
    case PmBackendCapability::kProductionCudaFft:
      return "production_cuda_fft";
  }
  return "unknown";
}

bool pmBackendProductionReady() noexcept {
  const PmBackendCapability capability = pmBackendCapability();
  return capability == PmBackendCapability::kProductionFftw ||
      capability == PmBackendCapability::kProductionCudaFft;
}

PmBackendArchitecture pmBackendArchitecture() noexcept {
  PmBackendArchitecture architecture;
  architecture.capability = pmBackendCapability();
#if COSMOSIM_ENABLE_FFTW
  architecture.host_fft = true;
#endif
#if COSMOSIM_ENABLE_CUDA
  // The present CUDA path performs assignment/interpolation on device but
  // returns through the host FFT solver. Keep the missing resident stages
  // explicit so provenance/capability code cannot over-promote it.
  architecture.device_assignment = true;
  architecture.device_interpolation = true;
  architecture.persistent_device_buffers = true;
  // The staging backend now retains its stream and device workspaces across
  // solves. Full persistent residency remains false because density/force still
  // cross the host boundary for FFT/Poisson/gradient.
  architecture.persistent_residency = false;
#endif
  return architecture;
}

std::string_view pmDecompositionTopologyName(PmDecompositionTopology topology) noexcept {
  switch (topology) {
    case PmDecompositionTopology::kSerial:
      return "serial";
    case PmDecompositionTopology::kXSlab:
      return "x_slab";
    case PmDecompositionTopology::kXSlabTransposedSpectral:
      return "x_slab_transposed_spectral";
    case PmDecompositionTopology::kPencil2D:
      return "pencil_2d";
  }
  return "unknown";
}

PmDecompositionView::PmDecompositionView(
    PmGridShape shape,
    parallel::PmSlabLayout layout,
    core::PmDecompositionMode mode)
    : m_shape(shape),
      m_layout(std::move(layout)),
      m_mode(mode),
      m_descriptor(describePmDecomposition(m_shape, m_layout, m_mode)),
      m_world_rank(m_layout.world_rank),
      m_world_size(m_layout.world_size) {}

PmDecompositionView::PmDecompositionView(
    PmGridShape shape,
    PmDecompositionDescriptor descriptor,
    int world_rank,
    int world_size)
    : m_shape(shape),
      m_descriptor(std::move(descriptor)),
      m_world_rank(world_rank),
      m_world_size(world_size) {
  if (!m_shape.isValid() || m_world_size <= 0 || m_world_rank < 0 ||
      m_world_rank >= m_world_size) {
    throw std::invalid_argument("PM decomposition view requires valid mesh/rank metadata");
  }
  if (m_descriptor.process_grid_x == 0U || m_descriptor.process_grid_y == 0U ||
      m_descriptor.process_grid_x * m_descriptor.process_grid_y !=
          static_cast<std::size_t>(m_world_size)) {
    throw std::invalid_argument("PM decomposition process grid does not match world size");
  }
}

const PmDecompositionDescriptor& PmDecompositionView::descriptor() const noexcept {
  return m_descriptor;
}

PmOwnershipExtent3D PmDecompositionView::realExtentForRank(int rank) const {
  if (rank < 0 || rank >= m_world_size) {
    throw std::out_of_range("PM real-space extent rank is outside decomposition world");
  }
  const auto block_range = [](std::size_t count, std::size_t parts, std::size_t coord) {
    const std::size_t base = count / parts;
    const std::size_t remainder = count % parts;
    const std::size_t begin = coord * base + std::min(coord, remainder);
    const std::size_t extent = base + (coord < remainder ? 1U : 0U);
    return std::pair<std::size_t, std::size_t>{begin, extent};
  };
  if (m_descriptor.topology == PmDecompositionTopology::kPencil2D) {
    const std::size_t rank_u = static_cast<std::size_t>(rank);
    const std::size_t px = rank_u % m_descriptor.process_grid_x;
    const std::size_t py = rank_u / m_descriptor.process_grid_x;
    const auto [x_begin, x_count] = block_range(
        m_shape.nx, m_descriptor.process_grid_x, px);
    const auto [y_begin, y_count] = block_range(
        m_shape.ny, m_descriptor.process_grid_y, py);
    return PmOwnershipExtent3D{
        .x_begin = x_begin,
        .x_count = x_count,
        .y_begin = y_begin,
        .y_count = y_count,
        .z_begin = 0U,
        .z_count = m_shape.nz,
    };
  }
  const parallel::PmSlabRange owned =
      parallel::pmOwnedXRangeForRank(m_shape.nx, m_world_size, rank);
  return PmOwnershipExtent3D{
      .x_begin = owned.begin_x,
      .x_count = owned.extentX(),
      .y_begin = 0U,
      .y_count = m_shape.ny,
      .z_begin = 0U,
      .z_count = m_shape.nz,
  };
}

int PmDecompositionView::realOwnerRank(PmGlobalCell global_cell) const {
  if (global_cell.x >= m_shape.nx || global_cell.y >= m_shape.ny ||
      global_cell.z >= m_shape.nz) {
    throw std::out_of_range("PM real-space global cell is outside the mesh");
  }
  if (m_descriptor.topology == PmDecompositionTopology::kPencil2D) {
    const auto owner_coord = [](std::size_t index, std::size_t count, std::size_t parts) {
      const std::size_t base = count / parts;
      const std::size_t remainder = count % parts;
      const std::size_t wide = (base + 1U) * remainder;
      if (index < wide) {
        return index / (base + 1U);
      }
      return remainder + (index - wide) / std::max<std::size_t>(base, 1U);
    };
    const std::size_t px = owner_coord(
        global_cell.x, m_shape.nx, m_descriptor.process_grid_x);
    const std::size_t py = owner_coord(
        global_cell.y, m_shape.ny, m_descriptor.process_grid_y);
    return static_cast<int>(py * m_descriptor.process_grid_x + px);
  }
  return parallel::pmOwnerRankForGlobalX(
      m_shape.nx, m_world_size, global_cell.x);
}

int PmDecompositionView::spectralOwnerRank(PmGlobalCell global_cell) const {
  if (global_cell.x >= m_shape.nx || global_cell.y >= m_shape.ny ||
      global_cell.z >= m_shape.nz) {
    throw std::out_of_range("PM spectral global cell is outside the mesh");
  }
  if (m_descriptor.topology == PmDecompositionTopology::kXSlabTransposedSpectral) {
    return parallel::pmOwnerRankForGlobalX(
        m_shape.ny, m_world_size, global_cell.y);
  }
  return realOwnerRank(global_cell);
}

bool PmDecompositionView::ownsRealCell(
    std::size_t global_x, std::size_t global_y, std::size_t global_z) const noexcept {
  if (global_x >= m_shape.nx || global_y >= m_shape.ny || global_z >= m_shape.nz) {
    return false;
  }
  try {
    return realOwnerRank(PmGlobalCell{global_x, global_y, global_z}) == m_world_rank;
  } catch (...) {
    return false;
  }
}

std::size_t PmDecompositionView::globalToLocalRealIndex(
    std::size_t global_x, std::size_t global_y, std::size_t global_z) const {
  if (!ownsRealCell(global_x, global_y, global_z)) {
    throw std::out_of_range("PM real-space global index is not owned by this backend view");
  }
  const PmOwnershipExtent3D extent = realExtentForRank(m_world_rank);
  return core::checkedSizeAdd(
      core::checkedSizeProduct3(
          global_x - extent.x_begin, extent.y_count, extent.z_count,
          "PM real local index"),
      core::checkedSizeAdd(
          core::checkedSizeProduct3(
              global_y - extent.y_begin, extent.z_count, 1U,
              "PM real local index"),
          global_z - extent.z_begin,
          "PM real local index"),
      "PM real local index");
}

bool PmDecompositionView::ownsSpectralCell(
    std::size_t global_x, std::size_t global_y, std::size_t global_z) const noexcept {
  const auto& e = m_descriptor.spectral_extent;
  return e.x_count != 0U && e.y_count != 0U && e.z_count != 0U &&
      global_x >= e.x_begin && global_x - e.x_begin < e.x_count &&
      global_y >= e.y_begin && global_y - e.y_begin < e.y_count &&
      global_z >= e.z_begin && global_z - e.z_begin < e.z_count;
}

std::size_t PmDecompositionView::globalToLocalSpectralIndex(
    std::size_t global_x, std::size_t global_y, std::size_t global_z) const {
  if (!ownsSpectralCell(global_x, global_y, global_z)) {
    throw std::out_of_range(
        "PM spectral global index is not owned by this backend view");
  }
  const auto& e = m_descriptor.spectral_extent;
  const std::size_t lx = global_x - e.x_begin;
  const std::size_t ly = global_y - e.y_begin;
  const std::size_t lz = global_z - e.z_begin;
  return core::checkedSizeAdd(
      core::checkedSizeProduct3(lx, e.y_count, e.z_count,
                                "PM spectral local index"),
      core::checkedSizeAdd(
          core::checkedSizeProduct3(ly, e.z_count, 1U, "PM spectral local index"),
          lz, "PM spectral local index"),
      "PM spectral local index");
}

PmDecompositionDescriptor describePmDecomposition(
    PmGridShape shape,
    const parallel::PmSlabLayout& layout,
    core::PmDecompositionMode mode) {
  if (!shape.isValid() || !layout.isValid() ||
      layout.global_nx != shape.nx || layout.global_ny != shape.ny || layout.global_nz != shape.nz) {
    throw std::invalid_argument("PM decomposition descriptor requires a valid layout matching the mesh shape");
  }
  PmDecompositionDescriptor descriptor;
  descriptor.process_grid_x = static_cast<std::size_t>(std::max(layout.world_size, 1));
  descriptor.process_grid_y = 1U;
  descriptor.real_extent = PmOwnershipExtent3D{
      .x_begin = layout.owned_x.begin_x,
      .x_count = layout.local_nx(),
      .y_begin = 0U,
      .y_count = shape.ny,
      .z_begin = 0U,
      .z_count = shape.nz,
  };
  // Spectral ownership is backend-specific. The current transposed-slab FFTW
  // implementation fills it internally; zero counts here deliberately mean
  // "owned by the selected backend contract", not a fabricated pencil range.
  if (layout.ownsFullDomain()) {
    descriptor.topology = PmDecompositionTopology::kSerial;
    descriptor.spectral_extent = descriptor.real_extent;
  } else if (mode == core::PmDecompositionMode::kPencil) {
    descriptor.topology = PmDecompositionTopology::kXSlabTransposedSpectral;
  } else {
    descriptor.topology = PmDecompositionTopology::kXSlab;
  }
  return descriptor;
}

void requireTreePmSupportOrThrow(
    core::GravitySolver gravity_solver,
    bool allow_diagnostic_naive_dft) {
  if (gravity_solver != core::GravitySolver::kTreePm) {
    return;
  }
  const PmBackendCapability capability = pmBackendCapability();
  if (pmBackendProductionReady()) {
    return;
  }
  if (capability == PmBackendCapability::kDiagnosticNaiveDft && allow_diagnostic_naive_dft) {
    return;
  }
  throw std::runtime_error(
      "TreePM production run requires a production FFT backend; this build exposes only " +
      std::string(pmBackendCapabilityName(capability)) +
      ". Enable FFTW or use the explicit diagnostic naive-DFT override for tiny validation cases.");
}

}  // namespace cosmosim::gravity
