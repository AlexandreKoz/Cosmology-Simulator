#pragma once

#include <cstdint>
#include <exception>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/parallel/distributed_memory.hpp"

namespace cosmosim::io::failure_protocol_internal {

#if COSMOSIM_ENABLE_HDF5 && COSMOSIM_ENABLE_MPI

[[nodiscard]] std::string collectiveFailureMessage(
    const parallel::MpiContext& mpi_context,
    const std::exception_ptr& local_failure,
    std::string_view phase);

void incrementCollectivePhaseCounter() noexcept;

class CollectivePhaseCounterScope {
 public:
  explicit CollectivePhaseCounterScope(std::uint64_t& counter) noexcept;
  CollectivePhaseCounterScope(const CollectivePhaseCounterScope&) = delete;
  CollectivePhaseCounterScope& operator=(const CollectivePhaseCounterScope&) =
      delete;
  ~CollectivePhaseCounterScope();

 private:
  std::uint64_t* m_previous = nullptr;
};

template <typename T, typename Function>
[[nodiscard]] T runCollectivePhase(
    const parallel::MpiContext& mpi_context,
    std::string_view phase,
    Function&& function) {
  incrementCollectivePhaseCounter();
  std::optional<T> value;
  std::exception_ptr local_failure;
  try {
    value.emplace(std::invoke(std::forward<Function>(function)));
  } catch (...) {
    local_failure = std::current_exception();
  }
  const std::string failure =
      collectiveFailureMessage(mpi_context, local_failure, phase);
  if (!failure.empty()) {
    throw std::runtime_error(failure);
  }
  return std::move(*value);
}

template <typename Function>
void runCollectivePhaseVoid(
    const parallel::MpiContext& mpi_context,
    std::string_view phase,
    Function&& function) {
  incrementCollectivePhaseCounter();
  std::exception_ptr local_failure;
  try {
    std::invoke(std::forward<Function>(function));
  } catch (...) {
    local_failure = std::current_exception();
  }
  const std::string failure =
      collectiveFailureMessage(mpi_context, local_failure, phase);
  if (!failure.empty()) {
    throw std::runtime_error(failure);
  }
}

void injectIcTestFault(
    const parallel::MpiContext& mpi_context,
    std::string_view phase);

void mutateIcTestRoute(
    const parallel::MpiContext& mpi_context,
    std::vector<std::vector<std::uint8_t>>& per_rank);

[[nodiscard]] std::string broadcastRootString(
    const parallel::MpiContext& mpi_context,
    std::string root_value);

[[nodiscard]] std::vector<std::uint8_t> alltoallBytes(
    const parallel::MpiContext& mpi_context,
    const std::vector<std::vector<std::uint8_t>>& per_rank,
    std::uint64_t& sent_bytes,
    std::uint64_t& received_bytes,
    std::uint64_t* exchange_peak_bytes = nullptr);

#endif

}  // namespace cosmosim::io::failure_protocol_internal
