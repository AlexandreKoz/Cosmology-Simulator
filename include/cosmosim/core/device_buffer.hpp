#pragma once

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/core/checked_arithmetic.hpp"

#if COSMOSIM_ENABLE_CUDA
#include <cuda_runtime.h>
#endif

namespace cosmosim::core {

class CudaStreamView {
 public:
  constexpr CudaStreamView() noexcept = default;
#if COSMOSIM_ENABLE_CUDA
  explicit constexpr CudaStreamView(cudaStream_t stream) noexcept : m_stream(stream) {}
  [[nodiscard]] constexpr cudaStream_t nativeHandle() const noexcept { return m_stream; }
  void synchronize() const {
    if (cudaStreamSynchronize(m_stream) != cudaSuccess) {
      throw std::runtime_error("cudaStreamSynchronize failed");
    }
  }
#else
  explicit constexpr CudaStreamView(std::nullptr_t) noexcept {}
  void synchronize() const noexcept {}
#endif

 private:
#if COSMOSIM_ENABLE_CUDA
  cudaStream_t m_stream = nullptr;
#endif
};

class DeviceBufferDouble {
 public:
  DeviceBufferDouble() = default;
  explicit DeviceBufferDouble(std::size_t count) { resize(count); }

  ~DeviceBufferDouble() { release(); }

  DeviceBufferDouble(const DeviceBufferDouble&) = delete;
  DeviceBufferDouble& operator=(const DeviceBufferDouble&) = delete;

  DeviceBufferDouble(DeviceBufferDouble&& other) noexcept { moveFrom(std::move(other)); }

  DeviceBufferDouble& operator=(DeviceBufferDouble&& other) noexcept {
    if (this != &other) {
      release();
      moveFrom(std::move(other));
    }
    return *this;
  }

  // Strong exception guarantee: allocation/overflow failure leaves the original
  // buffer and its externally visible size unchanged.
  void resize(std::size_t count) {
    if (count == m_count) {
      return;
    }
    const std::size_t bytes = checkedSizeMultiply(
        count, sizeof(double), "DeviceBufferDouble::resize");
#if COSMOSIM_ENABLE_CUDA
    double* replacement = nullptr;
    if (count > 0 && cudaMalloc(&replacement, bytes) != cudaSuccess) {
      throw std::runtime_error("cudaMalloc failed for DeviceBufferDouble");
    }
    double* old = m_device_ptr;
    m_device_ptr = replacement;
    m_count = count;
    m_high_water_count = std::max(m_high_water_count, count);
    if (old != nullptr) {
      (void)cudaFree(old);
    }
#else
    (void)bytes;
    std::vector<double> replacement(count, 0.0);
    m_host_shadow.swap(replacement);
    m_count = count;
    m_high_water_count = std::max(m_high_water_count, count);
#endif
  }

  [[nodiscard]] std::size_t size() const noexcept { return m_count; }
  [[nodiscard]] std::size_t highWaterSize() const noexcept { return m_high_water_count; }
  [[nodiscard]] std::size_t sizeBytes() const {
    return checkedSizeMultiply(m_count, sizeof(double), "DeviceBufferDouble::sizeBytes");
  }
  [[nodiscard]] std::size_t highWaterBytes() const {
    return checkedSizeMultiply(
        m_high_water_count, sizeof(double), "DeviceBufferDouble::highWaterBytes");
  }

#if COSMOSIM_ENABLE_CUDA
  [[nodiscard]] double* data() noexcept { return m_device_ptr; }
  [[nodiscard]] const double* data() const noexcept { return m_device_ptr; }
#else
  [[nodiscard]] double* data() noexcept { return m_host_shadow.data(); }
  [[nodiscard]] const double* data() const noexcept { return m_host_shadow.data(); }
#endif

  // CUDA mode is asynchronous. The host source must remain valid until the
  // supplied stream completes. CPU fallback completes before returning.
  std::size_t copyFromHostAsync(
      std::span<const double> host,
      CudaStreamView stream = {}) {
    if (host.size() != m_count) {
      throw std::invalid_argument("DeviceBufferDouble::copyFromHostAsync span size mismatch");
    }
    const std::size_t bytes = checkedSizeMultiply(
        m_count, sizeof(double), "DeviceBufferDouble::copyFromHostAsync");
#if COSMOSIM_ENABLE_CUDA
    const cudaError_t status = cudaMemcpyAsync(
        m_device_ptr, host.data(), bytes, cudaMemcpyHostToDevice, stream.nativeHandle());
    if (status != cudaSuccess) {
      throw std::runtime_error("cudaMemcpyAsync host->device failed");
    }
#else
    (void)stream;
    std::memcpy(m_host_shadow.data(), host.data(), bytes);
#endif
    return bytes;
  }

  // CUDA mode is asynchronous. The host destination must not be consumed until
  // the supplied stream completes. CPU fallback completes before returning.
  std::size_t copyToHostAsync(
      std::span<double> host,
      CudaStreamView stream = {}) const {
    if (host.size() != m_count) {
      throw std::invalid_argument("DeviceBufferDouble::copyToHostAsync span size mismatch");
    }
    const std::size_t bytes = checkedSizeMultiply(
        m_count, sizeof(double), "DeviceBufferDouble::copyToHostAsync");
#if COSMOSIM_ENABLE_CUDA
    const cudaError_t status = cudaMemcpyAsync(
        host.data(), m_device_ptr, bytes, cudaMemcpyDeviceToHost, stream.nativeHandle());
    if (status != cudaSuccess) {
      throw std::runtime_error("cudaMemcpyAsync device->host failed");
    }
#else
    (void)stream;
    std::memcpy(host.data(), m_host_shadow.data(), bytes);
#endif
    return bytes;
  }

  std::size_t copyFromHost(std::span<const double> host) {
#if COSMOSIM_ENABLE_CUDA
    const std::size_t bytes = checkedSizeMultiply(
        m_count, sizeof(double), "DeviceBufferDouble::copyFromHost");
    if (host.size() != m_count) {
      throw std::invalid_argument("DeviceBufferDouble::copyFromHost span size mismatch");
    }
    if (cudaMemcpy(m_device_ptr, host.data(), bytes, cudaMemcpyHostToDevice) != cudaSuccess) {
      throw std::runtime_error("cudaMemcpy host->device failed");
    }
    return bytes;
#else
    return copyFromHostAsync(host);
#endif
  }

  std::size_t copyToHost(std::span<double> host) const {
#if COSMOSIM_ENABLE_CUDA
    const std::size_t bytes = checkedSizeMultiply(
        m_count, sizeof(double), "DeviceBufferDouble::copyToHost");
    if (host.size() != m_count) {
      throw std::invalid_argument("DeviceBufferDouble::copyToHost span size mismatch");
    }
    if (cudaMemcpy(host.data(), m_device_ptr, bytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
      throw std::runtime_error("cudaMemcpy device->host failed");
    }
    return bytes;
#else
    return copyToHostAsync(host);
#endif
  }

 private:
  void release() noexcept {
#if COSMOSIM_ENABLE_CUDA
    if (m_device_ptr != nullptr) {
      (void)cudaFree(m_device_ptr);
      m_device_ptr = nullptr;
    }
#else
    m_host_shadow.clear();
#endif
    m_count = 0;
  }

  void moveFrom(DeviceBufferDouble&& other) noexcept {
    m_count = other.m_count;
    m_high_water_count = std::max(m_high_water_count, other.m_high_water_count);
#if COSMOSIM_ENABLE_CUDA
    m_device_ptr = other.m_device_ptr;
    other.m_device_ptr = nullptr;
#else
    m_host_shadow = std::move(other.m_host_shadow);
#endif
    other.m_count = 0;
    other.m_high_water_count = 0;
  }

  std::size_t m_count = 0;
  std::size_t m_high_water_count = 0;
#if COSMOSIM_ENABLE_CUDA
  double* m_device_ptr = nullptr;
#else
  std::vector<double> m_host_shadow;
#endif
};

}  // namespace cosmosim::core
