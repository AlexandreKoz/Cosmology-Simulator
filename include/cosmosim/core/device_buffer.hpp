#pragma once

#include <cstddef>
#include <cstring>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "cosmosim/core/build_config.hpp"

#if COSMOSIM_ENABLE_CUDA
#include <cuda_runtime.h>
#endif

namespace cosmosim::core {

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

  void resize(std::size_t count) {
    if (count == m_count) {
      return;
    }
    release();
    m_count = count;
#if COSMOSIM_ENABLE_CUDA
    if (m_count > 0 && cudaMalloc(&m_device_ptr, m_count * sizeof(double)) != cudaSuccess) {
      m_device_ptr = nullptr;
      m_count = 0;
      throw std::runtime_error("cudaMalloc failed for DeviceBufferDouble");
    }
#else
    m_host_shadow.resize(m_count, 0.0);
#endif
  }

  [[nodiscard]] std::size_t size() const { return m_count; }

#if COSMOSIM_ENABLE_CUDA
  [[nodiscard]] double* data() { return m_device_ptr; }
  [[nodiscard]] const double* data() const { return m_device_ptr; }
#else
  [[nodiscard]] double* data() { return m_host_shadow.data(); }
  [[nodiscard]] const double* data() const { return m_host_shadow.data(); }
#endif

  double copyFromHost(std::span<const double> host, void* stream_handle = nullptr) {
    if (host.size() != m_count) {
      throw std::invalid_argument("DeviceBufferDouble::copyFromHost span size mismatch");
    }
#if COSMOSIM_ENABLE_CUDA
    auto stream = reinterpret_cast<cudaStream_t>(stream_handle);
    cudaError_t status = cudaMemcpyAsync(
        m_device_ptr,
        host.data(),
        m_count * sizeof(double),
        cudaMemcpyHostToDevice,
        stream);
    if (status != cudaSuccess) {
      throw std::runtime_error("cudaMemcpyAsync host->device failed");
    }
    return static_cast<double>(m_count * sizeof(double));
#else
    std::memcpy(m_host_shadow.data(), host.data(), m_count * sizeof(double));
    return static_cast<double>(m_count * sizeof(double));
#endif
  }

  double copyToHost(std::span<double> host, void* stream_handle = nullptr) const {
    if (host.size() != m_count) {
      throw std::invalid_argument("DeviceBufferDouble::copyToHost span size mismatch");
    }
#if COSMOSIM_ENABLE_CUDA
    auto stream = reinterpret_cast<cudaStream_t>(stream_handle);
    cudaError_t status = cudaMemcpyAsync(
        host.data(),
        m_device_ptr,
        m_count * sizeof(double),
        cudaMemcpyDeviceToHost,
        stream);
    if (status != cudaSuccess) {
      throw std::runtime_error("cudaMemcpyAsync device->host failed");
    }
    return static_cast<double>(m_count * sizeof(double));
#else
    std::memcpy(host.data(), m_host_shadow.data(), m_count * sizeof(double));
    return static_cast<double>(m_count * sizeof(double));
#endif
  }

 private:
  void release() {
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

  void moveFrom(DeviceBufferDouble&& other) {
    m_count = other.m_count;
#if COSMOSIM_ENABLE_CUDA
    m_device_ptr = other.m_device_ptr;
    other.m_device_ptr = nullptr;
#else
    m_host_shadow = std::move(other.m_host_shadow);
#endif
    other.m_count = 0;
  }

  std::size_t m_count = 0;
#if COSMOSIM_ENABLE_CUDA
  double* m_device_ptr = nullptr;
#else
  std::vector<double> m_host_shadow;
#endif
};

}  // namespace cosmosim::core
