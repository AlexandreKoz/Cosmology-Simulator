#include "cosmosim/gravity/pm_cuda_kernels.hpp"

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>

#include <cuda_runtime.h>

namespace cosmosim::gravity {
namespace {

__device__ std::size_t wrapIndexDevice(long long i, std::size_t n) {
  const long long n_signed = static_cast<long long>(n);
  long long wrapped = i % n_signed;
  if (wrapped < 0) {
    wrapped += n_signed;
  }
  return static_cast<std::size_t>(wrapped);
}

__device__ double wrapPositionDevice(double x, double box_size) {
  const double wrapped = fmod(x, box_size);
  return wrapped < 0.0 ? wrapped + box_size : wrapped;
}

__global__ void assignDensityKernel(
    PmCudaAssignLaunch launch,
    const double* pos_x,
    const double* pos_y,
    const double* pos_z,
    const double* mass,
    double* density) {
  const std::size_t particle_index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (particle_index >= launch.particle_count) {
    return;
  }

  const double inv_dx = static_cast<double>(launch.nx) / launch.box_size_mpc_comoving;
  const double inv_dy = static_cast<double>(launch.ny) / launch.box_size_mpc_comoving;
  const double inv_dz = static_cast<double>(launch.nz) / launch.box_size_mpc_comoving;

  const double x = wrapPositionDevice(pos_x[particle_index], launch.box_size_mpc_comoving) * inv_dx;
  const double y = wrapPositionDevice(pos_y[particle_index], launch.box_size_mpc_comoving) * inv_dy;
  const double z = wrapPositionDevice(pos_z[particle_index], launch.box_size_mpc_comoving) * inv_dz;

  const long long ix0 = static_cast<long long>(floor(x));
  const long long iy0 = static_cast<long long>(floor(y));
  const long long iz0 = static_cast<long long>(floor(z));

  const double tx = x - floor(x);
  const double ty = y - floor(y);
  const double tz = z - floor(z);

  const double wx[2] = {1.0 - tx, tx};
  const double wy[2] = {1.0 - ty, ty};
  const double wz[2] = {1.0 - tz, tz};

  for (std::size_t dx = 0; dx < 2; ++dx) {
    const std::size_t ix = wrapIndexDevice(ix0 + static_cast<long long>(dx), launch.nx);
    for (std::size_t dy = 0; dy < 2; ++dy) {
      const std::size_t iy = wrapIndexDevice(iy0 + static_cast<long long>(dy), launch.ny);
      for (std::size_t dz = 0; dz < 2; ++dz) {
        const std::size_t iz = wrapIndexDevice(iz0 + static_cast<long long>(dz), launch.nz);
        const std::size_t linear_index = (ix * launch.ny + iy) * launch.nz + iz;
        atomicAdd(&density[linear_index], mass[particle_index] * wx[dx] * wy[dy] * wz[dz]);
      }
    }
  }
}

__global__ void interpolateForcesKernel(
    PmCudaInterpLaunch launch,
    const double* pos_x,
    const double* pos_y,
    const double* pos_z,
    const double* force_x,
    const double* force_y,
    const double* force_z,
    double* accel_x,
    double* accel_y,
    double* accel_z) {
  const std::size_t particle_index = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
  if (particle_index >= launch.particle_count) {
    return;
  }

  const double inv_dx = static_cast<double>(launch.nx) / launch.box_size_mpc_comoving;
  const double inv_dy = static_cast<double>(launch.ny) / launch.box_size_mpc_comoving;
  const double inv_dz = static_cast<double>(launch.nz) / launch.box_size_mpc_comoving;

  const double x = wrapPositionDevice(pos_x[particle_index], launch.box_size_mpc_comoving) * inv_dx;
  const double y = wrapPositionDevice(pos_y[particle_index], launch.box_size_mpc_comoving) * inv_dy;
  const double z = wrapPositionDevice(pos_z[particle_index], launch.box_size_mpc_comoving) * inv_dz;

  const long long ix0 = static_cast<long long>(floor(x));
  const long long iy0 = static_cast<long long>(floor(y));
  const long long iz0 = static_cast<long long>(floor(z));

  const double tx = x - floor(x);
  const double ty = y - floor(y);
  const double tz = z - floor(z);

  const double wx[2] = {1.0 - tx, tx};
  const double wy[2] = {1.0 - ty, ty};
  const double wz[2] = {1.0 - tz, tz};

  double gx = 0.0;
  double gy = 0.0;
  double gz = 0.0;

  for (std::size_t dx = 0; dx < 2; ++dx) {
    const std::size_t ix = wrapIndexDevice(ix0 + static_cast<long long>(dx), launch.nx);
    for (std::size_t dy = 0; dy < 2; ++dy) {
      const std::size_t iy = wrapIndexDevice(iy0 + static_cast<long long>(dy), launch.ny);
      for (std::size_t dz = 0; dz < 2; ++dz) {
        const std::size_t iz = wrapIndexDevice(iz0 + static_cast<long long>(dz), launch.nz);
        const std::size_t linear_index = (ix * launch.ny + iy) * launch.nz + iz;
        const double weight = wx[dx] * wy[dy] * wz[dz];
        gx += weight * force_x[linear_index];
        gy += weight * force_y[linear_index];
        gz += weight * force_z[linear_index];
      }
    }
  }

  accel_x[particle_index] = gx;
  accel_y[particle_index] = gy;
  accel_z[particle_index] = gz;
}

void throwOnCudaError(cudaError_t status, const char* context) {
  if (status != cudaSuccess) {
    throw std::runtime_error(std::string(context) + ": " + cudaGetErrorString(status));
  }
}

}  // namespace

void pmCudaAssignDensityCic(
    const PmCudaAssignLaunch& launch,
    const double* pos_x_device,
    const double* pos_y_device,
    const double* pos_z_device,
    const double* mass_device,
    double* density_device,
    void* stream_handle) {
  const cudaStream_t stream = static_cast<cudaStream_t>(stream_handle);
  throwOnCudaError(
      cudaMemsetAsync(density_device, 0, launch.nx * launch.ny * launch.nz * sizeof(double), stream),
      "cudaMemsetAsync density");

  constexpr int threads_per_block = 256;
  const int blocks = static_cast<int>((launch.particle_count + threads_per_block - 1U) / threads_per_block);
  assignDensityKernel<<<blocks, threads_per_block, 0, stream>>>(
      launch,
      pos_x_device,
      pos_y_device,
      pos_z_device,
      mass_device,
      density_device);
  throwOnCudaError(cudaGetLastError(), "assignDensityKernel launch");
}

void pmCudaInterpolateForcesCic(
    const PmCudaInterpLaunch& launch,
    const double* pos_x_device,
    const double* pos_y_device,
    const double* pos_z_device,
    const double* force_x_device,
    const double* force_y_device,
    const double* force_z_device,
    double* accel_x_device,
    double* accel_y_device,
    double* accel_z_device,
    void* stream_handle) {
  const cudaStream_t stream = static_cast<cudaStream_t>(stream_handle);
  constexpr int threads_per_block = 256;
  const int blocks = static_cast<int>((launch.particle_count + threads_per_block - 1U) / threads_per_block);
  interpolateForcesKernel<<<blocks, threads_per_block, 0, stream>>>(
      launch,
      pos_x_device,
      pos_y_device,
      pos_z_device,
      force_x_device,
      force_y_device,
      force_z_device,
      accel_x_device,
      accel_y_device,
      accel_z_device);
  throwOnCudaError(cudaGetLastError(), "interpolateForcesKernel launch");
}

}  // namespace cosmosim::gravity
