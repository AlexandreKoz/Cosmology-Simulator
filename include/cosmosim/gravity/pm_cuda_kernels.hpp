#pragma once

#include <cstddef>

namespace cosmosim::gravity {

struct PmCudaAssignLaunch {
  std::size_t particle_count = 0;
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
  double box_size_mpc_comoving = 0.0;
};

struct PmCudaInterpLaunch {
  std::size_t particle_count = 0;
  std::size_t nx = 0;
  std::size_t ny = 0;
  std::size_t nz = 0;
  double box_size_mpc_comoving = 0.0;
};

void pmCudaAssignDensityCic(
    const PmCudaAssignLaunch& launch,
    const double* pos_x_device,
    const double* pos_y_device,
    const double* pos_z_device,
    const double* mass_device,
    double* density_device,
    void* stream_handle);

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
    void* stream_handle);

}  // namespace cosmosim::gravity
