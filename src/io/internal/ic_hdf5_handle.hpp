#pragma once

#include "cosmosim/core/build_config.hpp"

#if COSMOSIM_ENABLE_HDF5
#include <hdf5.h>

namespace cosmosim::io::internal {

class Hdf5Handle {
 public:
  explicit Hdf5Handle(hid_t handle = -1) : m_handle(handle) {}
  Hdf5Handle(const Hdf5Handle&) = delete;
  Hdf5Handle& operator=(const Hdf5Handle&) = delete;
  Hdf5Handle(Hdf5Handle&& other) noexcept : m_handle(other.m_handle) {
    other.m_handle = -1;
  }
  Hdf5Handle& operator=(Hdf5Handle&& other) noexcept {
    if (this != &other) {
      close();
      m_handle = other.m_handle;
      other.m_handle = -1;
    }
    return *this;
  }
  ~Hdf5Handle() { close(); }

  [[nodiscard]] hid_t get() const noexcept { return m_handle; }
  [[nodiscard]] bool valid() const noexcept { return m_handle >= 0; }

 private:
  void close() noexcept {
    if (m_handle < 0) {
      return;
    }
    switch (H5Iget_type(m_handle)) {
      case H5I_FILE:
        H5Fclose(m_handle);
        break;
      case H5I_GROUP:
        H5Gclose(m_handle);
        break;
      case H5I_DATASET:
        H5Dclose(m_handle);
        break;
      case H5I_DATASPACE:
        H5Sclose(m_handle);
        break;
      case H5I_ATTR:
        H5Aclose(m_handle);
        break;
      case H5I_DATATYPE:
        H5Tclose(m_handle);
        break;
      default:
        break;
    }
    m_handle = -1;
  }

  hid_t m_handle = -1;
};

}  // namespace cosmosim::io::internal
#endif
