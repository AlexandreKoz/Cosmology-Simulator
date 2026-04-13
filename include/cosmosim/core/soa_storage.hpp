#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <new>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace cosmosim::core {

// Fixed-alignment allocator used by hot SoA arrays. Keeping alignment explicit
// avoids hidden assumptions in solver kernels and simplifies future host/device
// mirror contracts.
template <typename T, std::size_t k_alignment>
class AlignedAllocator {
 public:
  using value_type = T;

  AlignedAllocator() noexcept = default;

  template <typename U>
  constexpr AlignedAllocator(const AlignedAllocator<U, k_alignment>&) noexcept {}

  [[nodiscard]] T* allocate(std::size_t count) {
    if (count == 0) {
      return nullptr;
    }
    if (count > (static_cast<std::size_t>(-1) / sizeof(T))) {
      throw std::bad_alloc{};
    }
    void* memory = ::operator new(count * sizeof(T), std::align_val_t(k_alignment));
    return static_cast<T*>(memory);
  }

  void deallocate(T* ptr, std::size_t) noexcept { ::operator delete(ptr, std::align_val_t(k_alignment)); }

  template <typename U>
  struct rebind {
    using other = AlignedAllocator<U, k_alignment>;
  };

  using is_always_equal = std::true_type;
};

template <typename T, typename U, std::size_t k_alignment>
[[nodiscard]] constexpr bool operator==(
    const AlignedAllocator<T, k_alignment>&,
    const AlignedAllocator<U, k_alignment>&) noexcept {
  return true;
}

template <typename T, typename U, std::size_t k_alignment>
[[nodiscard]] constexpr bool operator!=(
    const AlignedAllocator<T, k_alignment>&,
    const AlignedAllocator<U, k_alignment>&) noexcept {
  return false;
}

template <typename T>
using AlignedVector = std::vector<T, AlignedAllocator<T, 64>>;

// Single field lane for SoA storage:
// - contiguous memory
// - explicit size/capacity controls
// - optional debug-only bounds checks
// This class intentionally does not add indirection or virtual dispatch.
template <typename T>
class SoaFieldArray {
 public:
  using value_type = T;

  void resize(std::size_t count) { m_values.resize(count); }
  void reserve(std::size_t count) { m_values.reserve(count); }
  [[nodiscard]] std::size_t size() const noexcept { return m_values.size(); }
  [[nodiscard]] std::size_t capacity() const noexcept { return m_values.capacity(); }
  [[nodiscard]] bool empty() const noexcept { return m_values.empty(); }

  [[nodiscard]] T* data() noexcept { return m_values.data(); }
  [[nodiscard]] const T* data() const noexcept { return m_values.data(); }

  [[nodiscard]] std::span<T> span() noexcept { return std::span<T>(m_values.data(), m_values.size()); }
  [[nodiscard]] std::span<const T> span() const noexcept {
    return std::span<const T>(m_values.data(), m_values.size());
  }

  [[nodiscard]] T& operator[](std::size_t index) noexcept {
#ifndef NDEBUG
    if (index >= m_values.size()) {
      throw std::out_of_range("SoaFieldArray.operator[]: index out of range");
    }
#endif
    return m_values[index];
  }

  [[nodiscard]] const T& operator[](std::size_t index) const noexcept {
#ifndef NDEBUG
    if (index >= m_values.size()) {
      throw std::out_of_range("SoaFieldArray.operator[]: index out of range");
    }
#endif
    return m_values[index];
  }

  void swapErase(std::size_t index) {
#ifndef NDEBUG
    if (index >= m_values.size()) {
      throw std::out_of_range("SoaFieldArray.swapErase: index out of range");
    }
#endif
    if (index + 1 < m_values.size()) {
      m_values[index] = m_values.back();
    }
    m_values.pop_back();
  }

  void clear() noexcept { m_values.clear(); }

  [[nodiscard]] std::size_t stableCompact(std::span<const std::uint8_t> keep_flags) {
    if (keep_flags.size() != m_values.size()) {
      throw std::invalid_argument("SoaFieldArray.stableCompact: keep_flags size mismatch");
    }

    // Stable in-place compaction:
    // preserve relative order of retained entries while rewriting into the
    // front of the buffer and truncating the tail once.
    std::size_t write_index = 0;
    for (std::size_t read_index = 0; read_index < m_values.size(); ++read_index) {
      if (keep_flags[read_index] != 0U) {
        if (write_index != read_index) {
          m_values[write_index] = m_values[read_index];
        }
        ++write_index;
      }
    }
    m_values.resize(write_index);
    return write_index;
  }

 private:
  AlignedVector<T> m_values;
};

template <typename T>
void gatherSpan(std::span<const T> source, std::span<const std::uint32_t> indices, std::span<T> destination) {
  // Gather into a compact destination span for active-set kernels.
  if (destination.size() != indices.size()) {
    throw std::invalid_argument("gatherSpan: destination size must match indices size");
  }
  for (std::size_t i = 0; i < indices.size(); ++i) {
    const std::uint32_t source_index = indices[i];
#ifndef NDEBUG
    if (source_index >= source.size()) {
      throw std::out_of_range("gatherSpan: source index out of range");
    }
#endif
    destination[i] = source[source_index];
  }
}

template <typename T>
void scatterSpan(std::span<const T> source, std::span<const std::uint32_t> indices, std::span<T> destination) {
  // Scatter compact computed values back into sparse destination indices.
  if (source.size() != indices.size()) {
    throw std::invalid_argument("scatterSpan: source size must match indices size");
  }
  for (std::size_t i = 0; i < indices.size(); ++i) {
    const std::uint32_t destination_index = indices[i];
#ifndef NDEBUG
    if (destination_index >= destination.size()) {
      throw std::out_of_range("scatterSpan: destination index out of range");
    }
#endif
    destination[destination_index] = source[i];
  }
}

enum class ParticleSoaField : std::uint8_t {
  kPosX,
  kPosY,
  kPosZ,
  kVelX,
  kVelY,
  kVelZ,
  kMass,
  kId,
  kRho,
  kUInt,
};

// Canonical particle-oriented field pack used by modules that need SoA without
// hand-rolling N separate vectors. Field names intentionally mirror common IC
// / runtime naming (`pos_*`, `vel_*`, `mass`, `id`, `rho`, `u_int`).
class ParticleSoaStorage {
 public:
  // Resize every particle lane to the same logical row count.
  void resize(std::size_t count);
  // Reserve capacity in every particle lane at once.
  void reserve(std::size_t count);
  // Logical row count shared by every lane.
  [[nodiscard]] std::size_t size() const noexcept;
  // Usable shared capacity (minimum lane capacity).
  [[nodiscard]] std::size_t capacity() const noexcept;

  template <ParticleSoaField field>
  [[nodiscard]] auto span() noexcept;

  template <ParticleSoaField field>
  [[nodiscard]] auto span() const noexcept;

  // Validate that all lanes still share identical logical size.
  [[nodiscard]] bool isConsistent() const noexcept;
  // Fast unordered erase mirrored over all lanes.
  void swapErase(std::size_t index);
  // Stable keep-mask compaction mirrored over all lanes.
  [[nodiscard]] std::size_t stableCompact(std::span<const std::uint8_t> keep_flags);

 private:
  SoaFieldArray<double> m_pos_x;
  SoaFieldArray<double> m_pos_y;
  SoaFieldArray<double> m_pos_z;
  SoaFieldArray<double> m_vel_x;
  SoaFieldArray<double> m_vel_y;
  SoaFieldArray<double> m_vel_z;
  SoaFieldArray<double> m_mass;
  SoaFieldArray<std::uint64_t> m_id;
  SoaFieldArray<double> m_rho;
  SoaFieldArray<double> m_u_int;
};

template <ParticleSoaField field>
[[nodiscard]] auto ParticleSoaStorage::span() noexcept {
  // Compile-time field selection keeps access strongly typed and branch-free
  // after optimization.
  if constexpr (field == ParticleSoaField::kPosX) {
    return m_pos_x.span();
  } else if constexpr (field == ParticleSoaField::kPosY) {
    return m_pos_y.span();
  } else if constexpr (field == ParticleSoaField::kPosZ) {
    return m_pos_z.span();
  } else if constexpr (field == ParticleSoaField::kVelX) {
    return m_vel_x.span();
  } else if constexpr (field == ParticleSoaField::kVelY) {
    return m_vel_y.span();
  } else if constexpr (field == ParticleSoaField::kVelZ) {
    return m_vel_z.span();
  } else if constexpr (field == ParticleSoaField::kMass) {
    return m_mass.span();
  } else if constexpr (field == ParticleSoaField::kId) {
    return m_id.span();
  } else if constexpr (field == ParticleSoaField::kRho) {
    return m_rho.span();
  } else {
    static_assert(field == ParticleSoaField::kUInt);
    return m_u_int.span();
  }
}

template <ParticleSoaField field>
[[nodiscard]] auto ParticleSoaStorage::span() const noexcept {
  // Const overload mirrors mutable accessor for read-only kernel paths.
  if constexpr (field == ParticleSoaField::kPosX) {
    return m_pos_x.span();
  } else if constexpr (field == ParticleSoaField::kPosY) {
    return m_pos_y.span();
  } else if constexpr (field == ParticleSoaField::kPosZ) {
    return m_pos_z.span();
  } else if constexpr (field == ParticleSoaField::kVelX) {
    return m_vel_x.span();
  } else if constexpr (field == ParticleSoaField::kVelY) {
    return m_vel_y.span();
  } else if constexpr (field == ParticleSoaField::kVelZ) {
    return m_vel_z.span();
  } else if constexpr (field == ParticleSoaField::kMass) {
    return m_mass.span();
  } else if constexpr (field == ParticleSoaField::kId) {
    return m_id.span();
  } else if constexpr (field == ParticleSoaField::kRho) {
    return m_rho.span();
  } else {
    static_assert(field == ParticleSoaField::kUInt);
    return m_u_int.span();
  }
}

}  // namespace cosmosim::core
