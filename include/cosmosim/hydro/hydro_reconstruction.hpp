#pragma once

#include <cstdint>
#include <string_view>

namespace cosmosim::hydro {

struct HydroFace;
struct HydroPrimitiveState;
class HydroConservedStateSoa;
class HydroPrimitiveCacheSoa;

enum class HydroSlopeLimiter {
  kMinmod,
  kMonotonizedCentral,
  kVanLeer,
};

[[nodiscard]] std::string_view hydroSlopeLimiterToString(HydroSlopeLimiter limiter);
[[nodiscard]] HydroSlopeLimiter hydroSlopeLimiterFromString(std::string_view name);
[[nodiscard]] double applyHydroSlopeLimiter(HydroSlopeLimiter limiter, double delta_minus, double delta_plus);

struct HydroReconstructionPolicy {
  HydroSlopeLimiter limiter = HydroSlopeLimiter::kMonotonizedCentral;
  double dt_over_dx_code = 0.0;
  double rho_floor = 1.0e-12;
  double pressure_floor = 1.0e-12;
  bool enable_muscl_hancock_predictor = true;
  double adiabatic_index = 5.0 / 3.0;
};

class HydroReconstruction {
 public:
  virtual ~HydroReconstruction() = default;

  [[nodiscard]] virtual bool reconstructFaceFromCache(
      const HydroPrimitiveCacheSoa& primitive_cache,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state) const;

  virtual void reconstructFace(
      const HydroConservedStateSoa& conserved,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state,
      double adiabatic_index) const = 0;

  [[nodiscard]] virtual std::uint64_t limiterClipCount() const { return 0; }
  [[nodiscard]] virtual std::uint64_t positivityFallbackCount() const { return 0; }
};

class PiecewiseConstantReconstruction final : public HydroReconstruction {
 public:
  [[nodiscard]] bool reconstructFaceFromCache(
      const HydroPrimitiveCacheSoa& primitive_cache,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state) const override;

  void reconstructFace(
      const HydroConservedStateSoa& conserved,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state,
      double adiabatic_index) const override;
};

// MUSCL-Hancock face reconstruction with runtime-selectable slope limiter.
// The predictor currently assumes a 1D contiguous cell ordering; otherwise it falls back to piecewise-constant.
class MusclHancockReconstruction final : public HydroReconstruction {
 public:
  explicit MusclHancockReconstruction(HydroReconstructionPolicy policy = {});

  [[nodiscard]] bool reconstructFaceFromCache(
      const HydroPrimitiveCacheSoa& primitive_cache,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state) const override;

  void reconstructFace(
      const HydroConservedStateSoa& conserved,
      const HydroFace& face,
      HydroPrimitiveState& left_state,
      HydroPrimitiveState& right_state,
      double adiabatic_index) const override;

  [[nodiscard]] HydroReconstructionPolicy policy() const;
  [[nodiscard]] std::uint64_t limiterClipCount() const override;
  [[nodiscard]] std::uint64_t positivityFallbackCount() const override;

 private:
  HydroReconstructionPolicy m_policy;
  mutable std::uint64_t m_limiter_clip_count = 0;
  mutable std::uint64_t m_positivity_fallback_count = 0;
};

}  // namespace cosmosim::hydro
