#pragma once

#include <cstdint>

namespace cosmosim::hydro {

struct HydroConservedState;
struct HydroFace;
struct HydroPrimitiveState;

class HydroRiemannSolver {
 public:
  virtual ~HydroRiemannSolver() = default;

  [[nodiscard]] virtual HydroConservedState computeFlux(
      const HydroPrimitiveState& left_state,
      const HydroPrimitiveState& right_state,
      const HydroFace& face,
      double adiabatic_index) const = 0;

  [[nodiscard]] virtual std::uint64_t fallbackCount() const { return 0; }
};

class HlleRiemannSolver final : public HydroRiemannSolver {
 public:
  [[nodiscard]] HydroConservedState computeFlux(
      const HydroPrimitiveState& left_state,
      const HydroPrimitiveState& right_state,
      const HydroFace& face,
      double adiabatic_index) const override;
};

// HLLC is the default approximate Riemann solver with HLLE fallback for positivity/speed degeneracy.
class HllcRiemannSolver final : public HydroRiemannSolver {
 public:
  [[nodiscard]] HydroConservedState computeFlux(
      const HydroPrimitiveState& left_state,
      const HydroPrimitiveState& right_state,
      const HydroFace& face,
      double adiabatic_index) const override;

  [[nodiscard]] std::uint64_t fallbackCount() const override;

 private:
  mutable std::uint64_t m_fallback_count = 0;
};

}  // namespace cosmosim::hydro
