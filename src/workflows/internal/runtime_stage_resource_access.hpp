#pragma once

#include <initializer_list>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>

#include "cosmosim/core/time_integration.hpp"
#include "cosmosim/workflows/runtime_resources.hpp"

namespace cosmosim::workflows::internal {

// Opaque source-private bridge between the core dispatcher and built-in runtime
// owners. Public headers never expose this bundle or a StepContext recovery API.
struct RuntimeStageResourceBundle {
  explicit RuntimeStageResourceBundle(core::StepContext& stage_context) noexcept
      : context(&stage_context) {}

  core::StepContext* context = nullptr;
};

class RuntimeStageAccess final {
 public:
  [[nodiscard]] static core::StepContext& driftContext(
      DriftParticleStageView& view,
      std::initializer_list<RuntimeResourceAccess> required) {
    return requireContext(view, required, "drift runtime");
  }

  [[nodiscard]] static core::StepContext& gravityContext(
      GravityStageView& view,
      std::initializer_list<RuntimeResourceAccess> required) {
    return requireContext(view, required, "gravity runtime");
  }

  [[nodiscard]] static core::StepContext& hydroAmrContext(
      HydroAmrStageView& view,
      std::initializer_list<RuntimeResourceAccess> required) {
    return requireContext(view, required, "hydro/AMR runtime");
  }

  [[nodiscard]] static core::StepContext& sourceContext(
      SourceMutationStageView& view,
      std::initializer_list<RuntimeResourceAccess> required) {
    return requireContext(view, required, "source runtime");
  }

  [[nodiscard]] static core::StepContext& analysisContext(
      AnalysisStageView& view,
      std::initializer_list<RuntimeResourceAccess> required) {
    return requireContext(view, required, "analysis runtime");
  }

  [[nodiscard]] static core::StepContext& outputRestartContext(
      OutputRestartStageView& view,
      std::initializer_list<RuntimeResourceAccess> required) {
    return requireContext(view, required, "output/restart runtime");
  }

 private:
  template <class View>
  [[nodiscard]] static core::StepContext& requireContext(
      View& view,
      std::initializer_list<RuntimeResourceAccess> required,
      std::string_view caller) {
    view.requireFresh();
    view.requireDeclaredResources(
        std::span<const RuntimeResourceAccess>(required.begin(), required.size()),
        caller);
    if (view.m_resources == nullptr || view.m_resources->context == nullptr) {
      throw std::logic_error(
          std::string(caller) + " received a view without runtime resources");
    }
    return *view.m_resources->context;
  }
};

}  // namespace cosmosim::workflows::internal
