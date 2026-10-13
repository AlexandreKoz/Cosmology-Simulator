// 20261012 regression repairs: NOT RUN — deferred to maintainer.
#ifdef NDEBUG
#error "CosmoSim test targets must keep standard assertions enabled"
#endif

// First preserve the existing target-level assertion-policy check above, then
// exercise the scientific helper with NDEBUG actually defined in this TU.
#define NDEBUG 1
#include <cassert>
#include "../support/test_require.hpp"

int main() {
  int standard_evaluations = 0;
  assert(++standard_evaluations == 1);
  int require_evaluations = 0;
  CHUI_REQUIRE(++require_evaluations == 1);
  bool rejected = false;
  try {
    CHUI_REQUIRE(++require_evaluations == 0);
  } catch (const std::runtime_error& error) {
    const std::string message = error.what();
    rejected = message.find("requirement failed") != std::string::npos &&
        message.find("++require_evaluations == 0") != std::string::npos;
  }
  // Independent return condition detects even a disabled CHUI_REQUIRE macro.
  return standard_evaluations == 0 && require_evaluations == 2 && rejected ? 0 : 1;
}
