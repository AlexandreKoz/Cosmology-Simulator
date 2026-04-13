#include <cassert>
#include <string>

#include "cosmosim/cosmosim.hpp"

int main() {
  const auto modules = cosmosim::core::moduleNames();
  assert(!modules.empty());
  assert(cosmosim::core::projectName() == "cosmosim");
  assert(cosmosim::gravity::GravityModule::name() == "gravity");
  assert(cosmosim::hydro::HydroModule::name() == "hydro");
  assert(cosmosim::amr::AmrModule::name() == "amr");
  assert(cosmosim::physics::PhysicsModule::name() == "physics");
  assert(cosmosim::io::IoModule::name() == "io");
  assert(cosmosim::analysis::AnalysisModule::name() == "analysis");
  assert(cosmosim::parallel::ParallelModule::name() == "parallel");
  assert(cosmosim::utils::UtilsModule::name() == "utils");
  assert(!cosmosim::architectureSummary().empty());
  assert(!cosmosim::core::buildProvenance().empty());
  assert(cosmosim::core::buildProvenance().find("preset=") != std::string::npos);
  return 0;
}
