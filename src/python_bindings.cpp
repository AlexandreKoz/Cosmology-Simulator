#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

#include "cosmosim/analysis/diagnostics.hpp"
#include "cosmosim/core/config.hpp"
#include "cosmosim/core/provenance.hpp"
#include "cosmosim/core/simulation_state.hpp"
#include "cosmosim/io/snapshot_hdf5.hpp"

namespace py = pybind11;

namespace cosmosim::python {
namespace {

py::array makeReadonlyDoubleView(const std::vector<double>& values, py::handle owner) {
  py::array array(
      py::dtype::of<double>(),
      {static_cast<py::ssize_t>(values.size())},
      {static_cast<py::ssize_t>(sizeof(double))},
      values.data(),
      owner);
  array.attr("setflags")(false);
  return array;
}

py::array makeReadonlyDoubleView(const core::AlignedVector<double>& values, py::handle owner) {
  py::array array(
      py::dtype::of<double>(),
      {static_cast<py::ssize_t>(values.size())},
      {static_cast<py::ssize_t>(sizeof(double))},
      values.data(),
      owner);
  array.attr("setflags")(false);
  return array;
}

py::array makeReadonlyU64View(const core::AlignedVector<std::uint64_t>& values, py::handle owner) {
  py::array array(
      py::dtype::of<std::uint64_t>(),
      {static_cast<py::ssize_t>(values.size())},
      {static_cast<py::ssize_t>(sizeof(std::uint64_t))},
      values.data(),
      owner);
  array.attr("setflags")(false);
  return array;
}

core::SimulationState makeUniformDarkMatterState(
    std::size_t particle_count,
    double mass_code,
    double box_size_mpc_comov) {
  core::SimulationState state;
  state.resizeParticles(particle_count);
  state.species.count_by_species = {static_cast<std::uint64_t>(particle_count), 0, 0, 0, 0};

  const double delta = particle_count > 0 ? box_size_mpc_comov / static_cast<double>(particle_count) : 0.0;
  for (std::size_t i = 0; i < particle_count; ++i) {
    state.particle_sidecar.particle_id[i] = static_cast<std::uint64_t>(i + 1);
    state.particle_sidecar.species_tag[i] = static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter);
    state.particle_sidecar.sfc_key[i] = static_cast<std::uint64_t>(i);
    state.particle_sidecar.particle_flags[i] = 0;
    state.particle_sidecar.owning_rank[i] = 0;

    state.particles.position_x_comoving[i] = delta * static_cast<double>(i);
    state.particles.position_y_comoving[i] = 0.5 * box_size_mpc_comov;
    state.particles.position_z_comoving[i] = 0.25 * box_size_mpc_comov;
    state.particles.velocity_x_peculiar[i] = 0.0;
    state.particles.velocity_y_peculiar[i] = 0.0;
    state.particles.velocity_z_peculiar[i] = 0.0;
    state.particles.mass_code[i] = mass_code;
    state.particles.time_bin[i] = 0;
  }

  state.rebuildSpeciesIndex();
  return state;
}

std::uint64_t diagnosticClassToInt(analysis::DiagnosticClass value) {
  return static_cast<std::uint64_t>(value);
}

}  // namespace
}  // namespace cosmosim::python

PYBIND11_MODULE(_cosmosim, module) {
  using cosmosim::analysis::DiagnosticsEngine;
  using cosmosim::analysis::RunHealthCounters;
  using cosmosim::core::FrozenConfig;
  using cosmosim::core::SimulationConfig;
  using cosmosim::core::SimulationMode;
  using cosmosim::core::SimulationState;
  using cosmosim::io::SnapshotReadResult;

  module.doc() = "CosmoSim Python bindings for snapshot analysis workflows.";

  py::enum_<SimulationMode>(module, "SimulationMode")
      .value("cosmo_cube", SimulationMode::kCosmoCube)
      .value("zoom_in", SimulationMode::kZoomIn)
      .value("isolated_galaxy", SimulationMode::kIsolatedGalaxy)
      .value("isolated_cluster", SimulationMode::kIsolatedCluster);

  py::class_<SimulationConfig>(module, "SimulationConfig")
      .def(py::init<>())
      .def_readwrite("schema_version", &SimulationConfig::schema_version)
      .def_property(
          "mode",
          [](const SimulationConfig& config) { return config.mode.mode; },
          [](SimulationConfig& config, SimulationMode mode) { config.mode.mode = mode; })
      .def_property(
          "box_size_mpc_comov",
          [](const SimulationConfig& config) { return config.cosmology.box_size_mpc_comoving; },
          [](SimulationConfig& config, double value) { config.cosmology.box_size_mpc_comoving = value; })
      .def_property(
          "output_directory",
          [](const SimulationConfig& config) { return config.output.output_directory; },
          [](SimulationConfig& config, const std::string& value) { config.output.output_directory = value; })
      .def_property(
          "run_name",
          [](const SimulationConfig& config) { return config.output.run_name; },
          [](SimulationConfig& config, const std::string& value) { config.output.run_name = value; });

  py::class_<FrozenConfig>(module, "FrozenConfig")
      .def_readonly("config", &FrozenConfig::config)
      .def_readonly("normalized_text", &FrozenConfig::normalized_text)
      .def_property_readonly(
          "config_hash_hex",
          [](const FrozenConfig& frozen) { return frozen.provenance.config_hash_hex; });

  py::class_<SimulationState>(module, "SimulationState")
      .def(py::init<>())
      .def("resize_particles", &SimulationState::resizeParticles)
      .def("rebuild_species_index", &SimulationState::rebuildSpeciesIndex)
      .def_property_readonly(
          "particle_count", [](const SimulationState& state) { return state.particles.size(); })
      .def("particle_id", [](SimulationState& state) {
        return cosmosim::python::makeReadonlyU64View(
            state.particle_sidecar.particle_id,
            py::cast(&state, py::return_value_policy::reference_internal));
      })
      .def("position_x_comoving", [](SimulationState& state) {
        return cosmosim::python::makeReadonlyDoubleView(
            state.particles.position_x_comoving,
            py::cast(&state, py::return_value_policy::reference_internal));
      })
      .def("position_y_comoving", [](SimulationState& state) {
        return cosmosim::python::makeReadonlyDoubleView(
            state.particles.position_y_comoving,
            py::cast(&state, py::return_value_policy::reference_internal));
      })
      .def("position_z_comoving", [](SimulationState& state) {
        return cosmosim::python::makeReadonlyDoubleView(
            state.particles.position_z_comoving,
            py::cast(&state, py::return_value_policy::reference_internal));
      })
      .def("velocity_x_peculiar", [](SimulationState& state) {
        return cosmosim::python::makeReadonlyDoubleView(
            state.particles.velocity_x_peculiar,
            py::cast(&state, py::return_value_policy::reference_internal));
      })
      .def("mass_code", [](SimulationState& state) {
        return cosmosim::python::makeReadonlyDoubleView(
            state.particles.mass_code,
            py::cast(&state, py::return_value_policy::reference_internal));
      });

  py::class_<RunHealthCounters>(module, "RunHealthCounters")
      .def_readonly("particle_count", &RunHealthCounters::particle_count)
      .def_readonly("cell_count", &RunHealthCounters::cell_count)
      .def_readonly("star_count", &RunHealthCounters::star_count)
      .def_readonly("ownership_invariants_ok", &RunHealthCounters::ownership_invariants_ok)
      .def_readonly("unique_particle_ids_ok", &RunHealthCounters::unique_particle_ids_ok)
      .def_readonly("non_finite_particles", &RunHealthCounters::non_finite_particles)
      .def_readonly("non_finite_cells", &RunHealthCounters::non_finite_cells);

  py::class_<SnapshotReadResult>(module, "SnapshotReadResult")
      .def_readonly("state", &SnapshotReadResult::state)
      .def_readonly("normalized_config_text", &SnapshotReadResult::normalized_config_text)
      .def_property_readonly("schema_name", [](const SnapshotReadResult& value) { return value.report.schema_name; })
      .def_property_readonly("schema_version", [](const SnapshotReadResult& value) { return value.report.schema_version; })
      .def_property_readonly(
          "present_aliases",
          [](const SnapshotReadResult& value) { return value.report.present_aliases; });

  py::class_<DiagnosticsEngine>(module, "DiagnosticsEngine")
      .def(py::init<SimulationConfig>())
      .def("compute_run_health", &DiagnosticsEngine::computeRunHealth)
      .def("compute_gas_xy_projection_density", [](const DiagnosticsEngine& engine, const SimulationState& state, std::size_t grid_n) {
        return engine.computeGasXyProjectionDensity(state, grid_n);
      })
      .def("generate_bundle_summary", [](const DiagnosticsEngine& engine, const SimulationState& state, std::uint64_t step_index, double scale_factor) {
        const auto bundle = engine.generateBundle(
            state,
            step_index,
            scale_factor,
            cosmosim::analysis::DiagnosticClass::kScienceLight);
        py::dict summary;
        summary["step_index"] = bundle.step_index;
        summary["scale_factor"] = bundle.scale_factor;
        summary["diagnostic_class"] = cosmosim::python::diagnosticClassToInt(bundle.diagnostic_class);
        summary["particle_count"] = bundle.health.particle_count;
        summary["power_bin_count"] = bundle.power_spectrum.size();
        summary["quicklook_grid_n"] = bundle.quicklook_grid_n;
        return summary;
      });

  module.def(
      "load_frozen_config",
      [](const std::filesystem::path& path) {
        return cosmosim::core::loadFrozenConfigFromFile(path, {});
      },
      py::arg("path"));

  module.def(
      "read_snapshot",
      [](const std::filesystem::path& path, const SimulationConfig& config) {
        return cosmosim::io::readGadgetArepoSnapshotHdf5(path, config, {});
      },
      py::arg("path"),
      py::arg("config"));

  module.def(
      "write_snapshot",
      [](const std::filesystem::path& path,
         const SimulationState& state,
         const SimulationConfig& config,
         const std::string& normalized_config_text,
         const std::string& git_sha) {
        cosmosim::io::SnapshotWritePayload payload;
        payload.state = &state;
        payload.config = &config;
        payload.normalized_config_text = normalized_config_text;
        payload.git_sha = git_sha;
        payload.provenance = cosmosim::core::makeProvenanceRecord(
            cosmosim::core::stableConfigHashHex(normalized_config_text),
            git_sha,
            0);
        cosmosim::io::writeGadgetArepoSnapshotHdf5(path, payload, {});
      },
      py::arg("path"),
      py::arg("state"),
      py::arg("config"),
      py::arg("normalized_config_text"),
      py::arg("git_sha") = "unknown");

  module.def(
      "make_uniform_dark_matter_state",
      &cosmosim::python::makeUniformDarkMatterState,
      py::arg("particle_count"),
      py::arg("mass_code"),
      py::arg("box_size_mpc_comov"));

  module.def("__version__", []() { return std::string("0.1.0"); });
}
