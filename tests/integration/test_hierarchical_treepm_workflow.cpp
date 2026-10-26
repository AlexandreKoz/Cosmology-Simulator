#include <algorithm>
#include <array>
#include <source_location>
#include <stdexcept>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>
#include <utility>
#include <sstream>
#include <functional>
#include <filesystem>
#include <string_view>
#include <unordered_map>

#include "cosmosim/cosmosim.hpp"
#include "cosmosim/core/build_config.hpp"
#include "cosmosim/io/restart_checkpoint.hpp"
#include "../support/test_temp_workspace.hpp"
#if COSMOSIM_ENABLE_MPI
#include <mpi.h>
#endif

namespace {
using namespace cosmosim;
void require(bool condition,const std::source_location location=std::source_location::current()) {
  if (!condition) throw std::runtime_error(std::string(location.file_name())+":"+
      std::to_string(location.line())+": hierarchical workflow invariant failed");
}
core::FrozenConfig config(int rung,int ranks) {
  std::ostringstream text;
  text << "[mode]\nmode = cosmo_cube\nic_file = generated\n"
       << "[cosmology]\nbox_size = 1 mpc\nomega_matter = 1\nomega_lambda = 0\nomega_baryon = 0\n"
       << "[numerics]\na_begin = 0.5\nintegrator_time_variable = code_time\n"
       << "t_code_begin = 0\nt_code_end = 0.00004\nmax_global_steps = 256\n"
       << "hierarchical_max_rung = " << rung << "\ntreepm_pm_grid = 16\n"
       << "gravity_softening = 1 kpc\ntreepm_asmth_cells = 1.25\ntreepm_rcut_cells = 6.25\n"
       << "[physics]\nenable_cooling = false\nenable_star_formation = false\nenable_feedback = false\n"
       << "enable_stellar_evolution = false\nenable_black_hole_agn = false\nenable_tracers = false\n"
       << "enable_metal_diffusion = false\n"
       << "[parallel]\nmpi_ranks_expected = " << ranks << "\ndecomposition_runtime_rebalance_enabled = false\n"
       << "[output]\nrun_name = hierarchical_treepm_qualification\nwrite_restarts = true\n"
       << "snapshot_interval_steps = 1\n";
  return core::loadFrozenConfigFromString(text.str(),"hierarchical_treepm_qualification");
}
core::SimulationState initialState(int rank,int ranks) {
  core::SimulationState state;
  std::vector<std::uint32_t> indices;
  for (std::uint32_t i=0;i<64U;++i) if (i%static_cast<unsigned>(ranks)==static_cast<unsigned>(rank)) indices.push_back(i);
  state.resizeParticles(indices.size());
  const double pi=std::acos(-1.0);
  for (std::size_t row=0;row<indices.size();++row) {
    const auto i=indices[row];
    const double x=(static_cast<double>(i%4U)+0.5)/4.0;
    state.particle_sidecar.particle_id[row]=1001U+i;
    state.particle_sidecar.species_tag[row]=static_cast<std::uint32_t>(core::ParticleSpecies::kDarkMatter);
    state.particle_sidecar.owning_rank[row]=static_cast<std::uint32_t>(rank);
    state.particles.position_x_comoving[row]=x+0.001*std::sin(2*pi*x);
    state.particles.position_y_comoving[row]=(static_cast<double>((i/4U)%4U)+0.5)/4.0;
    state.particles.position_z_comoving[row]=(static_cast<double>(i/16U)+0.5)/4.0;
    state.particles.velocity_x_peculiar[row]=100.0*std::sin(2*pi*x);
    state.particles.mass_code[row]=1e9;
    if (i<2U) { // Resolves a different short-force scale than the smooth mode.
      state.particles.position_x_comoving[row]=i==0U?0.4999:0.5001;
      state.particles.position_y_comoving[row]=0.5;state.particles.position_z_comoving[row]=0.5;
      state.particles.velocity_x_peculiar[row]=i==0U?1.0:-1.0;
    }
  }
  state.rebuildSpeciesIndex();
  return state;
}
double phaseError(const core::SimulationState& value,const core::SimulationState& reference) {
  require(value.particles.size()==reference.particles.size());
  std::unordered_map<std::uint64_t,std::size_t> row;
  for (std::size_t i=0;i<reference.particles.size();++i) row.emplace(reference.particle_sidecar.particle_id[i],i);
  double e2=0.0;
  for (std::size_t i=0;i<value.particles.size();++i) {
    const auto j=row.at(value.particle_sidecar.particle_id[i]);
    require(value.particles.mass_code[i]==reference.particles.mass_code[j]);
    const std::array v{value.particles.position_x_comoving[i],value.particles.position_y_comoving[i],value.particles.position_z_comoving[i]};
    const std::array r{reference.particles.position_x_comoving[j],reference.particles.position_y_comoving[j],reference.particles.position_z_comoving[j]};
    for (std::size_t axis=0;axis<3U;++axis) {
      double delta=v[axis]-r[axis];delta-=std::nearbyint(delta);e2+=delta*delta;
    }
    const std::array vu{value.particles.velocity_x_peculiar[i],value.particles.velocity_y_peculiar[i],value.particles.velocity_z_peculiar[i]};
    const std::array ru{reference.particles.velocity_x_peculiar[j],reference.particles.velocity_y_peculiar[j],reference.particles.velocity_z_peculiar[j]};
    for (std::size_t axis=0;axis<3U;++axis) e2+=std::pow((vu[axis]-ru[axis])/100.0,2);
  }
  return e2;
}
void assertSynchronized(const io::RestartReadResult& checkpoint,unsigned rung) {
  const auto& time=checkpoint.integrator_state;
  require(!time.inside_kdk_step && time.last_completed_restart_safe);
  // Global KDK tags each completed step by its output purpose
  // (kCheckpointPoint/kSnapshotPoint), all of which are synchronized and
  // restart-safe; a hierarchical block must close at its synchronization point.
  require(rung==0U ? core::isRestartSafeBoundary(time.last_completed_boundary_kind)
                   : time.last_completed_boundary_kind==core::StepBoundaryKind::kGlobalSynchronizationPoint);
  require(checkpoint.scheduler_state.current_tick%(1U<<rung)==0U);
  require(checkpoint.scheduler_state.current_tick==checkpoint.gas_cell_scheduler_state.current_tick);
  for (std::size_t i=0;i<checkpoint.state.particles.size();++i) {
    require(checkpoint.state.particleLastDriftTimeCode(i)==time.current_time_code);
    require(checkpoint.state.particleLastDriftScaleFactor(i)==time.current_scale_factor);
    require(checkpoint.scheduler_state.next_activation_tick[i]==checkpoint.scheduler_state.current_tick);
    require(checkpoint.scheduler_state.active_flag[i]==0U);
  }
}
std::string legacyTreePmDefaultsText(std::string text) {
  for (const std::string_view line : {"treepm_gaussian_pair_lookup_enabled = false\n",
                                     "treepm_full_mac_diagnostics = true\n",
                                     "treepm_full_residual_lookup_enabled = false\n",
                                     "treepm_leaf_simd_enabled = false\n",
                                     "treepm_radial_coefficients_enabled = false\n",
                                     "treepm_residual_block_size = 64\n",
                                     "treepm_tree_max_leaf_size = 16\n"}) {
    const auto position = text.find(line);
    require(position != std::string::npos);
    text.erase(position, line.size());
  }
  return text;
}

void testLegacyRestartConfigurationCompatibility(
    const core::FrozenConfig& frozen, const io::RestartReadResult& source,
    const std::filesystem::path& path, const std::function<void(const io::RestartReadResult&)>& resume) {
  core::HierarchicalTimeBinScheduler particles(source.scheduler_state.max_bin), cells(source.gas_cell_scheduler_state.max_bin);
  particles.importPersistentState(source.scheduler_state);
  cells.importPersistentState(source.gas_cell_scheduler_state);
  // Change only configuration serialization. All dynamics, IDs, scheduler,
  // source/PM ownership, stochastic and output epochs are the real dispatcher output.
  auto state = source.state;
  io::RestartWritePayload payload;
  payload.persistent_state.simulation_state = &state;
  payload.integrator_state = &source.integrator_state;
  payload.scheduler = &particles;
  payload.gas_cell_scheduler = &cells;
  payload.gravity_force_cache = &source.gravity_force_cache;
  payload.normalized_config_text = legacyTreePmDefaultsText(frozen.normalized_text);
  payload.normalized_config_hash_hex = core::stableConfigHashHex(payload.normalized_config_text);
  state.metadata.normalized_config_hash_hex = payload.normalized_config_hash_hex;
  state.metadata.normalized_config_hash = core::stableConfigHash(payload.normalized_config_text);
  payload.provenance = source.provenance;
  payload.provenance.config_hash_hex = payload.normalized_config_hash_hex;
  payload.provenance.normalized_config_hash_hex = payload.normalized_config_hash_hex;
  payload.provenance.normalized_config_sha256_hex = core::strongConfigHashSha256Hex(payload.normalized_config_text);
  payload.provenance.normalized_config = payload.normalized_config_text;
  payload.provenance.raw_input_config = payload.normalized_config_text;
  payload.distributed_gravity_state = source.distributed_gravity_state;
  payload.output_cadence_state = source.output_cadence_state;
  payload.stochastic_state = source.stochastic_state;
  io::writeRestartCheckpointHdf5(path, payload);
  auto legacy = io::readRestartCheckpointHdf5(path);
  const auto original_serialized = core::serializeProvenanceRecord(legacy.provenance);
  const auto record = io::restartConfigurationCompatibilityRecord(legacy, frozen);
  require(record.find("restart_config_compatibility_begin=treepm_normalized_defaults_v2") != std::string::npos);
  require(record.find("derived_normalized_hash_hex="+frozen.provenance.config_hash_hex) != std::string::npos);
  require(legacy.normalized_config_hash_hex == payload.normalized_config_hash_hex);
  require(core::serializeProvenanceRecord(legacy.provenance) == original_serialized);
  const auto refused = [&](const io::RestartReadResult& checkpoint, const core::FrozenConfig& config,
                           const std::string& context) {
    bool rejected = false;
    try { (void)io::restartConfigurationCompatibilityRecord(checkpoint, config); }
    catch (const std::exception& error) { rejected = std::string(error.what()).find(context) != std::string::npos; }
    require(rejected);
  };
  auto text = frozen.normalized_text;
  // Change all axes together: the incompatible candidate must itself satisfy
  // the spherical split/coarsest-cell guard, rather than fail in config parsing.
  for (const char axis : {'x', 'y', 'z'}) {
    const std::string key = std::string("treepm_pm_grid_n") + axis;
    const std::string mesh_line = key + " = 16\n";
    const auto mesh_position = text.find(mesh_line);
    require(mesh_position != std::string::npos);
    text.replace(mesh_position, mesh_line.size(), key + " = 32\n");
  }
  const auto different = core::loadFrozenConfigFromString(text, "compatibility_scientific_change");
  refused(legacy, different, "exact legacy projection");
  for (const std::string key : {"treepm_gaussian_pair_lookup_enabled = false", "treepm_full_mac_diagnostics = true"}) {
    text = frozen.normalized_text;
    const auto position = text.find(key);
    require(position != std::string::npos);
    text.replace(position, key.size(), key.find("lookup") != std::string::npos ?
        "treepm_gaussian_pair_lookup_enabled = true" : "treepm_full_mac_diagnostics = false");
    refused(legacy, core::loadFrozenConfigFromString(text, "compatibility_nonlegacy_flags"), "legacy defaults");
  }
  auto bad = legacy;
  bad.payload_integrity_verified = false;
  refused(bad, frozen, "integrity-verified v23");
  bad = legacy; bad.provenance.normalized_config.clear();
  refused(bad, frozen, "complete original normalized");
  bad = legacy; bad.normalized_config_text.clear();
  refused(bad, frozen, "non-empty normalized_config_text");
  bad = legacy; bad.integrator_state.inside_kdk_step = true;
  refused(bad, frozen, "closed global KDK");
  bad = legacy; bad.provenance.schema_version = "provenance_v6";
  refused(bad, frozen, "v7 provenance");
  bad = legacy; bad.diagnostics.restart_schema_version = 22U;
  refused(bad, frozen, "integrity-verified v23");
  // A partially introduced key is a different historical normalization revision.
  bad = legacy;
  bad.normalized_config_text += "\n[numerics]\ntreepm_full_mac_diagnostics = true\n";
  bad.normalized_config_hash_hex = core::stableConfigHashHex(bad.normalized_config_text);
  bad.provenance.config_hash_hex = bad.normalized_config_hash_hex;
  bad.provenance.normalized_config_hash_hex = bad.normalized_config_hash_hex;
  bad.provenance.normalized_config_sha256_hex = core::strongConfigHashSha256Hex(bad.normalized_config_text);
  bad.provenance.normalized_config = bad.normalized_config_text;
  bad.state.metadata.normalized_config_hash_hex = bad.normalized_config_hash_hex;
  bad.state.metadata.normalized_config_hash = core::stableConfigHash(bad.normalized_config_text);
  refused(bad, frozen, "exact legacy projection");
  resume(legacy);
}

void run(int rank,int ranks,std::uint64_t token) {
  auto scratch=test_support::TestTempWorkspace::createMpiShared("hierarchical_treepm_workflow",token,rank==0);
  const auto initial=initialState(rank,ranks);
  parallel::MpiContext mpi;
  const auto hierarchical_config=config(2,ranks),global_config=config(0,ranks);
  workflows::ReferenceWorkflowRunner hierarchical(hierarchical_config),global(global_config);
  auto execute=[&](const workflows::ReferenceWorkflowRunner& runner,const char* label,double dt,
                   const io::RestartReadResult* restart=nullptr,std::uint64_t max_steps=0U) {
    const auto report=runner.run(scratch.root()/label,workflows::ReferenceWorkflowOptions{
        .dt_time_code=dt,.write_outputs=true,.initial_state_override=restart?nullptr:&initial,
        .restart_state_override=restart,.max_steps_override=max_steps});
    require(report.restart_roundtrip_ok);
    return std::pair{report,io::readRestartCheckpointHdf5(report.restart_path)};
  };
  auto reference=execute(global,"global_fine",2.5e-6);
  auto coarse=execute(hierarchical,"hierarchical_coarse",2e-5);
  auto refined=execute(hierarchical,"hierarchical_refined",1e-5);
  auto replay=execute(hierarchical,"hierarchical_replay",2e-5);
  for (const auto* checkpoint:{&reference.second,&coarse.second,&refined.second,&replay.second}) {
    require(std::abs(checkpoint->integrator_state.current_time_code-4e-5)<1e-14);
    require(checkpoint->integrator_state.current_scale_factor>0.5);
  }
  assertSynchronized(reference.second,0U);assertSynchronized(coarse.second,2U);
  assertSynchronized(refined.second,2U);
  require(coarse.first.final_state_digest==replay.first.final_state_digest);
  const double coarse_error=mpi.allreduceSumDouble(phaseError(coarse.second.state,reference.second.state));
  const double refined_error=mpi.allreduceSumDouble(phaseError(refined.second.state,reference.second.state));
  // Short cosmological fixture: 1 Mpc periodic box, velocities normalized by
  // 100 km/s. This gate checks trajectory refinement, not long-time P(k).
  require(std::sqrt(refined_error/64.0)<1e-3);
  require(refined_error<=coarse_error+1e-18);
  auto optimized_text=hierarchical_config.normalized_text;
  for (const std::string key : {"treepm_full_residual_lookup_enabled", "treepm_leaf_simd_enabled",
                               "treepm_radial_coefficients_enabled", "treepm_identical_source_tree_reuse_enabled",
                               "treepm_topology_refit_enabled"}) {
    const auto position=optimized_text.find(key+" = false");
    require(position != std::string::npos);
    optimized_text.replace(position,(key+" = false").size(),key+" = true");
  }
  workflows::ReferenceWorkflowRunner optimized(core::loadFrozenConfigFromString(optimized_text,"optimized_hierarchical"));
  auto optimized_full=execute(optimized,"optimized_hierarchical_full",2e-5);
  assertSynchronized(optimized_full.second,2U);
  require(mpi.allreduceSumDouble(phaseError(optimized_full.second.state,coarse.second.state))<1e-18);
  auto optimized_first=execute(optimized,"optimized_hierarchical_first",2e-5,nullptr,1U);
  auto optimized_resume=execute(optimized,"optimized_hierarchical_resume",2e-5,&optimized_first.second);
  require(mpi.allreduceSumDouble(phaseError(optimized_resume.second.state,optimized_full.second.state))<1e-22);
  auto first=execute(hierarchical,"restart_first",2e-5,nullptr,1U);
  assertSynchronized(first.second,2U);
  auto resumed=execute(hierarchical,"restart_resumed",2e-5,&first.second);
  assertSynchronized(resumed.second,2U);
  require(mpi.allreduceSumDouble(phaseError(resumed.second.state,coarse.second.state))<1e-22);
  require(resumed.second.scheduler_state.current_tick==coarse.second.scheduler_state.current_tick);
  // Real global KDK checkpoint migration, then real resumed dynamics.
  auto global_first = execute(global, "global_legacy_first", 2.5e-6, nullptr, 1U);
  testLegacyRestartConfigurationCompatibility(global_config, global_first.second,
      scratch.root()/("legacy_rank_"+std::to_string(rank)+".hdf5"), [&](const io::RestartReadResult& legacy) {
        auto migrated = execute(global, "global_legacy_resumed", 2.5e-6, &legacy);
        require(!migrated.first.restart_config_compatibility_record.empty());
        require(migrated.second.provenance.derived_runtime_state.find("treepm_normalized_defaults_v2") != std::string::npos);
        require(migrated.second.normalized_config_hash_hex == global_config.provenance.config_hash_hex);
        require(mpi.allreduceSumDouble(phaseError(migrated.second.state, reference.second.state)) < 1e-22);
        require(io::restartConfigurationCompatibilityRecord(migrated.second, global_config) ==
                migrated.first.restart_config_compatibility_record);
      });
  // All ranks corrupt the same admission invariant and verify collective refusal.
  const auto refused_resume = [&](io::RestartReadResult invalid, const char* label, const char* expected) {
    bool rejected = false;
    try { (void)execute(hierarchical, label, 2e-5, &invalid); }
    catch (const std::exception& error) { rejected = std::string(error.what()).find(expected) != std::string::npos; }
    require(mpi.allreduceSumUint64(rejected ? 1U : 0U) == static_cast<std::uint64_t>(ranks));
  };
  auto invalid = first.second;
  invalid.integrator_state.inside_kdk_step = true;
  refused_resume(invalid, "reject_open_kdk", "closed collisionless synchronized KDK");
  invalid = first.second;
  if (!invalid.scheduler_state.active_flag.empty()) invalid.scheduler_state.active_flag.front() = 1U;
  refused_resume(invalid, "reject_active_interval", "open particle interval");
  bool saw_fine=false,saw_pm_endpoint=false;
  for (const auto& record:coarse.first.treepm_cadence_records) {
    saw_fine=saw_fine||record.short_range_only;
    saw_pm_endpoint=saw_pm_endpoint||record.refreshed_long_range_field;
    if (record.short_range_only) require(!record.refreshed_long_range_field);
  }
  require(saw_fine && saw_pm_endpoint);
  // Multirate outer PM level (2026-10-26, NOT RUN): treepm_pm_reuse_policy =
  // multirate_outer_kick with two coarse blocks per PM cycle. Provisional,
  // unqualified gates: the integrator must stay synchronized, converge toward
  // the fine global reference with refinement, refresh PM fewer times than
  // the exact hierarchical scheme, and resume bit-identically from a
  // checkpoint written at a cycle end.
  {
    auto multirate_text=hierarchical_config.normalized_text;
    for (const auto& [from,to] : {std::pair<std::string,std::string>{"treepm_pm_reuse_policy = exact","treepm_pm_reuse_policy = multirate_outer_kick"},
                                   std::pair<std::string,std::string>{"treepm_update_cadence_steps = 1","treepm_update_cadence_steps = 2"}}) {
      const auto position=multirate_text.find(from);
      require(position!=std::string::npos);
      multirate_text.replace(position,from.size(),to);
    }
    const auto multirate_config=core::loadFrozenConfigFromString(multirate_text,"multirate_hierarchical");
    workflows::ReferenceWorkflowRunner multirate(multirate_config);
    auto multirate_coarse=execute(multirate,"multirate_coarse",1e-5);   // 4 blocks = 2 PM cycles
    auto multirate_refined=execute(multirate,"multirate_refined",5e-6); // 8 blocks = 4 PM cycles
    for (const auto* checkpoint:{&multirate_coarse.second,&multirate_refined.second}) {
      require(std::abs(checkpoint->integrator_state.current_time_code-4e-5)<1e-14);
      assertSynchronized(*checkpoint,2U);
    }
    const double multirate_coarse_error=mpi.allreduceSumDouble(phaseError(multirate_coarse.second.state,reference.second.state));
    const double multirate_refined_error=mpi.allreduceSumDouble(phaseError(multirate_refined.second.state,reference.second.state));
    require(multirate_refined_error<=multirate_coarse_error+1e-18);
    require(std::sqrt(multirate_refined_error/64.0)<1e-3);
    std::uint64_t multirate_refreshes=0U, exact_refreshes=0U;
    for (const auto& record:multirate_refined.first.treepm_cadence_records) multirate_refreshes+=record.refreshed_long_range_field?1U:0U;
    auto exact_same_dt=execute(hierarchical,"multirate_exact_comparison",5e-6);
    for (const auto& record:exact_same_dt.first.treepm_cadence_records) exact_refreshes+=record.refreshed_long_range_field?1U:0U;
    require(multirate_refreshes<exact_refreshes);
    // Checkpoint at a cycle end (2 blocks), then resume: identical to the
    // uninterrupted run because cycles always open at restart-safe boundaries.
    auto multirate_first=execute(multirate,"multirate_first",1e-5,nullptr,2U);
    assertSynchronized(multirate_first.second,2U);
    auto multirate_resumed=execute(multirate,"multirate_resumed",1e-5,&multirate_first.second);
    require(mpi.allreduceSumDouble(phaseError(multirate_resumed.second.state,multirate_coarse.second.state))<1e-22);
    // A one-block budget truncates the cycle to one block: still synchronized.
    auto multirate_single=execute(multirate,"multirate_single_block",1e-5,nullptr,1U);
    assertSynchronized(multirate_single.second,2U);
    if (rank==0) std::cout<<"hierarchical_multirate_outer_kick ranks="<<ranks
        <<" coarse_phase_l2="<<std::sqrt(multirate_coarse_error/64.0)
        <<" refined_phase_l2="<<std::sqrt(multirate_refined_error/64.0)
        <<" pm_refreshes="<<multirate_refreshes<<" exact_pm_refreshes="<<exact_refreshes<<'\n';
  }
  double local_mass=0.0;
  for (const double mass:coarse.second.state.particles.mass_code) local_mass+=mass;
  require(mpi.allreduceSumDouble(local_mass)==64e9);
  if (rank==0) std::cout<<"hierarchical_treepm_qualification ranks="<<ranks
      <<" coarse_phase_l2="<<std::sqrt(coarse_error/64.0)<<" refined_phase_l2="<<std::sqrt(refined_error/64.0)<<'\n';
#if COSMOSIM_ENABLE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}
}  // namespace
int main() {
  int rank=0,ranks=1;
  auto token=cosmosim::test_support::TestTempWorkspace::uniqueRunToken();
#if COSMOSIM_ENABLE_MPI
  MPI_Init(nullptr,nullptr);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
  MPI_Bcast(&token,1,MPI_UINT64_T,0,MPI_COMM_WORLD);
#endif
  try {run(rank,ranks,token);}
  catch (const std::exception& error) {
    std::cerr<<"hierarchical TreePM qualification rank="<<rank<<": "<<error.what()<<'\n';
#if COSMOSIM_ENABLE_MPI
    MPI_Abort(MPI_COMM_WORLD,1);
#endif
    return 1;
  }
#if COSMOSIM_ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
