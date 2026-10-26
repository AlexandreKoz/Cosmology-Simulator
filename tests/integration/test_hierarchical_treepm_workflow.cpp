#include <algorithm>
#include <array>
#include <source_location>
#include <stdexcept>
#include <cmath>
#include <cstdint>
#include <iomanip>
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
// skewed: rank 0 owns almost every particle, ranks 1..P-2 own one each and
// rank P-1 owns none (empty rank) -- uneven-ownership MPI fixture.
core::SimulationState initialState(int rank,int ranks,bool skewed=false) {
  core::SimulationState state;
  std::vector<std::uint32_t> indices;
  const auto owner=[&](std::uint32_t i) {
    if (!skewed) return static_cast<int>(i%static_cast<unsigned>(ranks));
    return static_cast<int>(i)<ranks-2 ? static_cast<int>(i)+1 : 0;
  };
  for (std::uint32_t i=0;i<64U;++i) if (owner(i)==rank) indices.push_back(i);
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
                   const io::RestartReadResult* restart=nullptr,std::uint64_t max_steps=0U,
                   std::uint64_t forced_multirate_rollbacks=0U,
                   const std::function<void(workflows::ReferenceWorkflowOptions&)>& customize={}) {
    workflows::ReferenceWorkflowOptions options{
        .dt_time_code=dt,.write_outputs=true,.initial_state_override=restart?nullptr:&initial,
        .restart_state_override=restart,.max_steps_override=max_steps,
        .multirate_test_forced_cycle_rollbacks=forced_multirate_rollbacks};
    if (customize) customize(options);
    const auto report=runner.run(scratch.root()/label,options);
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
  // Multirate outer PM level (2026-10-26 round 3; AUTHORED, NOT RUN):
  // treepm_pm_reuse_policy = multirate_outer_kick. Owner gates for the
  // stabilized contract in docs/time_integration.md:
  //  * timestep invariants: no block quantum exceeds the global particle
  //    criterion, blocks tile each PM cycle exactly, cycles close exactly on
  //    their planned end and the run ends exactly on t_end;
  //  * run-segmentation invariance: one run equals any split by step budget
  //    (soft stop at cycle ends), including budgets that end inside a cycle;
  //  * transactional rollback: forced rejections leave the trajectory
  //    bit-identical to an uninterrupted run;
  //  * convergence toward the fine global reference and toward the exact
  //    hierarchical scheme at the same dt (splitting error) under refinement;
  //  * cross-rank agreement of the block schedule (np 1/2/4 registrations).
  {
    const auto multirate_config_with=[&](int blocks_per_cycle,const std::string& base,const char* label) {
      auto text=base;
      for (const auto& [from,to] : {std::pair<std::string,std::string>{"treepm_pm_reuse_policy = exact","treepm_pm_reuse_policy = multirate_outer_kick"},
                                     std::pair<std::string,std::string>{"treepm_update_cadence_steps = 1","treepm_update_cadence_steps = "+std::to_string(blocks_per_cycle)}}) {
        const auto position=text.find(from);
        require(position!=std::string::npos);
        text.replace(position,from.size(),to);
      }
      return core::loadFrozenConfigFromString(text,label);
    };
    workflows::ReferenceWorkflowRunner multirate(multirate_config_with(2,hierarchical_config.normalized_text,"multirate_m2"));
    workflows::ReferenceWorkflowRunner multirate4(multirate_config_with(4,hierarchical_config.normalized_text,"multirate_m4"));
    workflows::ReferenceWorkflowRunner multirate64(multirate_config_with(64,hierarchical_config.normalized_text,"multirate_m64"));
    workflows::ReferenceWorkflowRunner multirate4_optimized(multirate_config_with(4,optimized_text,"multirate_m4_optimized"));
    using Record=workflows::ReferenceWorkflowReport::HierarchicalBlockRecord;
    const auto same_on_all_ranks=[&](double value) {
      return mpi.allreduceMinDouble(value)==-mpi.allreduceMinDouble(-value);
    };
    // Per-segment invariants: criteria honoured, exact tiling, closed at the end.
    const auto check_blocks=[&](const workflows::ReferenceWorkflowReport& report) {
      const auto& records=report.hierarchical_block_records;
      require(!records.empty() && records.size()<512U);
      require(report.completed_steps==records.size());
      require(mpi.allreduceMinUint64(records.size())==mpi.allreduceMaxUint64(records.size()));
      std::uint64_t closed=0U;
      for (std::size_t i=0;i<records.size();++i) {
        const Record& block=records[i];
        require(block.multirate);
        require(block.quantum_time_code>0.0 &&
                block.quantum_time_code<=block.global_min_particle_criterion_time_code*(1.0+1e-12));
        require(block.time_end_code>block.time_begin_code);
        require(block.pm_cycle_opening_time_code<=block.time_begin_code && block.time_end_code<=block.pm_cycle_end_time_code);
        require(block.pm_cycle_end_time_code<=4e-5);
        // The accepted PM interval honours the displacement criterion at every block.
        require(block.pm_cycle_end_time_code-block.pm_cycle_opening_time_code<=
                block.pm_displacement_limit_time_code*(1.0+1e-12));
        require(same_on_all_ranks(block.time_end_code) && same_on_all_ranks(block.quantum_time_code));
        if (block.opens_pm_cycle) require(block.pm_cycle_block_index==0U && block.time_begin_code==block.pm_cycle_opening_time_code);
        if (block.closes_pm_cycle) { require(block.time_end_code==block.pm_cycle_end_time_code); ++closed; }
        if (i>0U) {
          const Record& previous=records[i-1U];
          require(block.time_begin_code==previous.time_end_code);
          require(block.opens_pm_cycle==previous.closes_pm_cycle);
          if (!block.opens_pm_cycle) {
            require(block.pm_cycle_index==previous.pm_cycle_index && block.pm_cycle_block_index==previous.pm_cycle_block_index+1U);
            require(block.pm_cycle_opening_time_code==previous.pm_cycle_opening_time_code &&
                    block.pm_cycle_end_time_code==previous.pm_cycle_end_time_code);
          }
        }
      }
      require(records.front().opens_pm_cycle && records.back().closes_pm_cycle);
      require(report.multirate_pm_cycles_completed==closed);
      // Soft stop: the budget is honoured or exceeded only to close a cycle.
      require(report.completed_steps<=report.requested_steps+report.multirate_soft_stop_extension_blocks);
      if (report.multirate_soft_stop_extension_blocks>0U) {
        require(report.completed_steps==report.requested_steps+report.multirate_soft_stop_extension_blocks);
      }
    };
    const auto refreshes=[](const workflows::ReferenceWorkflowReport& report) {
      std::uint64_t count=0U;
      for (const auto& record:report.treepm_cadence_records) count+=record.refreshed_long_range_field?1U:0U;
      return count;
    };
    struct Segmented {
      std::pair<workflows::ReferenceWorkflowReport,io::RestartReadResult> last;
      std::vector<Record> records;
      std::uint64_t refreshes=0U, extension_blocks=0U, segments=0U;
    };
    // Runs the listed budgets back to back from checkpoints, then the rest.
    const auto run_segments=[&](const workflows::ReferenceWorkflowRunner& runner,const std::string& label,double dt,
                                const std::vector<std::uint64_t>& budgets) {
      Segmented result;
      const io::RestartReadResult* restart=nullptr;
      std::size_t next=0U;
      while (true) {
        const std::uint64_t budget=next<budgets.size()?budgets[next]:0U;
        auto segment=execute(runner,(label+"_segment_"+std::to_string(result.segments)).c_str(),dt,restart,budget);
        check_blocks(segment.first);
        assertSynchronized(segment.second,2U);
        result.records.insert(result.records.end(),segment.first.hierarchical_block_records.begin(),
                              segment.first.hierarchical_block_records.end());
        result.refreshes+=refreshes(segment.first);
        result.extension_blocks+=segment.first.multirate_soft_stop_extension_blocks;
        ++result.segments;++next;
        result.last=std::move(segment);
        restart=&result.last.second;
        if (result.last.second.integrator_state.current_time_code>=4e-5) break;
        require(result.segments<64U);
      }
      return result;
    };
    const auto same_schedule=[&](const std::vector<Record>& a,const std::vector<Record>& b) {
      require(a.size()==b.size());
      for (std::size_t i=0;i<a.size();++i) {
        require(a[i].time_begin_code==b[i].time_begin_code && a[i].time_end_code==b[i].time_end_code);
        require(a[i].coarse_interval_time_code==b[i].coarse_interval_time_code && a[i].quantum_time_code==b[i].quantum_time_code);
        require(a[i].opens_pm_cycle==b[i].opens_pm_cycle && a[i].closes_pm_cycle==b[i].closes_pm_cycle);
        require(a[i].pm_cycle_opening_time_code==b[i].pm_cycle_opening_time_code &&
                a[i].pm_cycle_end_time_code==b[i].pm_cycle_end_time_code);
      }
    };
    const auto same_endpoint=[&](const io::RestartReadResult& value,const io::RestartReadResult& reference,double tolerance) {
      require(value.integrator_state.current_time_code==reference.integrator_state.current_time_code);
      require(value.integrator_state.current_scale_factor==reference.integrator_state.current_scale_factor);
      require(value.scheduler_state.current_tick==reference.scheduler_state.current_tick);
      require(mpi.allreduceSumDouble(phaseError(value.state,reference.state))<tolerance);
    };

    // (1) Uninterrupted reference runs: invariants, exact end time, fewer PM solves.
    auto multirate_coarse=execute(multirate,"multirate_coarse",1e-5);
    auto multirate_refined=execute(multirate,"multirate_refined",5e-6);
    for (const auto* outcome:{&multirate_coarse,&multirate_refined}) {
      check_blocks(outcome->first);
      // The final cycle is clipped to t_end and closes exactly on it.
      require(outcome->second.integrator_state.current_time_code==4e-5);
      assertSynchronized(outcome->second,2U);
    }
    const double multirate_coarse_error=mpi.allreduceSumDouble(phaseError(multirate_coarse.second.state,reference.second.state));
    const double multirate_refined_error=mpi.allreduceSumDouble(phaseError(multirate_refined.second.state,reference.second.state));
    require(multirate_refined_error<=multirate_coarse_error+1e-18);
    require(std::sqrt(multirate_refined_error/64.0)<1e-3);
    const std::uint64_t multirate_refreshes=refreshes(multirate_refined.first);
    auto exact_same_dt=execute(hierarchical,"multirate_exact_comparison",5e-6);
    const std::uint64_t exact_refreshes=refreshes(exact_same_dt.first);
    require(multirate_refreshes<exact_refreshes);

    // (2) Run-segmentation invariance (M = 4, nominally 8 blocks = 2 cycles).
    auto full4=run_segments(multirate4,"multirate4_full",5e-6,{});
    require(full4.segments==1U && full4.extension_blocks==0U);
    bool has_multi_block_cycle=false;
    for (const auto& block:full4.records) has_multi_block_cycle=has_multi_block_cycle||block.pm_cycle_block_index>0U;
    const std::vector<std::pair<std::string,std::vector<std::uint64_t>>> splits{{"multirate4_4x2",{4U,4U}},
                                        {"multirate4_2x4",{2U,2U,2U,2U}},
                                        {"multirate4_1x8",{1U,1U,1U,1U,1U,1U,1U,1U}},
                                        {"multirate4_3_rest",{3U}},
                                        {"multirate4_5_rest",{5U}}};
    for (const auto& [label,budgets] : splits) {
      auto split=run_segments(multirate4,label,5e-6,budgets);
      same_schedule(split.records,full4.records);
      same_endpoint(split.last.second,full4.last.second,1e-22);
      require(split.refreshes==full4.refreshes);
      // A budget that ends inside a multi-block cycle must soft-stop.
      if (has_multi_block_cycle && budgets.front()==1U) require(split.extension_blocks>0U);
    }
    // Resume at a cycle boundary with the M = 2 configuration (2 blocks/cycle).
    auto multirate_first=execute(multirate,"multirate_first",1e-5,nullptr,2U);
    assertSynchronized(multirate_first.second,2U);
    auto multirate_resumed=execute(multirate,"multirate_resumed",1e-5,&multirate_first.second);
    same_endpoint(multirate_resumed.second,multirate_coarse.second,1e-22);
    // A one-block budget no longer truncates the PM cycle: it soft-stops at
    // the cycle end and reports the extension.
    auto multirate_single=execute(multirate,"multirate_single_block",1e-5,nullptr,1U);
    check_blocks(multirate_single.first);
    assertSynchronized(multirate_single.second,2U);
    require(multirate_single.first.completed_steps>=1U);
    require(multirate_single.first.hierarchical_block_records.back().closes_pm_cycle);

    // (3) Transactional rollback: forced rejections at interior boundaries are
    // retried with the identical outer interval and must be invisible.
    if (has_multi_block_cycle) {
      auto forced=execute(multirate4,"multirate4_forced_rollback",5e-6,nullptr,0U,2U);
      check_blocks(forced.first);
      require(forced.first.multirate_forced_test_rollbacks>0U);
      require(forced.first.multirate_pm_cycle_rollbacks>=forced.first.multirate_forced_test_rollbacks);
      require(forced.first.multirate_rejected_blocks>=forced.first.multirate_forced_test_rollbacks);
      // Rejected attempts are excluded from the trajectory diagnostics.
      require(refreshes(forced.first)==full4.refreshes);
      require(forced.first.treepm_long_range_refresh_count==full4.last.first.treepm_long_range_refresh_count);
      same_schedule(forced.first.hierarchical_block_records,full4.records);
      same_endpoint(forced.second,full4.last.second,1e-22);
      // Generation-keyed tree reuse/refit may sum in a different order after
      // a rollback: roundoff-level agreement, identical schedule.
      auto optimized_plain=execute(multirate4_optimized,"multirate4_optimized_plain",5e-6);
      auto optimized_forced=execute(multirate4_optimized,"multirate4_optimized_forced",5e-6,nullptr,0U,2U);
      check_blocks(optimized_forced.first);
      require(optimized_forced.first.multirate_forced_test_rollbacks>0U);
      same_schedule(optimized_forced.first.hierarchical_block_records,optimized_plain.first.hierarchical_block_records);
      same_endpoint(optimized_forced.second,optimized_plain.second,1e-18);
    }

    // (4) Boundary-limited cycle: M = 64 exceeds the run; the cycle is clipped
    // to t_end (or shortened by PM displacement) and never truncated by budget.
    auto multirate_long=execute(multirate64,"multirate64",1e-5);
    check_blocks(multirate_long.first);
    require(multirate_long.second.integrator_state.current_time_code==4e-5);
    require(multirate_long.first.hierarchical_block_records.back().pm_cycle_end_time_code==4e-5);

    // (5) Convergence under refinement (M = 4): error against the fine global
    // reference and splitting error against exact hierarchical at the same dt.
    std::array<double,3> global_error{},splitting_error{};
    const std::array<double,3> ladder{1e-5,5e-6,2.5e-6};
    for (std::size_t level=0;level<ladder.size();++level) {
      const std::string suffix=std::to_string(level);
      auto mr=execute(multirate4,("multirate4_convergence_"+suffix).c_str(),ladder[level]);
      auto ex=execute(hierarchical,("exact_convergence_"+suffix).c_str(),ladder[level]);
      check_blocks(mr.first);
      global_error[level]=mpi.allreduceSumDouble(phaseError(mr.second.state,reference.second.state));
      splitting_error[level]=mpi.allreduceSumDouble(phaseError(mr.second.state,ex.second.state));
    }
    for (std::size_t level=1;level<ladder.size();++level) {
      require(global_error[level]<=global_error[level-1U]+1e-18);
      // Below roundoff the splitting error carries no order information.
      if (splitting_error[level-1U]>1e-24) require(splitting_error[level]<splitting_error[level-1U]);
    }

    // (6) Genuine PM-displacement rejection (no test seam). Particles start
    // at rest, so the RMS displacement criterion is unbounded when the first
    // cycle opens (T = 4 blocks) and finite once gravity has accelerated them.
    // The criterion is linear in its factor and the dynamics before the first
    // interior boundary t1 do not depend on it, so a calibration run at factor
    // 1 (bound D1 at t1) fixes factor = 1.5 t1 / D1: the interior re-check at
    // t1 then sees 1.5 t1 < T and must reject, roll back and retry with
    // T = 1.5 t1, which the same re-check accepts.
    {
      const auto set_key=[](std::string text,const std::string& key,const std::string& value) {
        const auto position=text.find("\n"+key+" = ");
        require(position!=std::string::npos);
        const auto line_end=text.find('\n',position+1U);
        text.replace(position+1U,line_end-(position+1U),key+" = "+value);
        return text;
      };
      const auto rejecting_runner=[&](double factor) {
        std::ostringstream value;
        value<<std::setprecision(17)<<factor;
        auto text=set_key(hierarchical_config.normalized_text,"hierarchical_coarse_displacement_criterion","rms_split_scale");
        text=set_key(text,"hierarchical_coarse_rms_displacement_factor",value.str());
        return workflows::ReferenceWorkflowRunner(multirate_config_with(4,text,"multirate_m4_genuine_rejection"));
      };
      auto at_rest=initial;
      for (auto* lane:{&at_rest.particles.velocity_x_peculiar,&at_rest.particles.velocity_y_peculiar,&at_rest.particles.velocity_z_peculiar}) {
        std::fill(lane->begin(),lane->end(),0.0);
      }
      const auto run_at_rest=[&](const workflows::ReferenceWorkflowRunner& runner,const char* label) {
        const auto report=runner.run(scratch.root()/label,workflows::ReferenceWorkflowOptions{
            .dt_time_code=1e-5,.write_outputs=true,.initial_state_override=&at_rest});
        require(report.restart_roundtrip_ok);
        return std::pair{report,io::readRestartCheckpointHdf5(report.restart_path)};
      };
      const auto calibration_runner=rejecting_runner(1.0);
      auto calibration=run_at_rest(calibration_runner,"multirate_rejection_calibration");
      check_blocks(calibration.first);
      const Record* first_interior=nullptr;
      for (const auto& block:calibration.first.hierarchical_block_records) {
        if (block.time_begin_code>0.0) { first_interior=&block; break; }
      }
      require(first_interior!=nullptr && std::isfinite(first_interior->pm_displacement_limit_time_code) &&
              first_interior->pm_displacement_limit_time_code>0.0);
      const double factor=1.5*first_interior->time_begin_code/first_interior->pm_displacement_limit_time_code;
      require(std::isfinite(factor) && factor>0.0);
      // factor >= 1: the calibration run itself already rejected at t1.
      const auto rejecting=rejecting_runner(std::min(factor,1.0));
      const auto run_rejecting=[&](const char* label) { return run_at_rest(rejecting,label); };
      auto genuine=run_rejecting("multirate_genuine_rejection");
      check_blocks(genuine.first);
      require(genuine.first.multirate_forced_test_rollbacks==0U);
      require(genuine.first.multirate_pm_cycle_rollbacks>0U);
      require(genuine.first.multirate_rejected_blocks>0U);
      // Interior blocks are short-range only: a rejected attempt never solved PM.
      require(genuine.first.multirate_rejected_pm_refreshes==0U);
      // Trajectory diagnostics exclude the rejected attempts.
      require(refreshes(genuine.first)==genuine.first.treepm_long_range_refresh_count);
      bool retried_cycle_published=false;
      for (const auto& block:genuine.first.hierarchical_block_records) retried_cycle_published=retried_cycle_published||block.pm_cycle_retries>0U;
      require(retried_cycle_published);
      require(genuine.second.integrator_state.current_time_code==4e-5);
      assertSynchronized(genuine.second,2U);
      // Rollback/retry is deterministic.
      auto replayed=run_rejecting("multirate_genuine_rejection_replay");
      require(replayed.first.final_state_digest==genuine.first.final_state_digest);
      same_schedule(replayed.first.hierarchical_block_records,genuine.first.hierarchical_block_records);
      require(replayed.first.multirate_pm_cycle_rollbacks==genuine.first.multirate_pm_cycle_rollbacks);
    }

    // (7) Memory-pressure fallback. Snapshot admission is reduced
    // collectively: a refusal (governor or std::bad_alloc) on ANY rank must
    // open the same single-block cycle on EVERY rank, every fallback block
    // must still honour the particle criteria and the PM displacement bound
    // (check_blocks), no PM half-kick may be duplicated (one PM refresh per
    // closed cycle plus at most the segment bootstrap), and no governed
    // snapshot storage may outlive the segment.
    if (has_multi_block_cycle) {
      const auto refusal_case=[&](const char* label,std::vector<int> refusing_ranks,bool allocation_failure,
                                  std::uint64_t refusals) {
        auto fallback=execute(multirate4,label,5e-6,nullptr,0U,0U,[&](workflows::ReferenceWorkflowOptions& options) {
          options.multirate_test_snapshot_refusals=refusals;
          options.multirate_test_snapshot_refusal_ranks=refusing_ranks;
          options.multirate_test_snapshot_allocation_failure=allocation_failure;
        });
        check_blocks(fallback.first);
        require(fallback.first.multirate_memory_limited_cycles>0U);
        require(mpi.allreduceMinUint64(fallback.first.multirate_memory_limited_cycles)==
                mpi.allreduceMaxUint64(fallback.first.multirate_memory_limited_cycles));
        bool memory_limited_block=false;
        for (const auto& block:fallback.first.hierarchical_block_records) {
          if (block.pm_cycle_length_limiter!="rollback_memory") continue;
          memory_limited_block=true;
          require(block.opens_pm_cycle && block.closes_pm_cycle && block.pm_cycle_block_index==0U);
        }
        require(memory_limited_block);
        const auto pm_refreshes=fallback.first.treepm_long_range_refresh_count;
        const auto cycles=fallback.first.multirate_pm_cycles_completed;
        require(pm_refreshes>=cycles && pm_refreshes<=cycles+1U);
        require(fallback.first.multirate_rollback_snapshot_retained_bytes==0U);
        require(fallback.second.integrator_state.current_time_code==4e-5);
        assertSynchronized(fallback.second,2U);
        return fallback;
      };
      // Every admission refused on every rank: only single-block cycles.
      auto all_refused=refusal_case("multirate4_refused_all",{},false,1000000U);
      for (const auto& block:all_refused.first.hierarchical_block_records) {
        require(block.opens_pm_cycle && block.closes_pm_cycle);
      }
      require(all_refused.first.multirate_rollback_snapshot_peak_bytes==0U);
      // Refusal on rank 0 only, on several ranks, and an allocation failure
      // on the last rank: identical fallback everywhere.
      auto rank0_refused=refusal_case("multirate4_refused_rank0",{0},false,1000000U);
      same_schedule(rank0_refused.first.hierarchical_block_records,all_refused.first.hierarchical_block_records);
      same_endpoint(rank0_refused.second,all_refused.second,1e-22);
      auto several_refused=refusal_case("multirate4_refused_several",{0,ranks-1},false,1000000U);
      same_schedule(several_refused.first.hierarchical_block_records,all_refused.first.hierarchical_block_records);
      auto allocation_failed=refusal_case("multirate4_allocation_failure",{ranks-1},true,1000000U);
      same_schedule(allocation_failed.first.hierarchical_block_records,all_refused.first.hierarchical_block_records);
      same_endpoint(allocation_failed.second,all_refused.second,1e-22);
      // Only the first admission refused: that cycle falls back, later cycles
      // are multi-block again; the fallback is visible in the block records.
      auto first_refused=refusal_case("multirate4_refused_first",{0},false,1U);
      require(first_refused.first.hierarchical_block_records.front().pm_cycle_length_limiter=="rollback_memory");
      require(first_refused.first.multirate_memory_limited_cycles==1U);
      // Reproducibility mode: the same refusal fails collectively instead.
      {
        auto strict_text=multirate_config_with(4,hierarchical_config.normalized_text,"multirate_m4_strict_source").normalized_text;
        const std::string from="treepm_multirate_rollback_memory_policy = single_block_fallback";
        const auto position=strict_text.find(from);
        require(position!=std::string::npos);
        strict_text.replace(position,from.size(),"treepm_multirate_rollback_memory_policy = fail");
        workflows::ReferenceWorkflowRunner strict(core::loadFrozenConfigFromString(strict_text,"multirate_m4_memory_fail"));
        bool failed=false;
        try {
          (void)execute(strict,"multirate4_memory_fail",5e-6,nullptr,0U,0U,[&](workflows::ReferenceWorkflowOptions& options) {
            options.multirate_test_snapshot_refusals=1U;
            options.multirate_test_snapshot_refusal_ranks={0};
          });
        } catch (const std::exception& error) {
          failed=std::string(error.what()).find("memory governor")!=std::string::npos;
        }
        require(mpi.allreduceSumUint64(failed?1U:0U)==static_cast<std::uint64_t>(ranks));
        // Without a refusal the strict policy is the uninterrupted schedule.
        auto strict_ok=execute(strict,"multirate4_memory_fail_admitted",5e-6);
        same_schedule(strict_ok.first.hierarchical_block_records,full4.records);
        same_endpoint(strict_ok.second,full4.last.second,1e-22);
      }
    }

    // (8) Consecutive rejections of the same cycle (forced, identical T):
    // bitwise equal to the uninterrupted run; the accepted-trajectory
    // diagnostics (cadence history, refresh counts, block records) equal the
    // uninterrupted run's although the attempts generated more records; the
    // governed snapshot storage is reused across retries and released.
    if (has_multi_block_cycle) {
      auto repeated=execute(multirate4,"multirate4_consecutive_rejections",5e-6,nullptr,0U,6U,
          [](workflows::ReferenceWorkflowOptions& options) { options.multirate_test_forced_rejections_per_cycle=3U; });
      check_blocks(repeated.first);
      require(repeated.first.multirate_forced_test_rollbacks>=3U);
      unsigned max_retries=0U;
      for (const auto& block:repeated.first.hierarchical_block_records) max_retries=std::max(max_retries,block.pm_cycle_retries);
      require(max_retries>=3U);
      same_schedule(repeated.first.hierarchical_block_records,full4.records);
      same_endpoint(repeated.second,full4.last.second,1e-22);
      const auto& accepted=repeated.first.treepm_cadence_records;
      const auto& uninterrupted=full4.last.first.treepm_cadence_records;
      require(accepted.size()==uninterrupted.size());
      for (std::size_t i=0;i<accepted.size();++i) {
        require(accepted[i].step_index==uninterrupted[i].step_index);
        require(accepted[i].refreshed_long_range_field==uninterrupted[i].refreshed_long_range_field);
        require(accepted[i].short_range_only==uninterrupted[i].short_range_only);
        require(accepted[i].field_version==uninterrupted[i].field_version);
      }
      require(repeated.first.multirate_rejected_cadence_records>0U);
      require(repeated.first.treepm_long_range_refresh_count==full4.last.first.treepm_long_range_refresh_count);
      require(repeated.first.multirate_rollback_snapshot_peak_bytes>0U);
      require(repeated.first.multirate_rollback_snapshot_retained_bytes==0U);
      // A rank-local failure inside the rollback (last rank) must fail every
      // rank through the collective failure protocol, never deadlock.
      {
        bool failed=false;
        try {
          (void)execute(multirate4,"multirate4_rollback_failure",5e-6,nullptr,0U,1U,
              [&](workflows::ReferenceWorkflowOptions& options) { options.multirate_test_rollback_failure_rank=ranks-1; });
        } catch (const std::exception&) {
          failed=true;
        }
        require(mpi.allreduceSumUint64(failed?1U:0U)==static_cast<std::uint64_t>(ranks));
      }
    }

    // (9) Rollback with every source-history cache enabled: identical-source
    // tree reuse, certified topology refit, residual lookups, SIMD leaves and
    // the (remote, under MPI) group walk; the sparse LET graph is reused
    // automatically. Rollback forces a full tree build, so a replay with
    // forced rejections keeps the uninterrupted schedule and agrees to
    // roundoff (summation order may differ from the refitted original).
    if (has_multi_block_cycle) {
      auto stress_text=optimized_text;
      const std::string from="treepm_group_walk_enabled = false";
      const auto position=stress_text.find(from);
      require(position!=std::string::npos);
      stress_text.replace(position,from.size(),"treepm_group_walk_enabled = true");
      workflows::ReferenceWorkflowRunner stress(multirate_config_with(4,stress_text,"multirate_m4_cache_stress"));
      auto plain=execute(stress,"multirate4_cache_stress_plain",5e-6);
      auto rejected=execute(stress,"multirate4_cache_stress_rejected",5e-6,nullptr,0U,4U,
          [](workflows::ReferenceWorkflowOptions& options) { options.multirate_test_forced_rejections_per_cycle=2U; });
      check_blocks(plain.first);
      check_blocks(rejected.first);
      require(rejected.first.multirate_forced_test_rollbacks>0U);
      same_schedule(rejected.first.hierarchical_block_records,plain.first.hierarchical_block_records);
      same_endpoint(rejected.second,plain.second,1e-18);
      // Deterministic replay of the same rejection pattern.
      auto replay=execute(stress,"multirate4_cache_stress_replay",5e-6,nullptr,0U,4U,
          [](workflows::ReferenceWorkflowOptions& options) { options.multirate_test_forced_rejections_per_cycle=2U; });
      require(replay.first.final_state_digest==rejected.first.final_state_digest);
    }

    // (10) Highly skewed ownership (rank 0 owns nearly everything, the last
    // rank owns nothing) with forced rollbacks: same physical trajectory as
    // the balanced run to roundoff, identical schedule on every rank, no
    // collective divergence (owner timeout on the _np2/_np4 registrations).
    if (has_multi_block_cycle) {
      const auto skewed=initialState(rank,ranks,true);
      auto skewed_run=execute(multirate4,"multirate4_skewed_rollback",5e-6,nullptr,0U,2U,
          [&](workflows::ReferenceWorkflowOptions& options) { options.initial_state_override=&skewed; });
      check_blocks(skewed_run.first);
      require(skewed_run.first.multirate_forced_test_rollbacks>0U);
      require(skewed_run.first.hierarchical_block_records.size()==full4.records.size());
      for (std::size_t i=0;i<full4.records.size();++i) {
        const auto& a=skewed_run.first.hierarchical_block_records[i];
        const auto& b=full4.records[i];
        require(std::abs(a.time_end_code-b.time_end_code)<=1e-12*std::abs(b.time_end_code));
        require(a.opens_pm_cycle==b.opens_pm_cycle && a.closes_pm_cycle==b.closes_pm_cycle);
      }
      require(skewed_run.second.integrator_state.current_time_code==4e-5);
      require(mpi.allreduceSumDouble(phaseError(skewed_run.second.state,full4.last.second.state))<1e-18);
    }
    if (rank==0) {
      std::cout<<"hierarchical_multirate_outer_kick ranks="<<ranks
               <<" coarse_phase_l2="<<std::sqrt(multirate_coarse_error/64.0)
               <<" refined_phase_l2="<<std::sqrt(multirate_refined_error/64.0)
               <<" pm_refreshes="<<multirate_refreshes<<" exact_pm_refreshes="<<exact_refreshes
               <<" m4_blocks="<<full4.records.size()<<'\n';
      // Observed splitting order (diagnostic only; not gated): log2 of the
      // RMS error ratio per halving of dt.
      for (std::size_t level=1;level<ladder.size();++level) {
        if (splitting_error[level]>0.0) {
          std::cout<<"hierarchical_multirate_splitting_order level="<<level
                   <<" observed="<<0.5*std::log2(splitting_error[level-1U]/splitting_error[level])<<'\n';
        }
      }
    }
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
