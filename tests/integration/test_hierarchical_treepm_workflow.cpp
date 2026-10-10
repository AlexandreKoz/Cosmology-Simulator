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
  require(time.last_completed_boundary_kind==core::StepBoundaryKind::kGlobalSynchronizationPoint);
  require(checkpoint.scheduler_state.current_tick%(1U<<rung)==0U);
  require(checkpoint.scheduler_state.current_tick==checkpoint.gas_cell_scheduler_state.current_tick);
  for (std::size_t i=0;i<checkpoint.state.particles.size();++i) {
    require(checkpoint.state.particleLastDriftTimeCode(i)==time.current_time_code);
    require(checkpoint.state.particleLastDriftScaleFactor(i)==time.current_scale_factor);
    require(checkpoint.scheduler_state.next_activation_tick[i]==checkpoint.scheduler_state.current_tick);
    require(checkpoint.scheduler_state.active_flag[i]==0U);
  }
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
  auto first=execute(hierarchical,"restart_first",2e-5,nullptr,1U);
  assertSynchronized(first.second,2U);
  auto resumed=execute(hierarchical,"restart_resumed",2e-5,&first.second);
  assertSynchronized(resumed.second,2U);
  require(mpi.allreduceSumDouble(phaseError(resumed.second.state,coarse.second.state))<1e-22);
  require(resumed.second.scheduler_state.current_tick==coarse.second.scheduler_state.current_tick);
  bool saw_fine=false,saw_pm_endpoint=false;
  for (const auto& record:coarse.first.treepm_cadence_records) {
    saw_fine=saw_fine||record.short_range_only;
    saw_pm_endpoint=saw_pm_endpoint||record.refreshed_long_range_field;
    if (record.short_range_only) require(!record.refreshed_long_range_field);
  }
  require(saw_fine && saw_pm_endpoint);
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
