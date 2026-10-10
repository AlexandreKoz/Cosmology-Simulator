// External diagnostic driver. Scientific production inputs remain .param.txt.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cosmosim/core/build_config.hpp"
#include "cosmosim/gravity/tree_pm_coupling.hpp"
#include "core_benchmark_process_memory.hpp"
#include "periodic_ewald_reference.hpp"

namespace {
struct Arguments {
  std::string particles, target_ids, forces;
  std::size_t count = 4096U, mesh = 64U, leaf = 16U, block = 64U;
  std::size_t targets = 4096U, repeats = 5U, warmups = 1U;
  double box_x = 1.0, box_y = 1.0, box_z = 1.0;
  double asmth = 1.25, rcut = 6.25, epsilon = 0.00035, g = 1.0, scale_factor = 1.0;
  bool adaptive = false, lookup = false, forensic = true, bootstrap = false;
  int ewald_level = 0;
};
std::size_t parseSizeValue(const std::string& value) {
  if (value.empty() || value.find_first_not_of("0123456789")!=std::string::npos) {
    throw std::invalid_argument("nonnegative integer required: "+value);
  }
  const auto number=std::stoull(value);
  if (number>std::numeric_limits<std::size_t>::max()) throw std::out_of_range("benchmark extent exceeds size_t");
  return static_cast<std::size_t>(number);
}
double parseRealValue(const std::string& value) {
  std::size_t consumed=0;
  const double number=std::stod(value,&consumed);
  if (consumed!=value.size()) throw std::invalid_argument("invalid real argument: "+value);
  return number;
}
Arguments parseArguments(int argc, char** argv) {
  Arguments a;
  for (int i=1; i<argc; ++i) {
    const std::string key(argv[i]);
    if (i+1 == argc) throw std::invalid_argument("option needs a value: "+key);
    const std::string value(argv[++i]);
    if (key == "--particles") a.particles=value;
    else if (key == "--target-ids") a.target_ids=value;
    else if (key == "--forces") a.forces=value;
    else if (key == "--count") a.count=parseSizeValue(value);
    else if (key == "--mesh") a.mesh=parseSizeValue(value);
    else if (key == "--leaf") a.leaf=parseSizeValue(value);
    else if (key == "--block") a.block=parseSizeValue(value);
    else if (key == "--targets") a.targets=parseSizeValue(value);
    else if (key == "--repeats") a.repeats=parseSizeValue(value);
    else if (key == "--warmups") a.warmups=parseSizeValue(value);
    else if (key == "--box-x") a.box_x=parseRealValue(value);
    else if (key == "--box-y") a.box_y=parseRealValue(value);
    else if (key == "--box-z") a.box_z=parseRealValue(value);
    else if (key == "--asmth") a.asmth=parseRealValue(value);
    else if (key == "--rcut") a.rcut=parseRealValue(value);
    else if (key == "--epsilon") a.epsilon=parseRealValue(value);
    else if (key == "--g") a.g=parseRealValue(value);
    else if (key == "--scale-factor") a.scale_factor=parseRealValue(value);
    else if (key == "--ewald-level") {
      const auto level=parseSizeValue(value);
      if (level>32U) throw std::invalid_argument("Ewald level must be <=32");
      a.ewald_level=static_cast<int>(level);
    }
    else if (key == "--policy" && (value=="strict" || value=="adaptive")) a.adaptive=value=="adaptive";
    else if (key == "--kernel" && (value=="analytic" || value=="lookup")) a.lookup=value=="lookup";
    else if (key == "--accounting" && (value=="full" || value=="fast")) a.forensic=value=="full";
    else if (key == "--history" && (value=="fallback" || value=="bootstrap")) a.bootstrap=value=="bootstrap";
    else throw std::invalid_argument("unknown option/value: "+key+" "+value);
  }
  for (const double x : {a.box_x,a.box_y,a.box_z,a.asmth,a.rcut,a.g,a.scale_factor}) {
    if (!std::isfinite(x) || x<=0.0) throw std::invalid_argument("lengths, split, G and scale must be positive finite");
  }
  if (!std::isfinite(a.epsilon) || a.epsilon<0.0 || a.mesh==0 || a.leaf==0 ||
      a.repeats==0 || a.targets==0 || a.count==0 || a.block==0 || a.block>4096U ||
      a.warmups>std::numeric_limits<std::size_t>::max()-a.repeats) {
    throw std::invalid_argument("invalid benchmark extent or softening");
  }
#if !COSMOSIM_ENABLE_FFTW
  if (a.mesh>16U) throw std::invalid_argument("mesh>16 requires FFTW; refusing a large diagnostic naive DFT");
#endif
  return a;
}
struct Sources {
  std::vector<std::uint64_t> id;
  std::vector<double> x,y,z,mass,epsilon;
};
Sources loadSources(const Arguments& a) {
  Sources s;
  if (a.particles.empty()) {
    for (std::size_t i=0; i<a.count; ++i) {
      s.id.push_back(i+1U);
      s.x.push_back(a.box_x*std::fmod((53.0*i+11.0)*9.31e-5,1.0));
      s.y.push_back(a.box_y*std::fmod((67.0*i+13.0)*7.73e-5,1.0));
      s.z.push_back(a.box_z*std::fmod((79.0*i+17.0)*6.29e-5,1.0));
      s.mass.push_back(1.0/static_cast<double>(a.count));
      s.epsilon.push_back(a.epsilon);
    }
  } else {
    std::ifstream input(a.particles);
    std::string line;
    if (!std::getline(input,line) || line!="id,x,y,z,mass,epsilon") {
      throw std::invalid_argument("particle CSV header must be id,x,y,z,mass,epsilon");
    }
    while (std::getline(input,line)) {
      std::replace(line.begin(),line.end(),',',' ');
      std::istringstream row(line);
      std::uint64_t id; double x,y,z,m,eps; std::string extra;
      if (!(row>>id>>x>>y>>z>>m>>eps) || (row>>extra) ||
          (!s.id.empty() && id<=s.id.back()) || !std::isfinite(x) || !std::isfinite(y) ||
          !std::isfinite(z) || !std::isfinite(m) || m<=0 || !std::isfinite(eps) || eps<0) {
        throw std::invalid_argument("particle rows must be finite, with increasing unique IDs, positive masses and nonnegative softenings");
      }
      s.id.push_back(id);s.x.push_back(x);s.y.push_back(y);s.z.push_back(z);
      s.mass.push_back(m);s.epsilon.push_back(eps);
    }
  }
  if (s.id.empty() || s.id.size()>=std::numeric_limits<std::uint32_t>::max()) {
    throw std::invalid_argument("particle population exceeds local tree index contract");
  }
  return s;
}
std::vector<std::uint32_t> selectTargets(const Sources& s,const Arguments& a) {
  std::vector<std::uint32_t> selected;
  if (a.target_ids.empty()) {
    const auto n=std::min(a.targets,s.id.size());
    for (std::size_t i=0;i<n;++i) selected.push_back(static_cast<std::uint32_t>((i*s.id.size())/n));
  } else {
    std::ifstream input(a.target_ids);std::uint64_t id;
    if (!input) throw std::invalid_argument("cannot open target-ID file");
    while (input>>id) {
      const auto it=std::lower_bound(s.id.begin(),s.id.end(),id);
      if (it==s.id.end() || *it!=id) throw std::invalid_argument("target ID absent from frozen sources");
      selected.push_back(static_cast<std::uint32_t>(it-s.id.begin()));
    }
    if (!input.eof()) throw std::invalid_argument("target-ID file requires integer IDs, one per line");
    auto sorted=selected;std::sort(sorted.begin(),sorted.end());
    if (std::adjacent_find(sorted.begin(),sorted.end())!=sorted.end()) throw std::invalid_argument("duplicate target ID");
  }
  if (selected.empty()) throw std::invalid_argument("empty numerical target set");
  return selected;
}
using Vector = std::array<double,3>;
Vector directResidual(const Sources& s,std::size_t target,const Arguments& a,double split,double cutoff) {
  std::array<long double,3> sum{};
  const long double pi=std::acos(-1.0L);
  for (std::size_t j=0;j<s.id.size();++j) {
    if (j==target) continue;
    // The double nearest-image/cutoff decisions match the solver; radial
    // algebra and accumulation are independent long-double reference work.
    double dx=s.x[j]-s.x[target],dy=s.y[j]-s.y[target],dz=s.z[j]-s.z[target];
    dx-=a.box_x*std::nearbyint(dx/a.box_x);dy-=a.box_y*std::nearbyint(dy/a.box_y);
    dz-=a.box_z*std::nearbyint(dz/a.box_z);
    const double r2=dx*dx+dy*dy+dz*dz;
    if (r2>cutoff*cutoff || r2==0.0) continue;
    const long double r=std::sqrt(static_cast<long double>(r2));
    const long double q=r/(2.0L*split), t=q*q;
    long double h;
    if (q<0.5L) {
      long double term=1.0L,series=0.0L;
      for (unsigned n=0;n<80U;++n) {
        const long double add=term/(2U*n+3U);series+=add;
        if (std::abs(add)<1e-30L) break;
        term*=-t/(n+1U);
      }
      h=4.0L/std::sqrt(pi)*series/(8.0L*split*split*split);
    } else h=(std::erf(q)-2.0L*q/std::sqrt(pi)*std::exp(-t))/(r*r*r);
    const long double eps=std::max(s.epsilon[j],s.epsilon[target]);
    const long double d=static_cast<long double>(r2)+eps*eps;
    const long double factor=a.g*s.mass[j]*(1.0L/(d*std::sqrt(d))-h);
    sum[0]+=factor*dx;sum[1]+=factor*dy;sum[2]+=factor*dz;
  }
  return {static_cast<double>(sum[0]),static_cast<double>(sum[1]),static_cast<double>(sum[2])};
}
void printError(const std::string& kind,const std::vector<Vector>& values,const std::vector<Vector>& reference) {
  long double err2=0,ref2=0;std::vector<double> norms,errors;
  for (std::size_t i=0;i<values.size();++i) {
    const double e=std::hypot(values[i][0]-reference[i][0],values[i][1]-reference[i][1],values[i][2]-reference[i][2]);
    const double r=std::hypot(reference[i][0],reference[i][1],reference[i][2]);
    err2+=static_cast<long double>(e)*e;ref2+=static_cast<long double>(r)*r;
    norms.push_back(r);errors.push_back(e);
  }
  // An explicit RMS floor makes small-force tails finite and interpretable.
  const double floor=1e-3*std::sqrt(static_cast<double>(ref2/reference.size()));
  std::vector<double> normalized;
  for (std::size_t i=0;i<norms.size();++i) normalized.push_back(errors[i]/std::max({norms[i],floor,1e-300}));
  double small_force_absolute_max=0.0;
  for (std::size_t i=0;i<norms.size();++i) if (norms[i]<=floor) {
    small_force_absolute_max=std::max(small_force_absolute_max,errors[i]);
  }
  auto absolute=errors;std::sort(absolute.begin(),absolute.end());
  std::sort(normalized.begin(),normalized.end());
  const auto quantile=[&](double p){return normalized[static_cast<std::size_t>(p*(normalized.size()-1U))];};
  std::cout<<"{\"record\":\"accuracy\",\"kind\":\""<<kind<<"\",\"targets\":"<<values.size()
    <<",\"relative_l2\":"<<std::sqrt(static_cast<double>(err2/std::max(ref2,1e-300L)))
    <<",\"absolute_max\":"<<*std::max_element(errors.begin(),errors.end())
    <<",\"absolute_p95\":"<<absolute[static_cast<std::size_t>(.95*(absolute.size()-1U))]
    <<",\"absolute_p99\":"<<absolute[static_cast<std::size_t>(.99*(absolute.size()-1U))]
    <<",\"small_force_absolute_max\":"<<small_force_absolute_max
    <<",\"normalization_floor\":"<<floor<<",\"p95\":"<<quantile(.95)<<",\"p99\":"<<quantile(.99)
    <<",\"maximum\":"<<normalized.back()<<"}\n";
}
int run(const Arguments& a) {
  const auto s=loadSources(a);const auto selected=selectTargets(s,a);
  const auto n=s.id.size();
  std::vector<std::uint32_t> active(n);std::iota(active.begin(),active.end(),0U);
  std::vector<double> ax(n),ay(n),az(n),tx,ty,tz,history;
  const bool accuracy=!a.forces.empty() || a.ewald_level>0;
  if (accuracy) {tx.resize(n);ty.resize(n);tz.resize(n);}
  const bool uniform=std::all_of(s.epsilon.begin(),s.epsilon.end(),[&](double e){return e==s.epsilon[0];});
  if (a.ewald_level>0 && (n>4096U || !uniform)) {
    throw std::invalid_argument("Ewald diagnostic is limited to <=4096 sources with uniform softening");
  }
  std::vector<std::uint8_t> overrides;
  if (!uniform) overrides.assign(n,1U);
  cosmosim::gravity::TreeSofteningView soften;
  if (!uniform) {soften.source_particle_epsilon_comoving=s.epsilon;soften.source_particle_epsilon_override_mask=overrides;}
  cosmosim::gravity::TreePmOptions options;
  options.pm_options.box_size_mpc_comoving=a.box_x;
  options.pm_options.box_size_x_mpc_comoving=a.box_x;options.pm_options.box_size_y_mpc_comoving=a.box_y;
  options.pm_options.box_size_z_mpc_comoving=a.box_z;options.pm_options.scale_factor=a.scale_factor;
  options.pm_options.gravitational_constant_code=a.g;
  options.pm_options.assignment_scheme=cosmosim::gravity::PmAssignmentScheme::kTsc;
  options.pm_options.enable_window_deconvolution=true;
  options.tree_options.gravitational_constant_code=a.g;options.tree_options.max_leaf_size=a.leaf;
  options.tree_options.multipole_order=cosmosim::gravity::TreeMultipoleOrder::kQuadrupole;
  options.tree_options.softening.epsilon_comoving=s.epsilon[0];
  options.tree_options.relative_force_tolerance=0.005;options.adaptive_maximum_opening_angle=0.25;
  options.split_policy=cosmosim::gravity::makeTreePmSplitPolicyFromMeshSpacing(a.asmth,a.rcut,
      std::cbrt((a.box_x/a.mesh)*(a.box_y/a.mesh)*(a.box_z/a.mesh)));
  options.gaussian_pair_lookup_enabled=a.lookup;options.full_mac_diagnostics=a.forensic;
  options.residual_block_size=a.block;
  cosmosim::gravity::TreePmForceAccumulatorView acc{active,ax,ay,az};
  if (accuracy) {acc.short_range_accel_x_comoving=tx;acc.short_range_accel_y_comoving=ty;acc.short_range_accel_z_comoving=tz;}
  cosmosim::gravity::TreePmCoordinator coordinator({a.mesh,a.mesh,a.mesh});
  if (a.adaptive && a.bootstrap) {
    // Controlled same-snapshot TOTAL history, explicitly not historical KDK.
    auto bootstrap_options=options;bootstrap_options.gaussian_pair_lookup_enabled=false;
    coordinator.solveActiveSet(s.x,s.y,s.z,s.mass,acc,bootstrap_options,nullptr,nullptr,soften);
    history.resize(n);
    for (std::size_t i=0;i<n;++i) history[i]=std::hypot(ax[i],ay[i],az[i]);
    acc.previous_acceleration_magnitude_code=history;
  }
  options.acceptance_policy=a.adaptive ? cosmosim::gravity::TreePmAcceptancePolicy::kAdaptiveRelative :
      cosmosim::gravity::TreePmAcceptancePolicy::kStrictReference;
  std::cout<<std::setprecision(17);
  std::cout<<"{\"record\":\"configuration\",\"particles\":"<<n<<",\"mesh\":"<<a.mesh<<",\"leaf\":"<<a.leaf
    <<",\"block\":"<<a.block<<",\"split\":"<<options.split_policy.split_scale_comoving
    <<",\"cutoff\":"<<options.split_policy.cutoff_radius_comoving<<",\"scale_factor\":"<<a.scale_factor
    <<",\"policy\":\""<<(a.adaptive?"adaptive":"strict")<<"\",\"kernel\":\""<<(a.lookup?"lookup":"analytic")
    <<"\",\"accounting\":\""<<(a.forensic?"full":"fast")<<"\",\"history\":\""<<(a.bootstrap?"same_snapshot_strict_bootstrap":"unavailable_fallback")
    <<"\",\"fftw\":"<<(COSMOSIM_ENABLE_FFTW?"true":"false")<<"}\n";
  for (std::size_t iteration=0;iteration<a.warmups+a.repeats;++iteration) {
    cosmosim::gravity::TreePmProfileEvent profile;cosmosim::gravity::TreePmDiagnostics d;
    const auto begin=std::chrono::steady_clock::now();
    coordinator.solveActiveSet(s.x,s.y,s.z,s.mass,acc,options,&profile,&d,soften);
    const double wall=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-begin).count();
    if (iteration<a.warmups) continue;
    const auto& p=profile.pm_profile;const auto memory=coordinator.memoryReport();
    const auto rss=cosmosim::bench::sampleProcessMemory();
    const auto& c=d.local_traversal;
    std::cout<<"{\"record\":\"timing\",\"iteration\":"<<iteration-a.warmups<<",\"wall_ms\":"<<wall
      <<",\"pm_ms\":"<<p.assign_ms+p.fft_forward_ms+p.poisson_ms+p.gradient_ms+p.fft_inverse_ms+p.interpolate_ms
      <<",\"tree_build_ms\":"<<profile.tree_profile.build_ms<<",\"multipole_ms\":"<<profile.tree_profile.multipole_ms
      <<",\"traversal_ms\":"<<profile.tree_profile.traversal_ms<<",\"pairs\":"<<c.direct_pair_evaluations
      <<",\"nodes\":"<<c.visited_nodes<<",\"opened\":"<<c.opened_nodes<<",\"multipoles\":"<<c.accepted_internal_multipoles
      <<",\"skipped_mac\":"<<c.skipped_mac_evaluations<<",\"workers\":"<<d.openmp_observed_workers
      <<",\"solver_retained_bytes\":"<<memory.totals.persistent_total_bytes+memory.totals.transient_total_bytes
      <<",\"rss_bytes\":"<<(rss.current_rss_bytes?std::to_string(*rss.current_rss_bytes):"null")
      <<",\"peak_rss_bytes\":"<<(rss.peak_rss_bytes?std::to_string(*rss.peak_rss_bytes):"null")<<"}\n";
  }
  if (accuracy) {
    std::vector<Vector> computed,total,reference;computed.reserve(selected.size());
    for (const auto i:selected) {
      computed.push_back({tx[i],ty[i],tz[i]});total.push_back({ax[i],ay[i],az[i]});
      reference.push_back(directResidual(s,i,a,options.split_policy.split_scale_comoving,options.split_policy.cutoff_radius_comoving));
    }
    printError("matching_split_direct_short_range",computed,reference);
    if (!a.forces.empty()) {
      std::ofstream output(a.forces);if (!output) throw std::runtime_error("cannot write force CSV");
      output<<std::setprecision(17)<<"id,total_x,total_y,total_z,tree_x,tree_y,tree_z,direct_x,direct_y,direct_z\n";
      for (std::size_t k=0;k<selected.size();++k) {
        output<<s.id[selected[k]];
        for (const auto& v:{total[k],computed[k],reference[k]}) for (const auto x:v) output<<','<<x;
        output<<'\n';
      }
      if (!output) throw std::runtime_error("force CSV write failed");
    }
    if (a.ewald_level>0) {
      std::vector<cosmosim::test_support::PeriodicEwaldSource> source;
      std::vector<cosmosim::test_support::PeriodicEwaldTarget> target;
      for (std::size_t i=0;i<n;++i) source.push_back({{s.x[i],s.y[i],s.z[i]},s.mass[i]});
      for (const auto i:selected) target.push_back({{s.x[i],s.y[i],s.z[i]},i});
      cosmosim::test_support::PeriodicEwaldOptions eo;
      eo.gravitational_constant=a.g;eo.alpha_inverse_length=2.0/std::min({a.box_x,a.box_y,a.box_z});
      eo.real_image_limits={a.ewald_level,a.ewald_level,a.ewald_level};
      eo.reciprocal_mode_limits={2*a.ewald_level,2*a.ewald_level,2*a.ewald_level};
      eo.plummer_softening_epsilon=s.epsilon[0];
      eo.softening_correction_image_limits={2*a.ewald_level,2*a.ewald_level,2*a.ewald_level};
      const auto ewald=cosmosim::test_support::periodicEwaldAccelerations(source,target,{a.box_x,a.box_y,a.box_z},eo);
      reference.clear();for (const auto& v:ewald) reference.push_back({v.x,v.y,v.z});
      printError("total_periodic_ewald_requires_image_convergence",total,reference);
    }
  }
  coordinator.shutdownMpiResources();
  return 0;
}
}  // namespace
int main(int argc,char** argv) {
  // Deliberately a serial-owner diagnostic even in an MPI-enabled build.
  // Distributed qualification uses the registered MPI workflow tests.
  try {return run(parseArguments(argc,argv));}
  catch (const std::exception& e) {std::cerr<<"bench_tree_pm_sweep: "<<e.what()<<'\n';return 1;}
}
