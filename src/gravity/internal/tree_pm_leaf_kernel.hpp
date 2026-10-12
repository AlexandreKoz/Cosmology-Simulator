#pragma once

#include <cmath>
#include "cosmosim/gravity/tree_pm_coupling.hpp"

#if (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#define COSMOSIM_TREEPM_AVX2 1
#else
#define COSMOSIM_TREEPM_AVX2 0
#endif

namespace cosmosim::gravity::internal {

[[nodiscard]] inline bool treePmLeafAvx2Available() noexcept {
#if COSMOSIM_TREEPM_AVX2
  return __builtin_cpu_supports("avx2");
#else
  return false;
#endif
}

// Full blocks are padded with valid finite inputs. Only the caller's count
// lanes contribute. No approximate reciprocal/sqrt, no FMA/reassociated sums.
#if COSMOSIM_TREEPM_AVX2
__attribute__((target("avx2"))) inline void treePmLeafGeometryAvx2(TreePmLeafScratch& b, bool softened_needed) noexcept {
  const auto one = _mm256_set1_pd(1.0);
  for (std::size_t i = 0; i < k_tree_pm_leaf_lanes; i += 4U) {
    const auto x = _mm256_loadu_pd(b.dx.data()+i);
    const auto y = _mm256_loadu_pd(b.dy.data()+i);
    const auto z = _mm256_loadu_pd(b.dz.data()+i);
    const auto r2 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(x,x), _mm256_mul_pd(y,y)), _mm256_mul_pd(z,z));
    const auto raw_d = _mm256_add_pd(r2, _mm256_loadu_pd(b.epsilon2.data()+i));
    const auto zero_mask = _mm256_cmp_pd(raw_d, _mm256_setzero_pd(), _CMP_EQ_OQ);
    const auto d = _mm256_blendv_pd(raw_d, one, zero_mask);
    _mm256_storeu_pd(b.r2.data()+i, r2);
    if (softened_needed) {
      const auto inverse = _mm256_div_pd(one, _mm256_mul_pd(d, _mm256_sqrt_pd(d)));
      _mm256_storeu_pd(b.softened.data()+i, _mm256_blendv_pd(inverse, _mm256_setzero_pd(), zero_mask));
    }
  }
}
__attribute__((target("avx2"))) inline void treePmLeafContributionsAvx2(TreePmLeafScratch& b, double g) noexcept {
  const auto gravity = _mm256_set1_pd(g);
  for (std::size_t i = 0; i < k_tree_pm_leaf_lanes; i += 4U) {
    const auto factor = _mm256_mul_pd(_mm256_mul_pd(_mm256_loadu_pd(b.coefficient.data()+i), gravity),
        _mm256_loadu_pd(b.mass.data()+i));
    _mm256_storeu_pd(b.ax.data()+i, _mm256_mul_pd(factor, _mm256_loadu_pd(b.dx.data()+i)));
    _mm256_storeu_pd(b.ay.data()+i, _mm256_mul_pd(factor, _mm256_loadu_pd(b.dy.data()+i)));
    _mm256_storeu_pd(b.az.data()+i, _mm256_mul_pd(factor, _mm256_loadu_pd(b.dz.data()+i)));
  }
}
#endif

inline void treePmLeafGeometry(TreePmLeafScratch& b, bool avx2, bool softened_needed = true) noexcept {
#if COSMOSIM_TREEPM_AVX2
  if (avx2) { treePmLeafGeometryAvx2(b,softened_needed); return; }
#else
  (void)avx2;
#endif
  for (std::size_t i = 0; i < k_tree_pm_leaf_lanes; ++i) {
    b.r2[i] = b.dx[i]*b.dx[i] + b.dy[i]*b.dy[i] + b.dz[i]*b.dz[i];
    const double d = b.r2[i]+b.epsilon2[i];
    if (softened_needed) b.softened[i] = d == 0.0 ? 0.0 : 1.0/(d*std::sqrt(d));
  }
}
inline void treePmLeafContributions(TreePmLeafScratch& b, double g, bool avx2) noexcept {
#if COSMOSIM_TREEPM_AVX2
  if (avx2) { treePmLeafContributionsAvx2(b,g); return; }
#else
  (void)avx2;
#endif
  for (std::size_t i = 0; i < k_tree_pm_leaf_lanes; ++i) {
    const double factor = (b.coefficient[i]*g)*b.mass[i];
    b.ax[i] = factor*b.dx[i]; b.ay[i] = factor*b.dy[i]; b.az[i] = factor*b.dz[i];
  }
}

} // namespace cosmosim::gravity::internal
