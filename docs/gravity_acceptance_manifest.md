# Gravity acceptance manifest

Gravity acceptance is **source-bound evidence**, not a statement hard-coded by the
simulator. `tools/gravity_acceptance_manifest.py` fingerprints every file that can
materially change the gravity algorithm, source ownership, PM/MPI decomposition,
or force/integration semantics and records external test evidence against that
fingerprint.

A material edit under `src/gravity/**`, `include/cosmosim/gravity/**`, the gravity
workflow/integration ownership surfaces, or the parallel/decomposition contracts
invalidates the previous manifest. Historical manifests may be archived, but a
manifest whose fingerprint does not match the current tree must never be called
current certification.

The manifest records the source revision, compiler/build profile, MPI/FFTW/HDF5/
CUDA versions where available, the TreePM profile (grid, assignment,
deconvolution, decomposition, MAC/theta, ASMTH/RCUT and softening), rank counts,
and machine-readable test evidence. `create` only records supplied evidence; it
does not invent PASS results. `verify` returns a non-zero status for stale source.

Example production evidence collection after the exact MPI+FFTW build has run:

```sh
python3 tools/gravity_acceptance_manifest.py create \
  --output validation/artifacts/gravity_acceptance_current.json \
  --source-revision "$(git rev-parse HEAD)" \
  --build-profile mpi-fftw-release \
  --profile-id poster_dmo_treepm_v1 \
  --pm-backend fftw --pm-backend-capability production_fftw \
  --pm-grid 256x256x256 --assignment cic --deconvolution true \
  --pm-decomposition transposed_slab --theta-mac 0.7 \
  --asmth-cells 1.25 --rcut-cells 4.5 \
  --softening-profile dmo_fixed_comoving \
  --rank-count 1 --rank-count 2 --rank-count 3 --rank-count 4 \
  --evidence validation/artifacts/treepm_ewald.json \
  --evidence validation/artifacts/zeldovich_rank_equivalence.json
python3 tools/gravity_acceptance_manifest.py verify \
  validation/artifacts/gravity_acceptance_current.json
```

Generated acceptance results are run artifacts and are not part of a clean source
package unless intentionally curated as immutable validation reference input.
