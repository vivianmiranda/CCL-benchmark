# Limber QAG false-success investigation, 2026-10-07

Permanent archive: `diagnostics/limber_false_success/`. Original work
folder `work/limber_false_success_20261007/` remains untouched.

## Results and interpretation

- Unmodified CCL 3.3.3 / GSL 2.7, PR #1313 frozen tomography at
  8978d57e02e07e6557428ff50e134f3f21d70ab9. Bins 4 × 5 one-based,
  ell 355.655882, density only, massless neutrinos.
- Rounded w0 -1.0389 did not trigger locally. A scan and initial-rule
  root search found -1.0389272675840322. Native CCL QAG equals the direct
  native-GSL result to roundoff: 3.2299469642377736e-7, status success,
  41 evaluations. CQUAD gives 3.3653154695468126e-7 at epsrel 1e-8.
  Relative QAG error -4.0224610897%. CQUAD epsrel 1e-8 → 1e-9 changes
  by 4.48e-11 fractionally. Split QAG into two intervals removes the spike.
- Geometry/kernel changes at fixed P move the crossing; changing P at
  fixed geometry retains false acceptance in neighboring controls.
  Endpoint shifts ±0.001 in log(k) remove it with <5e-10 fractional
  change in CQUAD. These isolate quadrature-node alignment, not unusual
  dark-energy physics. Do not claim it only happens at one cosmology:
  the 41-point local scan finds a narrow interval of false acceptance.
- Both embedded rules undersample the peak and make almost the same
  error. They are quadrature rules; avoid saying GSL explicitly fits
  interpolating polynomials as its implementation. A successful return
  is an error-estimate outcome, not a guarantee that the integral is right.

## CoCoA scope: keep the distinction explicit

- All-pairs covariance-module angular spectra were tested at three nearby
  cosmologies, 26 multipoles and 15 galaxy pairs, with RSD/IA off. The
  default uses 96 Gauss–Legendre nodes on seven prescribed a panels;
  levels use 96, 128, 256, 512, 1024. Default/reference discrepancy at
  the affected point is -0.0016053%; maximum over all tested galaxy
  spectra is 0.0069491%; 512 → 1024 changes at most 9.37e-8 fractionally.
- Ordinary data-vector C_gg_tomo_limber supports autos only in this
  binding and includes RSD. Default 128 vs 1024 nodes changes at most
  0.0196928%. Levels 3 and 4 both select 1024 DV nodes, so equality of
  those two is not an extra convergence step.
- These are within-code quadrature tests with fixed power interpolation.
  Do not interpret the absolute cross-code difference as quadrature error.
  No covariance matrix, derivatives or Fisher contour was calculated.
- Blue ticks in the figures are covariance-module spectrum nodes, not
  ordinary data-vector nodes. GL nodes crowd each panel's endpoints;
  here the existing z=1 boundary lies near the z=0.9944 integrand peak.
  There is no adaptive peak detection. Known support can inform a
  minimum grid, but fixed sampling still needs convergence checks.

## Packaging and reproduction

The original result JSON/NPZ files and figures are copied without changes.
`results/study_summary.json` preserves original file hashes; the separate
archive_manifest.json fingerprints the permanent package and records
the original scripts before their path adaptations. The three 7.8 MB
shared-input power bundles are regenerable; their original hashes are
retained. `plot_geometry.npz` copies only z_1D/chi from shared_input_1 for
figure-only regeneration.

The portable scripts use the pinned frozen input, the unmodified upstream
helper (with CCL license), native GSL calls and actual CoCoA bindings.
Numerical results write outside the archive. Fresh numerical reproduction
was not run during packaging because the parent owns the sequential
TJPCov numerical workload. No production source is modified.

Verification during packaging:

- Saved hashes and every reported CoCoA convergence row recomputed from
  the saved spectra; separate data-vector and covariance fields retained.
- All 672 actual CoCoA node positions agree with the configured fixed
  rule to 3e-16 in scale factor. Extra geometry export preserves both
  spectrum arrays bitwise.
- Three saved PNG figures inspected. Portable figure scripts regenerate
  all three byte-for-byte from archived arrays; no cosmology code is called.
- Portable Python scripts parsed, README file links/anchors checked,
  and one-command-per-Step reproduction instructions retained.

The default historical CCL lifecycle remains PR #1296. This diagnostic
uses the already installed TJPCov environment's released CCL 3.3.3;
do not silently replace or rebrand one environment as the other.


## Completed paired failure survey and native timings

The follow-up is archived in `diagnostics/limber_failure_survey/`.
`results/20261007/publication_manifest.json` verifies 428 saved result files
and 70 figure/metadata files. There are 71 paired cases: 60 uniform
input/cosmology cases, the original fluid control, seven new targeted
cosmologies and three targeted narrow-input cases. All 71 primary CCL and
CoCoA reference checks pass; individual failed split-QAG controls remain
unsuccessful in their records. Large regenerated power/background bundles
stay local, with hashes and source metadata preserved in the archive.

- Uniform design: 30 PCG64 draws, seed 20261007, w0 [-1.15,-0.90],
  wa [-0.2,0.2], Omega_m [0.27,0.33]. Original inputs and a 15% shared
  Gaussian at z=1, sigma_z=0.01 use the same cosmologies. The primary
  bins-4-by-5, ell=355.655882 test finds 0/30 errors above 0.01% in both
  codes for each family. Each Wilson 95% interval is 0–11.35%; do not
  pool paired families as 60 independent cosmologies or infer MCMC rates.
  Screened secondary pairs reach 0.02149% QAG error and are not incidence
  samples. New cosmologies use CAMB PPF; the original control uses fluid.
- Nine node-guided cosmology planes yield seven accepted additional
  false-success points with QAG errors from -3.96626% to -4.10755%.
  CoCoA's own 96/1024 quadrature changes stay below 0.001614% there.
  Two searches find no accepted crossing; a discontinuous-bracket search
  attempt is preserved separately, not mislabeled as a cosmology failure.
- Targeted 15% narrow populations retain the original broad supports.
  Widths 0.03, 0.01 and 0.003 yield QAG errors -30.703%, -58.518% and
  -78.751%, each with GSL_SUCCESS after 41 evaluations. References check
  full-domain CQUAD, physical subdivisions at the peak and its wings,
  and doubled input sampling. Tight split-QAG EROUND returns are retained;
  they are never relabeled as successful independent references.
- CoCoA also underresolves the width-0.003 input. The actual public
  integration_accuracy levels 0/1/2/3/4 use 96/128/256/512/1024 nodes
  per panel, with accuracy_boost=1 and fixed inputs/interpolation/panels.
  The primary errors relative to level 4 are -0.520360%, -0.015490%,
  -0.001516%, -0.000241%, and zero. Maximum changes over all 15 pairs
  are 0.705191%, 0.018437%, 0.002361%, 0.000326%, and zero.
  Three regression checks assert expected default rejection at rtol=1e-4,
  level-2 acceptance at that tolerance, and level-3/4 agreement at 1e-5,
  all with atol=0. They pass without changing scientific tolerances.
  This is a quadrature check, not accuracy_boost=1 versus 5, and not an
  ordinary data-vector check. The refined common-reference offset is
  -0.01762%; do not claim quadrature alone removes interpolation differences.

Thirty-three PNG/PDF figures show both codes: diagnostic-coordinate
triangles, distinct redshift-overlap plots/atlas pages, original failure,
new narrow-peak sampling and the integration_accuracy convergence ladder.
These are not posterior triangles. Show CoCoA's default failure, not only
its refined result. Figure metadata and source hashes are archived.

### Timing scope and values

`timing_results/20261007/summary.json` contains nine interleaved batches,
separate process controls and fingerprints. Actual CCL 3.3.3 native C
integrands are timed through public angular_cl. A process-local macOS
GSL shim selects native CQUAD at epsrel=1e-4 and reuses per-thread CQUAD
workspaces; CCL's outer QAG allocation remains. No Python integrand
callback, production source/binary modification or PR #1313 split-check
overhead is included. CAMB/setup/first calls/I/O are excluded. Parent
ran each job alone with BLAS one and verified the native 1/6 OpenMP
thread setting after setup; normal desktop background remained active.

- 15 pairs by 26 multipoles (30–3000, not the exact failing ell): native
  CQUAD/QAG ratios are 1.091 and 1.063 with one thread, 1.067 and 1.048
  with six, for the first uniform draw and known-failure cosmology.
- Isolated ell=355.655882, one thread: ratios 1.001 (ordinary), 2.933
  (original false success), 5.319 (width-0.003 false success). QAG's
  premature return makes the incorrect answer deceptively inexpensive.
- CQUAD point errors vs accepted accuracy references are 0.000087%,
  0.000221% and 0.000939%. Additional grid multipoles use stricter
  internal CQUAD controls; they were not all independently refined.
- Instrumented QAG matches unmodified QAG bitwise. One- and six-thread
  spectra are also bitwise equal. Separate unmodified timing controls
  differ by up to 4.8%, including process variation; paired ratios are
  the primary timing comparison. Original failed loader/workspace probes
  remain local and are excluded. No full-covariance or MCMC claim follows.

The primary production CoCoA/core/project source was not changed by this
study. The historical PR #1296 benchmark environment remains distinct.


Final publication review verified every archive hash, all 71 primary reference
flags, targeted coordinates/values, physical-support confirmations, timing
means/scatter and overlap-atlas ordering. It corrected the selected-secondary
maximum to 0.0214932% (u020 shared-peak bins 1 by 2); 0.0201407% was only
the next-largest case. README rendering and local links/anchors pass.
The independent review used the available review agent, not the unavailable
historical Fable model. Do not attribute this follow-up review to Fable.
