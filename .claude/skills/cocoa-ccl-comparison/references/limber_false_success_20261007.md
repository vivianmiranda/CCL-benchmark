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
