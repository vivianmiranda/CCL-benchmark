# Paired CCL / CoCoA quadrature stress survey

The completed single-case diagnostic remains in `../limber_false_success/`.
This document records the sampling and acceptance protocol; numerical results
belong to the saved manifests and summaries, not to this preparation document.

## Predeclared sampling and interpretation

- Draw 30 independent points with PCG64 seed 20261007, uniformly in
  w0 = [-1.15, -0.90], wa = [-0.2, 0.2], Omega_m = [0.27, 0.33].
  Keep Omega_b, h, n_s and sigma8 at the frozen reproducer values.
- Evaluate the same 30 points for two fixed input families: the original
  five-bin PR #1313 distributions, and a positive mixture of 85% of each
  original bin plus 15% of a common Gaussian population at z=1, sigma=0.01.
  Each mixture component is normalized on the original redshift interval.
  Broad support remains; the narrow family is an explicit synthetic stress
  input, not an observational claim. Sample it at 5001 uniform redshifts;
  repeat the primary CQUAD reference with 10001 input samples.
- The primary incidence statistic uses bins 4×5 (one-based) at
  ell=355.655882. Per family, retain the predeclared denominator of 30
  cosmologies, show exclusions separately, and report the number exceeding
  relative errors 0.01%, 0.1%, and 1% against a validated within-code
  reference. These are conditional frequencies under this artificial
  sampling measure, not the prevalence of failures in real analyses.
- Native QAG/spline screens cover all 15 galaxy pairs. Pairs found by that
  screen receive reference checks, but their selected subset is not a
  second unbiased incidence sample. The primary pair is always checked,
  including when QAG and spline agree.
- The separate node-guided search seeks initial-rule crossings on a small
  Omega_m/wa grid. Those deliberate failure searches and the known failing
  control are excluded from incidence denominators. A 4.02% integration
  error is never described as a 4.02% failure frequency.

CCL's CAMB fluid default rejects phantom-crossing wa cases. All new survey
points explicitly use PPF. The known-control case retains the original
fluid model (wa=0), and is outside the frequency calculation. Backend or
setup failures are recorded separately and never counted as quadrature
false successes. CoCoA receives the same generated power/background inputs,
including the actual Omega_m, bias factors, and redshift distributions.

## Numerical checks

CCL uses its native public QAG and spline methods. Direct GSL replay checks
the native QAG value/status, then evaluates CQUAD and independently split
2/4-interval QAG at epsrel=1e-7 using the same native CCL integrand. The
initial reference is accepted only if those controls agree within 1e-5
fractionally and return finite positive values with successful statuses.

The pilot returned GSL_EROUND from some tightly requested split-QAG controls,
even when their numerical values agreed. These are preserved as unsuccessful
controls. The supplementary acceptance protocol evaluates CQUAD at epsrel=1e-8
on the whole interval and on two independent halves, plus eight-way split QAG
at epsrel=1e-7. Both refined CQUAD evaluations must succeed and agree within
the unchanged 1e-5 comparison tolerance. Every available successful split-QAG
control must also agree. If all split-QAG controls return EROUND, a result can
therefore have a **two-CQUAD reference basis**, explicitly distinct from a
reference supported by successful QAG controls. Failed statuses are retained
and are never reported as passing. Original pilot files remain unchanged;
their supplements are separately saved with the original file hash.

Keep every warning, exception, status and exclusion. Input-grid refinement
must agree within 1e-5 before the stress-family point enters a denominator.
A QAG false success requires both the public call to return without a failure
warning and the replay to return GSL_SUCCESS. Explicit reported failures are
counted separately. Wilson 95% binomial intervals are given only after every
planned primary reference is valid; targeted folders have no rate denominator.

For the deliberately shifted narrow-population targets, the reference also
must resolve the known physical support. `validate_narrow_support.py` maps
the population's mean and mean ± 6 sigma to native log-wavenumber boundaries
and evaluates native CQUAD separately on those intervals, retaining the
entire original domain. A second partition additionally includes mean ± 1
and ± 3 sigma; the fine redshift-input grid is checked on that refined
partition. Every control must succeed and agree within the same 1e-5
fractional tolerance. Whole-domain CQUAD agreement alone is insufficient
for these targets. The separate support records gate their reported
reference validity; no failed or unresolved CoCoA check is hidden.

CoCoA uses its actual covariance-module all-pairs spectra at 96, 512 and
1024 nodes per fixed panel. All 15 galaxy pairs are retained. The reference
must agree between 512 and 1024 within 1e-5. Save any CoCoA failure too.
This is neither ordinary data-vector integration nor a covariance matrix.
Do not call either code universally immune from a finite survey.

The three-case pilot consists of the archived failing original case and
the first uniform cosmology in each input family. Its measured costs set
the size of subsequent bounded batches. Each process has a ten-minute
limit and runs sequentially with six OpenMP threads. Saved elapsed times
are resource-planning measurements, not published benchmarks.

## Commands for the parent runner

From the repository root, prepare the immutable design once:

```bash
python diagnostics/limber_failure_survey/scripts/prepare.py --output work/limber_failure_survey_20261007
```

In the CCL 3.3.3 runtime, run the first three cases:

```bash
python diagnostics/limber_failure_survey/scripts/run_batch.py --engine ccl --pilot --output work/limber_failure_survey_20261007
```

Then, in the activated CoCoA runtime with the LSST Y1 interface available:

```bash
python diagnostics/limber_failure_survey/scripts/run_batch.py --engine cocoa --pilot --output work/limber_failure_survey_20261007
```

Review both codes, their reference and source-grid checks, and the measured
pilot cost before replacing `--pilot` with `--uniform --resume` in the two
commands. Parent alone launches numerical jobs. No failed output is
overwritten; fresh runs use fresh output directories.

The bounded targeted search can then run on all nine predeclared planes,
with at most one candidate root per plane and 60 kernel evaluations per plane:

```bash
python diagnostics/limber_failure_survey/scripts/run_targeted.py --engine search --output work/limber_failure_survey_20261007
```

A sign-changing initial-rule difference is only a search bracket. A root
solver can approach a discontinuity instead of a zero. Individual failed
brackets are recorded and skipped; a candidate must also have a finite,
positive K41 estimate and both the GSL error estimate and |K41−G20| below
1e-4 times |K41|, with no detected local jump. Full native-spectrum and
reference checks remain mandatory. Every new sample records the native
integration bounds. Reaching the predeclared budget is explicitly labeled
`budget_exhausted`, not evidence that no failing cosmology exists.

An interrupted search may be preserved and supplied with `hunt.py --saved-search`.
The helper verifies its configuration, plane and input fingerprints before
reusing samples in a fresh output folder. Reused samples count against the
same total evaluation budget, and the original files remain unchanged.

Evaluate those candidates in the CCL runtime:

```bash
python diagnostics/limber_failure_survey/scripts/run_targeted.py --engine ccl --output work/limber_failure_survey_20261007
```

Then evaluate the same candidates in the activated CoCoA runtime:

```bash
python diagnostics/limber_failure_survey/scripts/run_targeted.py --engine cocoa --output work/limber_failure_survey_20261007
```

Collect their results without assigning an incidence estimate:

```bash
python diagnostics/limber_failure_survey/scripts/run_targeted.py --engine summary --output work/limber_failure_survey_20261007
```

## Figure design

The triangle has parameter-coordinate panels, not posterior contours or
KDE densities. Show every jointly evaluated case as a large point, using
paired CCL/CoCoA glyphs with a common CCL-CQUAD-reference fractional-offset
colour scale; give the
targeted-search arm a distinct marker and label. Diagonal panels show
tested parameter positions/error, never an inferred posterior. If a code's
reference is unresolved, use an explicit open/hatched glyph and retain it
in the count. Add per-family frequency tables with denominators, exclusions
and binomial uncertainty; do not pool the paired families as independent
draws. An all-zero sample does not establish zero failure probability.

Within-code quadrature convergence belongs in a separate, clearly labeled
panel. In the archived failing case, CoCoA default/refined differs by
0.001605%, but its default spectrum differs from CCL CQUAD by about 0.05677%.
Never put the former beside CCL's 4.0225% as if they share a denominator.
The shared-reference residual includes each code's interpolation and
kernel conventions, so it is not solely an integration-error statistic.
Frequency tables above use each code's own validated quadrature reference.

The redshift figure groups by input family and galaxy pair: changing only
cosmology leaves n_i(z)n_j(z) unchanged. Show those distributions/products
once, with a neighboring measured-error panel for the same cosmologies in
CCL and CoCoA. A distance-coordinate inset may show cosmology-dependent
mapping, clearly labeled. Every scientific figure must display actual
CoCoA convergence/error results alongside CCL, not merely blue node ticks.
