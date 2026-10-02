<!-- Fable 5 review C (referee_completeness) of the CoCoA vs DESC-CCL study, 2026-10-01/02.
     Saved so later sessions reuse it instead of re-running Fable. Its fixes were
     applied in CCL-benchmark commit 7d5dd04 and later; line numbers refer to the README
     and sources as they were then (cosmolike_core v5.03 + 4, PR #1296 647ad4a). -->

# Fable 5 review, part C: sufficiency, didactics, coverage (referee report)

Scope: CCL-benchmark/README.md from "# CoCoA vs DESC-CCL" down, the
cocoa_comparison scripts and tables, SKILL.md section 0, and the run outputs
in scratchpad/cmp. Builds on parts A (numbers/physics) and B (writing/
figures); the README under review already incorporates their fixes (verified:
the z_piv sentence in finding 2, the corrected caption maxima, the indexing
sentence, the Contents list, the Step blocks in Reproduce). Nothing from A or
B is re-reported here. Read-only; one new measurement was made from existing
npz files (the 13 left-out Roman points, gap 10).

## Verdict 1: is the study sufficient?

Sufficient for the claims it actually states, with one structural
reservation. The layered evidence — identical inputs by construction, Limber
$C_\ell$ agreeing to $\sim 3\times10^{-4}$ median, both-codes-in-Limber real
space at 0.16/0.031, per-code convergence rows, the separable diagnostic
collapsing 8.76 to 0.21 and 0.228 to 0.062, and the code-level mechanism
verified line-by-line in review A — establishes that (i) the two codes agree
below CoCoA's own test tolerance everywhere except DESC-CCL's FKEM
combination when the linear growth is $k$-dependent, and (ii) the mechanism
of that one difference is identified, measured at the $C_\ell$ level, and
reproduced by a diagnostic. What the study does not establish is the reading
most users will take away: that CoCoA's numbers are the ones to trust and
that 8.76 does or does not matter for an analysis. Both codes approximate
the exact non-Limber integral when growth is scale-dependent; the study
measures their difference and DESC-CCL's internal inconsistency, never
either code's distance from the exact answer (gap 1), it never translates
8.76 into a parameter statement (gap 2), the parameter axis that drives the
whole effect ($w$, $m_\nu$) is varied only across projects, confounded with
survey depth (gap 3), and the CCL reference's convergence is shown only
one-sidedly (gap 5). None of these threatens a stated number; they bound
what the study may be cited for.

## Verdict 2: is it didactic?

Largely yes, and unusually so for a README-hosted study: the FKEM equation
is written out, the growth scale-dependence and the $C_{gs}$ offset are
given as measured percentages with the responsible parameter identified, the
one-change ladder lets a reader rebuild 8.76 from parts, the failure table
quotes the exact error strings with the root cause, and every number in the
prose traces to results.md/json or a named diagnostic script, so trust is
checkable. The remaining holes are local and sentence-sized: the reader is
not told why FKEM's three-term decomposition exists (so the equation reads
as arbitrary), why non-cancellation produces an $\ell$-flat offset at every
$\ell < 150$ rather than an edge effect near the switch, how to weigh 8.76
against the 0.2 yardstick the Summary quotes (44 times a regression
tolerance is not an interpretation), why the diagnostic floor 0.21 nearly
equals the both-Limber floor 0.16 (the reconciliation that shows FKEM
accounts for essentially all of 8.76), what the benchmark-modeling row is
for, and why $\Omega_m$ and $n_s$ were the parameters varied. Each is one
or two sentences, proposed under the gaps below.

## Verdict 3: does it cover all the aspects such a comparison must cover?

It covers more layers than most published code comparisons: pinned inputs,
harmonic Limber $C_\ell$, both-Limber real space, each code's convergence
(one row each), per-point figures against $\sigma$, failure forensics, and
a reproduction path. Measured against the full checklist, the missing
aspects are: the non-Limber $C_\ell$ are never compared directly in
harmonic space (the FKEM layer is seen only through real space and through
CCL-internal diagnostics; a CoCoA non-Limber $C_{gg}$ binding exists and is
unused) (gap 4); the real-space transform is never isolated (the same
$C_\ell$ through both transforms), so the both-Limber row bundles two
layers (gap 6); the harmonic check stops at $\ell = 5000$ while the
transforms consume $C_\ell$ to 65000/100000 (gap 5); CCL's reference
settings have no two-sided convergence evidence (gap 5); the cosmology grid
does not move the parameters the mechanism depends on (gap 3); the model
subset (NLA, all systematics zero) is stated but its consequence for the
agreement claim is not, and no systematics-on case exists (gap 9); there is
no acceptance-threshold interpretation or parameter-bias translation
(gap 2) and no truth anchor for the non-Limber disagreement (gap 1); and
the environment is only partially pinned (gap 11). Five cosmologies are
enough for the claim the README makes (stability over the $\Omega_m$/$n_s$
range) and not enough for cosmology-independence, which the README —
correctly — never claims. The two left-out $\gamma_t$ pairs are handled
honestly; what is missing is only the one-line bound showing the omission
is negligible (gap 10, measured below).

## Ranked gaps

Format per entry: (a) why it matters; (b) smallest run or text that closes
it, with whether the existing scripts support it; (c) proposed README text
where it is only a writing gap.

### 1. No truth anchor: the study never argues which code's non-Limber C_l are closer to the exact integral

(a) The Summary juxtaposes 8.76 with "CoCoA's tests pass at 0.2" and
finding 2 shows DESC-CCL is internally inconsistent (the FKEM/Limber pair
does not cancel) while CoCoA is consistent — but self-consistency is not
correctness. Both codes approximate the exact non-Limber integral when the
growth is $k$-dependent: CoCoA replaces $P_{\rm lin}(k,z)$ by one separable
form anchored at the lens bin's mean redshift; DESC-CCL mixes two forms.
The study proves a difference and a mechanism, not which number is right,
yet every reader will leave believing CoCoA's is. A referee would require
either the argument or the disclaimer.
(b) Two steps. (i) Writing: the inference the study already supports but
never states — for well-separated lens-source pairs at $\ell$ near 150 the
Limber approximation is accurate, so an $\ell$-flat offset of the FKEM
combination against the Limber value there (what `diag_fkem_offset.py`
measures: the offset persists at $\ell = 140$, where ref and separable
Limber values coincide) is an error of the combination, not non-Limber
physics. This closes the question for the dominant, flat part of the offset
using an existing, already-run diagnostic. (ii) Optional run for the
genuinely non-Limber range ($\ell \lesssim 30$): a brute-force quadrature
of $C_{gs}$ and $C_{gg}$ at a few $\ell$ for one pair, with the
geometric-mean unequal-time spectrum
$\sqrt{P_{\rm lin}(k,z_1)P_{\rm lin}(k,z_2)}$ (the N5K-style reference),
placing both codes against one external number. No current script does
this; it is a new ~40-line diagnostic in the family of `diag_*.py`, allowed
under owner rule 2 (diagnose, do not fix), minutes of runtime, no
downloads.
(c) Proposed addition at the end of finding 2:

> Neither code evaluates the exact integral when the growth depends on
> $k$: CoCoA replaces $P_{\rm lin}(k,z)$ by one separable form, DESC-CCL by
> two. The study measures their difference, not either code's distance
> from the exact $C_\ell$. For well-separated pairs the offset persists at
> $\ell = 140$, where the Limber value is accurate and the two linear terms
> should cancel; there it measures the FKEM combination's error, not
> non-Limber physics.

### 2. 8.76 is never interpreted: no acceptance criterion for a cross-code comparison, no parameter impact

(a) The only yardstick the README offers is CoCoA's internal frozen-test
tolerance (0.2), which is a regression tolerance between runs of one code,
not an acceptance threshold between codes. The reader is left alone with
"is 44 times that bad?". The question an analyst actually asks — would
publishing with DESC-CCL's vector instead of CoCoA's move the contours? —
is answerable from data the study already has: the five cosmology exports
give finite-difference derivatives of the data vector with respect to
$\Omega_m$ and $n_s$, so the standard projection (bias
$= F^{-1}g$, $g_i = (\partial d/\partial\theta_i)^T C^{-1}\,\Delta d$,
$F$ the 2x2 Fisher matrix in that subspace) turns 8.76 into "$x\sigma$ on
$\Omega_m$, $y\sigma$ on $n_s$ (within this 2-parameter subspace)".
(b) A ~20-line script over the existing npz files (cocoa_*_omm_lo/hi,
ns_lo/hi, cov_*.npz, ccl_*_fid_ref.npz): no new CoCoA or CCL runs. Not
supported by any current script (`results.py` and `dchi2_files.py` compute
only quadratic forms); it is a natural sibling of `dchi2_files.py`. Also
one sentence of text regardless of the run. Do not propose a new numeric
acceptance threshold (the 0.2 is the owner's frozen yardstick and stays).
(c) Proposed, in Summary after the threshold parenthesis:

> 0.2 is CoCoA's regression tolerance between runs of one code, quoted as
> a scale, not as a cross-code acceptance criterion.

and, once the projection is run, in Results:

> Projected on the $(\Omega_m, n_s)$ subspace with the data vector's
> finite-difference derivatives, the fiducial LSST-Y1 difference biases
> $\Omega_m$ by $x\sigma$ and $n_s$ by $y\sigma$ (statistical errors of
> this subspace; the full-parameter bias needs the full Fisher matrix).

### 3. The mechanism's parameter axis (w, m_nu) is never varied within a project

(a) The five cosmologies vary $\Omega_m$ and $n_s$, to which the FKEM
offset is nearly insensitive (the five LSST rows sit within 25% of each
other), so the grid mostly measures one number five times. The attribution
of the LSST-vs-Roman contrast to $w$ (finding 2) is deductively sound
(review A) but rests on a cross-project comparison confounded with
different covariances, depths, binnings and $n(z)$. One within-project run
— LSST-Y1 at $w = -1$ — turns the deduction into a measurement: finding 2
predicts the fiducial 8.76 collapses toward the Roman-like level because
the growth scale-dependence drops from 0.5-1.2% to the 0.03-0.3%
neutrino floor.
(b) One export + one CCL ref run + `dchi2_files.py`. Nearly supported
today: $w$ is a sampled parameter in `EXAMPLE_EVALUATE2.yaml` (prior block
present), so `cocoa_export.py` needs exactly one entry in its `MODELS`
dict (`"w_m1": {"w": -1.0}`); `ccl_compute.py` reads $w$ from the export's
`pars` and needs nothing. An $m_\nu = 0$ companion (isolating the neutrino
part of the Roman residual) is not supported — `mnu` would need a params
patch in the export script. The $w$ run is the one to do.
(c) If the run is not done, one sentence in Method:

> The five cosmologies test stability over $\Omega_m$ and $n_s$; $w$ and
> $m_\nu$, which set the growth's scale dependence (finding 2), vary only
> between the two projects.

### 4. Non-Limber C_l are never compared directly in harmonic space

(a) The FKEM difference — the study's one real discrepancy — is measured
only after the real-space transform and the $\Delta\chi^2$ projection, or
within CCL (ref vs separable, `diag_fkem_offset.py`). CoCoA's own
non-Limber $C_\ell$ never appear next to DESC-CCL's at the same $\ell$, so
the layer where the codes actually disagree is the one layer with no
direct table. The harmonic table's Limber-only scope is easy to misread as
covering the comparison's full harmonic content.
(b) Partially supported today. CoCoA side: the lsst_y1 interface exposes a
non-Limber `C_gg_tomo` binding (interface.cpp line 581) — one line in
`cocoa_export.py` exports it at $\ell < 150$; there is no non-Limber
`C_gs_tomo` binding, so the $gs$ column would need a new interface binding
(a cosmolike edit, the owner's call, not this study's). CCL side: trivial —
`angular_cl` with `l_limber=150` already computes FKEM $C_\ell$; a 5-line
extension of the `harmonic` variant evaluates $C_{gg}$ and $C_{gs}$ at,
say, $\ell = 2$-149. Cheapest full closure without any cosmolike edit:
promote the already-measured CCL-internal numbers (`diag_fkem_offset.py`:
FKEM ref/separable $-1$ and Limber ref/separable $-1$ at
$\ell = 10$-140, per pair) from two prose percentages into a small table in
finding 2, plus the $C_{gg}$ row from the existing binding.
(c) Writing part (if only the diagnostic table is added): a four-row table
under finding 2, columns $\ell = 10, 20, 50, 100, 140$, rows "LSST-Y1
$C_{gs}$ FKEM ref/separable $-1$", same for Roman-Real, plus the two
Limber rows (all zeros), replacing the two percentage ranges in prose.

### 5. CCL's reference convergence is one-sided, and the harmonic check stops at l = 5000 while the transforms use C_l to 65000/100000

(a) "Default FKEM sampling vs reference = 0.34" shows the default is
unconverged; nothing shows the reference (fkem_Nchi 2000, 1500 log $\ell$,
N_ELL_CORR, ELL_MAX_CORR) is converged. The unexplained floors — 0.16
both-Limber, 0.21 diagnostic — are exactly the size that unconverged CCL
numerics would produce, so the attribution of the floors is open. Separately,
the harmonic table validates $C_\ell$ only to 5000, but the Legendre sums
consume $C_\ell$ to 65000 (LSST) and 100000 (Roman); $\xi_-$ at 2.5' weighs
precisely the unvalidated range. Shear's 0.0013 covers it empirically, but
the layer-validation story has a hole where the reader cannot see it.
(b) Three small pieces, all nearly supported. (i) `fkem_scan.py` already
scans fkem_Nchi up to 8000 and fkem_chi_min on one pair — run it on one
LSST pair and report one line. (ii) A full-vector two-sided check needs one
new entry in `ccl_compute.py`'s CHANGES dict (e.g. `"ccl_hi":
dict(fkem_nchi=4000, nell_log=3000)`) and one run per project, then
`dchi2_files.py` against ref. (iii) The high-$\ell$ harmonic extension is a
one-line change to `ell_h` in `cocoa_export.py` (geomspace to lmax instead
of 5000); the `harmonic` variant and `harmonic.py` follow automatically
(one more range tuple in RANGES).
(c) After the runs, one sentence in finding 3:

> The reference itself moves by $x$ against fkem_Nchi 4000 / 3000 log
> $\ell$, so the floors of the one-change table are (not) CCL sampling.

### 6. The real-space transform layer is never isolated (same C_l through both codes' transforms)

(a) The both-Limber row (0.158 / 0.031) bundles Limber $C_\ell$
differences with transform differences, and the transform is where the
newest, least-reviewed code sits (PR #1296's full-sky bin-averaged
Legendre). The study implicitly validates that PR; an isolated transform
test would make that validation explicit and attribute the both-Limber
floor.
(b) Feasible in one direction only without touching cosmolike: CoCoA's
Limber $C_\ell$ on a dense $\ell$ grid through `ccl.correlation` (the PR
build) vs CoCoA's own `fidlimber` real-space vector. Needs: a dense-grid
$C_\ell$ export (the `*_tomo_limber` bindings accept arbitrary $\ell$
arrays; the current export's 24 points are too few to drive a transform) —
a few lines in `cocoa_export.py` — plus a ~30-line script reusing
`ccl_compute.py`'s `_corr` call with externally supplied `C_ell` (the
`correlation(cosmo, ell=, C_ell=, ...)` signature already takes arrays, so
no CCL adaptation beyond the allowed ones). Not supported by any current
variant (all variants recompute $C_\ell$ inside CCL). The reverse direction
(CCL $C_\ell$ through CoCoA's transform) needs cosmolike edits — out of
scope here.
(c) Until run, one clause in finding 1: "the 0.16 both-Limber residual
bundles the Limber $C_\ell$ differences (table above) with the two
transforms; it is not attributed between them."

### 7. The one-change table mixes categories without labels, and the floors are never reconciled

(a) The table's rows are four different kinds of thing — CoCoA
convergence, CCL convergence, modeling conventions, CCL-internal behavior,
and one composite — and the reader must classify them alone. This is the
task's "can a reader tell CCL modeling choices from CCL numerics from
CoCoA" question, and today the answer is "only a careful one". The same
table contains the study's best unstated result: the separable diagnostic
(0.21) lands nearly on the both-Limber floor (0.16), i.e. the FKEM form
accounts for essentially all of 8.76, with the remainder shared by the
Limber-level numerics and the anchor difference ($z = 0$ in the diagnostic
table, per-lens-bin $z_{\rm piv}$ in CoCoA).
(b) Writing only; the numbers all exist in results.md.
(c) Add a class column to the one-change table (values: CoCoA numerics /
DESC-CCL numerics / convention / DESC-CCL non-Limber / composite), and
after the table:

> The separable diagnostic (0.21) lands near the both-Limber floor
> (0.16): the FKEM term accounts for essentially all of the fiducial
> 8.76. The remainder includes the Limber-level residual and the anchor
> difference — the diagnostic table is anchored at $z = 0$, CoCoA per
> lens bin.

### 8. Benchmark-modeling row: the P(k)-source choice is never isolated and the row's purpose is unstated

(a) 154/56 is a composite of five changes; the one-change ladder isolates
no-RSD (115), Limber $C_{gs}$ (9.2) and flat-sky/centers (12.7) but never
the Eisenstein-Hu + halofit $P(k)$ — the one modeling ingredient unique to
this row — so a reader cannot rank the benchmark scripts' choices, and
the row risks reading as an indictment of DESC-CCL rather than of the
scripts' settings. The row is fair in substance (the NLA-vs-TATT deviation
is disclosed), but its purpose — the NOTE block's timings came from those
scripts, and 154 warns that their outputs were never a physics cross-check
of CoCoA — is left for the reader to infer.
(b) One-line CHANGES addition in `ccl_compute.py` (`"eh":
dict(native=True)`), one run per project, `dchi2_files.py`. Text for the
purpose sentence either way.
(c) After the benchmark-modeling bullet list:

> The row bounds what the benchmark scripts' timing runs computed; it is a
> property of those scripts' settings, not of DESC-CCL, and their outputs
> are not a physics cross-check of CoCoA.

### 9. The model subset's consequence for the agreement claim is unstated

(a) Method states that photo-z shifts, shear calibration, magnification
and point masses are zero and IA is NLA, but not the consequence: the
agreement claim covers exactly that subset, and the untested ingredients
are the convention-sensitive ones — photo-z shifts stress the $n(z)$
convention row (the one place the codes read the same file differently),
magnification runs through the same Limber machinery with its own kernel,
TATT through FAST-PT on both sides. A reader citing "CoCoA and CCL agree"
for a systematics-on analysis would overreach.
(b) Writing now; runs are future work and not currently supported
(`ccl_compute.py` has no `mag_bias` or photo-z shift machinery; the
cheapest future run is a nonzero photo-z shift — drop the `_DZ_` zeroing
regex in `cocoa_export.py` line 57 for one model and apply the same shift
to the $n(z)$ passed to the tracers, ~5 lines).
(c) One sentence at the end of Method's Models bullet:

> The comparison therefore covers the NLA + linear-bias model with every
> systematic at zero; magnification, photo-z shifts, shear calibration,
> TATT and baryons are outside it.

### 10. The 13 left-out Roman points are not bounded in the text

(a) The Roman $\Delta\chi^2$ silently omits 13 kept points where CCL's
transform fails; the reader cannot tell how much the omission could hide.
Measured here from the existing npz (read-only): CoCoA's $\gamma_t$ on
those 13 points is $1.6\times10^{-15}$ to $2.8\times10^{-12}$ (lens behind
source: the signal is zero to roundoff against typical $\gamma_t$ values),
and the $\Delta\chi^2$ of CoCoA's values alone on those points is
$1.6\times10^{-11}$. The omission is negligible unless DESC-CCL's
FKEM-only $\xi$ there were large, and its $C_\ell$ are nonzero only below
$\ell = 150$.
(b) Writing (the number above), optionally backed by a ~10-line sibling of
`dchi2_files.py`; the existing npz files support it (this review ran it).
(c) After "13 points the mask keeps" in Results:

> On those points CoCoA's $\gamma_t$ is below $3\times10^{-12}$ (lens
> behind source) and contributes $\Delta\chi^2 < 10^{-10}$ on its own, so
> the omission does not move the table.

### 11. Environment pinning is partial

(a) The pyccl commit is pinned and numpy/FAST-PT appear anecdotally in
What-it-took, but CAMB's version, GSL/FFTW, python, the OS, and the thread
counts used for the exports (cosmolike $\chi^2$ reproducibility is
thread-sensitive in some configurations) are unstated. A code comparison's
numbers should be re-derivable to the quoted digits.
(b) Writing; the versions are read off the two conda environments (no new
runs). A six-row table in Reproduce: python, numpy, pyccl (commit, already
there), CAMB, GSL/FFTW, OMP threads.
(c) Table as in (b); no prose needed.

## Closing note

Gaps 1-3 shape what the study may be cited for and should land before the
study is advertised beyond the group; 4-6 are evidence gaps a journal
referee would request; 7-11 are sentence- or script-sized hardenings. The
cheapest high-value actions, in order: the two sentences of gap 1(i) and
gap 2(c) (pure writing), the gap 2 bias projection and the gap 3 $w = -1$
run (one MODELS line plus existing scripts), and the gap 4 table from the
already-run `diag_fkem_offset.py`.
