<!-- Fable 5 review A (physics_and_numbers) of the CoCoA vs DESC-CCL study, 2026-10-01/02.
     Saved so later sessions reuse it instead of re-running Fable. Its fixes were
     applied in CCL-benchmark commit 17ee00d and later; line numbers refer to the README
     and sources as they were then (cosmolike_core v5.03 + 4, PR #1296 647ad4a). -->

# Fable 5 review, part A: physics and numbers
CoCoA v5.02 vs DESC-CCL real-space 3x2pt study (CCL-benchmark/README.md).
Written incrementally; checks in the task's order.

## Check 1: Finding 2 (FKEM)

DESC-CCL side (scratchpad ccl_pr, `pyccl/nonlimber/_nonlimber_FKEM.py`): CONFIRMED.
- FKEM term: chi integrand multiplies each tracer kernel by the calculator's
  growth factor `growfac_arr = ccl.growth_factor(cosmo, a_arr)` (lines 98-99,
  129-135), and the k-sum uses `pk(k, 1.0, cosmo)` with `pk` = the linear
  spectrum (lines 399-402, 473-485, same in the helper at 251-265). So the
  FKEM term is C^FKEM[D^2(z) P_lin(k,0)] exactly as the README writes.
  (There is a `transfer_low/transfer_avg` factor at k_low = 1e-5, lines
  111-135 — a tracer-transfer correction, not a growth k-dependence
  correction; it does not change the README's statement.)
- Subtracted Limber term: `cl_limber_lin` via `angular_cl_vec_limber(...,
  psp_lin, ...)` (lines 418-428), psp_lin = the calculator's full P_lin(k,z)
  (line 386). Combination `cl_limber_nonlin - cl_limber_lin +
  cls_nonlimber_lin` at line 490. So DESC-CCL does mix the two forms:
  P_lin(k,z) in the subtracted Limber term, D^2(z) P_lin(k,0) in the FKEM
  term. The README's description of DESC-CCL is accurate.

CoCoA v5.02 side (`external_modules/code/cosmolike_core/cosmolike/cosmo2D.c`;
working tree is v5.02 + 4 commits whose only cosmo2D.c change is the
FPT k_min/k_max -> krange[] rename, so the cited lines are v5.02 content):
- Both linear terms use ONE identical separable spectrum, but it is NOT
  D^2(z) P_lin(k,0). It is (D(a)/D(a_piv))^2 P_lin(k, a_piv) anchored PER
  LENS BIN at a_piv = 1/(1 + zmean(lens bin)):
  - pivot setup: lines 4162-4176 (gs Limber work), 5423-5435 (gg Limber
    work), 9160-9176 (gg FFTLog assemble), 9816-9825 (gs FFTLog assemble).
    The comment at 9160-9165: "the separable linear spectrum is anchored per
    lens bin at a_piv = 1/(1 + zmean(bin)) ... COSMO2D_FKEM_PIVOT_Z0
    fallback reproduces the z = 0 anchor exactly."
  - `COSMO2D_FKEM_PIVOT_Z0` is referenced only inside cosmo2D.c and defined
    nowhere in cosmolike_core or the project interfaces, so the default (and
    the study's) build uses the zmean anchor. The commit that introduced the
    per-bin pivot (3c9ed69 "cosmo2D: per-lens-bin pivot for the FKEM
    separable spectrum") is contained in tag v5.00, hence in v5.02.
  - subtracted Limber linear term: gg KG[3] = gf^2 * invgf2w[zl] *
    P_lin(k, a_piv) (lines 5495-5506); gs KIA[9] = gf^2 * invgf2w[zl] *
    P_lin(k, a_piv) (lines 4246-4268). Comment at 4255-4258: "never
    p_lin(k,a): only an identical separable form on both sides lets the
    FFTLog/Limber pair cancel at high l."
  - FFTLog term: gg `p_lin(k1cH0, apiv[i])*invgf2piv[i]` (line 9261),
    combined as `Cl = tcl*dlnk*2/pi + CLnl - CLlin` (line 9274); gs
    `p_lin(k1cH0, apivL[zl])*invgf2L[zl]` (line 9885), combined at 9915.
- RSD gates (lines 59-64): include_RSD_GG = 1, include_RSD_GS = 0, as the
  conventions table says.

VERDICT on the README sentence "CoCoA v5.02 evaluates both linear terms with
D^2(z) P_lin(k,0), while DESC-CCL mixes the two forms" (Findings item 2; same
claim in the Summary, "CoCoA v5.02 uses D^2(z) P_lin(k,0) in both terms"):
- "both linear terms with ONE separable form, so the pair cancels at the
  Limber limit" — correct, and it is the point that matters.
- "with D^2(z) P_lin(k,0)" — INCORRECT for v5.02. The anchor is the lens
  bin's mean redshift: (D(z)/D(z_piv))^2 P_lin(k, z_piv), z_piv =
  zmean(lens bin). D^2(z) P_lin(k,0) is only the COSMO2D_FKEM_PIVOT_Z0
  fallback, which is not compiled in.
- Side effect of the fix: it explains two numbers the README leaves
  unexplained — the separable-table diagnostic floor (0.21 / 0.062: the
  diagnostic hands DESC-CCL a z = 0-anchored separable table, CoCoA anchors
  at zmean, and the two separable forms differ in k-shape when growth is
  scale-dependent) and the 1.98 (CoCoA) vs 2.06 (DESC-CCL separable)
  non-Limber effect on LSST-Y1 gamma_t.

Proposed replacement text (Summary bullet, lines 41-46): keep the DESC-CCL
sentence; replace "CoCoA `v5.02` uses $D^2(z)\,P_{\rm lin}(k,0)$ in both
terms." with:
  "CoCoA `v5.02` uses one separable spectrum,
  $(D(z)/D(z_{\rm piv}))^2\,P_{\rm lin}(k, z_{\rm piv})$ anchored at each
  lens bin's mean redshift, in both terms."
Proposed replacement (Finding 2, lines 181-182):
  "CoCoA `v5.02` evaluates both linear terms with the same separable
  spectrum, $(D(z)/D(z_{\rm piv}))^2\,P_{\rm lin}(k, z_{\rm piv})$ with
  $z_{\rm piv}$ the lens bin's mean redshift, so the pair cancels at the
  Limber limit; DESC-CCL mixes the two forms."

## Check 2: the "w accounts for the difference" attribution

- Verified the exports' cosmologies (npz meta `point`): the two projects are
  identical in As, ns, H0, omegab, omegam, mnu; only w differs (-0.9 vs -1).
- Re-ran `scripts/diag_growth.py` on the run folder:
  LSST-Y1 (w = -0.9): D(k,z)/D(k0,z) - 1 = +0.0052 to +0.0119 over
  k = 0.01-0.2, z = 0.5-2; Roman-Real (w = -1): +0.0003 to +0.0030.
  The README's "+0.5% to +1.2%" and "+0.03% to +0.3%" are exact to the
  stated precision.
- JUDGMENT: the attribution is justified as written, and deductively so:
  D(k,z)/D(k0,z) is a function of the cosmology alone, both projects run
  through the same exporter/CAMB path, and the cosmologies differ only in
  w. The sentence claims no mechanism, only the parameter responsible.
- Mechanism: no diagnostic run isolates CAMB's dark-energy perturbations
  (e.g. same cosmology with perturbations switched, or PPF vs fluid), so
  naming "DE perturbations change the growth at k0 = 5e-4/Mpc relative to
  k = 0.01-0.2/Mpc" would be an inference, not a measured fact. It is
  physically plausible (w = -1 has no DE perturbations; for w > -1 they
  exist only near/above the sound horizon, i.e. at k ~ k0 and not at
  0.01-0.2/Mpc; the residual +0.03-0.3% at w = -1 is the massive-neutrino
  scale dependence, with the expected positive sign). RECOMMENDATION: keep
  the current sentence, or add the mechanism explicitly marked as an
  inference, e.g. "consistent with dark-energy perturbations ($w \neq -1$)
  altering the growth near $k_0 = 5\times10^{-4}\,{\rm Mpc}^{-1}$; not
  isolated by a dedicated run." Do not state it as established.

(Also verified: the exporter's growth matches the README's Method bullet —
`_cosmolike_prototype_base.py` lines 437-440, G = sqrt(PKL.P(z, 5e-4)/
PKL.P(0, 5e-4)) (1+z), normalized so D(0) = 1 downstream.)

## Check 3: Finding 4 (transform failures) — CONFIRMED, one precision nit

- `ccl_f1d.c` (`ccl_f1d_t_new`): with high-end `ccl_f1d_extrap_logx_logy`,
  `y[n-1]*y[n-2] <= 0` sets CCL_ERROR_SPLINE (lines 98-103) and the
  function frees the spline and returns NULL (lines 107-114).
- `ccl_correlation.c`: the Legendre path builds the C_l spline with
  `ccl_f1d_extrap_logx_logy` at the high end (lines 528-530) and maps the
  NULL to CCL_ERROR_MEMORY "ran out of memory" (lines 531-537) — the
  message observed in results.md. The FFTLog helper (lines 70-79) reports
  "failed to create spline". So both observed errors are the same root
  cause, exactly as finding 4 says.
- Lens-behind-source zero Limber C_gs: verified from the run exports.
  CoCoA's Limber C_gs for LSST-Y1 l4-s0 is < 1.2e-17 and Roman-Real l7-s2
  < 5.7e-23 at all 24 ell (healthy pairs: ~1e-7/1e-8), i.e. zero to
  roundoff. (Roman pairs l6-s0, l7-s0, l7-s1 are identically zero because
  they are CoCoA's excluded pairs, consistent with the Layout bullet.)
- `ccl_cls.c` lines 108-150 also confirm the "What it took" table row: the
  Limber RSD path (der_bessel >= 1) evaluates
  `chi_lp = (l+1.5)/k = chi (l+3/2)/(l+1/2)` and
  `ccl_scale_factor_of_chi(chi_lp)` — 1.4 chi at l = 2.
- Precision nit (low priority): "which needs positive values at the end"
  (finding 4) — the code condition is that the last two C_l be nonzero with
  the same sign (`y[n-1]*y[n-2] <= 0` fails); here they are exactly zero.
  Proposed: "which needs the last two values nonzero and of one sign; here
  they are exactly zero."

## Check 4: conventions table

All five rows verified against code:
- NLA: pyccl `tracers.py` lines 971-978 with `use_A_ia=True` applies
  `-ia_bias * 5e-14 * RHO_CRITICAL * Omega_m / D`; `ccl_compute.py` lines
  141-144 scales `ia_bias` by `0.01389/(5e-14*RHO_CRITICAL)` and includes
  `A1 ((1+z)/1.62)^eta`, with D = the CoCoA table passed to the calculator.
  Net factor = cosmolike's `0.01389 * Omega_m / D`. Row correct.
- n(z): cosmolike `redshift_spline.c` lines 907-908, 1000-1001 — Z_LOW is
  the default convention, "the file z column holds left bin edges, so the
  tabulated value belongs at the cell center z + dz/2"; `ccl_compute.py`
  line 135 does `z = d[:,0] + 0.5*dz`. Row correct.
- RSD: `cosmo2D.c` lines 61-62 (`include_RSD_GS = 0`, `include_RSD_GG =
  1`); `ccl_compute.py` lines 145-149 (lens tracer has_rsd=True for w,
  separate gamma_t lens tracer has_rsd=False in ref). Row correct, and the
  sentence under the table about `has_rsd` applying to every correlation
  of a tracer is correct (0.358 / 0.076 effects match results.md/json).
- non-Limber: cosmolike `structs.c` line 42 `LMAX_NOLIMBER = 150`;
  `ccl_compute.py` `l_limber=150`, FKEM, `fkem_Nchi=2000` (ref). Row and
  Method bullet correct.
- lmax: `projects/lsst_y1/EXAMPLE_EVALUATE2.yaml` line 20 `lmax: 65000`,
  `projects/roman_real/EXAMPLE_EVALUATE2.yaml` line 16 `lmax: 100000`;
  `ccl_compute.py` line 66 (65000 / 100000) and line 82 (ELL_MAX_CORR).
  Matches "6.5e4 / 1e5, CoCoA's lmax, also used as ELL_MAX_CORR".
- covariance (SKILL row): CoCoA `generic_interface.cpp` `IP::set_inv_cov`,
  10-column case (~line 4219): C(j,k) = table(i,8) + table(i,9), i.e.
  1-indexed columns 9 + 10; `plots.py` line 83 `float(t[8]) + float(t[9])`.
  Consistent.
- Layout bullets verified from the run exports: LSST-Y1 26 bins 2.5'-900',
  mask 959 of 1560; Roman-Real 15 bins 2.5'-250', mask 1950 of 2115; the
  three Roman gamma_t exclusions (6,0), (7,0), (7,1) are the identically
  zero pairs in the export. Method table ranges verified: zg 0-6,
  kg 1e-5-200 /Mpc, chi to z = 49.9, growth to z = 49.
- Nit: "every integer to 400" — the ell grid is `arange(2, 400)` plus
  log-spaced from 400 (`ccl_compute.py` line 156), so "every integer from
  2 to 400" would be exact. Low priority.

## Check 5: numbers in Summary, Results, Findings

- Summary table and both Results tables: every entry matches results.md /
  results.json under correct rounding (checked all 40 reference-vs-CoCoA
  entries and all 13 one-change rows for both projects; e.g. 8.76 <-
  8.7567, 0.21 <- 0.2081, 0.062 <- 0.06244, 115 <- 115.402, 154 <-
  153.590, 56 <- 55.938). "w 6.26 to 0.16" <- 0.16490: correct.
- Finding 2 growth numbers: re-ran `diag_growth.py` — LSST-Y1 +0.52% to
  +1.19%, Roman-Real +0.03% to +0.30%; README's "+0.5% to +1.2%" and
  "+0.03% to +0.3%" exact.
- Finding 2 offset numbers: re-ran `diag_fkem_offset.py` (pyccl, 2
  threads) — LSST-Y1 FKEM ref/separable - 1 = -0.70% to -1.54%,
  Roman-Real -0.03% to -0.39%, Limber ratio 0 everywhere; README's
  "-0.7% to -1.5%", "-0.03% to -0.4%", "the Limber C_l do not change"
  all correct.
- Finding 2 non-Limber effects 1.98 / 11.2 / 2.06: match results.json
  (1.9787 / 11.2123 / 2.0626).
- Finding 3: 0.34 (in w(theta): w column 0.395, gamma_t 0.027) and
  0.0005 — match.
- Finding 5: Delta ln theta — ln(900/2.5)/26 = 0.226 -> 0.23,
  ln(250/2.5)/15 = 0.307 -> 0.31. Correct.
- Harmonic table: identical to harmonic.md (header "1000 <= l <= 5000" vs
  harmonic.md "< 5001" — equivalent).
- Figures prose, recomputed from the run outputs (dv, sigma_*.npy,
  cov mask):
  - LSST gamma_t "reaches -0.5 sigma near 100'": measured minimum -0.49
    (fiducial) / -0.53 (across the five cosmologies) at 105'. OK.
  - LSST w "-0.5 sigma to -1.3 sigma per point below 100'": measured
    per-point range below 100' is -0.35 to -1.35 (fiducial), -0.31 to
    -1.39 across the five cosmologies. BOTH BOUNDS ARE OFF: should read
    "-0.3 sigma to -1.4 sigma". "At most 0.84 sigma on the points the
    mask keeps": measured 0.840, correct.
  - LSST xi "at most 0.01 sigma": measured 0.008. OK.
  - Roman gamma_t 0.03 / w 0.09 / xi 0.04: measured 0.029 / 0.088 /
    0.039 on kept points. OK.
- Roman "13 kept points are left out": results.json left_out_kept = 13 for
  every DESC-CCL Roman row. Correct.

## Check 6: owner-rule wording

- No "bug", "fixed in X", or fault-dating language anywhere in the study
  section; finding 2 follows the prescribed "CoCoA v5.02 does X
  consistently, while DESC-CCL ..." pattern (modulo the check-1 anchor
  correction).
- Diagnostics labeled: "Diagnostic (not used in the figures)" in the
  Summary, "(diagnostic)" on the separable rows, "The diagnostics behind
  findings 2 and 4" in Reproduce. OK.
- No timing claims in the study section (the NOTE block above the heading
  carries the times and is out of scope). No adjectives flagged.

## Ranked list

1. Finding 2 / Summary — "CoCoA v5.02 ... D^2(z) P_lin(k,0) in both
   terms" is factually wrong for v5.02 (per-lens-bin zmean anchor;
   COSMO2D_FKEM_PIVOT_Z0 not compiled in). README lines 45-46 and
   181-182. Replacements in Check 1 above. Also makes the separable-
   diagnostic floor (0.21 / 0.062) and the 1.98-vs-2.06 pair explicable;
   an optional one-line addition after the diagnostic numbers: "The
   remaining 0.21 includes the anchor difference: the diagnostic table is
   anchored at z = 0, CoCoA per lens bin."
2. LSST-Y1 w(theta) figure caption (line 222): "-0.5 sigma to -1.3 sigma
   per point below 100'" -> "-0.3 sigma to -1.4 sigma per point below
   100'" (measured -0.31 to -1.39 across the five cosmologies).
3. Finding 4 (line 199): "which needs positive values at the end" ->
   "which needs the last two values nonzero and of one sign; here they
   are exactly zero" (ccl_f1d.c: y[n-1]*y[n-2] <= 0 fails).
4. Finding 2 attribution (line 174): "w accounts for the difference
   between the projects" — keep; it is deductively justified (identical
   cosmologies except w, verified from the exports). If a mechanism is
   wanted, add it as an inference only (wording in Check 2).
5. Method ell grid (line 89): "every integer to 400" -> "every integer
   from 2 to 400" (nit).


