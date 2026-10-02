---
name: cocoa-ccl-comparison
description: How to run, diagnose and write up a CoCoA vs DESC-CCL code comparison of real-space 3x2pt data vectors (LSST-Y1, Roman-Real) with the projects' covariances, ratio plots in the cosmolike notebook style, and the owner's rules for what counts as a fair comparison. Use for any CoCoA/CosmoLike vs CCL (or other code) comparison or benchmark write-up in this repository.
---

# CoCoA vs DESC-CCL code comparison

Scripts: `cocoa_comparison/scripts/`. Results, figures and the write-up
live in this repository (the main README carries the study; the main
Cocoa README only cites it).

## 0. Owner's rules (non-negotiable)

1. **Real CCL vs CoCoA.** The only adaptations allowed on the CCL side:
   - CoCoA's CAMB tables passed to `pyccl.CosmologyCalculator`
     (linear and nonlinear P(k), chi(z), H(z), growth D and f), so no
     Boltzmann or background difference contaminates the comparison;
   - pyccl built from LSSTDESC/CCL PR #1296 (full-sky xi+- and exact
     bin averaging; CCL 3.3 refuses full-sky xi+-, error 1040). Never
     code a transform yourself when CCL (or a CCL PR) has it.
   - CCL's documented accuracy settings (fkem_Nchi, l sampling,
     ELL_MAX_CORR = CoCoA's lmax), reported separately.
2. **Do not fix CCL.** When CCL fails or disagrees, find why (diagnostic
   runs and tricks are fine) and REPORT it; the plots and main tables stay
   real CCL vs CoCoA. A CCL transform that fails for a bin pair is
   recorded and the pair left out (NaN), never worked around.
3. Wording: never "bug fixed in X on date". Say "CoCoA v5.02 does this
   consistently, while CCL ...".
4. Never time anything (the owner times on the Intel test machine).
   Downloads (git clones, conda/pip packages) need explicit permission.
   <= 8 threads total. Never copy the cocoa environment dump (it holds a
   secret token) anywhere.
5. At the end: a Fable 5 review of the comparison (precision of every
   claim) and of the README (the cocoa-maintenance skill's README rules:
   less is more, tables over prose, LaTeX math, no adjectives).

## 0.5 Saved Fable 5 results (read before asking Fable again)

Fable 5 runs are expensive: every insight they produce is saved here and
reused. Read the relevant file before spawning a new review.

| file | what it verified (with file:line evidence) |
|---|---|
| `references/fable_review_A_physics_and_numbers.md` | CCL FKEM mixes P_lin(k,z) (Limber subtraction) with D^2 P_lin(k,0) (FKEM term); CoCoA's per-lens-bin pivot in cosmo2D.c; the transform-failure path (ccl_f1d log-log extrapolation); every convention row (NLA factor, Z_LOW n(z), RSD gates, lmax, covariance columns 9+10); all study numbers |
| `references/fable_review_B_writing_and_figures.md` | README rules (cocoa-maintenance Section 3) applied to the study; caption self-sufficiency; figure checks (alpha, indexing, excluded panels); Reproduce as Step blocks |
| `references/fable_review_C_referee_completeness.md` | what a code comparison must contain: truth anchor, parameter-shift interpretation, varying the parameter that drives a difference within one project, direct harmonic comparison of the disagreeing layer, two-sided convergence, transform isolation, scope statement, bounds on left-out points, version table |
| `references/tatt_conventions.md` | CoCoA (CFASTPT, IA_code 0) vs pyccl TATT convention map: C1/C2/C_delta normalizations and signs, FAST-PT kernel mapping, B-modes, the gamma_t split (C1 in FKEM, extras in Limber), PT settings, the recipe of the "tatt:" variant |
| `references/fable_review_D_tatt.md` | the TATT harness checked against that map; every TATT number recomputed; finding 3 per source pair |

## 1. Pipeline

1. CoCoA side, cocoa env (`start_cocoa.sh` sourced), from `Cocoa/`:
   `python cocoa_export.py <project> <model> <out.npz>` with
   `EXPORT_MODE=dv` (lsst_y1: a scratch copy of the data folder whose
   dataset uses `ones.mask`, full-length vector) or `EXPORT_MODE=cov`
   (the shipped mask; also the masked inverse covariance). Saves the full
   theory vector, block sizes, theta binning, nuisance values, CAMB
   P_lin/P_nl on (z 0-6, k 1e-5-200 /Mpc), chi(z) on the likelihood's
   z_interp_1D, growth D(z) to z = 49 (from P_lin at k = 5e-4 /Mpc), and
   CoCoA's Limber C_l (`C_ss/gs/gg_tomo_limber` bindings) at 24 ell.
   `COCOA_OVR='{json}'` overrides likelihood options (high accuracy,
   Limber-only runs). Photo-z shifts, shear calibration, point masses and
   magnification are set to zero.
2. CCL side, ccl conda env with the PR build first on PYTHONPATH:
   `ccl_compute.py <cocoa npz> <variant> <out.npz>` (variants in its
   docstring: ref, limgs, limber, norsd, rsd_gs, flat, points, separable
   (diagnostic), ccl_numerics, bench, harmonic). Needs `COCOA_ROOTDIR`.
   A failing CCL call (angular_cl or correlation) is recorded per pair in
   the output's `failures` and the pair is NaN.
3. `CMP_WORK=<run folder>`: `results.py` -> results.md/json (Delta chi2 =
   d^T C^-1 d, masked inverse covariance, per probe, failures table);
   `harmonic.py` (Limber C_l medians/maxima); `plots.py <outdir>` (cocoa
   env). `run_cocoa.sh` and `run_ccl.sh` run the whole campaign.

Models: fiducial, Omega_m 0.25/0.35, n_s 0.92/1.01.

## 2. Matching conventions (check each, state each in the README)

| item | CoCoA | CCL setting used |
|---|---|---|
| n(z) z column | read as Z_LOW (value at z + dz/2) | z + dz/2 passed to tracers |
| NLA amplitude | A1 ((1+z)/1.62)^eta Omega_m 0.01389/D | ia_bias scaled by 0.01389/(5e-14 RHO_CRITICAL) |
| galaxy bias | b1 per bin, constant | bias=(z, b1) |
| non-Limber switch | l < 150 (gg and gs) | l_limber = 150, FKEM |
| RSD | in gg only (include_RSD_GG = 1, include_RSD_GS = 0) | has_rsd=True in the w tracer, False in the gamma_t lens tracer (has_rsd applies to EVERY correlation of a tracer) |
| transform | full-sky, bin-averaged | PR #1296 legendre with theta_max |
| growth | D = sqrt(P_lin(k0,z)/P_lin(k0,0)), k0 = 5e-4/Mpc | the same table, to z = 49 |
| background | chi(z) on z_interp_1D to z = 50 | the full table; H from a cubic-spline derivative |
| covariance | C = cov_g + cov_ng (columns 9, 10) | sigma = sqrt(C_ii) for the figures |

## 3. Lessons (2026-10-01 study)

- Validate in layers: harmonic C_l first (agreed to <= 0.02% for
  l >= 66), then both codes fully Limber (LSST-Y1 Delta chi2 0.23), then
  non-Limber. Convergence of each code separately (CoCoA default vs
  high accuracy: 0.004; CCL default FKEM sampling vs converged: 0.34).
- CCL FKEM with CAMB's scale-dependent P_lin (massive neutrinos): the
  FFTLog term assumes separable growth D^2 P_lin(k,0) while the
  subtracted Limber term uses P_lin(k,z); a flat ~-0.85% offset in C_gs
  below l = 150 even for well-separated lens-source pairs. Diagnostic:
  pass P_lin(k,0) D(z)^2 -> offset gone (LSST-Y1 3x2pt 8.7 -> 0.68).
  CoCoA v5.02 uses one separable form on both terms, anchored per lens
  bin: (D(z)/D(z_piv))^2 P_lin(k, z_piv), z_piv = the bin's mean redshift
  (cosmo2D.c), not z = 0 (a first draft said z = 0; Fable caught it).
- CCL Legendre and FFTLog transforms fail ("ran out of memory" /
  "failed to create spline") on C_l with non-positive tails
  (lens-behind-source gamma_t pairs); real failure, pair left out.
- CCL's Limber RSD term is ~0 (same as CoCoA's Limber C_gg); RSD enters
  through the non-Limber range, so ignoring RSD costs Delta chi2 ~115 in
  LSST-Y1 w(theta).
- Input pitfalls that were MY bugs, not CCL's: growth table must cover
  the RSD kernel's range (to z ~ 49), H(z) from a derivative of chi(z) on
  a grid whose spacing jumps at z = 3; EH needs sigma8 (use ccl.sigma8 of
  the CAMB-table cosmology); relative data paths.
- CoCoA side: `init_data_real` cannot be called twice in one process;
  Roman-Real's covariance is not positive definite with every point kept
  (use the shipped mask; masked CoCoA entries are zero -> blank in plots).
- Local env: FAST-PT 4.0.0 calls np.trapz (removed in numpy 2.4): alias
  `np.trapz = np.trapezoid` before importing pyccl; building the PR needs
  swig (conda-forge, installed in the ccl env with the owner's OK).
- RSD in gamma_t: CCL's `has_rsd` puts RSD in gamma_t too; CoCoA does not
  (include_RSD_GS = 0). Use a separate lens tracer without RSD for
  gamma_t; report the size (LSST-Y1 0.36, Roman-Real 0.076).
- CCL's Limber RSD kernel evaluates the background at
  chi_{l+1} = chi (l + 3/2)/(l + 1/2): 1.4 chi at l = 2 (z ~ 14 for lenses
  at z = 4). Pass chi(z) to z = 50 or l = 2 fails ("integration error").
- The FKEM offset depends on the scale dependence of the linear growth:
  with w = -0.9 (LSST-Y1 fiducial) CAMB's D(k,z)/D(k0,z) - 1 is 0.5-1.2%
  (k = 0.01-0.2, z = 0.5-2), with w = -1 (Roman-Real) 0.03-0.3%. Offset in
  C_gs below l = 150: -0.7 to -1.5% (LSST-Y1), -0.03 to -0.4% (Roman).
- Failing gamma_t pairs: lens behind source, Limber C_gs exactly 0; FKEM
  gives C_l != 0 only below l = 150; CCL's correlation builds the C_l
  spline with log-log extrapolation beyond lmax (needs positive end
  values) and reports "ran out of memory" / "failed to create spline".
- Referee checks that closed the study (Fable 5, part C): swap w between
  the projects (LSST-Y1 at w = -1: 0.34; Roman-Real at w = -0.9: 7.57) to
  turn the w attribution into a measurement; push CoCoA's Limber C_l
  through CCL's transform (variant cocoa_cl: 0.009 / 0.0006) to isolate
  the transform layer; a finer CCL run (ccl_hi: 0.0014 / 0.0002) for
  two-sided convergence; the harmonic check to lmax and non-Limber C_gg
  directly (CCL 0.8-1.7% below CoCoA in LSST-Y1, flat in l); parameter
  shifts with a Fisher matrix marginalized over the lens biases and NLA
  (bias.py: LSST-Y1 biases absorb 6.74 of 8.76, ~1 sigma each; Omega_m,
  n_s < 0.11 sigma). State the scope (NLA, linear bias, no systematics)
  and that neither code is compared with the exact integral.
- TATT (2026-10-02): CoCoA IA_model 1 with CFASTPT (IA_code 0) vs pyccl
  EulerianPTCalculator; recipe and conventions in references/tatt_conventions.md
  (variant prefix "tatt:" in ccl_compute.py). Check first: TATT with
  A2 = b_TA = 0 must reproduce the NLA reference (8e-5). Results: xi_pm agree
  as with NLA (0.0018 / 0.058) against a TATT signal of 1677 / 1293; gamma_t
  and w keep the FKEM offset; PT tables converged (320 vs 160 per decade,
  < 1e-6). High-l C_EE (l > 5e4) differs where TATT dominates: CoCoA's TATT
  kernels end at k = 334 h/Mpc (1e6 H0/c), CCL extrapolates; no effect on
  xi_pm above 2.5'. FAST-PT needs an even number of k nodes. A variant name
  with a prefix must be compared through its base name (the harmonic branch
  once tested `variant == "harmonic"` and silently ran the real-space run).
  The repository's benchmark scripts had use_A_ia=True on the IA-only tracer
  (IA normalization twice, sign flipped), b_TA fixed, no B-modes: fixed in
  b8bfec7.
- Final numbers (2026-10-01): LSST-Y1 8.76 (6.6-10.5 over five models;
  gamma_t 7.29, w 6.26, shear 0.0013), separable diagnostic 0.21;
  Roman-Real 0.228 (0.17-0.29), diagnostic 0.062. Benchmark-script
  modeling: 154 and 56.

## 4. Plots (owner's preferences, learned the hard way)

The general figure style (notebook rcParams, plotter defaults, the rules
the owner applies by eye) is in the Cocoa skill, Section 8.8, "The
maintainer's figure style": read it first. This study's figures:

Use cosmolike_core `cosmolike_notebook_utils/plot_datavectors.py`
(`plot_xi`, `plot_gammat_tomo_limber`, `plot_wtheta_tomo`) in ratio mode,
with these adaptations (all in `plots.py`):

- Plot (CCL - CoCoA)/sigma, sigma = sqrt(C_ii) (1 = one sigma). Curves
  enter as 1 + that, reference = ones (0 for excluded pairs). Ratios are
  rejected: they spike where w(theta) and gamma_t cross zero.
- Only real CCL vs CoCoA in figures (one curve per cosmology). Modeling
  variants live in tables only ("way too many cases").
- One y-range PER ROW, fitted to the data, not symmetric, 0 included, 10%
  padding; outlier panels (spread > 2.5x the row median) get alpha and the
  panel prints "1/alpha = f" with f in 2, 3, 5, 10, 20, ... Pass the
  rows' bands through the plotters' per-row `ylim` (a list of one
  [1 + lo, 1 + hi] per row; cosmolike_core df58d62). Never monkey-patch
  matplotlib or the plotters (Cocoa skill, Section 8.5): a missing option
  is added to the plotter. If every panel of a row needs alpha, the range
  is wrong.
- Legend above the panels in ONE row; bin label and 1/alpha in the free
  left corner (top or bottom).
- Large fonts: y label 24, y ticks 19, x ticks 22, x label 24, bin text
  22, legend 22; dpi 180.
- Colors: twilight_shifted without its pale middle, luminance capped at
  0.5; fiducial solid; the lightest colors get the long dash and wide
  lines, the darkest the dots (a light dotted thin line is invisible).
