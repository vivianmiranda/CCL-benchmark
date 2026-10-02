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
   docstring: ref, limgs, limber, norsd, flat, points, separable
   (diagnostic), ccl_numerics, bench, harmonic). Needs `COCOA_ROOTDIR`.
3. `results.py` -> results.md/json (Delta chi2 = d^T C^-1 d, masked
   inverse covariance, per probe); `plots.py <outdir>` (cocoa env).

Models: fiducial, Omega_m 0.25/0.35, n_s 0.92/1.01.

## 2. Matching conventions (check each, state each in the README)

| item | CoCoA | CCL setting used |
|---|---|---|
| n(z) z column | read as Z_LOW (value at z + dz/2) | z + dz/2 passed to tracers |
| NLA amplitude | A1 ((1+z)/1.62)^eta Omega_m 0.01389/D | ia_bias scaled by 0.01389/(5e-14 RHO_CRITICAL) |
| galaxy bias | b1 per bin, constant | bias=(z, b1) |
| non-Limber switch | l < 150 (gg and gs) | l_limber = 150, FKEM |
| RSD | in gg (Limber and non-Limber) | has_rsd=True |
| transform | full-sky, bin-averaged | PR #1296 legendre with theta_max |
| spin-2 l factor | exact | exact (WL tracer prefactor) |

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
  CoCoA v5.02 uses the separable form consistently on both terms.
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

## 4. Plots

Use cosmolike_core `cosmolike_notebook_utils/plot_datavectors.py`
(`plot_xi`, `plot_gammat_tomo_limber`, `plot_wtheta_tomo`) in ratio mode:
curves = CCL/CoCoA (one per model), reference = ones (zero for excluded
pairs). Copy the notebooks' style exactly (rcParams block of
EXAMPLE_EVALUATE notebooks, twilight_shifted, linewidth/linestyle lists,
bintextsize 20, axis labels 17, legend 17, figsize (18, 13) grids and
(18, 13/5) for w, tight ylim from the data). The owner rejects plots with
white space, small fonts or small captions.
