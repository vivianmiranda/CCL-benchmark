# README audit, 2026-10-06

## Scope

Documentation-only audit of the saved October 1–2 comparison. No numerical
runs, compiler jobs, downloads or dependency installation were performed.
Saved results, figures and numerical scripts remain unchanged.

## Findings

- Saved comparison: CoCoA v5.02/core26191ec and CCL PR1296 commit647ad4a.
  The original top note was a different Intel/default benchmark; the
  README now leads with the saved study's scope and keeps the M2 timings
  together. The removed Intel note is preserved below.
- Current LSST/Roman likelihood prototypes read growth_k, default0.05/Mpc.
  LSST project commit a257e7e introduced the new reference. The comparison
  cocoa_export.py still fixes D0/DD at5e-4/Mpc. Consequently its current-code
  export does not supply CCL the same growth input as CoCoA. The README
  states this before the scientific conclusions and reproduction commands.
- The growth_sub variant alters CCL growth while retaining the old exported
  IA amplitude. It is a diagnostic of the old comparison, not a measured
  comparison with today's CoCoA growth convention.
- set_installation_options.sh assigns CCL_PATH unconditionally. Exporting a
  shell value before setup would be overwritten. README instructs editing
  the key in the options file, consistent with the actual script.
- Lifecycle scripts preserve the pinned PR and old CMake recipe. Their
  syntax/activation checks do not establish a successful clean native build.
- All twelve linked figures and saved result/timing files exist. Figures
  were not regenerated. The README no longer cites the Claude skill.

## Why the growth reference changed

The README retains the original Section4 growth/FKEM derivation, dark-energy
on/off tables, growth_sub diagnostics and per-bin pivot comparisons. A
follow-up subsection connects them to the current CoCoA default; a mere
mismatch warning does not explain the physical motivation.

Evidence:

- LSST project commit a257e7ec9c1f8c68df631afbb654b255929f2619 changes the
  fixed linear-power reference from5e-4 to0.05/Mpc in both table construction
  and normalization. It records the sub-horizon motivation. Current LSST
  and Roman prototype files retain that choice.
- The README's existing massless-neutrino CAMB dark-energy toggle isolates
  the horizon-scale step for w=-0.9. Do not recast this as proof of a gauge
  bug, a problem specific to all non-Lambda cosmologies, or a new numerical
  result from the current checkout.
- CAMB transfer-variable documentation specifies synchronous-gauge density
  inputs: https://camb.readthedocs.io/en/latest/transfer_variables.html
- Lesgourgues & Pastor, arXiv1212.6154v1, Section0.6.3 describes the separate
  massive-neutrino free-streaming effect and scale-dependent linear growth:
  https://arxiv.org/html/1212.6154v1
- A sub-horizon reference does not make the true P(k,z) separable. The old
  growth_sub experiment keeps the old IA amplitude, whereas changing
  CoCoA's growth table also changes its IA factors. Fresh current-code
  measurements are still required.

## Work required for a current-code numerical comparison

1. Record current Cocoa/core/project commits and all resolved settings.
2. Make the exported growth table follow the same growth_k and normalization
   convention used by each current likelihood. Check NLA/TATT amplitudes;
   do not relabel growth_sub outputs as a current-model comparison.
3. Confirm the current growth convention in the timing driver as well as
   the exporter. Both read live likelihood defaults, but only the exporter
   separately constructs CCL's growth input.
4. Regenerate five-cosmology NLA/TATT vectors, within-code refinements,
   Limber/transform isolation and growth diagnostics; preserve old results.
5. Recompute figures and tables from those new outputs, and time sequentially
   on a quiet machine. Do not reuse saved timing numbers as new measurements.
6. Test a fresh pinned Conda/private-environment build when dependency
   installation and compiler work are authorized.

## Verification

- Markdown-it commonmark with tables enabled rendered README.md to HTML:
  26 tables and 5 blockquotes.
- All34 internal anchors and23 local file/image links resolve;12 image links.
- All25 shell blocks contain one command each.
- Four timing rows match timing_macos.txt exactly at its recorded precision;
  the integer ratios match rounded CCL/CoCoA values.
- Headline fiducial statistics match results.json.
- No numerical source, saved result or figure is changed; git diff --check
  passes. This audit did not check remote URLs or execute numerical recipes.

Follow-up explanation checks: the full original table rows and image links
remain unchanged; Markdown-it renders26 tables and5 blockquotes, and all36
internal anchors and23 local links resolve. No numerical code was edited.
The CAMB and arXiv primary-source pages were opened during this follow-up.

## Archived Intel note

This note preceded the study in README.md through commit deda88b. It is
archived here because its modeling and sampling differ from the M2 study.

> [!NOTE]
> CoCoA `v5.02` benchmark (cosmolike only times)
> CPU: `Intel(R) Core(TM) i9-10940X CPU @ 3.30GHz` (`1/8 OpenMP cores`).
>
> Modeling: Full-sky on real functions except for DESC-CCL (unknown for CLOE-LIB).
>
> Modeling: **IA=TATT** in ($\xi_{\pm}, \gamma_t$) except in CLOE-LIB and DES-Y3-Real 6x2pt+N.
>
> Modeling: **Non-limber** $C_{gg}(l)$ in real space; **non-limber** $C_{gs}(l)$ in CoCoA LSST-Y1/Roman-Real.
>
> Modeling: Roman-Fourier and Roman-Real-KL CoCoA compute the exact (**non-Limber**) $C_{gg}(l)$ and $C_{gs}(l)$
> below $l = 150$ with RSD; DESC-CCL uses Limber without RSD.
>
> - **LSST-Y1-Real 3x2pt**: (CoCoA) `0.27/0.06s`, (DESC-CCL)`7.96/1.72s`, (CLOE-LIB) 0.23/0.23s.
> - **Roman-Real 3x2pt**: (CoCoA) `0.45/0.09s`, (DESC-CCL) `8.17/1.96s`, (CLOE-LIB) 0.27/0.27s.
> - **Roman-Fourier 3x2pt**:  (CoCoA) `0.27/0.08s`, (DESC-CCL) `0.65/0.36s`.
> - **Roman-Real-KL 3x2pt**: (CoCoA) `0.162/0.08s`.
> - **DES-Y3xPlanck 6x2pt (w/ CMB)**  (CoCoA) `0.34/0.07s`.
> - **DES-Y3-Real 3x2pt**  (CoCoA) `0.27/0.06s`.
> - **DES-YX-Real 6x2pt+N (clusters)**  (CoCoA) `0.71/0.14s` (YX = not yet production cov, n(z), dv, Y6 analysis).
>
