<!-- Fable 5 review B (writing_and_figures) of the CoCoA vs DESC-CCL study, 2026-10-01/02.
     Saved so later sessions reuse it instead of re-running Fable. Its fixes were
     applied in CCL-benchmark commit 6e78e37 and later; line numbers refer to the README
     and sources as they were then (cosmolike_core v5.03 + 4, PR #1296 647ad4a). -->

# Fable review, part B: writing and figures (CCL-benchmark/README.md)

Scope: README.md from "# CoCoA vs DESC-CCL: real-space 3x2pt code comparison"
down; the eight figures in cocoa_comparison/figures/ (all eight viewed, with
pixel-level crops of the panels behind each caption number); rules from
cocoa/.claude/skills/cocoa-maintenance/SKILL.md Section 3 and
CCL-benchmark/.claude/skills/cocoa-ccl-comparison/SKILL.md sections 0 and 4.
Read-only review; no repository file touched.

## 1. Caption numbers vs what the figures show (the one factual problem)

The per-figure "At most X sigma" numbers disagree with the figures once the
1/alpha panels are unscaled (the Figures intro itself tells the reader to
multiply). Checked panel by panel:

| caption | claim | figure reads | verdict |
|---|---|---|---|
| LSST-Y1 gamma_t | reaches -0.5 sigma near 100' | rows (.,4)/(.,5) dip to -0.45..-0.5; (1,2) 1/alpha=10 x -0.05 = -0.5 | consistent |
| LSST-Y1 w | -0.5 to -1.3 sigma below 100'; 0.84 sigma on kept points | panel (5) minimum approx -1.35 to -1.4 (below the -1.2 tick by ~0.4 tick) | low end understated; "to -1.3" should be approx -1.4 |
| LSST-Y1 xi+- | at most 0.01 sigma | xi+ (0,0)/(0,1) Omega_m=0.35 peak approx 0.013; xi- (0,0) approx 0.010 | 0.01 too low; approx 0.014 |
| Roman gamma_t | at most 0.03 sigma | panel (3,8): Omega_m=0.35 exceeds the 0.04 tick (crop confirms), approx 0.045 | 0.03 wrong; approx 0.05 |
| Roman w | at most 0.09 sigma | panels (7)/(8) dip to approx -0.085/-0.09 | consistent |
| Roman xi+- | at most 0.04 sigma | xi+ (7,7) 1/alpha=3 x approx 0.017 = 0.05; xi- (7,7) 1/alpha=2 x approx 0.027 = 0.054 (curve is above the 0.024 top tick, so >= 0.048 is certain) | 0.04 wrong; approx 0.06 |

results.json stores only Delta chi^2, not per-point (d_CCL - d_CoCoA)/sigma,
so the exact maxima must be recomputed from the run folder's npz files (the
arrays plots.py already builds). Recommendation: recompute max |Delta/sigma|
per figure (over kept points for Roman-Real, all points for LSST-Y1, mask
status stated) and put the recomputed value in each caption; my figure-read
estimates above say which captions change. A cheap guard: have plots.py print
the per-figure max next to its existing "1/alpha != 1 in N panels" line, so
the caption number has a printed source.

## 2. Bin-pair indexing (text vs figures, and figure vs figure)

Three conventions coexist:

- gamma_t and w figures: 1-indexed panels, (1,1)..(5,5), (1)..(8), from the
  gammat/wtheta plotters;
- xi+- figures: 0-indexed panels, (0,0)..(4,4), (0,0)..(7,7), from plot_xi;
- README text: 0-indexed pairs, stated once ("(6,0), (7,0) and (7,1),
  0-indexed" in Layout), then used bare in Results ("lens 7 - source 2") and
  Finding 4 ("lens 4 - source 0"), while the Roman gamma_t caption switches
  to 1-indexed ("lens 8 - source 3 (1-indexed)").

A reader going from Finding 4 to the figure must translate twice. Cheapest
repair in text only (no replot): one sentence in the Figures intro declaring
both conventions, plus the panel coordinate wherever the text names a pair.
Proposed intro sentence (after "One curve per cosmology."):

> gamma_t and w panels are labeled (lens, source) and (lens), 1-indexed;
> xi_pm panels are labeled (bin i, bin j), 0-indexed, CoCoA's order.

and in Finding 4:

> The pairs: LSST-Y1 lens 4 - source 0 and Roman-Real lens 7 - source 2
> (0-indexed; panel (8, 3) of the Roman-Real gamma_t figure).

The clean fix is regenerating all figures with one convention (0-indexed,
matching the text and the exclusion lists in scripts/), but that is a part-A
(code) change to plots.py's label handling.

## 3. Figures vs skill section 4 (compliance list)

Checked on all eight PNGs:

- Legend above the panels in one row: yes, all eight.
- Fonts and dpi: consistent with the skill's sizes (y label/x label large,
  readable tick labels); no violation seen.
- One y-range per row, zero included, not symmetric: yes, all grids.
- 1/alpha panels print the factor in the free left corner, f drawn from the
  2, 3, 5, 10, 20, ... ladder (observed 2, 3, 5, 10, 20, 200, 300): yes.
- No row has alpha on every panel: holds everywhere I checked (worst is
  LSST gamma_t row 3 with 2 of 5).
- Line styles: fiducial solid dark; darkest variation dash-dot-dot (dense);
  lightest colors carry the wide long dashes; no light dotted thin line.
  Compliant.
- Violation: Roman gamma_t panels (7,1), (8,1), (8,2) print "excluded" AND
  "1/alpha = 300" / "1/alpha = 200" on empty panels. The alpha fit ran on a
  pair whose curve is zeroed by construction (skill: reference = ones, 0 for
  excluded pairs), so the printed factor is an artifact and the number is
  meaningless to a reader. Either suppress the alpha label on excluded
  panels (plots.py change, part A) or say in the caption what "excluded"
  means. The caption currently explains only the blank (8,3) panel, not the
  three "excluded" ones. Proposed caption sentence (Roman gamma_t):
  > Panels marked "excluded" are CoCoA's gamma_t exclusions (Layout above);
  > pair (8, 3) is blank because DESC-CCL's transform fails (finding 4).

## 4. Section-3 rule pass, sentence by sentence

Clean on: hard-zero words and phrases (grep: zero hits), adjectives and
judgment adverbs (none found), em dashes (none), LaTeX math (consistent:
$\Delta\chi^2$, $\xi_\pm$, $\sigma$, settings in code spans), no raw | inside
table-cell math, tables over prose (the Results and What-it-took sections are
exactly the house shape), no measured-value staleness problem (this page IS
the measurement record of a dated study against v5.02; the numbers are its
content, not test documentation).

Findings that remain:

1. No Contents list. Seven sections; Section 3.1 requires a numbered
   Contents list with `<a name>` anchors for any README over three sections.
2. Summary duplication. The FKEM mechanism paragraph (Summary bullet 2)
   restates Finding 2; the diagnostic 0.21/0.062 appears four times (Summary
   table column, Summary bullet 3, one-change table, Finding 2); RSD 115 and
   benchmark 154/56 each appear three times (Summary bullet 4, one-change
   table, Finding 5). A summary may repeat the headline, but bullets 2-4
   carry full mechanism text. Proposed: compress Summary bullets 2-4 to one
   sentence each with a pointer ("finding 2", "finding 5"), keep the table.
3. Vague column header "after the diagnostic below" (Summary table): points
   at a bullet instead of naming the thing. Proposed: "separable
   $P_{\rm lin}$ diagnostic", with bullet 3 (or its one-line survivor)
   defining it.
4. Reproduce section shape: three bare command blocks with comment lines;
   Section 3.1 demands `**Step :one:**:` blocks, an assumption paragraph in
   the flow, and self-containment. Missing facts a runner needs:
   - commands use relative paths, so the working directory (the
     CCL-benchmark repository root) is never stated;
   - `<run folder>` is created by step 1 (run_cocoa.sh does `mkdir -p`) and
     must be the same folder in steps 2-3 - not stated;
   - `CCL_PR=<PR build>` is a PYTHONPATH entry pointing at the PR #1296
     pyccl build (run_ccl.sh's own comment says so) - "<PR build>" alone is
     not runnable;
   - step 3's plots.py needs `start_cocoa.sh` sourced (ROOTDIR set), per its
     docstring; the README says only "in the cocoa environment";
   - results.py / harmonic.py environment unstated (any Python with NumPy).
5. "the harness" (What-it-took rows 3-4) is a coined label never defined.
   Replace with the file the reader can open: `cocoa_export.py` (or "the
   export script" defined at first use).
6. One name per object: "the tables cosmolike receives" (Method bullet 1)
   vs "CoCoA" everywhere else. Use CoCoA.
7. Internal names: `init_data_real` (CoCoA-side bullet) is internal
   bookkeeping - "cosmolike reads the data files once per process" says the
   same without it. `ccl_f1d_extrap_logx_logy` in Finding 4 is the precise
   locus of the reported failure (evidence, owner's rule 2) - keep it.
8. Finding 2's closing sentence lists 1.98 / 11.2 / 2.06 without naming the
   quantity; prepend "$\Delta\chi^2$ of the non-Limber effect ...".
9. Minor wording: "13 kept Roman-Real points are left out" (kept ... left
   out) - "13 points the mask keeps are left out of the comparison".
10. Method's DESC-CCL-build bullet and What-it-took row 1 state the PR #1296
    / error 1040 fact twice; the row adds the build detail (CMake, SWIG), so
    trim the Method bullet to the pointer and the commit.

Not flagged (deliberate judgment): e-notation (2.7e-4) in the harmonic
table matches results.md and the NOTE block's number style; "both codes in
Limber" names a two-member set already named in the title; the Figures
intro's "1 is a one-sigma difference" is exactly the didactic line a caption
needs.

## 5. Caption self-sufficiency check (task item 2)

The shared Figures intro carries: what is plotted, sigma's definition, the
1/alpha rule, per-row y-ranges, and the masking asymmetry (LSST shows every
theta, Roman blanks masked points). With the indexing sentence of item 2 and
the "excluded" sentence of item 3 added, each caption plus the intro reads
the figure. One optional clause: why Roman blanks masked points ("its
covariance is defined only on the mask", already stated in the CoCoA-side
bullets) - a cross-reference would close the "why the asymmetry" question a
reader asks at the Figures section.

## 6. Structure (task item 4)

Summary -> Method -> Results -> Findings -> Figures -> What it took ->
Reproduce is right for a reader arriving from the main Cocoa README: verdict
first, the two tables a citer quotes, mechanisms, pictures, then operational
detail. Captions reference findings backward (Figures after Findings), which
works. The only forward reference is inside Summary ("the diagnostic below",
item 4.3). Things said twice are listed in item 4.2 and 4.10. No section
needs moving.

## 7. Ranked list

1. HIGH - caption maxima wrong (Figures section, three captions):
   Roman gamma_t "At most $0.03\sigma$", Roman xi "At most $0.04\sigma$",
   LSST xi "At most $0.01\sigma$", LSST w "to $-1.3\sigma$". Recompute from
   the npz and replace; figure-read values: 0.05, 0.06, 0.014, -1.4.
2. HIGH - indexing: add the convention sentence to the Figures intro and the
   panel coordinate in Finding 4 (exact text in item 2); or regenerate the
   figures 0-indexed (part A).
3. MEDIUM - Roman gamma_t "excluded" panels print meaningless 1/alpha; add
   the caption sentence (item 3) or suppress the label in plots.py (part A).
4. MEDIUM - Reproduce: Step :one: blocks + assumption paragraph + the five
   missing facts (item 4.4).
5. MEDIUM - add the numbered Contents list with anchors (item 4.1).
6. LOW - Summary compression and "after the diagnostic below" header
   (items 4.2, 4.3); "the harness" -> `cocoa_export.py` (4.5);
   "cosmolike" -> CoCoA (4.6); `init_data_real` rephrase (4.7); name
   $\Delta\chi^2$ in Finding 2's last sentence (4.8); wording 4.9; Method
   build-bullet trim (4.10).
