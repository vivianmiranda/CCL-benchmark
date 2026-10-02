<!-- Fable 5 review D (TATT section precision), 2026-10-02. Saved for reuse; its four
     wording fixes were applied in the commit after 1f3402b. -->

# Fable 5 precision review D: TATT (code + claims)

Scope: `cocoa_comparison/scripts/ccl_compute.py` "tatt:" mode vs
`.claude/skills/cocoa-ccl-comparison/references/tatt_conventions.md`; the
README "## TATT intrinsic alignments" section and the Summary TATT bullet vs
`cocoa_comparison/results.md`, `harmonic.md`, the cmp/ run outputs, the
cosmolike C code, and commit b8bfec7. Read-only review; every number below
was recomputed from the npz outputs (compare.py / harmonic.py / plots.py
logic reproduced), GG for the finding-3 check via the study's pyccl build.

## Ranked issues

### 1. Finding 3: "Roman-Real sources 1-4" reads as a bin range (medium, wording)

README, finding 3: "In pairs with the first source bin, where TATT dominates
$C_{EE}$ (Roman-Real sources 1-4: TATT is $-2$ times NLA), the two differ by
8% at $\ell = 6.9\times10^4$ and 49% at $\ell = 10^5$".

"sources 1-4" means the single source pair (bin 1, bin 4), but it reads as
"source bins 1 through 4". Verified per pair (CoCoA harmonic, TATT/NLA
$C_{EE}$ ratio at $\ell = 6.9\times10^4$ / $10^5$): s1-s4 is $-1.93$/$-2.01$
("$-2$ times NLA" holds for exactly this pair); s1-s1 is $+6.1$, s1-s3
$-1.4$/$-1.5$, s2-s2 $+1.4$ — not $-2$. The 8% and 49% are also this one
pair (|CCL/CoCoA - 1| for s1-s4: 0.082 and 0.491); the worst
signal-carrying pair at $\ell = 6.9\times10^4$ is the s1 auto at 10%
(0.101), so "the two differ by 8%" understates the band if read as a
maximum. Replacement text (keeps every number, names the pair):

> In pairs with the first source bin, where TATT dominates $C_{EE}$ (in the
> Roman-Real source pair 1-4, $C_{EE}$ with TATT is $-2$ times its NLA
> value), the two differ: the 1-4 pair by 8% at $\ell = 6.9\times10^4$ and
> 49% at $\ell = 10^5$, the source-1 auto by 10% at $\ell = 6.9\times10^4$
> (2.5% at $\ell = 6.5\times10^4$ in LSST-Y1).

### 2. Finding 2: "the fiducial 8.89 and 0.237 ARE the NLA FKEM offset" (minor)

The NLA reference rows are 8.757 and 0.228 (results.md); the TATT rows are
8.890 and 0.237 — 1.5% and 4% higher (the TATT-extra residuals ride on
top). "are" overstates. Replacement:

> the fiducial 8.89 and 0.237 are essentially the NLA FKEM offset
> (NLA: 8.76 and 0.23)

### 3. Finding 1 / Summary bullet: "$\xi_\pm$ differ by 0.0018 ..., the NLA values" (minor)

Roman matches exactly (0.058 = 0.058) but LSST-Y1 is 0.0018 (TATT) vs
0.0013 (NLA): "the NLA values" / "as with NLA" is loose for LSST. Suggested
in finding 1: "..., at the NLA level (NLA: 0.0013 and 0.058)"; the Summary
bullet's "as with NLA" can stay. Same nit later in finding 3: "The
real-space $\xi_\pm$ difference equals the NLA one" — "matches the NLA one"
is accurate for both projects.

### 4. LSST-Y1 $\gamma_t$ figure caption: "at most $0.53\sigma$" (cosmetic)

Recomputed max over the five cosmologies (plots.py definition):
$0.534\sigma$, so "at most $0.53\sigma$" is literally exceeded. Write
"$0.54\sigma$" (or "$0.534\sigma$"). All other caption maxima are exact to
their quoted digits: LSST $\xi_+$ 0.020, $\xi_-$ 0.010; Roman $\xi_+$
0.049, $\xi_-$ 0.059, $\gamma_t$ 0.045.

## Verified correct (no action)

Code, `ccl_compute.py` "tatt:" mode vs the convention map:

- Amplitudes (lines 177-184): c1_ratio = 0.01389/(5e-14 RHO_CRITICAL),
  pivot (1+z)/1.62, a1delta = b_TA*A1(z), Om_m2_for_c2=False; z grid
  0..6 (601 pts) covers all kernels (sources end below z = 4; PTTracer and
  the Pk2D sampling hold beyond). Matches map Sec 2/6.
- EulerianPTCalculator (186-191): log10k_min=-5, log10k_max=2.3, a_arr =
  the export's pk a grid, with_matter_1loop=False, low_extrap=-6,
  high_extrap=3, pad_factor=1.0, P_window=[0.2,0.2], C_window=0.65,
  b1_pk_kind="nonlinear", no k_cutoff. nk_per_decade=160 instead of the
  map's 150 is documented (FAST-PT needs an even count; int(7.3*150)=1095
  is odd; 7.3*160=1168) and convergence-checked by tatt:pt_hi (7.9e-9 /
  4.0e-7).
- m:cdelta / m:c2 templates (196-201): verified in the PR build's ept.py
  that get_pk2d_template returns g4*(a00e+c00e) and g4*(a0e2+b0e2), i.e.
  pk_mi_extra = P_GI minus its c1*P_nl term, exactly _get_pgi's bracket at
  b1=1 (exp_cutoff = 1 with k_cutoff=None). The interp of cdt, c2t onto the
  template a grid is correct.
- xi_pm (213-217, 284-287): C_EE = GG(nonlinear table) + GI + IG + II,
  C_BB = II_BB (return_ia_bb=True); xi+ from EE+BB (GG+), xi- from EE-BB
  (GG-). Matches CoCoA's assembly.
- gamma_t (218-220, 296-298): the reference FKEM call with the NLA
  composite source tracer (C1 only, as CoCoA's non-Limber C_gs), TATT
  extras Limber-only at every ell. Harmonic tatt mode (232-249) is the same
  assembly at the exported multipoles, cssbb saved.

Claims, all confirmed against results.md / harmonic.md / npz / C code / git:

- Every number in the two README TATT chi2 tables and the harmonic TATT
  table matches results.md / harmonic.md (README rounds 1741.599 -> 1742,
  1292.799 -> 1293, 4.8e-3 -> 0.005, 7.7e-3 -> 0.008, 0.222 -> 0.22).
  The harmonic TATT table reproduces from the npz with harmonic.py's
  stated 1%-signal cut (top-band maxima 2.5e-2, 4.9e-1 recovered).
- TATT point A1=0.7, eta1=-1.7, A2=-1.36, eta2=-2.5, b_TA=1: in the
  exports' model string. FAST-PT 4.0.0: confirmed in the ccl env.
- "With A2 = b_TA = 0 the tatt: variant reproduces the NLA reference to
  Delta chi2 = 8e-5": recomputed 8.096e-5 (tattcheck vs ref, masked icov).
- Convention-table PT row: FPTIA.krange = [0.05, 1e6] (k in H0/c units) =
  [1.67e-5, 334] h/Mpc at pt_cfastpt.c:368-369, N = 1100 at default
  FPTboost (line 372); kernels zeroed outside the table (cosmo2D.c: KIA
  zero3d + the lnk-in-[limTATT] gate); FPTIA.k_cutoff = 1e4 is assigned
  (line 370) and never applied (grep: assignments only). CCL side: 160 per
  decade, 1e-5..200 Mpc^-1, Pk2D extrapolates (extrap_order_hik=2) — as
  finding 3 states.
- Finding 4: commit b8bfec7 exists in CCL-benchmark and its message states
  exactly the three listed changes (use_A_ia=False with the double
  normalization + sign flip, b_TA as a parameter, C_BB in xi_pm) and that
  the timings predate it.
- "these multipoles do not reach the data vector above 2.5'": theta_min =
  2.5 arcmin in both projects' meta; the TATT xi_pm rows equal (Roman) or
  nearly equal (LSST, issue 3) the NLA ones.
- Summary bullet numbers (0.0018 / 0.058; 1677 and 1293) match results.md.
