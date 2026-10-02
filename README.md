> [!NOTE]
> (From CoCoA notes)
> `v4.11.1` benchmark: do not include CAMB (or the Hybrid Emulator); **Includes TATT in** ($\xi_{\pm}, \gamma_t$) **and non-limber in** $w_{gg}(\theta)$.
> CLOE-LIB caveat: We were not able to make TATT work on CloeLib. We were also not able to speed-up cloelib with `OMP_NUM_THREADS` flag
>
> CPU: `Intel(R) Core(TM) i9-10940X CPU @ 3.30GHz` (`1/8 OpenMP cores`). *Times are approximate*.
>
> - **LSST-Y1-Real 3x2pt**: (CoCoA) `0.29/0.06s`, (DESC-CCL)`7.96/1.72s`, (CLOE-LIB) 0.23/0.23s. **CoCoA speed-up (CCL)**: `27/28x`
> - **Roman-Real 3x2pt**: (CoCoA) `0.45/0.095s`, (DESC-CCL) `8.17/1.96s`, (CLOE-LIB) 0.27/0.27s. **CoCoA speed-up (CCL)**: `18/20x`
> - **Roman-Fourier 3x2pt**:  (CoCoA) `0.08/0.03s`, (DESC-CCL) `0.65/0.36s`. **CoCoA speed-up**: `7.5/21x`
> - **DES-Y3xPlanck 6x2pt**  (CoCoA) `0.40/0.075s`
> - **DES-Y3-Real 3x2pt (des_y3 repo)**  (CoCoA)~`0.25/0.05s`
> 

# CoCoA vs DESC-CCL: real-space 3x2pt code comparison

CoCoA `v5.02` and DESC-CCL compute the same real-space 3x2pt data vectors
($\xi_\pm$, $\gamma_t$, $w(\theta)$) of LSST-Y1 and Roman-Real at five
cosmologies. The differences are measured with each project's covariance.
Scripts, tables and figures: [`cocoa_comparison/`](cocoa_comparison).

Contents:
1. [Summary](#ccl_summary)
2. [Method](#ccl_method)
3. [Results](#ccl_results)
4. [Findings](#ccl_findings)
5. [Figures](#ccl_figures)
6. [TATT intrinsic alignments](#ccl_tatt)
7. [What it took to run DESC-CCL](#ccl_run)
8. [Reproduce](#ccl_reproduce)

## Summary <a name="ccl_summary"></a>

$\Delta\chi^2 = (d_{\rm CCL} - d_{\rm CoCoA})^{T} C^{-1} (d_{\rm CCL} - d_{\rm CoCoA})$,
with the projects' masks (CoCoA's tests pass at $\Delta\chi^2 < 0.2$):

| project | five cosmologies | fiducial: $\xi_\pm$ / $\gamma_t$ / $w(\theta)$ | separable $P_{\rm lin}$ diagnostic |
|---|---|---|---|
| LSST-Y1 | 6.6 to 10.5 | 0.0013 / 7.29 / 6.26 | 0.21 |
| Roman-Real | 0.17 to 0.29 | 0.058 / 0.27 / 0.096 | 0.062 |

- Cosmic shear agrees in both projects.
- $\gamma_t$ and $w(\theta)$ differ through DESC-CCL's non-Limber (FKEM)
  term when the linear growth depends on $k$ ([finding 2](#ccl_findings)).
  The difference follows $w$: LSST-Y1 at $w = -1$ gives 0.34, Roman-Real at
  $w = -0.9$ gives 7.6.
- In LSST-Y1 the difference moves the lens biases by $-0.85\sigma$ to
  $-1.17\sigma$; $\Omega_m$ and $n_s$ move by $0.07\sigma$ and $0.11\sigma$
  ([parameter shifts](#ccl_results)).
- With TATT intrinsic alignments the codes agree on the TATT terms:
  $\xi_\pm$ differ by 0.0018 (LSST-Y1) and 0.058 (Roman-Real), as with NLA,
  while TATT moves $\xi_\pm$ by 1677 and 1293 ([TATT](#ccl_tatt)).
- 0.2 is CoCoA's tolerance between runs of one code, quoted as a scale, not
  as a cross-code acceptance criterion.
- The modeling of this repository's DESC-CCL benchmark scripts differs from
  CoCoA by $\Delta\chi^2 = 154$ (LSST-Y1) and $56$ (Roman-Real)
  ([finding 5](#ccl_findings)).

## Method <a name="ccl_method"></a>

- **Same inputs.** DESC-CCL receives, through `pyccl.CosmologyCalculator`,
  the tables CoCoA receives at the same point:
  - CAMB's linear and nonlinear $P(k,z)$ ($z \le 6$,
    $10^{-5} \le k \le 200\ {\rm Mpc}^{-1}$);
  - $\chi(z)$ and $H(z)$ (from $d\chi/dz$) to $z = 50$;
  - the growth factor CoCoA uses,
    $D(z) = \sqrt{P_{\rm lin}(k_0, z)/P_{\rm lin}(k_0, 0)}$ with
    $k_0 = 5\times 10^{-4}\ {\rm Mpc}^{-1}$, to $z = 49$, and
    $f = d\ln D/d\ln a$.

  No Boltzmann or background difference enters.
- **DESC-CCL build.** pyccl from
  [LSSTDESC/CCL PR #1296](https://github.com/LSSTDESC/CCL/pull/1296) (commit
  `647ad4a`): full-sky $\xi_\pm$ and exact bin averaging.
- **Layout.**
  - Bins: CoCoA's log $\theta$ bins (LSST-Y1: 26 bins from 2.5' to 900';
    Roman-Real: 15 bins from 2.5' to 250').
  - Pairs: CoCoA's bin pairs, ordering and $\gamma_t$ exclusions
    (Roman-Real: lens-source pairs (7, 1), (8, 1) and (8, 2)). Bins count
    from 1 in the text and the figures.
  - Masks: LSST-Y1 `lsst_y1_M1_GGLOLAP0.05.mask` (959 of 1560 points);
    Roman-Real `example1.mask` (1950 of 2115 points).
- **Models.**
  - Likelihood: each project's `EXAMPLE_EVALUATE2.yaml` (NLA, linear bias),
    with photo-z shifts, shear calibration, magnification and point masses
    set to zero.
  - Fiducial: $A_s = 2.1\times10^{-9}$, $n_s = 0.96605$, $H_0 = 67.32$,
    $\Omega_m = 0.3$, $\Omega_b = 0.04$, $m_\nu = 0.06$ eV;
    $w = -0.9$ (LSST-Y1) and $w = -1$ (Roman-Real).
  - Variations: $\Omega_m = 0.25, 0.35$ and $n_s = 0.92, 1.01$, at fixed $A_s$;
    each project also at the other project's $w$ (LSST-Y1 at $w = -1$,
    Roman-Real at $w = -0.9$). $m_\nu$ is not varied.
  - Scope: intrinsic alignment NLA in both codes in every section except
    [TATT](#ccl_tatt), which repeats the comparison with TATT; linear bias;
    every systematic at zero. Magnification, photo-z shifts, shear
    calibration and baryons are outside this comparison.
- **DESC-CCL settings** (reference run):
  - `angular_cl(..., l_limber=150, non_limber_integration_method="FKEM", fkem_Nchi=2000)`
    for $C_{gg}$ and $C_{gs}$, the same switch as CoCoA;
  - $\ell$: every integer from 2 to 400, then 1500 log-spaced values to
    $\ell_{\max} = 6.5\times10^4$ (LSST-Y1) or $10^5$ (Roman-Real), CoCoA's
    `lmax`, also used as `ELL_MAX_CORR`;
  - `correlation(..., theta=lower edges, theta_max=upper edges, method="legendre")`.

| convention | CoCoA | DESC-CCL setting |
|---|---|---|
| $n(z)$ | first column read as the bin's lower edge | `dndz` at $z + \Delta z/2$ |
| NLA | $A_1[(1+z)/1.62]^{\eta}\,C_1\rho_{\rm crit}\,\Omega_m/D$, $C_1\rho_{\rm crit} = 0.01389$ | `ia_bias` $\times\ 0.01389/(5\times10^{-14}\rho_{\rm crit})$, `use_A_ia=True` |
| bias | $b_1$ per lens bin | `bias=(z, b1)` |
| non-Limber | $\ell < 150$ in $C_{gg}$ and $C_{gs}$ | `l_limber=150`, FKEM |
| RSD | in $C_{gg}$ only (`include_RSD_GS = 0`) | `has_rsd=True` in the $w(\theta)$ tracer, `False` in the $\gamma_t$ tracer |
| real space | full-sky Legendre sums, bin-averaged | PR #1296, `method="legendre"` with `theta_max` |

DESC-CCL's `has_rsd` applies RSD to every correlation of a tracer, so the
$\gamma_t$ lens tracer has it off to match CoCoA's model. Switching it on
moves DESC-CCL's $\gamma_t$ by $\Delta\chi^2 = 0.36$ (LSST-Y1) and $0.076$
(Roman-Real).

## Results <a name="ccl_results"></a>

DESC-CCL (reference settings) vs CoCoA:

| model | LSST-Y1 3x2pt | $\xi_\pm$ | $\gamma_t$ | $w(\theta)$ | Roman-Real 3x2pt | $\xi_\pm$ | $\gamma_t$ | $w(\theta)$ |
|---|---|---|---|---|---|---|---|---|
| fiducial | 8.76 | 0.0013 | 7.29 | 6.26 | 0.228 | 0.058 | 0.269 | 0.096 |
| $\Omega_m = 0.25$ | 6.56 | 0.0005 | 4.87 | 5.48 | 0.166 | 0.023 | 0.193 | 0.087 |
| $\Omega_m = 0.35$ | 10.50 | 0.0032 | 9.44 | 6.46 | 0.291 | 0.123 | 0.360 | 0.098 |
| $n_s = 0.92$ | 9.07 | 0.0013 | 7.56 | 6.55 | 0.236 | 0.055 | 0.280 | 0.099 |
| $n_s = 1.01$ | 8.45 | 0.0014 | 7.06 | 5.97 | 0.215 | 0.055 | 0.259 | 0.091 |

Per-probe columns zero the other probes' entries of $d$; the cross-covariance
makes them not add up to the 3x2pt value. In Roman-Real, the pair lens 8 -
source 3 is left out: 13 points the mask keeps
([finding 4](#ccl_findings)). CoCoA's $\gamma_t$ on them is below
$3\times10^{-12}$ (lens behind source) and gives $\Delta\chi^2 = 2\times10^{-11}$
on its own.

Fiducial cosmology, one change at a time (3x2pt $\Delta\chi^2$):

| layer | comparison | LSST-Y1 | Roman-Real |
|---|---|---|---|
| CoCoA numerics | CoCoA default vs CoCoA high accuracy (`accuracyboost: 2`, `integration_accuracy: 1`) | 0.004 | 0.007 |
| DESC-CCL numerics | DESC-CCL default FKEM sampling, 500 log $\ell$ vs reference | 0.34 | 0.0005 |
| DESC-CCL numerics | DESC-CCL finer sampling (`fkem_Nchi=4000`, 3000 log $\ell$) vs reference | 0.0014 | 0.0002 |
| transform | DESC-CCL transform of CoCoA's Limber $C_\ell$ vs CoCoA in Limber | 0.009 | 0.0006 |
| Limber $C_\ell$ | both codes in Limber | 0.16 | 0.031 |
| non-Limber | DESC-CCL separable $P_{\rm lin}$ (diagnostic) vs CoCoA | 0.21 | 0.062 |
| non-Limber | DESC-CCL vs CoCoA at the other project's $w$ | 0.34 ($w = -1$) | 7.57 ($w = -0.9$) |
| convention | DESC-CCL, RSD also in $\gamma_t$, vs CoCoA | 8.70 | 0.258 |
| modeling | DESC-CCL $C_{gs}$ in Limber vs CoCoA | 9.23 | 0.71 |
| modeling | DESC-CCL without RSD vs CoCoA | 115.4 | 6.98 |
| modeling | DESC-CCL full-sky at bin centers (no bin averaging) vs CoCoA | 10.3 ($\xi_\pm$: 2.0) | 22.8 ($\xi_\pm$: 28.3) |
| modeling | DESC-CCL flat-sky FFTLog at bin centers vs CoCoA | 12.7 ($\xi_\pm$: 1.8) | 22.6 ($\xi_\pm$: 19.4) |
| $P(k)$ source | DESC-CCL Eisenstein-Hu + halofit vs DESC-CCL reference | 21.3 | 58.8 |
| composite | DESC-CCL with this repository's benchmark modeling vs CoCoA | 153.6 | 55.9 |

The first rows read as a ladder: the transforms agree (0.009); the Limber
$C_\ell$ leave 0.16, almost all in $w(\theta)$; with the separable table
DESC-CCL's non-Limber term stays near that level (0.21); with CAMB's table
it brings the difference to 8.76.

The benchmark-modeling row applies the modeling choices of
`ccl_test_lsst.py` and `ccl_test_roman.py` to CoCoA's layout:
- DESC-CCL's own Eisenstein-Hu linear and halofit nonlinear $P(k)$
  ($\sigma_8$ from CoCoA's linear table);
- $C_{gs}$ in Limber, no RSD, `l_limber=100`, `fkem_Nchi=500`;
- flat-sky FFTLog at the bin centers.

The row copies every modeling choice of those scripts except the intrinsic
alignment: the scripts use TATT, the row uses NLA, as every other row. It
describes the settings of the benchmark scripts behind the timings in the
note above, not DESC-CCL itself; the $P(k)$-source row isolates the one
ingredient unique to it.

Limber $C_\ell$ at 32 $\ell$ from 20 to $\ell_{\max}$,
$|C_\ell^{\rm CCL}/C_\ell^{\rm CoCoA} - 1|$ (median over bin pairs / maximum
over pairs with $|C_\ell|$ above 1% of the largest):

| project | probe | $20 \le \ell < 66$ | $66 \le \ell < 1000$ | $1000 \le \ell \le 5000$ | $5000 < \ell \le \ell_{\max}$ |
|---|---|---|---|---|---|
| LSST-Y1 | $C_{ss}$ | 2.7e-4 / 1.7e-3 | 8.4e-5 / 1.1e-3 | 7.4e-5 / 2.2e-3 | 1.8e-4 / 4.0e-3 |
| LSST-Y1 | $C_{gs}$ | 2.9e-4 / 2.2e-3 | 3.4e-4 / 3.8e-3 | 3.2e-4 / 5.7e-3 | 3.9e-4 / 3.7e-3 |
| LSST-Y1 | $C_{gg}$ | 5.4e-4 / 2.2e-3 | 3.0e-5 / 2.0e-4 | 2.5e-5 / 5.1e-5 | 3.3e-5 / 1.0e-4 |
| Roman-Real | $C_{ss}$ | 5.5e-5 / 6.5e-4 | 7.3e-5 / 5.8e-4 | 1.2e-4 / 1.0e-3 | 2.1e-4 / 9.3e-4 |
| Roman-Real | $C_{gs}$ | 3.0e-4 / 1.9e-3 | 3.6e-4 / 1.9e-3 | 3.3e-4 / 1.8e-3 | 5.8e-4 / 1.8e-3 |
| Roman-Real | $C_{gg}$ | 9.3e-4 / 4.3e-3 | 4.2e-5 / 1.0e-3 | 2.7e-5 / 3.0e-3 | 3.7e-5 / 3.5e-3 |

Non-Limber $C_{gg}$ (CoCoA's `C_gg_tomo`, DESC-CCL's FKEM at the reference
settings), $C_\ell^{\rm CCL}/C_\ell^{\rm CoCoA} - 1$, lowest to highest over
the lens bins:

| project | $\ell = 2$ | $\ell = 10$ | $\ell = 50$ | $\ell = 140$ |
|---|---|---|---|---|
| LSST-Y1 | $-1.60\%$ to $-0.93\%$ | $-1.67\%$ to $-0.89\%$ | $-1.74\%$ to $-0.80\%$ | $-1.67\%$ to $-0.93\%$ |
| Roman-Real | $-0.58\%$ to $-0.01\%$ | $-0.45\%$ to $-0.05\%$ | $-0.42\%$ to $+0.30\%$ | $-0.41\%$ to $-0.08\%$ |

Parameter shifts the fiducial difference would cause
(`scripts/bias.py`): Fisher matrix of CoCoA's data vector with derivatives
in $\Omega_m$, $n_s$, every lens $b_1$, $A_1$ and $\eta$ (finite
differences of the exports), shift $F^{-1} g$ with
$g_a = \partial_a d^{T} C^{-1} (d_{\rm CCL} - d_{\rm CoCoA})$, in units of
the marginalized error of this subspace (every other parameter fixed, so
the errors are smaller than in a full analysis):

| project | $\Delta\Omega_m/\sigma$ | $\Delta n_s/\sigma$ | lens $b_1$ | $\Delta\chi^2$ absorbed / left |
|---|---|---|---|---|
| LSST-Y1 | $+0.07$ | $-0.11$ | $-0.75\%$ to $-1.66\%$ ($-0.85\sigma$ to $-1.17\sigma$) | 6.74 / 2.02 |
| Roman-Real | $+0.08$ | $-0.04$ | $-0.03\%$ to $-0.06\%$ ($-0.07\sigma$ to $-0.13\sigma$) | 0.05 / 0.18 |

## Findings <a name="ccl_findings"></a>

1. **Limber layer and transforms.** The codes agree in harmonic space to
   $\ell_{\max}$ (table above) and in real space with both codes in Limber
   ($\Delta\chi^2$ 0.16 and 0.031). DESC-CCL's transform applied to CoCoA's
   Limber $C_\ell$ reproduces CoCoA's Limber real space to 0.009 and
   0.0006, so the both-Limber residual is in the Limber $C_\ell$ (mostly
   $C_{gg}$ below $\ell = 66$, in $w(\theta)$). CoCoA at high accuracy moves
   by 0.004 and 0.007.
2. **Non-Limber $C_{gs}$ and $C_{gg}$ (DESC-CCL FKEM).** DESC-CCL computes
   $$C_\ell = C_\ell^{\rm Limber}[P_{\rm nl}] - C_\ell^{\rm Limber}[P_{\rm lin}(k,z)] + C_\ell^{\rm FKEM}[D^2(z)\,P_{\rm lin}(k,0)].$$
   The FKEM term multiplies the kernels by DESC-CCL's growth factor and uses
   the linear table at $z = 0$ (`pyccl/nonlimber/_nonlimber_FKEM.py`).
   With CAMB tables the linear growth depends on $k$:
   $D(k,z)/D(k_0,z) - 1$ for $k = 0.01$ to $0.2\ {\rm Mpc}^{-1}$ and
   $z = 0.5$ to $2$ is $+0.5\%$ to $+1.2\%$ in LSST-Y1 ($w = -0.9$) and
   $+0.03\%$ to $+0.3\%$ in Roman-Real ($w = -1$). Both have
   $m_\nu = 0.06$ eV. Swapping $w$ moves the difference with it: LSST-Y1 at
   $w = -1$ gives 0.34, Roman-Real at $w = -0.9$ gives 7.57.
   The last two terms then do not cancel at the Limber limit:
   - DESC-CCL's $C_{gs}$ below $\ell = 150$ changes by $-0.7\%$ to $-1.5\%$
     (LSST-Y1) and $-0.03\%$ to $-0.4\%$ (Roman-Real) between the CAMB table
     and the separable table;
   - the Limber $C_\ell$ do not change;
   - compared directly, DESC-CCL's non-Limber $C_{gg}$ is below CoCoA's by
     $0.8\%$ to $1.7\%$ in LSST-Y1 at every $\ell$ from 2 to 140 (table
     above). The offset is flat in $\ell$, the signature of two terms that do
     not cancel rather than of a non-Limber effect, which changes with $\ell$.

   CoCoA `v5.02` evaluates both linear terms with the same separable
   spectrum, $(D(z)/D(z_{\rm piv}))^2\,P_{\rm lin}(k, z_{\rm piv})$ with
   $z_{\rm piv}$ the lens bin's mean redshift (`cosmo2D.c`), so the pair
   cancels at the Limber limit; DESC-CCL mixes the two forms. Neither form
   keeps the $k$ dependence of the growth, but CoCoA's anchor limits that
   error to the bin's redshift width instead of the whole range from
   $z = 0$: the largest $|P_{\rm sep}/P_{\rm lin} - 1|$ over
   $k = 0.01$ to $0.2\ {\rm Mpc}^{-1}$ and the central 90% of each lens
   bin's $n(z)$ is

   | project | anchor at $z = 0$ (DESC-CCL's FKEM term) | anchor at the bin's mean redshift (CoCoA) |
   |---|---|---|
   | LSST-Y1 | $1.0\%$ to $1.9\%$ | $0.15\%$ to $0.28\%$ |
   | Roman-Real | $0.16\%$ to $0.96\%$ | $0.05\%$ to $0.21\%$ |

   (`scripts/diag_anchor.py`). Diagnostic: DESC-CCL given the
   separable table lowers the fiducial $\Delta\chi^2$ from 8.76 to 0.21
   (LSST-Y1, $\gamma_t$ 7.29 to 0.064, $w$ 6.26 to 0.16) and from 0.228 to
   0.062 (Roman-Real). The $\Delta\chi^2$ of the non-Limber effect on
   LSST-Y1 $\gamma_t$ (non-Limber vs Limber in the same code) is 1.98 in
   CoCoA, 11.2 in DESC-CCL and 2.06 in DESC-CCL with the separable table.
   Neither code evaluates the exact integral when the growth depends on $k$:
   CoCoA uses one separable form per lens bin, DESC-CCL two forms. The study
   measures their difference, not either code's distance from the exact
   $C_\ell$. In parameters (table above), the LSST-Y1 difference goes into
   the lens biases: 6.74 of the 8.76 is absorbed by lowering every $b_1$ by
   about $1\sigma$, while $\Omega_m$ and $n_s$ move by $0.07\sigma$ and
   $0.11\sigma$.
3. **DESC-CCL numerics.** DESC-CCL's default FKEM sampling (`fkem_Nchi`
   default, 500 log $\ell$) is 0.34 from the reference in LSST-Y1 (in
   $w(\theta)$) and 0.0005 in Roman-Real. The reference moves by 0.0014 and
   0.0002 against a finer run (`fkem_Nchi=4000`, 3000 log $\ell$).
4. **DESC-CCL transform failures.**
   - The pairs: LSST-Y1 lens 5 - source 1 and Roman-Real lens 8 - source 3
     (panel (8, 3) of the Roman-Real $\gamma_t$ figure).
     In both, the lens bin lies behind the source bin, so the Limber
     $C_{gs}$ is zero at every $\ell$. FKEM gives a nonzero $C_{gs}$ below
     $\ell = 150$, and the Limber part above is zero.
   - `correlation` fails on these $C_\ell$: Legendre reports
     `ran out of memory` and FFTLog `failed to create spline`. The failing
     call builds the $C_\ell$ spline with log-log extrapolation beyond
     $\ell_{\max}$ (`ccl_f1d_extrap_logx_logy` in `ccl_correlation.c`),
     which needs the last two values nonzero and of one sign; here they are
     exactly zero.
   - Both pairs are left out of the comparison. CoCoA's mask removes the
     LSST-Y1 pair; in Roman-Real, 13 points the mask keeps are left out.
5. **Modeling choices** (one-change table): RSD is 115 (LSST-Y1, in
   $w(\theta)$); $C_{gs}$ in Limber 9.2; values at bin centers instead of
   bin averages 10.3 (LSST-Y1) and 22.8 (Roman-Real, whose bins are wider:
   $\Delta\ln\theta = 0.31$ against 0.23); DESC-CCL's own Eisenstein-Hu and
   halofit $P(k)$ instead of CAMB's tables 21.3 and 58.8.

## Figures <a name="ccl_figures"></a>

Each panel shows $(d_{\rm CCL} - d_{\rm CoCoA})/\sigma$ at the same
cosmology, with $\sigma = \sqrt{C_{ii}}$ from the project's covariance: 1 is a
one-sigma difference. One curve per cosmology. Each row of panels has its own
y-range. A panel marked $1/\alpha = f$ shows the difference divided by $f$:
the difference is $f$ times what the axis reads. LSST-Y1 shows every
$\theta$; the $\Delta\chi^2$ above uses the mask. Roman-Real leaves the
masked points blank. Panel labels: $\gamma_t$ (lens, source),
$w(\theta)$ (lens), $\xi_\pm$ (bin $i$, bin $j$).
The maxima below are over the five cosmologies and every plotted point.

**LSST-Y1, $\gamma_t$.** DESC-CCL's FKEM offset from finding 2 reaches
$-0.53\sigma$ near 100'. The same shape appears at every cosmology.

![LSST-Y1 gamma_t](cocoa_comparison/figures/lsst_y1_cosmologies_gammat.png)

**LSST-Y1, $w(\theta)$.** The same FKEM offset: $-0.3\sigma$ to $-1.39\sigma$
per point below 100' (at most $0.89\sigma$ on the points the mask keeps).

![LSST-Y1 w](cocoa_comparison/figures/lsst_y1_cosmologies_w.png)

**LSST-Y1, $\xi_+$ and $\xi_-$.** At most $0.014\sigma$ ($\xi_+$) and
$0.010\sigma$ ($\xi_-$).

![LSST-Y1 xi+](cocoa_comparison/figures/lsst_y1_cosmologies_xip.png)
![LSST-Y1 xi-](cocoa_comparison/figures/lsst_y1_cosmologies_xim.png)

**Roman-Real, $\gamma_t$.** At most $0.045\sigma$. Panels marked
"excluded" are CoCoA's $\gamma_t$ exclusions (Method, Layout); panel (8, 3)
is blank because DESC-CCL's transform fails ([finding 4](#ccl_findings)).

![Roman-Real gamma_t](cocoa_comparison/figures/roman_real_cosmologies_gammat.png)

**Roman-Real, $w(\theta)$.** At most $0.088\sigma$.

![Roman-Real w](cocoa_comparison/figures/roman_real_cosmologies_w.png)

**Roman-Real, $\xi_+$ and $\xi_-$.** At most $0.050\sigma$ ($\xi_+$) and
$0.059\sigma$ ($\xi_-$).

![Roman-Real xi+](cocoa_comparison/figures/roman_real_cosmologies_xip.png)
![Roman-Real xi-](cocoa_comparison/figures/roman_real_cosmologies_xim.png)

## TATT intrinsic alignments <a name="ccl_tatt"></a>

The comparison repeated with TATT: CoCoA with `IA_model: 1` and its TATT
terms from CFASTPT (`IA_code: 0`, the C port of FAST-PT in cosmolike);
DESC-CCL with `EulerianPTCalculator` (FAST-PT 4.0.0), variant `tatt:` of
`ccl_compute.py`. TATT point of this repository's benchmark scripts:
$A_1 = 0.7$, $\eta_1 = -1.7$, $A_2 = -1.36$, $\eta_2 = -2.5$, $b_{\rm TA} = 1$,
pivot $z = 0.62$ (CoCoA's $(1+z)/1.62$), at the five cosmologies.

| convention | CoCoA | DESC-CCL setting |
|---|---|---|
| $C_1$ | $A_1(z)\,C_1\rho_{\rm crit}\,\Omega_m/D$, $C_1\rho_{\rm crit} = 0.01389$ (sign in the integrand) | `translate_IA_norm(a1=A1(z) 0.01389/(5e-14 RHO_CRITICAL))`, so $c_1 = -C_1$ |
| $C_{1\delta}$ | $b_{\rm TA} C_1$ | `a1delta = b_TA a1` |
| $C_2$ | $A_2(z)\,C_1\rho_{\rm crit}\,\Omega_m/D^2$, times 5 in the kernels | `a2` with the same factor, `Om_m2_for_c2=False`, so $c_2 = 5 C_2$ |
| one-loop terms | FAST-PT kernels (CFASTPT) of $P_{\rm lin}(k,0)$ times $D^4$ | FAST-PT kernels of $P_{\rm lin}(k,0)$ times $D^4$ |
| tree-level term | $C_1 P_{\rm nl}$ | `b1_pk_kind="nonlinear"` |
| $\xi_\pm$ | $C_{EE}$ (GG + GI + IG + II) $\pm$ $C_{BB}$ (II) | `correlation` of $C_{EE} \pm C_{BB}$ (`return_ia_bb=True`) |
| $\gamma_t$ | non-Limber term with the $C_1$ part only (as NLA); the $b_{\rm TA}$ and $A_2$ terms in Limber | the reference FKEM call, plus Limber $C_{gI}$ from the `m:cdelta` and `m:c2` templates |
| PT tables | 1100 $k$, $1.7\times10^{-5}$ to $334\ h\,{\rm Mpc}^{-1}$ | 160 per decade, $10^{-5}$ to $200\ {\rm Mpc}^{-1}$, FAST-PT windows as CoCoA |

The mapping, with file:line evidence in both codes, is in
`.claude/skills/cocoa-ccl-comparison/references/tatt_conventions.md`. With
$A_2 = b_{\rm TA} = 0$ the `tatt:` variant reproduces the NLA reference to
$\Delta\chi^2 = 8\times10^{-5}$.

DESC-CCL (reference settings) vs CoCoA, TATT:

| model | LSST-Y1 3x2pt | $\xi_\pm$ | $\gamma_t$ | $w(\theta)$ | Roman-Real 3x2pt | $\xi_\pm$ | $\gamma_t$ | $w(\theta)$ |
|---|---|---|---|---|---|---|---|---|
| fiducial | 8.89 | 0.0018 | 7.46 | 6.26 | 0.237 | 0.058 | 0.280 | 0.096 |
| $\Omega_m = 0.25$ | 6.62 | 0.0006 | 4.96 | 5.48 | 0.168 | 0.022 | 0.196 | 0.087 |
| $\Omega_m = 0.35$ | 10.76 | 0.0059 | 9.74 | 6.46 | 0.306 | 0.123 | 0.378 | 0.098 |
| $n_s = 0.92$ | 9.21 | 0.0018 | 7.74 | 6.55 | 0.244 | 0.055 | 0.291 | 0.099 |
| $n_s = 1.01$ | 8.59 | 0.0020 | 7.22 | 5.97 | 0.222 | 0.056 | 0.269 | 0.091 |

| comparison (fiducial, TATT, 3x2pt) | LSST-Y1 | Roman-Real |
|---|---|---|
| TATT vs NLA in CoCoA (size of the TATT terms; $\xi_\pm$ alone) | 1742 (1677) | 1307 (1293) |
| CoCoA default vs CoCoA high accuracy | 0.005 | 0.008 |
| DESC-CCL PT tables at 320 vs 160 $k$ per decade | $8\times10^{-9}$ | $4\times10^{-7}$ |
| both codes in Limber | 0.159 | 0.031 |
| DESC-CCL separable $P_{\rm lin}$ (diagnostic) vs CoCoA | 0.22 | 0.066 |

Limber $C_\ell$ with TATT, $|C_\ell^{\rm CCL}/C_\ell^{\rm CoCoA} - 1|$ (median / maximum, as above):

| project | probe | $20 \le \ell < 66$ | $66 \le \ell < 1000$ | $1000 \le \ell \le 5000$ | $5000 < \ell \le \ell_{\max}$ |
|---|---|---|---|---|---|
| LSST-Y1 | $C_{ss}^{EE}$ | 2.7e-4 / 1.2e-3 | 6.1e-5 / 7.5e-4 | 6.4e-5 / 5.7e-3 | 2.2e-4 / 2.5e-2 |
| LSST-Y1 | $C_{ss}^{BB}$ | 7.1e-4 / 5.5e-3 | 7.0e-4 / 5.4e-3 | 1.4e-3 / 5.2e-3 | 1.4e-3 / 5.1e-3 |
| LSST-Y1 | $C_{gs}$ | 3.6e-5 / 2.5e-3 | 4.5e-5 / 3.3e-3 | 5.4e-5 / 5.1e-3 | 9.8e-5 / 2.6e-2 |
| Roman-Real | $C_{ss}^{EE}$ | 5.3e-5 / 6.2e-4 | 7.3e-5 / 1.2e-3 | 8.4e-5 / 1.8e-3 | 1.7e-4 / 4.9e-1 |
| Roman-Real | $C_{ss}^{BB}$ | 3.3e-5 / 1.2e-4 | 2.2e-5 / 1.1e-4 | 3.8e-5 / 1.2e-4 | 1.2e-4 / 4.3e-3 |
| Roman-Real | $C_{gs}$ | 3.1e-5 / 1.7e-3 | 5.4e-5 / 1.9e-3 | 6.8e-5 / 4.9e-3 | 1.2e-4 / 1.2e-1 |

1. **The TATT terms agree.** $\xi_\pm$ differ by 0.0018 (LSST-Y1) and
   0.058 (Roman-Real), the NLA values, while TATT moves $\xi_\pm$ by 1677
   and 1293 from NLA; both codes in Limber give 0.159 and 0.031, as with
   NLA. Both codes are converged (CoCoA high accuracy, DESC-CCL PT tables).
2. **The rest is finding 2.** The TATT $b_{\rm TA}$ and $A_2$ terms enter
   $\gamma_t$ in Limber in both codes; the fiducial 8.89 and 0.237 are the
   NLA FKEM offset, and the separable diagnostic lowers them to 0.22 and
   0.066.
3. **$k$ support above $\ell \approx 5\times10^4$.** CoCoA's TATT kernels end at
   $k = 334\ h\,{\rm Mpc}^{-1}$; DESC-CCL extrapolates the PT spectra. In
   pairs with the first source bin, where TATT dominates $C_{EE}$ (Roman-Real
   sources 1-4: TATT is $-2$ times NLA), the two differ by 8% at
   $\ell = 6.9\times10^4$ and 49% at $\ell = 10^5$ (2.5% at $\ell = 6.5\times10^4$
   in LSST-Y1). The real-space $\xi_\pm$ difference equals the NLA one, so
   these multipoles do not reach the data vector above 2.5'.
4. **The benchmark scripts.** This repository's DESC-CCL scripts set up TATT
   as above since commit `b8bfec7`: IA-only tracer with `use_A_ia=False`
   (with `use_A_ia=True`, as before, the IA normalization is applied a
   second time and its sign flips), $b_{\rm TA}$ as a parameter, and $C_{BB}$
   in $\xi_\pm$. The timings in the note above were measured before that
   commit.

Figures, as above ($(d_{\rm CCL} - d_{\rm CoCoA})/\sigma$, five
cosmologies). LSST-Y1 $\xi_+$ (at most $0.020\sigma$; $\xi_-$ $0.010\sigma$)
and $\gamma_t$ (the FKEM offset, at most $0.53\sigma$; $w(\theta)$ is the NLA
figure):

![LSST-Y1 xi+ TATT](cocoa_comparison/figures/lsst_y1_tatt_cosmologies_xip.png)
![LSST-Y1 gamma_t TATT](cocoa_comparison/figures/lsst_y1_tatt_cosmologies_gammat.png)

Roman-Real $\xi_+$ (at most $0.049\sigma$; $\xi_-$ $0.059\sigma$) and
$\gamma_t$ (at most $0.045\sigma$):

![Roman-Real xi+ TATT](cocoa_comparison/figures/roman_real_tatt_cosmologies_xip.png)
![Roman-Real gamma_t TATT](cocoa_comparison/figures/roman_real_tatt_cosmologies_gammat.png)

## What it took to run DESC-CCL <a name="ccl_run"></a>

| symptom | cause | resolution |
|---|---|---|
| `CCLError` 1040 for full-sky $\xi_\pm$ | CCL 3.3 has no full-sky $\xi_\pm$ | pyccl built from PR #1296 (CMake against the conda env's GSL and FFTW; SWIG from conda-forge) |
| `AttributeError: np.trapz` in FAST-PT 4.0.0 | numpy 2.4 removed `np.trapz` | `np.trapz = np.trapezoid` before importing pyccl (no package changed) |
| `ccl_angular_cls_limber(): integration error`, Roman-Real lens bins 7 and 8 with RSD | the growth table of `cocoa_export.py` ended at $z = 6$ | growth from CAMB to $z = 49$ |
| same error, lens bin 8, $\ell = 2$ only | the Limber RSD kernel evaluates the background at $\chi_{\ell+1} = \chi\,(\ell + 3/2)/(\ell + 1/2)$, $1.4\chi$ at $\ell = 2$ ($z \approx 14$ for $z = 4$); `ccl_compute.py` cut CoCoA's $\chi(z)$ table at $z = 10$ | CoCoA's $\chi(z)$ to $z = 50$ |
| $H(z)$ off by $3\times10^{-4}$ above $z = 3$ | `np.gradient` of $\chi(z)$ where the $z$ grid coarsens | cubic-spline derivative (agrees with CCL's own $H(z)$ to $10^{-5}$) |
| Eisenstein-Hu cosmology needs $\sigma_8$ | the benchmark modeling uses DESC-CCL's own $P(k)$ | `ccl.sigma8` of the CAMB-table cosmology |
| `correlation` fails for two $\gamma_t$ pairs | finding 4 | not worked around: pairs left out |

CoCoA side:
- cosmolike reads the data files once per process, so the full vector (a
  scratch dataset with `ones.mask`) and the masked covariance come from
  separate runs;
- the Roman-Real covariance is not positive definite with every point, so
  Roman-Real uses its mask;
- Cobaya returns $\chi(z)$ only at requested redshifts, so it is read on the
  likelihood's $z$ grid.

## Reproduce <a name="ccl_reproduce"></a>

We assume users run the commands from the CCL-benchmark root folder, have a
Cocoa installation with the lsst_y1 and roman_real projects and their data,
a conda environment `ccl` with pyccl's dependencies, and a pyccl build of
[LSSTDESC/CCL PR #1296](https://github.com/LSSTDESC/CCL/pull/1296) (the
folder that contains its `pyccl` package). `<run folder>` is the same path in
every step.

**Step :one:**: Export CoCoA's data vectors, tables and covariances (Conda
cocoa environment, from the `Cocoa/` folder: `source start_cocoa.sh`
first, then come back to the CCL-benchmark root)

```bash
bash cocoa_comparison/run_cocoa.sh <run folder>
```

**Step :two:**: Compute the DESC-CCL data vectors (`CCL_PYTHON` is the python
of the `ccl` environment; `CCL_PR` the PR #1296 build folder; `COCOA_ROOTDIR`
the `Cocoa/` folder)

```bash
CCL_PYTHON=<ccl env python> CCL_PR=<PR build folder> COCOA_ROOTDIR=<Cocoa folder> bash cocoa_comparison/run_ccl.sh <run folder>
```

**Step :three:**: Write the tables (any python with numpy)

```bash
CMP_WORK=<run folder> python cocoa_comparison/scripts/results.py
CMP_WORK=<run folder> python cocoa_comparison/scripts/harmonic.py
```

**Step :four:**: Compute the parameter shifts (any python with numpy)

```bash
CMP_WORK=<run folder> python cocoa_comparison/scripts/bias.py
```

**Step :five:**: Draw the figures (Conda cocoa environment with
`start_cocoa.sh` sourced: `plots.py` uses cosmolike_core's
`plot_datavectors.py`)

```bash
CMP_WORK=<run folder> python cocoa_comparison/scripts/plots.py cocoa_comparison/figures
```

The diagnostics behind findings 2 and 4 are `scripts/diag_fkem_offset.py`,
`scripts/diag_growth.py`, `scripts/diag_anchor.py` and `scripts/diag_l7s2.py`
(usage in each file).

Versions of this study's runs (4 OpenMP threads for every run):

| component | version |
|---|---|
| Cocoa | `v5.02`, cosmolike_core `v5.03` + 5 commits (`26191ec`): nested $N_a$ above `accuracyboost: 1` (enters only the high-accuracy row), a testing fix, plotting (per-row `ylim`, bin labels from 1) |
| CAMB (Cocoa environment) | 1.6.7, Python 3.11.12, NumPy 1.26.3 |
| pyccl | PR #1296 build, commit `647ad4a` (2026-09-18), first on `PYTHONPATH` |
| `ccl` conda environment | Python 3.12.13, NumPy 2.4.3, SciPy 1.17.1, FAST-PT 4.0.0, GSL 2.7, FFTW 3.3.11, SWIG 4.5.1 |
