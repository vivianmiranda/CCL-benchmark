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
6. [What it took to run DESC-CCL](#ccl_run)
7. [Reproduce](#ccl_reproduce)

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
    (Roman-Real: lens-source pairs (6,0), (7,0) and (7,1), 0-indexed).
  - Masks: LSST-Y1 `lsst_y1_M1_GGLOLAP0.05.mask` (959 of 1560 points);
    Roman-Real `example1.mask` (1950 of 2115 points).
- **Models.**
  - Likelihood: each project's `EXAMPLE_EVALUATE2.yaml` (NLA, linear bias),
    with photo-z shifts, shear calibration, magnification and point masses
    set to zero.
  - Fiducial: $A_s = 2.1\times10^{-9}$, $n_s = 0.96605$, $H_0 = 67.32$,
    $\Omega_m = 0.3$, $\Omega_b = 0.04$, $m_\nu = 0.06$ eV;
    $w = -0.9$ (LSST-Y1) and $w = -1$ (Roman-Real).
  - Variations: $\Omega_m = 0.25, 0.35$ and $n_s = 0.92, 1.01$, at fixed $A_s$.
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
makes them not add up to the 3x2pt value. In Roman-Real, the pair lens 7 -
source 2 (0-indexed) is left out: 13 points the mask keeps
([finding 4](#ccl_findings)).

Fiducial cosmology, one change at a time (3x2pt $\Delta\chi^2$):

| comparison | LSST-Y1 | Roman-Real |
|---|---|---|
| CoCoA default vs CoCoA high accuracy (`accuracyboost: 2`, `integration_accuracy: 1`) | 0.004 | 0.007 |
| DESC-CCL default FKEM sampling, 500 log $\ell$ vs reference | 0.34 | 0.0005 |
| both codes in Limber | 0.16 | 0.031 |
| DESC-CCL separable $P_{\rm lin}$ (diagnostic) vs CoCoA | 0.21 | 0.062 |
| DESC-CCL, RSD also in $\gamma_t$, vs CoCoA | 8.70 | 0.258 |
| DESC-CCL $C_{gs}$ in Limber vs CoCoA | 9.23 | 0.71 |
| DESC-CCL without RSD vs CoCoA | 115.4 | 6.98 |
| DESC-CCL full-sky at bin centers (no bin averaging) vs CoCoA | 10.3 ($\xi_\pm$: 2.0) | 22.8 ($\xi_\pm$: 28.3) |
| DESC-CCL flat-sky FFTLog at bin centers vs CoCoA | 12.7 ($\xi_\pm$: 1.8) | 22.6 ($\xi_\pm$: 19.4) |
| DESC-CCL with this repository's benchmark modeling vs CoCoA | 153.6 | 55.9 |

The benchmark-modeling row applies the modeling choices of
`ccl_test_lsst.py` and `ccl_test_roman.py` to CoCoA's layout:
- DESC-CCL's own Eisenstein-Hu linear and halofit nonlinear $P(k)$
  ($\sigma_8$ from CoCoA's linear table);
- $C_{gs}$ in Limber, no RSD, `l_limber=100`, `fkem_Nchi=500`;
- flat-sky FFTLog at the bin centers.

The intrinsic alignment stays NLA; the scripts use TATT.

Limber $C_\ell$ at 24 $\ell$ from 20 to 5000, $|C_\ell^{\rm CCL}/C_\ell^{\rm CoCoA} - 1|$
(median over bin pairs / maximum over pairs with $|C_\ell|$ above 1% of the
largest):

| project | probe | $20 \le \ell < 66$ | $66 \le \ell < 1000$ | $1000 \le \ell \le 5000$ |
|---|---|---|---|---|
| LSST-Y1 | $C_{ss}$ | 2.7e-4 / 1.7e-3 | 8.4e-5 / 1.1e-3 | 7.4e-5 / 2.2e-3 |
| LSST-Y1 | $C_{gs}$ | 2.9e-4 / 2.2e-3 | 3.4e-4 / 3.8e-3 | 3.2e-4 / 5.7e-3 |
| LSST-Y1 | $C_{gg}$ | 5.4e-4 / 2.2e-3 | 3.0e-5 / 2.0e-4 | 2.5e-5 / 5.1e-5 |
| Roman-Real | $C_{ss}$ | 5.5e-5 / 6.5e-4 | 7.3e-5 / 5.8e-4 | 1.2e-4 / 1.0e-3 |
| Roman-Real | $C_{gs}$ | 3.0e-4 / 1.9e-3 | 3.6e-4 / 1.9e-3 | 3.3e-4 / 1.8e-3 |
| Roman-Real | $C_{gg}$ | 9.3e-4 / 4.3e-3 | 4.2e-5 / 1.0e-3 | 2.7e-5 / 3.0e-3 |

## Findings <a name="ccl_findings"></a>

1. **Limber layer.** The codes agree in harmonic space (table above) and in
   real space with both codes in Limber ($\Delta\chi^2$ 0.16 and 0.031).
   CoCoA at high accuracy moves by 0.004 and 0.007.
2. **Non-Limber $C_{gs}$ and $C_{gg}$ (DESC-CCL FKEM).** DESC-CCL computes
   $$C_\ell = C_\ell^{\rm Limber}[P_{\rm nl}] - C_\ell^{\rm Limber}[P_{\rm lin}(k,z)] + C_\ell^{\rm FKEM}[D^2(z)\,P_{\rm lin}(k,0)].$$
   The FKEM term multiplies the kernels by DESC-CCL's growth factor and uses
   the linear table at $z = 0$ (`pyccl/nonlimber/_nonlimber_FKEM.py`).
   With CAMB tables the linear growth depends on $k$:
   $D(k,z)/D(k_0,z) - 1$ for $k = 0.01$ to $0.2\ {\rm Mpc}^{-1}$ and
   $z = 0.5$ to $2$ is $+0.5\%$ to $+1.2\%$ in LSST-Y1 ($w = -0.9$) and
   $+0.03\%$ to $+0.3\%$ in Roman-Real ($w = -1$). Both have
   $m_\nu = 0.06$ eV; $w$ accounts for the difference between the projects.
   The last two terms then do not cancel at the Limber limit:
   - DESC-CCL's $C_{gs}$ below $\ell = 150$ changes by $-0.7\%$ to $-1.5\%$
     (LSST-Y1) and $-0.03\%$ to $-0.4\%$ (Roman-Real) between the CAMB table
     and the separable table;
   - the Limber $C_\ell$ do not change.

   CoCoA `v5.02` evaluates both linear terms with the same separable
   spectrum, $(D(z)/D(z_{\rm piv}))^2\,P_{\rm lin}(k, z_{\rm piv})$ with
   $z_{\rm piv}$ the lens bin's mean redshift (`cosmo2D.c`), so the pair
   cancels at the Limber limit; DESC-CCL mixes the two forms. Diagnostic: DESC-CCL given the
   separable table lowers the fiducial $\Delta\chi^2$ from 8.76 to 0.21
   (LSST-Y1, $\gamma_t$ 7.29 to 0.064, $w$ 6.26 to 0.16) and from 0.228 to
   0.062 (Roman-Real). The $\Delta\chi^2$ of the non-Limber effect on
   LSST-Y1 $\gamma_t$ (non-Limber vs Limber in the same code) is 1.98 in
   CoCoA, 11.2 in DESC-CCL and 2.06 in DESC-CCL with the separable table.
3. **DESC-CCL numerics.** DESC-CCL's default FKEM sampling (`fkem_Nchi`
   default, 500 log $\ell$) is 0.34 from the converged run in LSST-Y1 (in
   $w(\theta)$) and 0.0005 in Roman-Real.
4. **DESC-CCL transform failures.**
   - The pairs (0-indexed): LSST-Y1 lens 4 - source 0 and Roman-Real lens
     7 - source 2 (panel (8, 3) of the Roman-Real $\gamma_t$ figure).
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
   $\Delta\ln\theta = 0.31$ against 0.23).

## Figures <a name="ccl_figures"></a>

Each panel shows $(d_{\rm CCL} - d_{\rm CoCoA})/\sigma$ at the same
cosmology, with $\sigma = \sqrt{C_{ii}}$ from the project's covariance: 1 is a
one-sigma difference. One curve per cosmology. Each row of panels has its own
y-range. A panel marked $1/\alpha = f$ shows the difference divided by $f$:
the difference is $f$ times what the axis reads. LSST-Y1 shows every
$\theta$; the $\Delta\chi^2$ above uses the mask. Roman-Real leaves the
masked points blank. Panel labels: $\gamma_t$ (lens, source) and
$w(\theta)$ (lens) count from 1; $\xi_\pm$ (bin $i$, bin $j$) count from 0.
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

## What it took to run DESC-CCL <a name="ccl_run"></a>

| symptom | cause | resolution |
|---|---|---|
| `CCLError` 1040 for full-sky $\xi_\pm$ | CCL 3.3 has no full-sky $\xi_\pm$ | pyccl built from PR #1296 (CMake against the conda env's GSL and FFTW; SWIG from conda-forge) |
| `AttributeError: np.trapz` in FAST-PT 4.0.0 | numpy 2.4 removed `np.trapz` | `np.trapz = np.trapezoid` before importing pyccl (no package changed) |
| `ccl_angular_cls_limber(): integration error`, Roman-Real lens bins 6 and 7 with RSD | the growth table of `cocoa_export.py` ended at $z = 6$ | growth from CAMB to $z = 49$ |
| same error, lens bin 7, $\ell = 2$ only | the Limber RSD kernel evaluates the background at $\chi_{\ell+1} = \chi\,(\ell + 3/2)/(\ell + 1/2)$, $1.4\chi$ at $\ell = 2$ ($z \approx 14$ for $z = 4$); `ccl_compute.py` cut CoCoA's $\chi(z)$ table at $z = 10$ | CoCoA's $\chi(z)$ to $z = 50$ |
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

**Step :four:**: Draw the figures (Conda cocoa environment with
`start_cocoa.sh` sourced: `plots.py` uses cosmolike_core's
`plot_datavectors.py`)

```bash
CMP_WORK=<run folder> python cocoa_comparison/scripts/plots.py cocoa_comparison/figures
```

The diagnostics behind findings 2 and 4 are `scripts/diag_fkem_offset.py`,
`scripts/diag_growth.py` and `scripts/diag_l7s2.py` (usage in each file).
