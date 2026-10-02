<!-- Fable 5 convention map, 2026-10-02: CoCoA (cosmolike) TATT vs pyccl (PR #1296 647ad4a)
     TATT. Saved so later sessions reuse it. Line numbers: cosmolike_core v5.03 + 5,
     pyccl PR build 647ad4a. Note (checked in session): both projects' EXAMPLE_EVALUATE2.yaml
     set IA_code: 0 (CFASTPT); the comparison sets IA_code: 0 explicitly. -->

# CoCoA (cosmolike) TATT <-> DESC-CCL (pyccl PR #1296) convention map

Purpose: configure CCL to reproduce CoCoA's real-space 3x2pt data vector with
TATT intrinsic alignments on (lsst_y1, roman_real), as a "tatt" variant of
`CCL-benchmark/cocoa_comparison/scripts/ccl_compute.py`. All CoCoA paths are
under `/Users/vivianmiranda/data/COCOA/september2026/test/cocoa/Cocoa/`
(cosmolike = `external_modules/code/cosmolike_core/cosmolike/`); all pyccl
paths under the PR #1296 build
`/private/tmp/claude-501/.../scratchpad/ccl_pr/pyccl/`.

---

## 1. How CoCoA selects TATT, and its parameters

Options (likelihood yaml, e.g. `projects/lsst_y1/likelihood/combo_3x2pt.yaml:36-45`):

- `IA_model`: 0 = NLA, 1 = TATT (yaml line 42; both projects ship with 0, the
  TATT comparison run sets `IA_model: 1`). Mapped to `nuisance.IA_MODEL`
  (`IA_MODEL_TATT = 1`, `IA.h:7-8`) by `init_IA_fastpt`
  (`generic_interface.cpp:1529-1545`), called from the likelihood at
  `_cosmolike_prototype_base.py:212-214`.
- `IA_redshift_evolution`: 3 = `IA_REDSHIFT_EVOLUTION` (`IA.h:10-13`,
  yaml line 40). Both projects use 3: ONE power law in z shared by all source
  bins, not per-bin amplitudes.
- `IA_code`: 0 = C cfastpt, 1 = python FAST-PT (`fastpt` package) via
  `set_IA_PS` (`combo_3x2pt.yaml:44-45`; lsst_y1 ships 1, roman_real ships 0).
  Both fill the same `FPTIA.tab` rows; the python path
  (`external_modules/code/PyFAST-PT/fastpt.py:184-195,229-234`) literally
  calls the same `fastpt.FASTPT.IA_tt/IA_ta/IA_mix` that CCL's
  `EulerianPTCalculator` calls (`ept.py:255-260`), so the kernel set is
  identical by construction.

Sampled parameters (`_cosmolike_prototype_base.py:619-623`, prefix `LSST_`
for lsst_y1, `roman_` for roman_real; slot mapping
`generic_interface.cpp:3531-3547`, mode `IA_REDSHIFT_EVOLUTION`):

| parameter      | meaning                               | slot            |
|----------------|---------------------------------------|-----------------|
| `LSST_A1_1`    | A1 (tidal-alignment amplitude)        | `nuisance.ia[0][0]` |
| `LSST_A1_2`    | eta_1 (A1 redshift power-law index)   | `nuisance.ia[0][1]` |
| `LSST_A2_1`    | A2 (tidal-torquing amplitude)         | `nuisance.ia[1][0]` |
| `LSST_A2_2`    | eta_2 (A2 redshift power-law index)   | `nuisance.ia[1][1]` |
| `LSST_BTA_1`   | b_TA (density weighting), z-constant  | `nuisance.ia[2][0]` |

Pivot: `nuisance.oneplusz0_ia = 1.62` (z0 = 0.62,
`generic_interface.cpp:3533`). `nuisance.c1rhocrit_ia = 0.01389`
(`generic_interface.cpp:3512`, default also `structs.c:223`).

Redshift dependence (`IA.c`):
- A1(z) = A1 * ((1+z)/1.62)^eta_1            (`IA.c:278-284`)
- A2(z) = A2 * ((1+z)/1.62)^eta_2            (`IA.c:342-349`)
- b_TA(z) = b_TA (constant)                  (`IA.c:398-402`)

(With `IA_REDSHIFT_BINNING` = 2 the same slots hold per-bin A1_i, A2_i,
b_TA_i, `generic_interface.cpp:3514-3529`; the comparison projects do not use
it, and the recipe below assumes the evolution mode — per-bin amplitudes would
need per-pair Pk2Ds in CCL.)

---

## 2. C1(z), C2(z), C_delta(z) and the exact pyccl arguments

### CoCoA amplitudes

- `IA_A1_Z1Z2` (`IA.c:238-296`): returns
  **C1_cocoa(z) = A1(z) * Omega_m * 0.01389 / D(z)**, POSITIVE sign
  (`IA.c:293-295`; the explicit warning `IA.c:246-247` says the minus sign of
  original cosmolike's C1_TA was dropped — the minus is instead carried by the
  `-WS*IA` structure of the integrand cores, see Sec. 3).
- `IA_A2_Z1Z2` (`IA.c:305-361`): returns
  **C2_cocoa(z) = A2(z) * Omega_m * 0.01389 / D(z)^2** WITHOUT the Blazek
  et al. 2019 factor 5 (`IA.c:312-313, 358-360`). The 5 is applied in the
  integrand cores, one per power of the quadratic field: 5*C2 in terms linear
  in C2, 25*C2^2 in the tt term (`cosmo2D.c:2379-2385`). So the effective
  Blazek amplitude is C2_eff = 5 * A2(z) * Omega_m * 0.01389 / D(z)^2.
- `IA_BTA_Z1Z2` (`IA.c:370-415`): returns b_TA; the density-weighted
  amplitude is C1delta = b_TA * C1 (it always multiplies C1 in the cores).

D(z) is `growfac(a)`: CAMB's growth table G(z) with D = G*a / G(0), i.e.
**normalized to D(z=0) = 1** (`cosmo3D.c:512-525, 568-603`), fed from the
same CAMB run as the P(k) tables. Omega_m is total matter (including
neutrinos) as CoCoA defines `cosmology.Omega_m`.

### pyccl mapping (PR #1296)

`translate_IA_norm` (`nl_pt/tracers.py:41-58`), with
`gz = cosmo.growth_factor(1/(1+z))`:

- `c1      = -1 * a1      * 5e-14 * RHO_CRITICAL * Om_m / gz`
- `c1delta = -1 * a1delta * 5e-14 * RHO_CRITICAL * Om_m / gz`
- `c2      = +a2 * 5 * 5e-14 * RHO_CRITICAL * Om_m / gz^2`   (DES convention,
  `Om_m2_for_c2=False`, line 55-56; `Om_m2_for_c2=True` would use
  Om_m^2/Om_m_fid — do NOT use it, CoCoA scales with one power of Omega_m)

`5e-14 * RHO_CRITICAL = 0.0138775...`, not 0.01389. Reuse the NLA harness's
ratio (`ccl_compute.py:151`):

```python
c1_ratio = 0.01389/(5e-14*ccl.physical_constants.RHO_CRITICAL)
A1z = c1_ratio * A1 * ((1+z)/1.62)**eta1
A2z = c1_ratio * A2 * ((1+z)/1.62)**eta2
c1, cdelta, c2 = pt.translate_IA_norm(cosmo, z=z, a1=A1z,
                                      a1delta=A1z*bta, a2=A2z,
                                      Om_m2_for_c2=False)
```

This gives exactly (with the CosmologyCalculator cosmology, so
`cosmo.growth_factor` is CoCoA's CAMB D and `cosmo['Omega_m']` the same
Omega_m):

- `c1(z)     = -C1_cocoa(z)`
- `cdelta(z) = -b_TA * C1_cocoa(z)`
- `c2(z)     = +5 * C2_cocoa(z)`  ( = C2_eff )

These are the `c1=(z,c1)`, `c2=(z,c2)`, `cdelta=(z,cdelta)` arguments of
`PTIntrinsicAlignmentTracer` (`nl_pt/tracers.py:196-224`). With this mapping
every term of CCL's EPT combinations equals CoCoA's integrand term for term
(verified sign-by-sign below).

### Kernel/normalization equivalence (cfastpt vs CCL's FAST-PT)

CoCoA `get_FPT_IA` (`pt_cfastpt.c:352-1066`) computes, from the z=0 LINEAR
power `p_lin(k, 1.0)` (`pt_cfastpt.c:415`), the ten spectra
(`pt_cfastpt.c:427-435, 806-818`); CCL `EulerianPTCalculator` computes the
same from `cosmo.linear_matter_power(k_s, 1.0)` (`ept.py:230, 255-260`):

| FPTIA.tab | cosmolike name   | fastpt / CCL name (`ept.py:363-364,440-442`) |
|-----------|------------------|-----------------------------------------------|
| tab[0]    | tt E  (`tt`)     | `ia_tt[0]` = ae2e2                            |
| tab[1]    | tt B  (`ttbb`)   | `ia_tt[1]` = ab2b2                            |
| tab[2]    | ta deltaE1       | `ia_ta[0]` = a00e                             |
| tab[3]    | ta deltaE2       | `ia_ta[1]` = c00e                             |
| tab[4]    | ta 0E0E (`ta`)   | `ia_ta[2]` = a0e0e                            |
| tab[5]    | ta 0B0B (`tabb`) | `ia_ta[3]` = a0b0b                            |
| tab[6]    | mix A            | `ia_mix[0]` = a0e2                            |
| tab[7]    | mix B (x4 folded in, `pt_cfastpt.c:1037-1043`) | `ia_mix[1]` = b0e2 |
| tab[8]    | mix D_EE         | `ia_mix[2]` = d0ee2                           |
| tab[9]    | mix D_BB         | `ia_mix[3]` = d0bb2                           |

The C tables are a port of the python `fastpt` outputs (the `IA_code = 1`
path installs the python outputs into the same rows,
`PyFAST-PT/fastpt.py:229-234`, `_cosmolike_prototype_base.py:497-506`), so
cfastpt vs CCL's FAST-PT is a numerics difference (grids/windows), not a
convention difference.

Growth scaling: identical scheme. CoCoA multiplies each kernel by
g4 = D(a)^4 at the Limber node (`cosmo2D.c:2675, 2708`); CCL multiplies by
`self._g4 = growth_factor(a_arr)**4` (`ept.py:231-233` and the `_g4` factors
in `_get_pim/_get_pii/_get_pgi`). Tree-level C1 terms: CoCoA uses the full
NONLINEAR P_delta(k,a) (`PK` in the cores); CCL `b1_pk_kind='nonlinear'`
(default, `ept.py:141, 264-266, 279`) uses `nonlin_matter_power` — with the
CosmologyCalculator fed CoCoA's CAMB halofit table, the same object.

### Term-by-term sign verification

CoCoA EE core (`cosmo2D.c:2411-2441`):

```
EE = WK1*WK2*PK
   - WS1*WK2*(C11*PK + C11*bta1*(ta_dE1+ta_dE2) - 5*C21*(mixA+mixB))
   - WS2*WK1*(C12*PK + C12*bta2*(ta_dE1+ta_dE2) - 5*C22*(mixA+mixB))
   + WS1*WS2*(C11*C12*PK
              + C11*C12*(bta1*bta2*ta + (bta1+bta2)*(ta_dE1+ta_dE2))
              - 5*(C11*C22 + C12*C21)*(mixA+mixB)
              - 5*(C11*bta1*C22 + C12*bta2*C21)*mixEE
              + 25*C21*C22*tt)
```

CCL (`ept.py:_get_pim` 490-493, `_get_pii` 457-462, `_get_pgi` 381-384):

```
pim    = c1*Pd1d1 + g4*cd*(a00e+c00e) + g4*c2*(a0e2+b0e2)
pii_ee = c1c1*Pd1d1 + (c1*cd' + c1'*cd)*g4*(a00e+c00e) + cd*cd'*g4*a0e0e
         + c2*c2'*g4*ae2e2 + (c1*c2' + c2*c1')*g4*(a0e2+b0e2)
         + (cd*c2' + cd'*c2)*g4*d0ee2
pii_bb = cd*cd'*g4*a0b0b + c2*c2'*g4*ab2b2 + (cd*c2' + c2*cd')*g4*d0bb2
```

With c1 = -C1, cd = -bta*C1, c2 = +5*C2 every single term matches,
including the CoCoA BB core (`cosmo2D.c:2476-2497`):
`BB = WS1*WS2*(C11*C12*bta1*bta2*ta_BB - 5*(C11*bta1*C22+C12*bta2*C21)*mix_BB
+ 25*C21*C22*tt_BB)` = `pii_bb` with the mapping. The CoCoA GI bracket (the
quantity multiplying -WS*WK) equals `pim` exactly. CCL's IA tracer
(`WeakLensingTracer(has_shear=False, ia_bias=(z, ones), use_A_ia=False)`,
`tracers.py:916-988`) has a POSITIVE density kernel with unit transfer
(`tracers.py:979-987`), so `angular_cl(t_wl, t_ia, p_of_k_a=pim)` lands the
GI term with CoCoA's sign: negative for A1 > 0. Do NOT also use
`use_A_ia=True` on the IA leg when the amplitudes are in the PT spectra
(that is the bug in the repo's `CCL-benchmark/ccl_test_lsst.py:117-122`).

---

## 3. Which terms CoCoA includes per probe, and the matching CCL recipe

### xi_+/- (shear-shear) — Limber at ALL multipoles (no non-Limber ss path)

Per pair (`cosmo2D.c:2411-2497, 2753-2829`): EE = GG + GI + IG + II_EE (one
nonlinear-P tree level + one-loop TATT kernels as above); BB = II_BB only (no
tree-level BB, `cosmo2D.c:2447-2452`). Real space
(`legendre_sums_xipm`, `cosmo2D.c:646-710`):

```
xi+ = sum_l Glp[theta][l] * (Cl_EE[l] + Cl_BB[l])
xi- = sum_l Glm[theta][l] * (Cl_EE[l] - Cl_BB[l])
```

CCL recipe (tracers per source bin i: `t_wl[i]` shear-only; `t_ia[i]` =
`WeakLensingTracer(has_shear=False, ia_bias=(z, ones), use_A_ia=False)`):

```
Cl_EE(i,j) = angular_cl(t_wl[i], t_wl[j], p_of_k_a='delta_matter:delta_matter')  # CAMB NL table
           + angular_cl(t_wl[i], t_ia[j], p_of_k_a=pk_mi)
           + angular_cl(t_ia[i], t_wl[j], p_of_k_a=pk_mi)
           + angular_cl(t_ia[i], t_ia[j], p_of_k_a=pk_ii)
Cl_BB(i,j) = angular_cl(t_ia[i], t_ia[j], p_of_k_a=pk_ii_bb)
xip = ccl.correlation(..., C_ell=Cl_EE + Cl_BB, type='GG+', ...)
xim = ccl.correlation(..., C_ell=Cl_EE - Cl_BB, type='GG-', ...)
```

with `pk_mi = ptc.get_biased_pk2d(PTMatterTracer(), tracer2=ia_ptt)`,
`pk_ii = ptc.get_biased_pk2d(ia_ptt, tracer2=ia_ptt)`,
`pk_ii_bb = ptc.get_biased_pk2d(ia_ptt, tracer2=ia_ptt, return_ia_bb=True)`
(`ept.py:511-589`). CCL's 'GG+'/'GG-' correlation takes a single C_ell; the
EE+BB / EE-BB sums are exactly CoCoA's assembly, so B-modes enter through the
input C_ell, not through a CCL option. Keep the bin-averaged Legendre method
of the ref variant. Note `get_biased_pk2d(m, ia)` and `(ia, m)` both route to
`_get_pim` (`ept.py:567-581`) — symmetric, as in CoCoA.

### gamma_t (galaxy-galaxy lensing)

Limber integrand core (`cosmo2D.c:3825-3859`), with IA = C1*PK +
C1*BTA*(ta_dE1+ta_dE2) - 5*C2*(mixA+mixB) (the same bracket as pim):

```
ft = WGAL*b1*(WK*PK - WS*IA) + WGAL*oneloop*(WK - WS*C1)
st = WRSD*(WK*PK - WS*IA)                   (include_RSD_GS = 0 in both projects)
tt = WMAG*ell_pf*bmag*(WK*PK - WS*IA)
```

So the gI term is linear galaxy bias times the FULL P_{m,IA} (`pim`), while
one-loop galaxy bias (b2 etc.) multiplies only the LINEAR IA (C1) — a
two-loop-consistency choice (`cosmo2D.c:3790-3796`). At linear bias (the
comparison configuration, B2 = 0) this distinction vanishes and CCL's
`_get_pgi` (b1 times pim, `ept.py:381-384`) is the same model.

NON-LIMBER gamma_t (`C_gs_tomo`, `cosmo2D.c:9428-9496, 9567-9944`): FKEM
split `C_l = C^fftlog(P_lin) + C^limber(full model) - C^limber(linear
model)`. The FFTLog (exact) term's source kernel is
`fx_src = ((W_kappa - W_source*C1)/fK)*D` (`cosmo2D.c:9466, 9742`): **the
FFTLog term carries ONLY the linear C1 (NLA-like) IA**; the linear Limber
subtraction term likewise has "IA through C1 only, no one-loop bias"
(`cosmo2D.c:3973-3976, 4002-4004`). The TATT extras (bta, C2 one-loop
kernels) enter ONLY through the Limber full-model term — i.e. additively, in
Limber, at every ell (the gs integrand is linear in the IA extras). The
FFTLog linear spectrum is anchored per lens bin at a_piv = 1/(1+zmean(lens))
(`cosmo2D.c:9820-9827`), P_lin(k, a_piv)/D(a_piv)^2 with D(chi) on the
kernels.

CCL recipe that reproduces this exactly (linear bias, bmag = 0):

```
# (a) NLA-equivalent composite source tracer: the SAME construction the
#     validated NLA harness uses (ccl_compute.py:150-153)
src_nla[j] = WeakLensingTracer(cosmo, dndz=(zs, nzs[j]),
                               ia_bias=(z, c1_ratio*A1*((1+z)/1.62)**eta1),
                               use_A_ia=True)
# (b) non-Limber part, identical to the ref variant's gamma_t call:
Cl_gs = angular_cl(lens_t[l], src_nla[j], ell, l_limber=150,
                   non_limber_integration_method='FKEM', fkem_Nchi=2000)
# (c) the TATT extras, pure Limber, at every ell:
pk_mi_extra = Pk2D(a_arr=ptc.a_s, lk_arr=np.log(ptc.k_s), is_logp=False,
                   pk_arr=(ptc._g4*cd)[:,None]*(ptc.ia_ta[0]+ptc.ia_ta[1])
                        + (ptc._g4*c2)[:,None]*(ptc.ia_mix[0]+ptc.ia_mix[1]))
Cl_gs += angular_cl(lens_t[l], t_ia[j], ell, p_of_k_a=pk_mi_extra)   # Limber
gt = ccl.correlation(..., C_ell=Cl_gs, type='NG', ...)
```

(`pk_mi_extra` = pim minus its `c1*Pd1d1` term; equivalently
`cd(z)*(m:cdelta template) + c2(z)*(m:c2 template)`, cf.
`ept.py:653-656`.) This is exactly CoCoA's model: FKEM handles gG + gI(C1)
non-Limber (CCL's FKEM carries a z-dependent, k-independent tracer transfer
such as the IA amplitude exactly, via the transfer_low/transfer_avg ratio,
`nonlimber/_nonlimber_FKEM.py:111-156`), and the TATT extras ride on top in
Limber — the same split CoCoA makes.

Can CCL's FKEM take the full TATT IA piece directly (a non-separable `pim`
Pk2D as `p_of_k_a` with a matching `p_of_k_a_lin`)? Mechanically yes
(`cells.py:128-146` passes any Pk2D pair), but its FFTLog term evaluates
`pk(k, 1.0)` with D(chi) growth on the kernels
(`_nonlimber_FKEM.py:473-485`), i.e. it would push the one-loop TATT terms
through the exact projection with a D^2 scaling instead of their D^4 — NOT
CoCoA's model, and not self-consistent. Do not do it; use the split above,
which is both CoCoA's model and the correct use of FKEM.

RSD: clustering tracer has_rsd=True, gamma_t lens tracer has_rsd=False
(CoCoA include_RSD_GG=1, include_RSD_GS=0), unchanged from the ref variant.

### w(theta) (clustering)

No IA enters C_gg; the TATT variant changes nothing (same calls as ref,
including FKEM below l=150). One-loop galaxy bias is off (B2=0) in the
comparison; CoCoA's one-loop bias machinery (`cosmo2D.c:3755-3768`) is out of
scope here.

---

## 4. P(k) inside the TATT terms; EulerianPTCalculator settings

CoCoA:
- One-loop kernels: from the CAMB LINEAR z=0 P(k)
  (`pt_cfastpt.c:413-416` for cfastpt; `PyFAST-PT/fastpt.py:168-182` for the
  python path: Cobaya `Pk_interpolator`, nonlinear=False, z=0,
  `extrap_kmax=250` 1/Mpc beyond CAMB's `kmax_boltzmann` 5-7.5 1/Mpc), scaled
  by D(a)^4 at each Limber node (`cosmo2D.c:2675,2708`).
- Grid: k in [0.05, 1e6] c/H0 = [1.668e-5, 333.6] h/Mpc, N = 1100
  (log-spaced, endpoint=False; ~150 pts/decade), cfastpt
  `pt_cfastpt.c:368-376`, python `fastpt.py:86-99`. Python FAST-PT settings:
  `low_extrap=-5, high_extrap=3, n_pad=0.5*N, P_window=[0.2,0.2],
  C_window=0.65` (`fastpt.py:110-117`). Outside the table's k range the
  kernels are ZERO (`cosmo2D.c:2643, 2694-2711` — KIA is zeroed and only
  filled inside [ln kmin, ln kmax]); `FPTIA.k_cutoff` (1e4 c/H0 = 3.34 h/Mpc)
  is stored but never applied — no exponential damping anywhere.
- Tree-level C1 terms: full CAMB nonlinear P_delta(k, a) (same table the NLA
  comparison already feeds CCL).

CCL `EulerianPTCalculator` (`ept.py:137-154, 221-284`): same architecture —
`pklz0 = linear_matter_power(k_s, 1.0)`, one-loop terms times
`growth_factor(a_arr)**4`, tree level from `nonlin_matter_power` (default
`b1_pk_kind='nonlinear'`). With the CosmologyCalculator of the harness
(CAMB linear + nonlinear tables, CAMB growth), the ingredients are the SAME
functions CoCoA uses.

Converged settings for the harness (defaults in parentheses are CCL's and are
too coarse for this comparison):

```python
ptc = pt.EulerianPTCalculator(
    with_NC=False, with_IA=True, with_matter_1loop=False,
    log10k_min=-5, log10k_max=2.3,      # ~1e-5 .. 200 1/Mpc; CCL default -4..2
    nk_per_decade=150,                  # CoCoA ~150/decade; CCL default 20
    a_arr=a_pk,                         # the export's a grid (same nodes as the
                                        # pk tables; CCL default get_pk_spline_a())
    low_extrap=-6, high_extrap=3,       # fastpt power-law extensions
    pad_factor=1.0,
    P_window=np.array([0.2, 0.2]),      # CoCoA's python FAST-PT uses [.2,.2]
    C_window=0.65,                      # CoCoA uses 0.65 (CCL default 0.75)
    b1_pk_kind='nonlinear', bk2_pk_kind='nonlinear',
    k_cutoff=None,                      # CoCoA applies no damping
    cosmo=cosmo)
```

Notes:
- `a_arr` must cover the full Limber range (z up to ~4 for roman); the
  export's `a_pk` grid does. c1/c2/cdelta `(z, array)` inputs must cover the
  same z range (PTTracer holds the last value beyond,
  `nl_pt/tracers.py:113-116`).
- CoCoA reads its tables with LINEAR interpolation in ln k at each Limber
  node; CCL builds a Pk2D spline over (a_arr, log k) of the already-summed
  spectrum. Pure numerics; at nk_per_decade ~150 and a_arr = a_pk it is far
  inside the 0.2 chi2 tolerance class of differences.
- `P_window`/`C_window` only matter near the grid edges; setting them to
  CoCoA's values removes a gratuitous difference since both codes call the
  same `fastpt` package (lsst_y1's shipped `IA_code: 1` literally; roman's
  cfastpt is a validated port of it).

---

## 5. What cannot be matched exactly (findings, not fixes)

1. **k-support of the one-loop kernels.** CoCoA zeroes the TATT kernels
   outside [1.67e-5, 334] h/Mpc (`cosmo2D.c:2694-2711`) and its table stops
   there; CCL's Pk2D extrapolates the summed IA spectra beyond `k_s` with
   `extrap_order_hik=2` in log k (`ept.py:583-588`). Affects only
   ell ~ tens of thousands at the lowest-chi nodes; a tail-behavior
   difference of the comparison, not configurable away exactly (log10k_max
   2.3 keeps it tiny).
2. **FKEM pivot with massive neutrinos.** CoCoA's non-Limber gs anchors the
   separable linear spectrum per lens bin at z_piv = zmean(bin)
   (`cosmo2D.c:9820-9827`); CCL's FKEM anchors at z=0 with scale-independent
   D(chi) on the kernels (`_nonlimber_FKEM.py:473-485`). Identical only for
   separable growth; with massive neutrinos CAMB's P_lin(k,z) is not exactly
   separable. This is the SAME known finding as the NLA comparison (the
   harness's `separable` variant quantifies it); TATT adds nothing new since
   the TATT extras bypass the FFTLog term in both codes.
3. **One-loop galaxy bias x IA.** CoCoA multiplies one-loop bias by the
   linear-IA-only factor (WK - WS*C1) (`cosmo2D.c:3855`); CCL's `_get_pgi`
   assumes linear galaxy bias for g-IA altogether (warns at `ept.py:372-376`).
   Identical at B2 = B3 = BK = 0 (the comparison configuration); a model
   difference if nonlinear bias is ever switched on with TATT.
4. **Amplitude z-interpolation.** CoCoA evaluates A1(z), A2(z), b_TA and
   D(a) exactly at every quadrature node (`cosmo2D.c:2679-2681`); CCL bakes
   c1(z), c2(z), cdelta(z), g4(a) into the Pk2D at the a_arr nodes and
   interpolates in a. Convergence knob (dense a_arr), not a model difference;
   listed because it is the one place the TATT amplitudes pass through a
   spline in CCL but not in CoCoA.
5. **cfastpt vs python FAST-PT (roman only).** roman_real ships
   `IA_code: 0` (cfastpt: internal grid 1100, `N_pad=1500`,
   `N_extrap=500/500`, `c_window=0.65`, plus two direct-convolution terms,
   `pt_cfastpt.c:773-786, 865-1047`); CCL uses python fastpt. For lsst_y1
   (`IA_code: 1`) both sides call the same python package. The cfastpt/python
   agreement is a validated property of CoCoA, but for roman it is a (tiny)
   numerics layer of the comparison, worth one sentence in the write-up.
6. **Repo benchmark script is NOT a reference.** `ccl_test_lsst.py` deviates
   from the CoCoA conventions in three ways (do not copy): `use_A_ia=True`
   on the IA-only tracer while the amplitudes are already in c1/c2/cdelta
   (double normalization + wrong sign, line 117-122); `a1delta=A1_z` (forces
   b_TA = 1, line 84); xi+- built from EE only (no BB, lines 172-177); and it
   omits the 0.01389/0.0138775 ratio.

---

## 6. Ordered recipe: "tatt" variant of ccl_compute.py

Assumes a CoCoA export made at a TATT point (`IA_model: 1`; the export's
`point` then contains `<pre>A2_1`, `<pre>A2_2`, `<pre>BTA_1` next to the
existing `<pre>A1_1`, `<pre>A1_2`; pre = `LSST_` / `roman_`). Everything not
listed is unchanged from the `ref` variant.

```python
# 0. cosmology: the existing CosmologyCalculator (CAMB linear+nonlinear pk,
#    chi, CAMB growth).  cosmo.growth_factor == CoCoA's growfac by table.

# 1. amplitudes (z grid covering the full n(z)/Limber range, e.g.
#    z = np.linspace(0, 4, 400)):
A1, eta1 = pt_[pre+"A1_1"], pt_[pre+"A1_2"]
A2, eta2 = pt_[pre+"A2_1"], pt_[pre+"A2_2"]
bta      = pt_[pre+"BTA_1"]
c1_ratio = 0.01389/(5e-14*ccl.physical_constants.RHO_CRITICAL)
A1z = c1_ratio*A1*((1+z)/1.62)**eta1
A2z = c1_ratio*A2*((1+z)/1.62)**eta2
c1, cdelta, c2 = ccl.nl_pt.translate_IA_norm(
    cosmo, z=z, a1=A1z, a1delta=A1z*bta, a2=A2z, Om_m2_for_c2=False)

# 2. PT calculator (Sec. 4 settings) + IA PT tracer + spectra:
ptc  = ccl.nl_pt.EulerianPTCalculator(with_IA=True, cosmo=cosmo, a_arr=a_pk,
         log10k_min=-5, log10k_max=2.3, nk_per_decade=150,
         P_window=np.array([0.2, 0.2]), C_window=0.65)
tia  = ccl.nl_pt.PTIntrinsicAlignmentTracer(c1=(z,c1), c2=(z,c2), cdelta=(z,cdelta))
tm   = ccl.nl_pt.PTMatterTracer()
pk_mi    = ptc.get_biased_pk2d(tm,  tracer2=tia)
pk_ii    = ptc.get_biased_pk2d(tia, tracer2=tia)
pk_ii_bb = ptc.get_biased_pk2d(tia, tracer2=tia, return_ia_bb=True)
pk_mi_extra = ccl.Pk2D(a_arr=ptc.a_s, lk_arr=np.log(ptc.k_s), is_logp=False,
    pk_arr=(ptc._g4*np.interp(ptc.z_s, z, cdelta))[:,None]*(ptc.ia_ta[0]+ptc.ia_ta[1])
         + (ptc._g4*np.interp(ptc.z_s, z, c2))[:,None]*(ptc.ia_mix[0]+ptc.ia_mix[1]))
# (pk_mi_extra == pk_mi minus its c1*P_nl term: the TATT-beyond-C1 part of P_gI)

# 3. tracers per source bin j (dndz as in ref):
t_wl[j]  = ccl.WeakLensingTracer(cosmo, dndz=..., ia_bias=None)           # shear only
t_ia[j]  = ccl.WeakLensingTracer(cosmo, dndz=..., has_shear=False,
                                 ia_bias=(zs, np.ones_like(zs)), use_A_ia=False)
src_nla[j] = ccl.WeakLensingTracer(cosmo, dndz=...,                        # as ref:
                                 ia_bias=(zs, c1_ratio*A1*((1+zs)/1.62)**eta1),
                                 use_A_ia=True)                            # NLA composite
# lens tracers unchanged from ref (linear bias; has_rsd per ref/rsd_gs rules)

# 4. shear-shear (Limber everywhere, as CoCoA):
Cl_EE = cl(t_wl[i], t_wl[j]) \
      + cl(t_wl[i], t_ia[j], p_of_k_a=pk_mi) + cl(t_ia[i], t_wl[j], p_of_k_a=pk_mi) \
      + cl(t_ia[i], t_ia[j], p_of_k_a=pk_ii)
Cl_BB = cl(t_ia[i], t_ia[j], p_of_k_a=pk_ii_bb)
xip[i,j] = corr(Cl_EE + Cl_BB, 'GG+')      # bin-averaged Legendre, as ref
xim[i,j] = corr(Cl_EE - Cl_BB, 'GG-')

# 5. gamma_t: the ref FKEM call with the NLA composite tracer, plus the
#    Limber-only TATT extras:
Cl_gs = cl(lens_tgs[l], src_nla[s], nonlimber=True)            # ref machinery
Cl_gs += cl(lens_tgs[l], t_ia[s], p_of_k_a=pk_mi_extra)        # Limber
gt[l,s] = corr(Cl_gs, 'NG')

# 6. w(theta): unchanged from ref (no IA in gg).
```

Checks to run first (cheap): (a) with A2 = bta = 0 the tatt variant must
reproduce the ref NLA data vector to transform accuracy (CoCoA's TATT cores
reduce identically to NLA, `cosmo2D.c:2740-2741`); (b) Cl_BB = 0 when
A2 = bta = 0; (c) the GI sign: Cl_EE(i,j) < Cl_GG(i,j) for A1 > 0 at large
scales.
