"""DESC-CCL side of the CoCoA vs DESC-CCL comparison (ccl env, with pyccl
built from LSSTDESC/CCL PR #1296, commit 647ad4a, first on PYTHONPATH: it
adds full-sky xi+- and bin-averaged correlations, absent from CCL 3.3).

python ccl_compute.py <cocoa export .npz> <variant> <out.npz>

The real-space 3x2pt data vector of a CoCoA project (lsst_y1, roman_real)
computed with CCL, on CoCoA's layout: xi+ and xi- for source pairs i <= j,
gamma_t for lens-major source-minor pairs (CoCoA's exclusions), w(theta)
for lens autos, each over CoCoA's log theta bins.

CCL is given CoCoA's inputs through pyccl.CosmologyCalculator: the CAMB
linear and nonlinear P(k), comoving distance and growth that cosmolike
received at the same point, so no Boltzmann or background difference
enters. The NLA amplitude uses cosmolike's convention (A1 ((1+z)/1.62)^eta,
C1 rho_crit = 0.01389), galaxy bias is linear, n(z) are read from CoCoA's
files with cosmolike's Z_LOW convention (values at z + dz/2).

CCL is used as it is: the only adaptations are CoCoA's CAMB tables passed
through pyccl.CosmologyCalculator, the PR #1296 build (full-sky xi+- and bin
averaging), and CCL's documented accuracy settings. When a CCL transform
fails for a bin pair, the failure is recorded and that pair is left out
(NaN) rather than worked around.

Variants (each a single change from "ref" unless stated):
  ref        non-Limber C_gg and C_gs below l = 150 (FKEM, the same switch as
             CoCoA), RSD in clustering only (CoCoA's model: include_RSD_GG = 1,
             include_RSD_GS = 0), full-sky (Legendre) transforms
             averaged over each theta bin, CAMB's linear and nonlinear tables,
             converged FKEM sampling (fkem_Nchi 2000, 1500 log-spaced l above
             400)
  limgs      C_gs in Limber (the CCL-benchmark scripts)
  limber     C_gs and C_gg in Limber
  norsd      no RSD (the CCL-benchmark scripts)
  rsd_gs     RSD also in the lens tracer of gamma_t (CCL's has_rsd applies
             to every correlation of the tracer; CoCoA omits it in C_gs)
  flat       flat-sky FFTLog transform at the bin centers (the CCL-benchmark
             scripts; FFTLog has no bin averaging)
  points     full-sky correlations at the bin centers, no bin averaging
  separable  diagnostic: the linear table in the separable form
             P_lin(k,0) D(z)^2 that FKEM's FFTLog term assumes
  growth_sub diagnostic: the growth factor at k = 0.05/Mpc instead of
             CoCoA's k0 = 5e-4/Mpc (NLA amplitude kept CoCoA's)
  ccl_numerics  CCL's default FKEM sampling and l sampling
  ccl_hi     finer than ref: fkem_Nchi 4000, 3000 log-spaced l above 400
  eh         CCL's own Eisenstein-Hu linear and halofit nonlinear P(k)
             (sigma8 of CoCoA's linear table), nothing else changed
  cocoa_cl   transform-only test: CoCoA's Limber C_l (the export's cld_*,
             on this script's l grid) through CCL's real-space transform;
             compare with CoCoA's Limber run (cocoa_<project>_fidlimber)
  bench      the CCL-benchmark real-space scripts on CoCoA's layout: CCL's own
             Eisenstein-Hu linear and halofit nonlinear P(k), C_gs in Limber,
             no RSD, flat-sky FFTLog at the bin centers, l_limber = 100,
             fkem_Nchi = 500
  harmonic   Limber C_l at CoCoA's exported multipoles (no real-space step)

TATT: a variant name prefixed "tatt:" (tatt:ref, tatt:limber, tatt:separable,
tatt:harmonic, tatt:pt_hi) computes the same layout with TATT intrinsic
alignments, on a CoCoA export made with IA_model 1 (CFASTPT). Conventions
(.claude/skills/cocoa-ccl-comparison/references/tatt_conventions.md): the
amplitudes enter pyccl's PTIntrinsicAlignmentTracer through translate_IA_norm
(a1 = A1(z), a1delta = b_TA A1(z), a2 = A2(z), each times 0.01389/(5e-14
rho_crit), Om_m2_for_c2=False); EulerianPTCalculator gives P_GI, P_II^EE and
P_II^BB; xi+ takes C_EE + C_BB, xi- C_EE - C_BB. gamma_t: the NLA part (C1)
goes through the same FKEM call as the reference and the TATT extras (b_TA,
A2 terms) are added in Limber, as CoCoA's non-Limber C_gs carries only C1.
  tatt:pt_hi  the PT tables at 320 k per decade instead of 160 (FAST-PT needs an
             even number of k points: 160 x 7.3 decades = 1168)

A CCL call that fails for a bin pair (angular_cl or correlation) is recorded
in the output's "failures" list and that pair is left out (NaN).
"""
import json, os, sys, time
import numpy as np
from scipy.interpolate import CubicSpline
import pyccl as ccl

src, variant, out = sys.argv[1:4]
H = os.path.dirname(os.path.abspath(__file__))
ex = np.load(src)
meta = json.loads(str(ex["meta"]))
project = meta["project"]; par = meta["pars"]; pt = meta["point"]
V = dict(nonlimber_gs=True, nonlimber_gg=True, rsd=True, rsd_gs=False, cocoa_cl=False, tatt=False, pt_nk=160, method="legendre",
         binavg=True, native=False, separable=False, growth_k=None, l_limber=150, fkem_nchi=2000,
         nell_log=1500, ell_max=65000 if project == "lsst_y1" else 100000)
CHANGES = {
  "ref": {}, "harmonic": {},
  "limgs": dict(nonlimber_gs=False),
  "limber": dict(nonlimber_gs=False, nonlimber_gg=False),
  "norsd": dict(rsd=False),
  "rsd_gs": dict(rsd_gs=True),
  "flat": dict(method="fftlog", binavg=False),
  "points": dict(binavg=False),
  "separable": dict(separable=True),
  "growth_sub": dict(growth_k=0.05),
  "ccl_numerics": dict(fkem_nchi=None, nell_log=500),
  "ccl_hi": dict(fkem_nchi=4000, nell_log=3000),
  "eh": dict(native=True),
  "cocoa_cl": dict(nonlimber_gs=False, nonlimber_gg=False, cocoa_cl=True),
  "pt_hi": dict(pt_nk=320),
  "bench": dict(native=True, nonlimber_gs=False, rsd=False, method="fftlog",
                binavg=False, l_limber=100, fkem_nchi=500, nell_log=500),
}
base = variant.split(":", 1)[1] if variant.startswith("tatt:") else variant
if base not in CHANGES: sys.exit("unknown variant " + variant)
V.update(CHANGES[base]); V["tatt"] = variant.startswith("tatt:")
ccl.spline_params.ELL_MAX_CORR = float(V["ell_max"])
ccl.spline_params.N_ELL_CORR = int(max(5000, V["nell_log"]*10))

# ---- cosmology: CoCoA's tables, or CCL's own (bench) ----------------------
h = par["H0"]/100.0
Onu = par["omnuh2"]/h**2
common = dict(Omega_c=par["omegam"] - par["omegab"] - Onu, Omega_b=par["omegab"], h=h,
              n_s=par["ns"], A_s=par["As"], m_nu=par["mnu"], mass_split="single",
              w0=par["w"], wa=par.get("wa", 0.0))
GROWTH_RATIO = [lambda z: np.ones_like(z)]   # growth_sub: D(k, z)/D(k0, z)
def calculator():
  # CoCoA's chi(z) over its whole table (z <= 50): CCL's Limber RSD kernel
  # at l = 2 reads the background at 1.4 chi (z ~ 14 for lenses at z = 4).
  # H(z) from a cubic-spline derivative: np.gradient is off by 3e-4 where
  # the z grid coarsens (z > 3); the spline matches CCL's own H(z) to 1e-5
  z1, chi1 = ex["z1"], ex["chi1"]
  dchidz = CubicSpline(z1, chi1).derivative()(z1)
  hoh0 = (ccl.physical_constants.CLIGHT/1e3/par["H0"])/dchidz
  zg, kg = ex["zg"], ex["kg"]
  a_bg = 1.0/(1.0 + z1[::-1])
  a_pk = 1.0/(1.0 + zg[::-1])
  D = ex["D0"][::-1]
  if "DD" in ex.files:   # growth over CAMB's full z range (to z = 49)
    a_D = 1.0/(1.0 + ex["zD"][::-1]); DDa = ex["DD"][::-1]
  else:
    a_D, DDa = a_pk, D
  if V["growth_k"] is not None:
    # diagnostic: the growth factor at a sub-horizon k instead of k0 = 5e-4/Mpc
    # (for w != -1, CAMB's dark-energy perturbations change the growth on
    # horizon scales): D -> D r(z), r = D(k, z)/D(k0, z) from CAMB's linear
    # table (z <= 6; held at r(6) above)
    lk = np.log(ex["kg"]); lp = ex["lnPL"]
    lnP = lambda kk: np.array([np.interp(np.log(kk), lk, row) for row in lp])
    r = np.exp(0.5*(lnP(V["growth_k"]) - lnP(V["growth_k"])[0]) - 0.5*(lnP(5e-4) - lnP(5e-4)[0]))
    GROWTH_RATIO[0] = lambda z: np.interp(z, ex["zg"], r)
    D = D*GROWTH_RATIO[0](ex["zg"])[::-1]
    DDa = DDa*GROWTH_RATIO[0](1.0/a_D - 1.0)
  f = np.gradient(np.log(DDa), np.log(a_D))
  # the linear table CCL receives: CAMB's P_lin(k, z), or (separable
  # variants) P_lin(k, 0) D(z)^2, the separable form FKEM's FFTLog term
  # assumes; with massive neutrinos CAMB's growth depends on k
  PLIN = np.exp(ex["lnPL"][::-1])
  if V["separable"]:
    PLIN = np.exp(ex["lnPL"][0])[None, :]*(D**2)[:, None]
  cosmo = ccl.CosmologyCalculator(
    **common,
    background={"a": a_bg, "chi": chi1[::-1], "h_over_h0": hoh0[::-1]},
    growth={"a": a_D, "growth_factor": DDa, "growth_rate": f},
    pk_linear={"a": a_pk, "k": kg, "delta_matter:delta_matter": PLIN},
    pk_nonlin={"a": a_pk, "k": kg, "delta_matter:delta_matter": np.exp(ex["lnPN"][::-1])})
  return cosmo
cosmo = calculator()
if V["native"]:
  # CCL's own Eisenstein-Hu + halofit; the analytic linear spectrum takes
  # sigma8, set to the value of CoCoA's CAMB linear table
  s8 = ccl.sigma8(cosmo)
  c2 = dict(common); c2.pop("A_s")
  cosmo = ccl.Cosmology(**c2, sigma8=s8, transfer_function="eisenstein_hu",
                        matter_power_spectrum="halofit")

# ---- tracers --------------------------------------------------------------
ds = meta["dataset"]
def read_nz(fn, n):
  base = meta["path"] if os.path.isabs(meta["path"]) else os.path.join(meta.get("rootdir") or os.environ["COCOA_ROOTDIR"], meta["path"])
  d = np.loadtxt(os.path.join(base, fn))
  z = d[:, 0] + 0.5*(d[1, 0] - d[0, 0])     # cosmolike: the column is Z_LOW
  return z, [d[:, i + 1] for i in range(n)]
L, S = meta["lens_ntomo"], meta["source_ntomo"]
zl, nzl = read_nz(ds["nz_lens_file"], L)
zs, nzs = read_nz(ds["nz_source_file"], S)
pre = "LSST_" if project == "lsst_y1" else "roman_"
A1, eta = pt[pre + "A1_1"], pt[pre + "A1_2"]
c1_ratio = 0.01389/(5e-14*ccl.physical_constants.RHO_CRITICAL)   # cosmolike / CCL C1 rho_crit
# CCL's NLA divides by its growth factor: growth_sub multiplies by the same
# ratio, so the IA amplitude stays CoCoA's
ia = lambda z: (z, c1_ratio*A1*((1.0 + z)/1.62)**eta*GROWTH_RATIO[0](z))
src_t = [ccl.WeakLensingTracer(cosmo, dndz=(zs, nzs[i]), ia_bias=ia(zs), use_A_ia=True) for i in range(S)]
lens_t = [ccl.NumberCountsTracer(cosmo, dndz=(zl, nzl[i]), has_rsd=V["rsd"],
          bias=(zl, pt["%sB1_%d" % (pre, i + 1)]*np.ones_like(zl))) for i in range(L)]
# the lens tracer of gamma_t: RSD only on request (CoCoA: include_RSD_GS = 0)
lens_tgs = [ccl.NumberCountsTracer(cosmo, dndz=(zl, nzl[i]), has_rsd=V["rsd"] and V["rsd_gs"],
            bias=(zl, pt["%sB1_%d" % (pre, i + 1)]*np.ones_like(zl))) for i in range(L)]
if V["tatt"]:
  # TATT (references/tatt_conventions.md): CoCoA's A1, eta1, A2, eta2, b_TA
  # with its pivot (1+z)/1.62 and C1 rho_crit = 0.01389
  from pyccl import nl_pt
  bta, A2, eta2 = pt[pre + "BTA_1"], pt[pre + "A2_1"], pt[pre + "A2_2"]
  zt = np.linspace(0.0, 6.0, 601)
  A1z = c1_ratio*A1*((1.0 + zt)/1.62)**eta
  A2z = c1_ratio*A2*((1.0 + zt)/1.62)**eta2
  c1t, cdt, c2t = nl_pt.translate_IA_norm(cosmo, z=zt, a1=A1z, a1delta=bta*A1z, a2=A2z,
                                          Om_m2_for_c2=False)
  tia = nl_pt.PTIntrinsicAlignmentTracer(c1=(zt, c1t), c2=(zt, c2t), cdelta=(zt, cdt))
  tm = nl_pt.PTMatterTracer()
  # FAST-PT settings as CoCoA's tables (references/tatt_conventions.md, 4)
  ptc = nl_pt.EulerianPTCalculator(with_IA=True, with_matter_1loop=False, cosmo=cosmo,
                                   log10k_min=-5, log10k_max=2.3, nk_per_decade=V["pt_nk"],
                                   a_arr=1.0/(1.0 + ex["zg"][::-1]), low_extrap=-6, high_extrap=3,
                                   pad_factor=1.0, P_window=np.array([0.2, 0.2]), C_window=0.65,
                                   b1_pk_kind="nonlinear")
  pk_mi = ptc.get_biased_pk2d(tm, tracer2=tia)
  pk_ii = ptc.get_biased_pk2d(tia, tracer2=tia)
  pk_bb = ptc.get_biased_pk2d(tia, tracer2=tia, return_ia_bb=True)
  # the TATT extras of P_GI (P_GI minus its C1 P_nl term), for gamma_t
  a_t, lk_t, pk_cd = ptc.get_pk2d_template("m:cdelta").get_spline_arrays()
  _, _, pk_c2 = ptc.get_pk2d_template("m:c2").get_spline_arrays()
  z_t = 1.0/a_t - 1.0
  pk_mi_extra = ccl.Pk2D(a_arr=a_t, lk_arr=lk_t, is_logp=False,
                         pk_arr=np.interp(z_t, zt, cdt)[:, None]*pk_cd + np.interp(z_t, zt, c2t)[:, None]*pk_c2)
  # shear-only and IA-only tracers (the IA amplitudes sit in c1, c2, cdelta)
  wl_t = [ccl.WeakLensingTracer(cosmo, dndz=(zs, nzs[i]), has_shear=True, ia_bias=None) for i in range(S)]
  ia_t = [ccl.WeakLensingTracer(cosmo, dndz=(zs, nzs[i]), has_shear=False,
                                ia_bias=(zs, np.ones_like(zs)), use_A_ia=False) for i in range(S)]
  lim_pk = lambda t1, t2, pk, ells: ccl.angular_cl(cosmo, t1, t2, ells, p_of_k_a=pk)
excl = set()
ef = os.path.join(H, "ggl_exclude_%s.txt" % project)
if os.path.isfile(ef) and os.path.getsize(ef) > 0:
  excl = {tuple(int(x) for x in l.split()) for l in open(ef) if l.strip()}

# ---- spectra --------------------------------------------------------------
ell = np.unique(np.concatenate([np.arange(2, 400), np.geomspace(400, V["ell_max"], V["nell_log"]).astype(int)])).astype(float)
def css_tatt(i, j, ells):
  # C_EE = GG + GI + IG + II and C_BB = II_BB, Limber (as CoCoA's C_ss)
  ee = (lim_pk(wl_t[i], wl_t[j], "delta_matter:delta_matter", ells) + lim_pk(wl_t[i], ia_t[j], pk_mi, ells)
        + lim_pk(ia_t[i], wl_t[j], pk_mi, ells) + lim_pk(ia_t[i], ia_t[j], pk_ii, ells))
  return ee, lim_pk(ia_t[i], ia_t[j], pk_bb, ells)
def cgs_tatt_extra(l, s, ells):
  # the TATT extras of C_gI (b_TA and A2 terms), Limber
  return lim_pk(lens_tgs[l], ia_t[s], pk_mi_extra, ells)
def cl(t1, t2, nonlimber):
  kw = {}
  if nonlimber:
    kw = dict(l_limber=V["l_limber"], non_limber_integration_method="FKEM")
    if V["fkem_nchi"]: kw["fkem_Nchi"] = V["fkem_nchi"]
  return ccl.angular_cl(cosmo, t1, t2, ell, **kw)

# ---- harmonic mode: Limber C_l at CoCoA's exported multipoles ----------
if base == "harmonic":
  eh = ex["ell_h"]
  lim = lambda t1, t2: ccl.angular_cl(cosmo, t1, t2, eh)
  if V["tatt"]:   # EE, BB and C_gs with TATT (Limber)
    EB = [[css_tatt(i, j, eh) for j in range(S)] for i in range(S)]
    css = np.array([[EB[i][j][0] for j in range(S)] for i in range(S)])
    cssbb = np.array([[EB[i][j][1] for j in range(S)] for i in range(S)])
    cgs = np.array([[lim(lens_tgs[l], src_t[s]) + cgs_tatt_extra(l, s, eh) for s in range(S)] for l in range(L)])
  else:
    css = np.array([[lim(src_t[i], src_t[j]) for j in range(S)] for i in range(S)])
    cgs = np.array([[lim(lens_tgs[l], src_t[s]) for s in range(S)] for l in range(L)])
  cgg = np.array([lim(lens_t[l], lens_t[l]) for l in range(L)])
  lens_nr = [ccl.NumberCountsTracer(cosmo, dndz=(zl, nzl[i]), has_rsd=False,
             bias=(zl, pt["%sB1_%d" % (pre, i + 1)]*np.ones_like(zl))) for i in range(L)]
  cgg_nr = np.array([lim(lens_nr[l], lens_nr[l]) for l in range(L)])
  extra = {}
  if "ell_nl" in ex.files:   # non-Limber C_gg (FKEM, reference settings)
    enl = ex["ell_nl"]
    kw = dict(l_limber=V["l_limber"], non_limber_integration_method="FKEM", fkem_Nchi=V["fkem_nchi"])
    extra = dict(ell_nl=enl, cgg_nl=np.array([ccl.angular_cl(cosmo, lens_t[l], lens_t[l], enl, **kw) for l in range(L)]))
  if V["tatt"]: extra["cssbb"] = cssbb
  np.savez(out, ell=eh, css=css, cgs=cgs, cgg=cgg, cgg_nr=cgg_nr, pyccl=ccl.__file__, **extra)
  print("CCL %s %s harmonic done" % (project, meta["model"])); sys.exit(0)

# ---- theta bins: CoCoA's log bins ----------------------------------------
edges = np.logspace(np.log10(meta["theta_min"]), np.log10(meta["theta_max"]), meta["ntheta"] + 1)/60.0
def _corr(c, typ):
  if V["binavg"]:   # exact bin average (integrated polynomials, PR #1296)
    return ccl.correlation(cosmo, ell=ell, C_ell=c, theta=edges[:-1], theta_max=edges[1:], type=typ, method=V["method"])
  return ccl.correlation(cosmo, ell=ell, C_ell=c, theta=np.sqrt(edges[:-1]*edges[1:]), type=typ, method=V["method"])
FAILURES = []
def safe_cl(t1, t2, nonlimber, typ, label):
  try:
    return cl(t1, t2, nonlimber)
  except Exception as e:   # a CCL limitation: recorded, the pair left out
    FAILURES.append(dict(pair=label, type=typ, stage="angular_cl",
                         error=str(e).strip().splitlines()[-1][:160]))
    return None
def corr(c, typ, label=""):
  if c is None: return np.full(meta["ntheta"], np.nan)
  try:
    return _corr(c, typ)
  except Exception as e:   # a CCL limitation: recorded, the pair left out
    FAILURES.append(dict(pair=label, type=typ, stage="correlation", error=str(e).strip().splitlines()[-1][:160],
                         nonpositive_last=bool(np.any(c[-2:] <= 0))))
    return np.full(meta["ntheta"], np.nan)

t0 = time.perf_counter()
xip, xim, gt, w = [], [], [], []
if V["cocoa_cl"]:   # CoCoA's Limber C_l, on this script's l grid
  assert np.array_equal(ex["ell_d"], ell), "the export's ell_d is not this script's l grid"
  ccs = lambda i, j: ex["cld_ss"][0][:, i, j]
  cgs_ = lambda l, s: ex["cld_gs"][:, l, s]
  cgg_ = lambda l: ex["cld_gg"][:, l, l]
pairs_ss = [(i, j) for i in range(S) for j in range(i, S)]
if V["tatt"]:     # xi+ from C_EE + C_BB, xi- from C_EE - C_BB
  EB = [css_tatt(i, j, ell) for i, j in pairs_ss]
  xip = [corr(e + b, "GG+", "s%d-s%d" % p) for (e, b), p in zip(EB, pairs_ss)]
  xim = [corr(e - b, "GG-", "s%d-s%d" % p) for (e, b), p in zip(EB, pairs_ss)]
else:
  css = ([ccs(i, j) for i in range(S) for j in range(i, S)] if V["cocoa_cl"] else
         [cl(src_t[i], src_t[j], False) for i in range(S) for j in range(i, S)])
  xip = [corr(c, "GG+", "s%d-s%d" % p) for c, p in zip(css, pairs_ss)]
  xim = [corr(c, "GG-", "s%d-s%d" % p) for c, p in zip(css, pairs_ss)]
for l in range(L):
  for s in range(S):
    if (l, s) in excl: continue
    c = cgs_(l, s) if V["cocoa_cl"] else safe_cl(lens_tgs[l], src_t[s], V["nonlimber_gs"], "NG", "l%d-s%d" % (l, s))
    if V["tatt"] and c is not None:   # TATT extras in Limber (CoCoA: C1 only in FKEM)
      c = c + cgs_tatt_extra(l, s, ell)
    gt.append(corr(c, "NG", "l%d-s%d" % (l, s)))
for l in range(L):
  c = cgg_(l) if V["cocoa_cl"] else safe_cl(lens_t[l], lens_t[l], V["nonlimber_gg"], "NN", "l%d-l%d" % (l, l))
  w.append(corr(c, "NN", "l%d-l%d" % (l, l)))
dv = np.concatenate(xip + xim + gt + w)
dt = time.perf_counter() - t0
sizes = list(ex["sizes"])
ok = len(dv) == len(ex["dv"])
np.savez(out, dv=dv, variant=variant, V=json.dumps(V), time=dt, pyccl=ccl.__file__,
         failures=json.dumps(FAILURES))
print("CCL %s %s %s len %d (CoCoA %d, %s) %.1f s, failed pairs: %s" % (project, meta["model"], variant, len(dv), len(ex["dv"]), "ok" if ok else "MISMATCH", dt, [f["pair"] + ":" + f["type"] + ":" + f["stage"] for f in FAILURES] or "none"), flush=True)
