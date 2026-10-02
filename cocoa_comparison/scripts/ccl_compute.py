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
             CoCoA), RSD in clustering, full-sky (Legendre) transforms
             averaged over each theta bin, CAMB's linear and nonlinear tables,
             converged FKEM sampling (fkem_Nchi 2000, 1500 log-spaced l above
             400)
  limgs      C_gs in Limber (the CCL-benchmark scripts)
  limber     C_gs and C_gg in Limber
  norsd      no RSD (the CCL-benchmark scripts)
  flat       flat-sky FFTLog transform at the bin centers (the CCL-benchmark
             scripts; FFTLog has no bin averaging)
  points     full-sky correlations at the bin centers, no bin averaging
  separable  diagnostic: the linear table in the separable form
             P_lin(k,0) D(z)^2 that FKEM's FFTLog term assumes
  ccl_numerics  CCL's default FKEM sampling and l sampling
  bench      the CCL-benchmark real-space scripts on CoCoA's layout: CCL's own
             Eisenstein-Hu linear and halofit nonlinear P(k), C_gs in Limber,
             no RSD, flat-sky FFTLog at the bin centers, l_limber = 100,
             fkem_Nchi = 500
  harmonic   Limber C_l at CoCoA's exported multipoles (no real-space step)
"""
import json, os, sys, time
import numpy as np
if not hasattr(np, "trapz"):     # FAST-PT/numpy 2.4 pairing of the local ccl env
  np.trapz = np.trapezoid
import pyccl as ccl

src, variant, out = sys.argv[1:4]
H = os.path.dirname(os.path.abspath(__file__))
ex = np.load(src)
meta = json.loads(str(ex["meta"]))
project = meta["project"]; par = meta["pars"]; pt = meta["point"]
V = dict(nonlimber_gs=True, nonlimber_gg=True, rsd=True, method="legendre",
         binavg=True, native=False, separable=False, l_limber=150, fkem_nchi=2000,
         nell_log=1500, ell_max=65000 if project == "lsst_y1" else 100000)
CHANGES = {
  "ref": {}, "harmonic": {},
  "limgs": dict(nonlimber_gs=False),
  "limber": dict(nonlimber_gs=False, nonlimber_gg=False),
  "norsd": dict(rsd=False),
  "flat": dict(method="fftlog", binavg=False),
  "points": dict(binavg=False),
  "separable": dict(separable=True),
  "ccl_numerics": dict(fkem_nchi=None, nell_log=500),
  "bench": dict(native=True, nonlimber_gs=False, rsd=False, method="fftlog",
                binavg=False, l_limber=100, fkem_nchi=500, nell_log=500),
}
if variant not in CHANGES: sys.exit("unknown variant " + variant)
V.update(CHANGES[variant])
ccl.spline_params.ELL_MAX_CORR = float(V["ell_max"])
ccl.spline_params.N_ELL_CORR = int(max(5000, V["nell_log"]*10))

# ---- cosmology: CoCoA's tables, or CCL's own (bench) ----------------------
h = par["H0"]/100.0
Onu = par["omnuh2"]/h**2
common = dict(Omega_c=par["omegam"] - par["omegab"] - Onu, Omega_b=par["omegab"], h=h,
              n_s=par["ns"], A_s=par["As"], m_nu=par["mnu"], mass_split="single",
              w0=par["w"], wa=par.get("wa", 0.0))
def calculator():
  z1, chi1 = ex["z1"], ex["chi1"]
  sel = z1 <= 10.0; z1, chi1 = z1[sel], chi1[sel]
  dchidz = np.gradient(chi1, z1)
  hoh0 = (ccl.physical_constants.CLIGHT/1e3/par["H0"])/dchidz
  zg, kg = ex["zg"], ex["kg"]
  a_bg = 1.0/(1.0 + z1[::-1])
  a_pk = 1.0/(1.0 + zg[::-1])
  D = ex["D0"][::-1]
  if "DD" in ex.files:   # growth over CAMB's full z range (to z = 49)
    a_D = 1.0/(1.0 + ex["zD"][::-1]); DDa = ex["DD"][::-1]
  else:
    a_D, DDa = a_pk, D
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
ia = lambda z: (z, c1_ratio*A1*((1.0 + z)/1.62)**eta)
src_t = [ccl.WeakLensingTracer(cosmo, dndz=(zs, nzs[i]), ia_bias=ia(zs), use_A_ia=True) for i in range(S)]
lens_t = [ccl.NumberCountsTracer(cosmo, dndz=(zl, nzl[i]), has_rsd=V["rsd"],
          bias=(zl, pt["%sB1_%d" % (pre, i + 1)]*np.ones_like(zl))) for i in range(L)]
excl = set()
ef = os.path.join(H, "ggl_exclude_%s.txt" % project)
if os.path.isfile(ef) and os.path.getsize(ef) > 0:
  excl = {tuple(int(x) for x in l.split()) for l in open(ef) if l.strip()}

# ---- spectra --------------------------------------------------------------
ell = np.unique(np.concatenate([np.arange(2, 400), np.geomspace(400, V["ell_max"], V["nell_log"]).astype(int)])).astype(float)
def cl(t1, t2, nonlimber):
  kw = {}
  if nonlimber:
    kw = dict(l_limber=V["l_limber"], non_limber_integration_method="FKEM")
    if V["fkem_nchi"]: kw["fkem_Nchi"] = V["fkem_nchi"]
  return ccl.angular_cl(cosmo, t1, t2, ell, **kw)

# ---- harmonic mode: Limber C_l at CoCoA's exported multipoles ----------
if variant == "harmonic":
  eh = ex["ell_h"]
  lim = lambda t1, t2: ccl.angular_cl(cosmo, t1, t2, eh)
  css = np.array([[lim(src_t[i], src_t[j]) for j in range(S)] for i in range(S)])
  cgs = np.array([[lim(lens_t[l], src_t[s]) for s in range(S)] for l in range(L)])
  cgg = np.array([lim(lens_t[l], lens_t[l]) for l in range(L)])
  lens_nr = [ccl.NumberCountsTracer(cosmo, dndz=(zl, nzl[i]), has_rsd=False,
             bias=(zl, pt["%sB1_%d" % (pre, i + 1)]*np.ones_like(zl))) for i in range(L)]
  cgg_nr = np.array([lim(lens_nr[l], lens_nr[l]) for l in range(L)])
  np.savez(out, ell=eh, css=css, cgs=cgs, cgg=cgg, cgg_nr=cgg_nr, pyccl=ccl.__file__)
  print("CCL %s %s harmonic done" % (project, meta["model"])); sys.exit(0)

# ---- theta bins: CoCoA's log bins ----------------------------------------
edges = np.logspace(np.log10(meta["theta_min"]), np.log10(meta["theta_max"]), meta["ntheta"] + 1)/60.0
def _corr(c, typ):
  if V["binavg"]:   # exact bin average (integrated polynomials, PR #1296)
    return ccl.correlation(cosmo, ell=ell, C_ell=c, theta=edges[:-1], theta_max=edges[1:], type=typ, method=V["method"])
  return ccl.correlation(cosmo, ell=ell, C_ell=c, theta=np.sqrt(edges[:-1]*edges[1:]), type=typ, method=V["method"])
FAILURES = []
def corr(c, typ, label=""):
  try:
    return _corr(c, typ)
  except Exception as e:   # a CCL limitation: recorded, the pair left out
    FAILURES.append(dict(pair=label, type=typ, error=str(e).strip().splitlines()[-1][:160],
                         nonpositive_last=bool(np.any(c[-2:] <= 0))))
    return np.full(meta["ntheta"], np.nan)

t0 = time.perf_counter()
xip, xim, gt, w = [], [], [], []
css = [cl(src_t[i], src_t[j], False) for i in range(S) for j in range(i, S)]
pairs_ss = [(i, j) for i in range(S) for j in range(i, S)]
xip = [corr(c, "GG+", "s%d-s%d" % p) for c, p in zip(css, pairs_ss)]
xim = [corr(c, "GG-", "s%d-s%d" % p) for c, p in zip(css, pairs_ss)]
for l in range(L):
  for s in range(S):
    if (l, s) in excl: continue
    gt.append(corr(cl(lens_t[l], src_t[s], V["nonlimber_gs"]), "NG", "l%d-s%d" % (l, s)))
for l in range(L):
  w.append(corr(cl(lens_t[l], lens_t[l], V["nonlimber_gg"]), "NN", "l%d-l%d" % (l, l)))
dv = np.concatenate(xip + xim + gt + w)
dt = time.perf_counter() - t0
sizes = list(ex["sizes"])
ok = len(dv) == len(ex["dv"])
np.savez(out, dv=dv, variant=variant, V=json.dumps(V), time=dt, pyccl=ccl.__file__,
         failures=json.dumps(FAILURES))
print("CCL %s %s %s len %d (CoCoA %d, %s) %.1f s, failed pairs: %s" % (project, meta["model"], variant, len(dv), len(ex["dv"]), "ok" if ok else "MISMATCH", dt, [f["pair"] + ":" + f["type"] for f in FAILURES] or "none"), flush=True)
