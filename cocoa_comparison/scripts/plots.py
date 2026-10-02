"""Ratio figures of the CoCoA vs DESC-CCL study, drawn with the
data-vector plotters of cosmolike_core (cosmolike_notebook_utils/
plot_datavectors.py: plot_xi, plot_gammat_tomo_limber, plot_wtheta_tomo)
in their ratio mode. Every curve is DESC-CCL / CoCoA - 1 at the same
point: the curves enter as the ratio and the reference is 1 (0 for pairs
a project excludes). Points removed by the project's mask, and pairs DESC-CCL
could not transform, are left blank.
python plots.py <outdir>   (cocoa env: matplotlib, ROOTDIR set)"""
import importlib.util, json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
# the plot options of the projects' EXAMPLE_EVALUATE notebooks
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['xtick.bottom'] = True
matplotlib.rcParams['xtick.top'] = False
matplotlib.rcParams['ytick.right'] = False
matplotlib.rcParams['axes.edgecolor'] = 'black'
matplotlib.rcParams['axes.linewidth'] = '1.0'
matplotlib.rcParams['axes.labelsize'] = 'medium'
matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['grid.linewidth'] = '0.0'
matplotlib.rcParams['grid.alpha'] = '0.18'
matplotlib.rcParams['grid.color'] = 'lightgray'
matplotlib.rcParams['legend.labelspacing'] = 0.77
matplotlib.rcParams['legend.fontsize'] = 17
matplotlib.rcParams['savefig.bbox'] = 'tight'
H = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
pd_path = os.path.join(os.environ["ROOTDIR"], "external_modules/code/cosmolike_core/cosmolike_notebook_utils/plot_datavectors.py")
spec = importlib.util.spec_from_file_location("plot_datavectors", pd_path)
pdv = importlib.util.module_from_spec(spec); spec.loader.exec_module(pdv)

PROJ = {"lsst_y1": "LSST-Y1", "roman_real": "Roman-Real"}
MODELS = [("fid", "fiducial"), ("omm_lo", r"$\Omega_m=0.25$"), ("omm_hi", r"$\Omega_m=0.35$"),
          ("ns_lo", r"$n_s=0.92$"), ("ns_hi", r"$n_s=1.01$")]
VARIANTS = [("ref", "reference"), ("limgs", r"$C_{gs}$ Limber"), ("norsd", "no RSD"),
            ("flat", "flat-sky, bin centers"), ("points", "full-sky, bin centers"),
            ("bench", "CCL-benchmark modeling")]

def unpack(p, dv):
  """CoCoA's data-vector layout -> xi+ / xi- (ntheta, S, S), gamma_t
  (ntheta, L, S), w (ntheta, L, L)."""
  meta = json.loads(str(np.load(os.path.join(H, "cocoa_%s_fid.npz" % p))["meta"]))
  S, L, nt = meta["source_ntomo"], meta["lens_ntomo"], meta["ntheta"]
  ef = os.path.join(H, "ggl_exclude_%s.txt" % p)
  excl = {tuple(int(x) for x in l.split()) for l in open(ef) if l.strip()} if os.path.isfile(ef) else set()
  xp = np.zeros((nt, S, S)); xm = np.zeros((nt, S, S)); gt = np.zeros((nt, L, S)); w = np.zeros((nt, L, L))
  k = 0
  for i in range(S):
    for j in range(i, S):
      xp[:, i, j] = xp[:, j, i] = dv[k:k+nt]; k += nt
  for i in range(S):
    for j in range(i, S):
      xm[:, i, j] = xm[:, j, i] = dv[k:k+nt]; k += nt
  for l in range(L):
    for s in range(S):
      if (l, s) in excl: continue
      gt[:, l, s] = dv[k:k+nt]; k += nt
  for l in range(L):
    w[:, l, l] = dv[k:k+nt]; k += nt
  assert k == len(dv)
  return xp, xm, gt, w, meta

def ratio(p, a, b):
  # every theta both codes computed (lsst_y1: the full vector; roman_real:
  # CoCoA's masked points are zero and stay blank), NaN where DESC-CCL failed
  return np.where(np.isfinite(a) & (b != 0), a/np.where(b == 0, 1.0, b), np.nan)

def band(curves, lo_cut=0.98):
  """Symmetric ratio band around 1 that holds 98% of the plotted points
  (w crosses zero at large theta, where single ratios run off any scale)."""
  r = np.abs(np.concatenate([c[np.isfinite(c)] for c in curves]) - 1.0)
  h = np.quantile(r, lo_cut)*1.15 if r.size else 0.02
  step = 0.002 if h < 0.02 else (0.005 if h < 0.05 else 0.01)
  h = max(step, np.ceil(h/step)*step)
  return [1.0 - h, 1.0 + h]

STYLE = dict(cmap="twilight_shifted", linewidth=[1.0, 1.3, 1.6, 1.9, 2.2, 2.5],
             linestyle=["solid", "dashed", "dashdot", "dotted", (0, (5, 1)), (0, (3, 1, 1, 1))],
             colorbar=None, show=None)

def figures(p, curves, labels, tag):
  ones = unpack(p, np.ones_like(curves[0]))
  meta = ones[4]
  e = np.logspace(np.log10(meta["theta_min"]), np.log10(meta["theta_max"]), meta["ntheta"] + 1)
  th = np.sqrt(e[:-1]*e[1:]); show = [e[0], e[-1]]
  U = [unpack(p, c) for c in curves]
  par = list(range(len(curves)))
  S, L = meta["source_ntomo"], meta["lens_ntomo"]
  big = (16 + 1.2*S, 12 + 1.1*S)
  for pm, nm, k in ((1, "xip", 0), (-1, "xim", 1)):
    yl = band([u[k] for u in U])
    fig, _ = pdv.plot_xi(pm, [(th, u[0], u[1]) for u in U], xi_ref=(th, ones[0], ones[1]), param=par,
                         legend=labels, legendloc=(0.62, 0.62), ylim=yl, thetashow=show, figsize=big,
                         bintextpos=[[0.2, 0.875], [0.2, 0.875]], bintextsize=20, yaxislabelsize=17,
                         yaxisticklabelsize=14, xaxisticklabelsize=17, wspace=0.4, **STYLE)
    fig.savefig(os.path.join(OUT, "%s_%s_%s.png" % (p, tag, nm)), dpi=110); plt.close(fig)
  yl = band([u[2] for u in U])
  fig, _ = pdv.plot_gammat_tomo_limber([(th, u[2]) for u in U], gammat_ref=(th, ones[2]), param=par,
                                       legend=labels, legendloc=(0.9, 0.55), ylim=yl, thetashow=show,
                                       figsize=(16 + 1.2*S, 12 + 1.1*L), bintextpos=[0.85, 0.2], bintextsize=20,
                                       yaxislabelsize=17, yaxisticklabelsize=14, xaxisticklabelsize=17, **STYLE)
  fig.savefig(os.path.join(OUT, "%s_%s_gammat.png" % (p, tag)), dpi=110); plt.close(fig)
  yl = band([u[3] for u in U])
  fig, _ = pdv.plot_wtheta_tomo([(th, u[3]) for u in U], theta_wtheta_ref=(th, ones[3]), param=par,
                                legend=labels, legendloc=(0.9, 0.05), ylim=yl, thetashow=show,
                                figsize=(18, 13./5), bintextpos=[0.85, 0.2], bintextsize=20,
                                yaxislabelsize=17, yaxisticklabelsize=14, xaxisticklabelsize=17, **STYLE)
  fig.savefig(os.path.join(OUT, "%s_%s_w.png" % (p, tag)), dpi=110); plt.close(fig)

for p in PROJ:
  ld = lambda f: np.load(os.path.join(H, f))["dv"]
  cur, lab = [], []
  for m, l in MODELS:
    f = "ccl_%s_%s_ref.npz" % (p, m)
    if os.path.isfile(os.path.join(H, f)):
      cur.append(ratio(p, ld(f), ld("cocoa_%s_%s.npz" % (p, m)))); lab.append(l)
  if cur: figures(p, cur, lab, "cosmologies")
  cur, lab = [], []
  for v, l in VARIANTS:
    f = "ccl_%s_fid_%s.npz" % (p, v)
    if os.path.isfile(os.path.join(H, f)):
      cur.append(ratio(p, ld(f), ld("cocoa_%s_fid.npz" % p))); lab.append(l)
  if cur: figures(p, cur, lab, "choices")
  print("figures for", p)
