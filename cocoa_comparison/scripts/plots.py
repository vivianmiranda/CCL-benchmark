"""Figures of the CoCoA vs DESC-CCL study, drawn with the data-vector
plotters of cosmolike_core (cosmolike_notebook_utils/plot_datavectors.py:
plot_xi, plot_gammat_tomo_limber, plot_wtheta_tomo) in their ratio mode.
Every curve is (DESC-CCL - CoCoA)/sigma at the same point, sigma = sqrt(C_ii)
of the project's covariance (Gaussian + non-Gaussian columns, as CoCoA reads
them): 1 means a one-sigma difference. The curves enter as 1 + that and the
reference is 1 (0 for pairs a project excludes). Points removed by the
project's mask (Roman-Real), and pairs DESC-CCL could not transform, are
left blank. Each row of panels has its own y-range fitted to its data
(row_ranges); a panel whose curves leave it shows alpha times the
difference, with 1/alpha printed in the panel (1/alpha = 3: the difference
is three times what the axis reads).
CMP_WORK=<run folder> python plots.py <outdir>   (cocoa env: matplotlib, ROOTDIR set)"""
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
matplotlib.rcParams['legend.fontsize'] = 22
matplotlib.rcParams['savefig.bbox'] = 'tight'
H = os.path.dirname(os.path.abspath(__file__))
# the run outputs (cocoa_*.npz, ccl_*.npz, cov_*.npz): CMP_WORK, or this folder
W = os.environ.get("CMP_WORK", H)
OUT = sys.argv[1]; os.makedirs(OUT, exist_ok=True)
pd_path = os.path.join(os.environ["ROOTDIR"], "external_modules/code/cosmolike_core/cosmolike_notebook_utils/plot_datavectors.py")
spec = importlib.util.spec_from_file_location("plot_datavectors", pd_path)
pdv = importlib.util.module_from_spec(spec); spec.loader.exec_module(pdv)

PROJ = {"lsst_y1": "LSST-Y1", "roman_real": "Roman-Real"}
MODELS = [("fid", "fiducial"), ("omm_lo", r"$\Omega_m=0.25$"), ("omm_hi", r"$\Omega_m=0.35$"),
          ("ns_lo", r"$n_s=0.92$"), ("ns_hi", r"$n_s=1.01$")]
def unpack(p, dv):
  """CoCoA's data-vector layout -> xi+ / xi- (ntheta, S, S), gamma_t
  (ntheta, L, S), w (ntheta, L, L)."""
  meta = json.loads(str(np.load(os.path.join(W, "cocoa_%s_fid.npz" % p))["meta"]))
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

def sigma(p):
  """sqrt(C_ii) for every point of the project's data vector, from its
  covariance file (10 columns: i j ... cov_g cov_ng; CoCoA sums the last
  two), cached as sigma_<project>.npy in the run folder."""
  f = os.path.join(W, "sigma_%s.npy" % p)
  if os.path.isfile(f): return np.load(f)
  ex = np.load(os.path.join(W, "cocoa_%s_fid.npz" % p)); meta = json.loads(str(ex["meta"]))
  base = meta["path"] if os.path.isabs(meta["path"]) else os.path.join(meta["rootdir"], meta["path"])
  c = np.full(len(ex["dv"]), np.nan)
  with open(os.path.join(base, meta["dataset"]["cov_file"])) as fin:
    for line in fin:
      t = line.split()
      if len(t) >= 3 and not t[0].startswith("#") and t[0] == t[1]:
        c[int(t[0])] = float(t[8]) + float(t[9]) if len(t) == 10 else float(t[2])
  np.save(f, np.sqrt(c)); return np.sqrt(c)

def nsigma(p, a, b):
  # 1 + (DESC-CCL - CoCoA)/sigma on every theta both codes computed
  # (lsst_y1: the full vector; roman_real: CoCoA's masked points are zero
  # and stay blank), NaN where DESC-CCL failed
  return np.where(np.isfinite(a) & (b != 0), 1.0 + (a - b)/sigma(p), np.nan)

# alpha = 1/f with f from this ladder (the panel prints f = 1/alpha)
ALPHAS = [1.0/f for f in (2, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500, 1000, 2000, 5000, 10000)]

def row_ranges(R, rows):
  """One y-range per row of panels, fitted to the data (not symmetric):
  from the lowest to the highest value of the row's typical panels, 0
  included, 10% padding. Typical = the panels whose spread is within 2.5
  times the median spread of the panels holding at least 10% of the row's
  largest. Any panel that leaves its row's range gets
  alpha = 1/f, f the first of the ladder that brings it inside.
  R: (ncurve, ntheta, n1, n2), 1 + (DESC-CCL - CoCoA)/sigma; rows: list of
  lists of (i, j)."""
  D = np.where(np.isfinite(R), R - 1.0, np.nan)
  # a pair the project excludes arrives as zeros (unpack): no range, no alpha
  D[:, :, np.all(R == 0, axis=(0, 1))] = np.nan
  with np.errstate(all="ignore"):
    lo_p, hi_p = np.nanmin(D, axis=(0, 1)), np.nanmax(D, axis=(0, 1))
  a = np.ones(lo_p.shape); lims = []
  for row in rows:
    core = [ij for ij in row if np.isfinite(lo_p[ij])]
    if not core:
      lims.append(None); continue
    spread = np.array([max(-lo_p[ij], hi_p[ij], 0.0) for ij in core])
    # the typical spread: median over the panels with at least 10% of the
    # row's largest (panels near zero, such as lens-behind-source gamma_t,
    # do not set the range)
    typical = np.median(spread[spread >= 0.1*spread.max()])
    keep = [ij for ij, x in zip(core, spread) if x <= 2.5*typical] or core
    lo = min(0.0, min(lo_p[ij] for ij in keep)); hi = max(0.0, max(hi_p[ij] for ij in keep))
    if hi - lo <= 0: lo, hi = -0.1, 0.1
    pad = 0.10*(hi - lo); lo, hi = lo - pad, hi + pad
    lims.append((lo, hi))
    for ij in row:
      if not np.isfinite(lo_p[ij]) or (lo_p[ij] >= lo and hi_p[ij] <= hi): continue
      ok = [x for x in ALPHAS if lo_p[ij]*x >= lo and hi_p[ij]*x <= hi]
      a[ij] = ok[0] if ok else ALPHAS[-1]
  return a, lims

def bands(lims):
  # the plotters' ratio-mode ylim: one [1 + lo, 1 + hi] band per row
  return [[1.0 + lo, 1.0 + hi] for lo, hi in (x if x is not None else (-0.1, 0.1) for x in lims)]

def set_rows(axrows, lims, ylabel):
  from matplotlib.ticker import MaxNLocator
  for axs, lim in zip(axrows, lims):
    axs[0].set_ylabel(ylabel, fontsize=24)
    if lim is None: continue
    # ticks at least 8% inside the row's range, so the labels of adjacent
    # (glued) rows never touch
    span = lim[1] - lim[0]
    ticks = [t for t in MaxNLocator(nbins=5).tick_values(lim[0], lim[1])
             if lim[0] + 0.08*span <= t <= lim[1] - 0.08*span]
    for ax in axs:
      ax.set_yticks(ticks)
      ax.tick_params(axis="y", labelsize=19)
      ax.tick_params(axis="x", labelsize=22)
      if ax.get_xlabel(): ax.set_xlabel(r"$\theta$ [arcmin]", fontsize=24)

def pct(lim):
  return lim

def scaled(r, a):
  # the plotters draw curve/ref - 1 = alpha (DESC-CCL - CoCoA)/sigma
  return 1.0 + (r - 1.0)*a[None, :, :]

def top_legend(fig, axes):
  """Move the plotter's legend above the panels, all entries in one row."""
  leg = fig.legends[0]
  hs = getattr(leg, "legend_handles", None) or leg.legendHandles
  labs = [t.get_text() for t in leg.get_texts()]
  fig.legends.remove(leg)
  boxes = [x.get_position() for x in np.ravel(axes) if x.axison]
  x0 = min(b.x0 for b in boxes); x1 = max(b.x1 for b in boxes); y1 = max(b.y1 for b in boxes)
  fig.legend(hs, labs, loc="lower center", bbox_to_anchor=(0.5*(x0 + x1), y1), ncols=len(labs),
             frameon=False, handlelength=2.2, columnspacing=1.4, borderaxespad=0.3)

def mark(ax, a, vals, lim, th, show):
  """Bin label (the plotter's) and 1/alpha in the left corner, top or bottom,
  that the curves leave free over the first 35% of the log-theta axis.
  vals: (ncurve, ntheta) plotted values; lim: the row's y-range."""
  x = np.log(th/show[0])/np.log(show[1]/show[0])
  v = vals[:, x < 0.35]; v = v[np.isfinite(v)]
  f = (v - lim[0])/(lim[1] - lim[0]) if (v.size and lim is not None) else np.array([0.5])
  top = f.max() < 0.6 or (f.min() <= 0.4 and 1.0 - f.max() >= f.min())
  ybin, ya = (0.86, 0.70) if top else (0.12, 0.28)
  for t in ax.texts:
    if t.get_text().strip().startswith("$(") or t.get_text().strip().startswith("("):
      t.set_position((0.04, ybin)); t.set_horizontalalignment("left"); t.set_verticalalignment("center")
  if a != 1.0:
    ax.text(0.04, ya, r"$1/\alpha = %d$" % round(1.0/a), fontsize=20, transform=ax.transAxes,
            horizontalalignment="left", verticalalignment="center")

def style(n):
  """Colors from twilight_shifted without its pale middle (positions
  0.05-0.35 and 0.65-0.95, cool and warm alternating), none lighter than
  luminance 0.5, as a ListedColormap the plotters index exactly (curve x
  gets cm(x/n)). Curve 0 (the fiducial) is solid; the others get line
  styles by luminance, the lightest the long dash and the darkest the
  dots; lighter colors and dotted lines get wider lines."""
  from matplotlib.colors import ListedColormap
  cool, warm = np.linspace(0.05, 0.35, (n + 1)//2), np.linspace(0.95, 0.65, n//2)
  pos = np.array([cool[x//2] if x % 2 == 0 else warm[x//2] for x in range(n)])
  cols = matplotlib.colormaps["twilight_shifted"](pos)
  lum = cols[:, :3] @ np.array([0.299, 0.587, 0.114])
  cols[:, :3] *= np.minimum(1.0, 0.5/lum)[:, None]; lum = np.minimum(lum, 0.5)
  ls, lw = ["solid"]*n, [2.4]*n
  for r, x in enumerate(sorted(range(1, n), key=lambda x: -lum[x])):   # lightest first
    ls[x] = VISIBLE[r]
    lw[x] = 2.2 + 2.0*lum[x] + (0.8 if VISIBLE[r] == "dotted" else 0.0)
  return dict(cmap=ListedColormap(cols, N=n), linewidth=lw, linestyle=ls, colorbar=None, show=None)

# from the most to the least visible
VISIBLE = [(0, (9, 3)), "dashed", "dashdot", (0, (3, 1, 1, 1)), "dotted", (0, (1, 1))]

def figures(p, curves, labels, tag):
  ones = unpack(p, np.ones_like(curves[0]))
  meta = ones[4]
  e = np.logspace(np.log10(meta["theta_min"]), np.log10(meta["theta_max"]), meta["ntheta"] + 1)
  th = np.sqrt(e[:-1]*e[1:]); show = [e[0], e[-1]]
  U = [unpack(p, c) for c in curves]
  par = list(range(len(curves)))
  S, L = meta["source_ntomo"], meta["lens_ntomo"]
  def prep(k, rows):
    R = np.array([u[k] for u in U])
    a, lims = row_ranges(R, rows)
    D = np.where(np.isfinite(R) & ~np.all(R == 0, axis=(0, 1))[None, None], np.abs(R - 1.0), np.nan)
    print("  %s %s: max |Delta/sigma| %.3f; 1/alpha != 1 in %d panels"
          % (p, "xip xim gammat w".split()[k], np.nanmax(D), int((a != 1).sum())))
    return [scaled(r, a) for r in R], a, lims
  for pm, nm, k in ((1, "xip", 0), (-1, "xim", 1)):
    rows = [[(i, j) for i in range(j + 1)] for j in range(S)]
    R, a, lims = prep(k, rows)
    fig, axes = pdv.plot_xi(pm, [(th, r, r) for r in R], xi_ref=(th, ones[0], ones[1]), param=par,
                            legend=labels, legendloc=(0.62, 0.62), ylim=bands(lims), thetashow=show,
                            figsize=(16 + 1.2*S, 12 + 1.1*S), bintextpos=[[0.1, 0.85], [0.1, 0.85]],
                            bintextsize=22, yaxislabelsize=24, yaxisticklabelsize=19, xaxisticklabelsize=22,
                            wspace=0.4, **style(len(curves)))
    set_rows([[axes[j, i] for i in range(j + 1)] for j in range(S)], lims,
             r"$\Delta\xi_{%s}/\sigma$" % ("+" if pm > 0 else "-"))
    for i in range(S):
      for j in range(i, S):
        mark(axes[j, i], a[i, j], np.array([r[:, i, j] for r in R]) - 1.0, pct(lims[j]), th, show)
    top_legend(fig, axes)
    fig.savefig(os.path.join(OUT, "%s_%s_%s.png" % (p, tag, nm)), dpi=180, bbox_inches="tight", bbox_extra_artists=fig.legends); plt.close(fig)
  R, a, lims = prep(2, [[(i, j) for i in range(L)] for j in range(S)])
  fig, axes = pdv.plot_gammat_tomo_limber([(th, r) for r in R], gammat_ref=(th, ones[2]), param=par,
                                          legend=labels, legendloc=(0.92, 0.40), ylim=bands(lims), thetashow=show,
                                          figsize=(16 + 1.2*S, 12 + 1.1*L), bintextpos=[0.1, 0.85], bintextsize=22,
                                          yaxislabelsize=24, yaxisticklabelsize=19, xaxisticklabelsize=22, **style(len(curves)))
  set_rows([[axes[j, i] for i in range(L)] for j in range(S)], lims, r"$\Delta\gamma_t/\sigma$")
  for i in range(L):
    for j in range(S):
      mark(axes[j, i], a[i, j], np.array([r[:, i, j] for r in R]) - 1.0, pct(lims[j]), th, show)
  top_legend(fig, axes)
  fig.savefig(os.path.join(OUT, "%s_%s_gammat.png" % (p, tag)), dpi=180, bbox_inches="tight", bbox_extra_artists=fig.legends); plt.close(fig)
  R, a, lims = prep(3, [[(i, i) for i in range(L)]])
  fig, axes = pdv.plot_wtheta_tomo([(th, r) for r in R], theta_wtheta_ref=(th, ones[3]), param=par,
                                   legend=labels, legendloc=(0.915, 0.05), ylim=bands(lims)[0], thetashow=show,
                                   figsize=(18, 13./5), bintextpos=[0.1, 0.85], bintextsize=22,
                                   yaxislabelsize=24, yaxisticklabelsize=19, xaxisticklabelsize=22, **style(len(curves)))
  set_rows([[axes[i] for i in range(L)]], lims, r"$\Delta w/\sigma$")
  for i in range(L): mark(axes[i], a[i, i], np.array([r[:, i, i] for r in R]) - 1.0, pct(lims[0]), th, show)
  top_legend(fig, axes)
  fig.savefig(os.path.join(OUT, "%s_%s_w.png" % (p, tag)), dpi=180, bbox_inches="tight", bbox_extra_artists=fig.legends); plt.close(fig)

for p in PROJ:
  ld = lambda f: np.load(os.path.join(W, f))["dv"]
  cur, lab = [], []
  for m, l in MODELS:
    f = "ccl_%s_%s_ref.npz" % (p, m)
    if os.path.isfile(os.path.join(W, f)):
      cur.append(nsigma(p, ld(f), ld("cocoa_%s_%s.npz" % (p, m)))); lab.append(l)
  if cur: figures(p, cur, lab, "cosmologies")
  print("figures for", p)
