"""Every Delta chi2 of the CoCoA vs DESC-CCL study, as Markdown tables
(results.md) and JSON (results.json). Delta chi2 = d^T C^-1 d with each
project's masked inverse covariance (shipped mask), per probe and in total.
python results.py"""
import json, os
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
PROJ = {"lsst_y1": "LSST-Y1", "roman_real": "Roman-Real"}
MODELS = [("fid", "fiducial"), ("omm_lo", "$\\Omega_m = 0.25$"), ("omm_hi", "$\\Omega_m = 0.35$"),
          ("ns_lo", "$n_s = 0.92$"), ("ns_hi", "$n_s = 1.01$")]
def load(f):
  p = os.path.join(H, f); return np.load(p)["dv"] if os.path.isfile(p) else None
def dchi2(p, a, b):
  cov = np.load(os.path.join(H, "cov_%s.npz" % p)); ic, m = cov["icov"], cov["mask"].astype(bool)
  s = [int(x) for x in cov["sizes"]]; e = np.cumsum([0, s[0], s[1], s[2]])
  bad = ~(np.isfinite(a) & np.isfinite(b))    # pairs a code could not compute: left out
  d = np.where(bad, 0.0, a - b)*m
  q = lambda x: float(x @ ic @ x); out = {"3x2pt": q(d), "left_out_kept": int((bad & m).sum())}
  for k, nm in enumerate(("shear", "gammat", "w")):
    x = np.zeros_like(d); x[e[k]:e[k+1]] = d[e[k]:e[k+1]]; out[nm] = q(x)
  return out
R = {}; md = []
row = lambda lab, r: "| %s | %.3f | %.4f | %.3f | %.3f |" % (lab, r["3x2pt"], r["shear"], r["gammat"], r["w"]) + (" (%d kept points left out)" % r["left_out_kept"] if r.get("left_out_kept") else "")
head = "| %s | 3x2pt | shear | $\\gamma_t$ | $w(\\theta)$ |\n|---|---|---|---|---|"
for p, P in PROJ.items():
  R[p] = {}
  co = {m: load("cocoa_%s_%s.npz" % (p, m)) for m, _ in MODELS}
  md.append("### %s\n" % P)
  md.append(head % "DESC-CCL (reference settings) vs CoCoA")
  for m, lab in MODELS:
    a = load("ccl_%s_%s_ref.npz" % (p, m))
    if a is not None and co[m] is not None:
      R[p]["ref_" + m] = dchi2(p, a, co[m]); md.append(row(lab, R[p]["ref_" + m]))
  md.append("")
  md.append(head % "Fiducial, one choice at a time")
  ref = load("ccl_%s_fid_ref.npz" % p)
  pairs = [("CoCoA default vs CoCoA high accuracy", "cocoa_%s_fid.npz" % p, "cocoa_%s_fidhi.npz" % p),
           ("DESC-CCL default FKEM/$\\ell$ sampling vs reference", "ccl_%s_fid_ccl_numerics.npz" % p, "ccl_%s_fid_ref.npz" % p),
           ("both in Limber: DESC-CCL vs CoCoA", "ccl_%s_fid_limber.npz" % p, "cocoa_%s_fidlimber.npz" % p),
           ("non-Limber effect in CoCoA", "cocoa_%s_fid.npz" % p, "cocoa_%s_fidlimber.npz" % p),
           ("non-Limber effect in DESC-CCL", "ccl_%s_fid_ref.npz" % p, "ccl_%s_fid_limber.npz" % p),
           ("non-Limber effect in DESC-CCL, separable $P_{\\rm lin}$ (diagnostic)", "ccl_%s_fid_separable.npz" % p, "ccl_%s_fid_limber.npz" % p),
           ("DESC-CCL, separable $P_{\\rm lin}$ (diagnostic), vs CoCoA", "ccl_%s_fid_separable.npz" % p, "cocoa_%s_fid.npz" % p)]
  for v, lab in (("limgs", "$C_{gs}$ in Limber"), ("norsd", "no RSD"), ("flat", "flat-sky FFTLog, bin centers"),
                 ("points", "full-sky, bin centers"), ("bench", "the CCL-benchmark modeling")):
    pairs.append(("DESC-CCL %s vs CoCoA" % lab, "ccl_%s_fid_%s.npz" % (p, v), "cocoa_%s_fid.npz" % p))
  for lab, fa, fb in pairs:
    a, b = load(fa), load(fb)
    if a is not None and b is not None:
      R[p][lab] = dchi2(p, a, b); md.append(row(lab, R[p][lab]))
  md.append("")
open(os.path.join(H, "results.md"), "w").write("\n".join(md) + "\n")
json.dump(R, open(os.path.join(H, "results.json"), "w"), indent=1)
print("\n".join(md))
