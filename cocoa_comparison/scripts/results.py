"""Every Delta chi2 of the CoCoA vs DESC-CCL study, as Markdown tables
(results.md) and JSON (results.json). Delta chi2 = d^T C^-1 d with each
project's masked inverse covariance (shipped mask), per probe and in total.
CMP_WORK=<run folder> python results.py"""
import json, os
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
# the run outputs (cocoa_*.npz, ccl_*.npz, cov_*.npz): CMP_WORK, or this folder
W = os.environ.get("CMP_WORK", H)
PROJ = {"lsst_y1": "LSST-Y1", "roman_real": "Roman-Real"}
MODELS = [("fid", "fiducial"), ("omm_lo", "$\\Omega_m = 0.25$"), ("omm_hi", "$\\Omega_m = 0.35$"),
          ("ns_lo", "$n_s = 0.92$"), ("ns_hi", "$n_s = 1.01$")]
def load(f):
  p = os.path.join(W, f); return np.load(p)["dv"] if os.path.isfile(p) else None
def dchi2(p, a, b):
  cov = np.load(os.path.join(W, "cov_%s.npz" % p)); ic, m = cov["icov"], cov["mask"].astype(bool)
  s = [int(x) for x in cov["sizes"]]; e = np.cumsum([0, s[0], s[1], s[2]])
  bad = ~(np.isfinite(a) & np.isfinite(b))    # pairs a code could not compute: left out
  d = np.where(bad, 0.0, a - b)*m
  q = lambda x: float(x @ ic @ x); out = {"3x2pt": q(d), "left_out_kept": int((bad & m).sum())}
  for k, nm in enumerate(("shear", "gammat", "w")):
    x = np.zeros_like(d); x[e[k]:e[k+1]] = d[e[k]:e[k+1]]; out[nm] = q(x)
  return out
R = {}; md = []
fm = lambda x: "%.3f" % x if abs(x) >= 0.01 else ("%.1e" % x if x != 0 else "0")
row = lambda lab, r: "| %s | %s | %s | %s | %s |" % (lab, fm(r["3x2pt"]), fm(r["shear"]), fm(r["gammat"]), fm(r["w"])) + " %d |" % r.get("left_out_kept", 0)
head = "| %s | 3x2pt | shear | $\\gamma_t$ | $w(\\theta)$ | points left out |\n|---|---|---|---|---|---|"
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
           ("DESC-CCL, separable $P_{\\rm lin}$ (diagnostic), vs CoCoA", "ccl_%s_fid_separable.npz" % p, "cocoa_%s_fid.npz" % p),
           ("RSD in $\\gamma_t$: effect in DESC-CCL", "ccl_%s_fid_rsd_gs.npz" % p, "ccl_%s_fid_ref.npz" % p),
           ("DESC-CCL finer sampling (ccl_hi) vs reference", "ccl_%s_fid_ccl_hi.npz" % p, "ccl_%s_fid_ref.npz" % p),
           ("DESC-CCL transform of CoCoA's Limber $C_\\ell$ vs CoCoA in Limber", "ccl_%s_fid_cocoa_cl.npz" % p, "cocoa_%s_fidlimber.npz" % p),
           ("DESC-CCL Eisenstein-Hu + halofit vs reference", "ccl_%s_fid_eh.npz" % p, "ccl_%s_fid_ref.npz" % p),
           ("DESC-CCL, growth factor at k = 0.05/Mpc (diagnostic), vs CoCoA", "ccl_%s_fid_growth_sub.npz" % p, "cocoa_%s_fid.npz" % p)]
  # each project also at the other project's w (LSST-Y1 fiducial w = -0.9, Roman-Real w = -1)
  ws, wl = ("w_m1", "$w = -1$") if p == "lsst_y1" else ("w_m09", "$w = -0.9$")
  pairs += [("DESC-CCL vs CoCoA, %s" % wl, "ccl_%s_%s_ref.npz" % (p, ws), "cocoa_%s_%s.npz" % (p, ws)),
            ("DESC-CCL, growth factor at k = 0.05/Mpc (diagnostic), vs CoCoA, %s" % wl,
             "ccl_%s_%s_growth_sub.npz" % (p, ws), "cocoa_%s_%s.npz" % (p, ws))]
  for v, lab in (("limgs", "$C_{gs}$ in Limber"), ("norsd", "no RSD"), ("rsd_gs", "RSD also in $\\gamma_t$"), ("flat", "flat-sky FFTLog, bin centers"),
                 ("points", "full-sky, bin centers"), ("bench", "the CCL-benchmark modeling")):
    pairs.append(("DESC-CCL %s vs CoCoA" % lab, "ccl_%s_fid_%s.npz" % (p, v), "cocoa_%s_fid.npz" % p))
  for lab, fa, fb in pairs:
    a, b = load(fa), load(fb)
    if a is not None and b is not None:
      R[p][lab] = dchi2(p, a, b); md.append(row(lab, R[p][lab]))
  md.append("")
  # TATT (IA_model 1, CFASTPT in CoCoA; pyccl EulerianPTCalculator in DESC-CCL)
  ta = {m: load("cocoa_%s_tatt_%s.npz" % (p, m)) for m, _ in MODELS}
  if any(x is not None for x in ta.values()):
    md.append(head % "TATT: DESC-CCL (reference settings) vs CoCoA")
    for m, lab in MODELS:
      a = load("ccl_%s_tatt_%s_ref.npz" % (p, m))
      if a is not None and ta[m] is not None:
        R[p]["tatt_ref_" + m] = dchi2(p, a, ta[m]); md.append(row(lab, R[p]["tatt_ref_" + m]))
    md.append("")
    md.append(head % "TATT, fiducial, one choice at a time")
    for lab, fa, fb in (
        ("TATT vs NLA in CoCoA (size of the TATT terms)", "cocoa_%s_tatt_fid.npz" % p, "cocoa_%s_fid.npz" % p),
        ("CoCoA default vs CoCoA high accuracy (TATT)", "cocoa_%s_tatt_fid.npz" % p, "cocoa_%s_tatt_fidhi.npz" % p),
        ("DESC-CCL PT tables at 320 k per decade vs 160 (TATT)", "ccl_%s_tatt_fid_pt_hi.npz" % p, "ccl_%s_tatt_fid_ref.npz" % p),
        ("both in Limber: DESC-CCL vs CoCoA (TATT)", "ccl_%s_tatt_fid_limber.npz" % p, "cocoa_%s_tatt_fidlimber.npz" % p),
        ("DESC-CCL, separable $P_{\\rm lin}$ (diagnostic), vs CoCoA (TATT)", "ccl_%s_tatt_fid_separable.npz" % p, "cocoa_%s_tatt_fid.npz" % p)):
      a, b = load(fa), load(fb)
      if a is not None and b is not None:
        R[p][lab] = dchi2(p, a, b); md.append(row(lab, R[p][lab]))
    md.append("")
  # CCL calls that failed (pair left out), per run
  fl = []
  for f in sorted(os.listdir(W)):
    if f.startswith("ccl_%s_" % p) and f.endswith(".npz"):
      z = np.load(os.path.join(W, f))
      for x in (json.loads(str(z["failures"])) if "failures" in z.files else []):
        fl.append("| %s | %s | %s | %s | %s |" % (f[len("ccl_%s_" % p):-4], x["pair"], x["type"], x.get("stage", "correlation"), x["error"]))
  R[p]["failures"] = fl
  if fl:
    md.append("| DESC-CCL run | pair | type | call | error |\n|---|---|---|---|---|")
    md += fl; md.append("")
open(os.path.join(W, "results.md"), "w").write("\n".join(md) + "\n")
json.dump(R, open(os.path.join(W, "results.json"), "w"), indent=1)
print("\n".join(md))
