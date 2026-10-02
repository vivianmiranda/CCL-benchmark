"""Limber C_l of both codes at CoCoA's exported multipoles (20 to the
real-space lmax):
|DESC-CCL / CoCoA - 1| per probe, median over bin pairs and maximum over
the pairs that carry signal (|C_l| above 1% of the probe's largest |C_l| at
that l; this drops lens-behind-source gamma_t pairs, whose C_l are 1e-4 to
1e-10 of the others), in four l ranges. Then the non-Limber C_gg of both
codes (CoCoA's C_gg_tomo, DESC-CCL's FKEM at the reference settings) below
l = 150: DESC-CCL / CoCoA - 1, lowest and highest over the lens bins.
CMP_WORK=<run folder> python harmonic.py"""
import os
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
W = os.environ.get("CMP_WORK", H)
RANGES = ((20, 66), (66, 1000), (1000, 5001), (5001, 100001))
print("| project | probe | " + " | ".join("$%d \\le \\ell < %d$" % r for r in RANGES) + " |")
print("|---|---|" + "---|"*len(RANGES))
for p, P in (("lsst_y1", "LSST-Y1"), ("roman_real", "Roman-Real")):
  co = np.load(os.path.join(W, "cocoa_%s_fid.npz" % p)); cc = np.load(os.path.join(W, "ccl_%s_fid_harmonic.npz" % p))
  ell = co["ell_h"]
  ef = os.path.join(H, "ggl_exclude_%s.txt" % p)
  excl = {tuple(int(x) for x in l.split()) for l in open(ef) if l.strip()} if os.path.isfile(ef) else set()
  S = cc["css"].shape[0]; L = cc["cgs"].shape[0]
  probes = {
    "$C_{ss}$": (np.array([cc["css"][i, j] for i in range(S) for j in range(i, S)]),
                 np.array([co["cl_ss"][0][:, i, j] for i in range(S) for j in range(i, S)])),
    "$C_{gs}$": (np.array([cc["cgs"][l, k] for l in range(L) for k in range(S) if (l, k) not in excl]),
                 np.array([co["cl_gs"][:, l, k] for l in range(L) for k in range(S) if (l, k) not in excl])),
    "$C_{gg}$": (np.array([cc["cgg"][l] for l in range(L)]), np.array([co["cl_gg"][:, l, l] for l in range(L)])),
  }
  for nm, (a, b) in probes.items():
    r = np.abs(a/b - 1.0)
    sig = np.abs(b) > 0.01*np.abs(b).max(axis=0)[None, :]
    cells = []
    for lo, hi in RANGES:
      s = (ell >= lo) & (ell < hi)
      cells.append("%.1e / %.1e" % (np.median(r[:, s]), r[:, s][sig[:, s]].max()))
    print("| %s | %s | %s |" % (P, nm, " | ".join(cells)))

print()
NL = (2.0, 10.0, 50.0, 140.0)
print("| project | " + " | ".join("$\\ell = %d$" % l for l in NL) + " |")
print("|---|" + "---|"*len(NL))
for p, P in (("lsst_y1", "LSST-Y1"), ("roman_real", "Roman-Real")):
  co = np.load(os.path.join(W, "cocoa_%s_fid.npz" % p)); cc = np.load(os.path.join(W, "ccl_%s_fid_harmonic.npz" % p))
  if "cgg_nl" not in cc.files: continue
  enl = co["ell_nl"]; L = cc["cgg_nl"].shape[0]
  r = np.array([cc["cgg_nl"][l]/co["cl_gg_nl"][:, l, l] - 1.0 for l in range(L)])   # (L, n_ell)
  cells = ["%+.2f%% to %+.2f%%" % (100*r[:, list(enl).index(l)].min(), 100*r[:, list(enl).index(l)].max()) for l in NL]
  print("| %s | %s |" % (P, " | ".join(cells)))
