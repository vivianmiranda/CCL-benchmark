"""Diagnostic (finding 2): error of the separable linear spectrum
P_lin(k,z) ~ (D(z)/D(z_a))^2 P_lin(k,z_a) against CAMB's P_lin(k,z), for two
anchors: z_a = 0 (DESC-CCL's FKEM term) and z_a = each lens bin's mean
redshift (CoCoA). D is the growth CoCoA uses (k0 = 5e-4/Mpc). For each lens
bin: the largest |P_sep/P_lin - 1| over k = 0.01-0.2/Mpc and the redshifts
holding the central 90% of the bin's n(z).
CMP_WORK=<run folder> python diag_anchor.py"""
import json, os
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
W = os.environ.get("CMP_WORK", H)
for p in ("lsst_y1", "roman_real"):
  ex = np.load(os.path.join(W, "cocoa_%s_fid.npz" % p)); meta = json.loads(str(ex["meta"]))
  zg, kg, lnPL = ex["zg"], ex["kg"], ex["lnPL"]
  base = meta["path"] if os.path.isabs(meta["path"]) else os.path.join(meta["rootdir"], meta["path"])
  d = np.loadtxt(os.path.join(base, meta["dataset"]["nz_lens_file"]))
  zn = d[:, 0] + 0.5*(d[1, 0] - d[0, 0])
  k0 = np.argmin(abs(kg - 5e-4)); ks = (kg >= 0.01) & (kg <= 0.2)
  lnD = 0.5*(lnPL[:, k0] - lnPL[0, k0])
  P = lambda z: np.array([np.interp(z, zg, lnPL[:, i]) for i in range(len(kg))])   # ln P_lin(k, z)
  D = lambda z: np.interp(z, zg, lnD)
  rows = []
  for b in range(meta["lens_ntomo"]):
    nz = d[:, b + 1]; c = np.cumsum(nz)/nz.sum()
    zlo, zhi = np.interp(0.05, c, zn), np.interp(0.95, c, zn); zm = np.sum(zn*nz)/nz.sum()
    zz = np.linspace(zlo, zhi, 25)
    err = {}
    for name, za in (("z=0", 0.0), ("z=zmean", zm)):
      e = [np.max(np.abs(np.exp(P(za)[ks] + 2*(D(z) - D(za)) - P(z)[ks]) - 1)) for z in zz]
      err[name] = max(e)
    rows.append((b + 1, zm, zlo, zhi, err["z=0"], err["z=zmean"]))
  print(p)
  for r in rows:
    print("  lens %d  zmean %.2f  (90%%: %.2f-%.2f)  max|P_sep/P_lin-1|: anchor z=0 %.2e, anchor zmean %.2e" % r)
