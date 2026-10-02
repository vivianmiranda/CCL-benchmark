"""Diagnostic (finding 2): scale dependence of CAMB's linear growth in the
CoCoA exports, D(k,z)/D(k0,z) - 1 with D(k,z) = sqrt(P_lin(k,z)/P_lin(k,0))
and k0 = 5e-4 /Mpc (the growth factor CoCoA uses and DESC-CCL receives).
CMP_WORK=<run folder> python diag_growth.py"""
import json, os
import numpy as np
W = os.environ.get("CMP_WORK", os.path.dirname(os.path.abspath(__file__)))
for p in ("lsst_y1", "roman_real"):
  d = np.load(os.path.join(W, "cocoa_%s_fid.npz" % p)); zg, kg, lnPL = d["zg"], d["kg"], d["lnPL"]
  w = json.loads(str(d["meta"]))["pars"]["w"]
  k0 = np.argmin(abs(kg - 5e-4)); out = []
  for kk in (0.01, 0.05, 0.2):
    i = np.argmin(abs(kg - kk)); row = []
    for zz in (0.5, 1.0, 2.0):
      j = np.argmin(abs(zg - zz))
      row.append(np.exp(0.5*(lnPL[j, i] - lnPL[0, i]))/np.exp(0.5*(lnPL[j, k0] - lnPL[0, k0])) - 1)
    out.append("k=%.2f: %s" % (kg[i], " ".join("%+.4f" % x for x in row)))
  print(p, "w = %.1f" % w, "D(k,z)/D(k0,z) - 1 at z = 0.5, 1, 2 |", " | ".join(out))
