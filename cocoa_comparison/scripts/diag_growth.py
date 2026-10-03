"""Diagnostic (Section 4 of the README): scale dependence of CAMB's linear
growth in the CoCoA exports, D(k,z)/D(k0,z) - 1 with D(k,z) =
sqrt(P_lin(k,z)/P_lin(k,0)) and k0 = 5e-4 /Mpc (where CoCoA's likelihood
measures its growth factor, the one DESC-CCL receives). For w != -1, CAMB's
dark-energy perturbations change the growth on horizon scales: the step
between k0 and k ~ 3e-3 /Mpc. Above it only the neutrino slope remains.
Prints every project at both w (the fiducial export and the w swap).
CMP_WORK=<run folder> python diag_growth.py"""
import json, os
import numpy as np
W = os.environ.get("CMP_WORK", os.path.dirname(os.path.abspath(__file__)))
KS = (1e-3, 3e-3, 0.01, 0.05, 0.2)
for f in ("cocoa_lsst_y1_fid.npz", "cocoa_lsst_y1_w_m1.npz", "cocoa_roman_real_w_m09.npz", "cocoa_roman_real_fid.npz"):
  if not os.path.isfile(os.path.join(W, f)): continue
  d = np.load(os.path.join(W, f)); zg, kg, lnPL = d["zg"], d["kg"], d["lnPL"]
  w = json.loads(str(d["meta"]))["pars"]["w"]
  lnP = lambda kk: np.array([np.interp(np.log(kk), np.log(kg), row) for row in lnPL])
  D = lambda kk: np.exp(0.5*(lnP(kk) - lnP(kk)[0]))
  print("%s (w = %.1f): D(k,z)/D(k0,z) - 1 [%%], k = %s /Mpc" % (f, w, ", ".join("%g" % k for k in KS)))
  for zz in (0.5, 1.0, 2.0):
    j = np.argmin(abs(zg - zz))
    print("  z = %.1f: %s" % (zg[j], " ".join("%+.3f" % (100*(D(k)[j]/D(5e-4)[j] - 1)) for k in KS)))
