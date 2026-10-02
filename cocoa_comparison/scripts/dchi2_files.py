"""Delta chi2 between any two data-vector files of one project, per probe.
python dchi2_files.py <project> <a.npz> <b.npz>   (key "dv" in both)"""
import os, sys
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
project, fa, fb = sys.argv[1:4]
cov = np.load(os.path.join(H, "cov_%s.npz" % project)); ic, mask = cov["icov"], cov["mask"].astype(bool)
s = [int(x) for x in cov["sizes"]]
a, b = np.load(os.path.join(H, fa))["dv"], np.load(os.path.join(H, fb))["dv"]
d = np.where(np.isfinite(a) & np.isfinite(b), a - b, 0.0)*mask
e = np.cumsum([0, s[0], s[1], s[2]])
q = lambda x: float(x @ ic @ x)
p = []
for i in range(3):
  x = np.zeros_like(d); x[e[i]:e[i+1]] = d[e[i]:e[i+1]]; p.append(q(x))
print("%-34s vs %-34s total %8.4f | shear %7.4f gammat %7.4f w %7.4f" % (fa, fb, q(d), *p))
