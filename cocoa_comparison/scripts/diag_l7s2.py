# diagnostic (finding 4): sign structure of C_gs for the failing pairs, Limber vs FKEM
# ccl env: CMP_WORK=<run folder> python diag_l7s2.py <project> <lens> <source>
import os, sys, json, numpy as np
H = os.path.dirname(os.path.abspath(__file__))
W = os.environ.get("CMP_WORK", H)
P = sys.argv[1]; L_, S_ = int(sys.argv[2]), int(sys.argv[3])
sys.argv = ["ccl_compute.py", os.path.join(W, "cocoa_%s_fid.npz" % P), "ref", "/dev/null"]
exec(open(os.path.join(H, "ccl_compute.py")).read().split("# ---- harmonic mode")[0])
for nl in (False, True):
  c = cl(lens_tgs[L_], src_t[S_], nl)
  neg = ell[c <= 0]
  print(P, "l%d-s%d" % (L_, S_), "FKEM" if nl else "Limber", "C(2..5)=", np.array2string(c[:4], precision=2),
        "C(lmax)=%.2e" % c[-1], "non-positive at %d of %d l; first/last:" % (len(neg), len(ell)), neg[:3], neg[-3:])
