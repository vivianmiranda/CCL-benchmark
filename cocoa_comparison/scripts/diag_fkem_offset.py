# diagnostic: C_gs below l = 150, CCL FKEM with CAMB's P_lin(k,z) (ref) vs
# the separable table P_lin(k,0) D^2 (diagnostic), and the Limber C_gs
# ccl env: CMP_WORK=<run folder> python diag_fkem_offset.py <project>
import os, sys, numpy as np
H = os.path.dirname(os.path.abspath(__file__))
W = os.environ.get("CMP_WORK", H)
P = sys.argv[1]
res = {}
for var in ("ref", "separable"):
  sys.argv = ["ccl_compute.py", os.path.join(W, "cocoa_%s_fid.npz" % P), var, "/dev/null"]
  exec(open(os.path.join(H, "ccl_compute.py")).read().split("# ---- harmonic mode")[0])
  e = np.array([10., 20., 50., 100., 140.])
  pairs = [(0, S - 1), (1, S - 1), (L - 2, S - 1)]
  res[var] = {pq: ccl.angular_cl(cosmo, lens_tgs[pq[0]], src_t[pq[1]], e, l_limber=150, non_limber_integration_method="FKEM", fkem_Nchi=2000) for pq in pairs}
  res[var + "_lim"] = {pq: ccl.angular_cl(cosmo, lens_tgs[pq[0]], src_t[pq[1]], e) for pq in pairs}
for pq in res["ref"]:
  print(P, "lens %d source %d" % pq, "l =", e.astype(int),
        "FKEM ref/separable - 1:", np.round(res["ref"][pq]/res["separable"][pq] - 1, 4),
        "| Limber ref/separable - 1:", np.round(res["ref_lim"][pq]/res["separable_lim"][pq] - 1, 4))
