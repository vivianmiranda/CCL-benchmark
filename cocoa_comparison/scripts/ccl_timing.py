"""DESC-CCL execution time per data vector from CAMB's tables (ccl env,
pyccl of PR #1296 first on PYTHONPATH).

CMP_WORK=<run folder> python ccl_timing.py <project> <nla|tatt> <rounds> [<variant>]

Runs ccl_compute.py (variant ref, or tatt:ref; or the variant given, e.g.
ccl_numerics, CCL's default FKEM and l sampling) in this one process on the
exported CoCoA tables of the five cosmologies, rounds times after two
warm-up evaluations, and reports the wall time per evaluation: the
CosmologyCalculator built from the CAMB tables, the tracers, the TATT PT
step, every C_l and the full-sky bin-averaged transforms (reading the
tables from disk is timed apart and subtracted)."""
import os, sys, time, contextlib, io
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
W = os.environ.get("CMP_WORK", H)
project, ia, rounds = sys.argv[1], sys.argv[2], int(sys.argv[3])
src = compile(open(os.path.join(H, "ccl_compute.py")).read(), os.path.join(H, "ccl_compute.py"), "exec")
tag = "tatt_" if ia == "tatt" else ""
files = [os.path.join(W, "cocoa_%s_%s%s.npz" % (project, tag, m)) for m in ("fid", "omm_lo", "omm_hi", "ns_lo", "ns_hi")]
base = sys.argv[4] if len(sys.argv) > 4 else "ref"
variant = ("tatt:" if ia == "tatt" else "") + base
tmp_out = os.path.join(W, "_ccl_timing_scratch.npz")   # the short data vector, overwritten
def run(f):
  sys.argv = ["ccl_compute.py", f, variant, tmp_out]
  with contextlib.redirect_stdout(io.StringIO()):
    exec(src, {"__name__": "__main__", "__file__": os.path.join(H, "ccl_compute.py")})
run(files[1]); run(files[2])                    # two warm-up evaluations, not timed
t, tload = [], []
for r in range(rounds):
  for f in files:
    t0 = time.perf_counter(); dict(np.load(f)); tload.append(time.perf_counter() - t0)
    t0 = time.perf_counter(); run(f); t.append(time.perf_counter() - t0)
t = np.array(t) - np.mean(tload)
print("CCL %s %s %s threads %s: mean %.3f median %.3f min %.3f max %.3f s/eval (n = %d after 2 warm-up calls; reading the tables %.3f s, subtracted)"
      % (project, ia, base, os.environ.get("OMP_NUM_THREADS"), t.mean(), np.median(t), t.min(), t.max(), len(t), np.mean(tload)), flush=True)
