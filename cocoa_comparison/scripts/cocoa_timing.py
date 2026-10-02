"""CoCoA execution time per likelihood evaluation, CAMB excluded (cocoa env).

python cocoa_timing.py <project> <nla|tatt> <rounds>

Cobaya times every component separately: the likelihood component's
time per evaluation, after two warm-up evaluations that are not timed,
is CoCoA's time (the interpolation of CAMB's tables and
cosmolike's data vector and chi2); the CAMB theory component is timed apart
and excluded. The model is the project's EXAMPLE_EVALUATE2 (NLA, or TATT
with IA_model 1 and CFASTPT, IA_code 0, at the TATT point of the study);
every evaluation is a new cosmology, cycling through the study's five
(fiducial, Omega_m 0.25 / 0.35, n_s 0.92 / 1.01), rounds times."""
import os, sys, re
from cobaya.model import get_model
from cobaya.yaml import yaml_load_file

project, ia, rounds = sys.argv[1], sys.argv[2], int(sys.argv[3])
R = os.environ["ROOTDIR"]; os.chdir(R)
info = yaml_load_file("./projects/%s/EXAMPLE_EVALUATE2.yaml" % project)
ov = dict(info["sampler"]["evaluate"]["override"]); info.pop("sampler"); info.pop("output", None)
info["timing"] = True
name = list(info["likelihood"])[0]
lb = info["likelihood"][name]
lb["IA_model"] = 1 if ia == "tatt" else 0
lb["IA_code"] = 0
model = get_model(info)
point = {p: float(ov[p]) for p in model.parameterization.sampled_params() if p in ov}
for p in model.parameterization.sampled_params():
  if p not in point:
    ref = model.info()["params"][p].get("ref")
    point[p] = float(ref["loc"] if isinstance(ref, dict) else ref)
for p in point:
  if re.search(r"_DZ_|_M[0-9]+$|_PM[0-9]+$|_BMAG_", p): point[p] = 0.0
if ia == "tatt":
  pre = "LSST" if project == "lsst_y1" else "roman"
  point.update({pre + "_A1_1": 0.7, pre + "_A1_2": -1.7, pre + "_A2_1": -1.36,
                pre + "_A2_2": -2.5, pre + "_BTA_1": 1.0})
models = [{}, {"omegam": 0.25}, {"omegam": 0.35}, {"ns": 0.92}, {"ns": 1.01}]
import numpy as np
tl = model.likelihood[name].timer; tc = model.theory["camb"].timer
# two warm-up evaluations, not timed: the once-per-run work (tables,
# FFTW plans, caches) happens in the first calls
for m in ({"ns": 0.95}, {"ns": 0.96}):
  model.logposterior(dict(point, **m), cached=False)
lik, cmb = [], []
for r in range(rounds):
  for m in models:
    s0, c0 = tl.time_sum, tc.time_sum
    model.logposterior(dict(point, **m), cached=False)
    lik.append(tl.time_sum - s0); cmb.append(tc.time_sum - c0)
lik = np.array(lik)
print("COCOA %s %s threads %s: likelihood mean %.4f median %.4f min %.4f max %.4f s/eval (n = %d after 2 warm-up calls); CAMB %.3f s/eval (excluded)"
      % (project, ia, os.environ.get("OMP_NUM_THREADS"), lik.mean(), np.median(lik), lik.min(), lik.max(), len(lik), np.mean(cmb)), flush=True)
