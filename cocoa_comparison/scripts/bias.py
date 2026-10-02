"""Parameter shift the DESC-CCL - CoCoA difference would cause at the
fiducial point: Fisher matrix F_ab = d_a^T C^-1 d_b (the project's masked
inverse covariance) with derivatives of CoCoA's data vector, shift
F^-1 g with g_a = d_a^T C^-1 (d_CCL - d_CoCoA), reported in units of the
parameter's error sqrt((F^-1)_aa).

Parameters: Omega_m and n_s (central differences of the omm_lo/omm_hi and
ns_lo/ns_hi exports), then, when the run folder has the cocoa_<p>_d_*.npz
exports of run_fisher.sh, the linear bias of every lens bin and the NLA
amplitude and redshift slope (one-sided differences; marginalized over).
Everything else (A_s, h, w, Omega_b, m_nu, photo-z, shear calibration)
stays fixed, so the errors are those of this subspace.
CMP_WORK=<run folder> python bias.py"""
import json, os
import numpy as np
H = os.path.dirname(os.path.abspath(__file__))
W = os.environ.get("CMP_WORK", H)
ld = lambda f: np.load(os.path.join(W, f))
NOTES = []
print("| project | comparison | parameters | $\\Delta\\Omega_m/\\sigma$ | $\\Delta n_s/\\sigma$ | $\\sigma(\\Omega_m)$ | $\\sigma(n_s)$ |")
print("|---|---|---|---|---|---|---|")
for p, P, pre in (("lsst_y1", "LSST-Y1", "LSST_"), ("roman_real", "Roman-Real", "roman_")):
  cov = ld("cov_%s.npz" % p); ic, m = cov["icov"], cov["mask"].astype(bool)
  ex = lambda mod: ld("cocoa_%s_%s.npz" % (p, mod))
  dv = lambda mod: ex(mod)["dv"]
  meta = lambda mod: json.loads(str(ex(mod)["meta"]))
  D = [(dv("omm_hi") - dv("omm_lo"))/(meta("omm_hi")["pars"]["omegam"] - meta("omm_lo")["pars"]["omegam"]),
       (dv("ns_hi") - dv("ns_lo"))/(meta("ns_hi")["pars"]["ns"] - meta("ns_lo")["pars"]["ns"])]
  sets = [("$\\Omega_m$, $n_s$", list(D))]
  L = meta("fid")["lens_ntomo"]; pt = meta("fid")["point"]
  tags = ["b%d" % (i + 1) for i in range(L)] + ["A1", "eta"]
  keys = ["%sB1_%d" % (pre, i + 1) for i in range(L)] + [pre + "A1_1", pre + "A1_2"]
  if all(os.path.isfile(os.path.join(W, "cocoa_%s_d_%s.npz" % (p, t))) for t in tags):
    Dn = []
    for t, k in zip(tags, keys):
      step = meta("d_" + t)["point"][k] - pt[k]
      Dn.append((dv("d_" + t) - dv("fid"))/step)
    sets.append(("+ %d $b_1$, $A_1$, $\\eta$ (marginalized)" % L, D + Dn))
  for lab, f in (("DESC-CCL reference", "ccl_%s_fid_ref.npz" % p),
                 ("DESC-CCL separable $P_{\\rm lin}$ (diagnostic)", "ccl_%s_fid_separable.npz" % p)):
    a, b = ld(f)["dv"], dv("fid")
    d = np.where(np.isfinite(a) & np.isfinite(b), a - b, 0.0)*m
    for plab, DD in sets:
      DD = [x*m for x in DD]
      F = np.array([[x @ ic @ y for y in DD] for x in DD]); Fi = np.linalg.inv(F)
      sig = np.sqrt(np.diag(Fi)); shift = Fi @ np.array([x @ ic @ d for x in DD])
      print("| %s | %s | %s | %+.2f | %+.2f | %.4f | %.4f |" % (P, lab, plab, shift[0]/sig[0], shift[1]/sig[1], sig[0], sig[1]))
      if len(DD) > 2 and lab == "DESC-CCL reference":
        chi2, absorbed = d @ ic @ d, shift @ np.array([x @ ic @ d for x in DD])
        rel = [shift[2 + i]/pt[keys[i]] for i in range(L)]
        NOTES.append("%s: Delta chi2 %.2f, absorbed by the parameter shifts %.2f, left %.2f; lens biases move by %+.2f%% to %+.2f%% (%+.2f to %+.2f sigma)"
                     % (P, chi2, absorbed, chi2 - absorbed, 100*min(rel), 100*max(rel), min(shift[2:2+L]/sig[2:2+L]), max(shift[2:2+L]/sig[2:2+L])))
print()
for n in NOTES: print("- " + n)
