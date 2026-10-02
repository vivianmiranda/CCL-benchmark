#!/bin/bash
# CoCoA side of the comparison: every export the tables and figures need.
# Run from a shell with the cocoa environment active and start_cocoa.sh
# sourced (ROOTDIR set). Usage: bash run_cocoa.sh <run folder>
[ -n "${ROOTDIR:-}" ] || { echo "source start_cocoa.sh first"; return 1 2>/dev/null || exit 1; }
mkdir -p "${1:?run folder}"; W="$(cd "${1}" && pwd)"
S="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"
ex() { python "${S:?}/cocoa_export.py" "$@" > "${3%.npz}.log" 2>&1 || echo "FAILED: $*"; }
for m in fid omm_lo omm_hi ns_lo ns_hi; do
  # lsst_y1: the full vector (a scratch dataset with ones.mask);
  # roman_real: its covariance is not positive definite with every point,
  # so the shipped mask is used and masked points stay blank
  EXPORT_MODE=dv  ex lsst_y1    ${m} "${W}/cocoa_lsst_y1_${m}.npz"
  EXPORT_MODE=cov ex roman_real ${m} "${W}/cocoa_roman_real_${m}.npz"
done
# the masked inverse covariances (shipped masks)
EXPORT_MODE=cov ex lsst_y1    fid "${W}/cov_lsst_y1.npz"
EXPORT_MODE=cov ex roman_real fid "${W}/cov_roman_real.npz"
# CoCoA in Limber, and CoCoA at high accuracy
L='{"adopt_limber_gs": 1, "adopt_limber_gg": 1}'
A='{"integration_accuracy": 1, "accuracyboost": 2.0}'
EXPORT_MODE=dv  COCOA_OVR="${L}" ex lsst_y1    fid "${W}/cocoa_lsst_y1_fidlimber.npz"
EXPORT_MODE=cov COCOA_OVR="${L}" ex roman_real fid "${W}/cocoa_roman_real_fidlimber.npz"
EXPORT_MODE=dv  COCOA_OVR="${A}" ex lsst_y1    fid "${W}/cocoa_lsst_y1_fidhi.npz"
EXPORT_MODE=cov COCOA_OVR="${A}" ex roman_real fid "${W}/cocoa_roman_real_fidhi.npz"
# each project at the other project's w (LSST-Y1: -0.9 -> -1; Roman-Real: -1 -> -0.9)
EXPORT_MODE=dv  ex lsst_y1    w_m1  "${W}/cocoa_lsst_y1_w_m1.npz"
EXPORT_MODE=cov ex roman_real w_m09 "${W}/cocoa_roman_real_w_m09.npz"
# Fisher derivatives for bias.py: every lens bias x 1.02, the NLA amplitude
# + 0.1 and its redshift slope + 0.2 (one export each)
python - "${W}" > "${W}/fisher_models.txt" <<'PY'
import json, sys
import numpy as np
for p, pre, mode in (("lsst_y1", "LSST_", "dv"), ("roman_real", "roman_", "cov")):
  meta = json.loads(str(np.load("%s/cocoa_%s_fid.npz" % (sys.argv[1], p))["meta"])); pt = meta["point"]
  for i in range(meta["lens_ntomo"]):
    k = "%sB1_%d" % (pre, i + 1); print(p, mode, "%s=%.8g" % (k, pt[k]*1.02), "b%d" % (i + 1))
  print(p, mode, "%sA1_1=%.8g" % (pre, pt[pre + "A1_1"] + 0.1), "A1")
  print(p, mode, "%sA1_2=%.8g" % (pre, pt[pre + "A1_2"] + 0.2), "eta")
PY
while read -r p mode set tag; do
  EXPORT_MODE=${mode} ex ${p} "${set}" "${W}/cocoa_${p}_d_${tag}.npz"
done < "${W}/fisher_models.txt"
