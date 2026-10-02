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
# TATT (IA_model 1) from CFASTPT (IA_code 0): the CCL-benchmark scripts' TATT
# point (A1 0.7, eta1 -1.7, A2 -1.36, eta2 -2.5, b_TA 1; pivot z = 0.62) at
# the five cosmologies, plus Limber and high-accuracy fiducial runs
for p in lsst_y1 roman_real; do
  pre=$([ ${p} = lsst_y1 ] && echo LSST || echo roman); mode=$([ ${p} = lsst_y1 ] && echo dv || echo cov)
  T="${pre}_A1_1=0.7,${pre}_A1_2=-1.7,${pre}_A2_1=-1.36,${pre}_A2_2=-2.5,${pre}_BTA_1=1"
  for m in "fid:" "omm_lo:omegam=0.25," "omm_hi:omegam=0.35," "ns_lo:ns=0.92," "ns_hi:ns=1.01,"; do
    EXPORT_MODE=${mode} COCOA_OVR='{"IA_model": 1, "IA_code": 0}' ex ${p} "${m#*:}${T}" "${W}/cocoa_${p}_tatt_${m%%:*}.npz"
  done
  EXPORT_MODE=${mode} COCOA_OVR='{"IA_model": 1, "IA_code": 0, "adopt_limber_gs": 1, "adopt_limber_gg": 1}' ex ${p} "${T}" "${W}/cocoa_${p}_tatt_fidlimber.npz"
  EXPORT_MODE=${mode} COCOA_OVR='{"IA_model": 1, "IA_code": 0, "integration_accuracy": 1, "accuracyboost": 2.0}' ex ${p} "${T}" "${W}/cocoa_${p}_tatt_fidhi.npz"
done
while read -r p mode set tag; do
  EXPORT_MODE=${mode} ex ${p} "${set}" "${W}/cocoa_${p}_d_${tag}.npz"
done < "${W}/fisher_models.txt"
