#!/bin/bash
# CoCoA side of the comparison: every export the tables and figures need.
# Run from a shell with the cocoa environment active and start_cocoa.sh
# sourced (ROOTDIR set). Usage: bash run_cocoa.sh <run folder>
[ -n "${ROOTDIR:-}" ] || { echo "source start_cocoa.sh first"; return 1 2>/dev/null || exit 1; }
mkdir -p "${1:?run folder}"; W="$(cd "${1}" && pwd)"
S="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"
ex() { python "${S:?}/cocoa_export.py" "$@" > "${W:?}/export_${1}_${2}.log" 2>&1 || echo "FAILED: $*"; }
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
