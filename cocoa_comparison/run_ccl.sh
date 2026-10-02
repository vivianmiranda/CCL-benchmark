#!/bin/bash
# DESC-CCL side of the comparison: every CCL run on the CoCoA exports.
# Usage: CCL_PYTHON=<python of the ccl env> CCL_PR=<pyccl build of PR #1296>
#        COCOA_ROOTDIR=<Cocoa folder> bash run_ccl.sh <run folder>
W="$(cd "${1:?run folder}" && pwd)"
S="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"
cc() { ( cd "${W:?}"; PYTHONPATH="${CCL_PR:?}" COCOA_ROOTDIR="${COCOA_ROOTDIR:?}" \
         "${CCL_PYTHON:?}" "${S:?}/ccl_compute.py" "$@" 2>&1 | grep -E "^CCL|Error" ) ; }
cc cocoa_lsst_y1_w_m1.npz ref ccl_lsst_y1_w_m1_ref.npz
cc cocoa_roman_real_w_m09.npz ref ccl_roman_real_w_m09_ref.npz
for p in lsst_y1 roman_real; do
  for m in fid omm_lo omm_hi ns_lo ns_hi; do cc cocoa_${p}_${m}.npz ref ccl_${p}_${m}_ref.npz; done
  for v in limgs limber norsd rsd_gs flat points separable ccl_numerics ccl_hi eh cocoa_cl bench harmonic; do
    cc cocoa_${p}_fid.npz ${v} ccl_${p}_fid_${v}.npz
  done
done
