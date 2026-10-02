#!/bin/bash
# Execution time per evaluation, CoCoA vs DESC-CCL, 8 OpenMP threads, one
# job at a time (run nothing else meanwhile). CoCoA: 20 evaluations per
# case; DESC-CCL (reference sampling): 10; each after two warm-up calls.
# Needs the run folder of run_cocoa.sh (the five cosmologies, NLA and TATT).
# Usage, in a shell with the cocoa environment active and start_cocoa.sh
# sourced: CCL_PYTHON=<ccl env python> CCL_PR=<PR build folder>
#          bash run_timing.sh <run folder>
W="$(cd "${1:?run folder}" && pwd)"
S="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts"
for p in lsst_y1 roman_real; do
  for ia in nla tatt; do
    OMP_NUM_THREADS=8 python "${S:?}/cocoa_timing.py" ${p} ${ia} 4 2>&1 | grep -E "^COCOA|Error"
  done
done
for p in lsst_y1 roman_real; do
  for ia in nla tatt; do
    ( cd "${W:?}"; env -i HOME="${HOME}" PATH=/usr/bin:/bin PYTHONPATH="${CCL_PR:?}" \
        COCOA_ROOTDIR="${ROOTDIR:?}" CMP_WORK="${W}" OMP_NUM_THREADS=8 \
        "${CCL_PYTHON:?}" "${S}/ccl_timing.py" ${p} ${ia} 2 ref 2>&1 | grep -E "^CCL|Error" )
  done
done
