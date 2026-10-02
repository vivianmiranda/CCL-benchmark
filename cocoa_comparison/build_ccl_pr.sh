#!/bin/bash
# pyccl of LSSTDESC/CCL PR #1296 (full-sky xi+-, bin-averaged correlations),
# commit 647ad4a, built in place the way CCL's setup.py builds it (Linux and
# macOS). Run it with the ccldev conda environment active (ccldev.yml at the
# repository root: compilers, CMake, SWIG, GSL, FFTW, NumPy, FAST-PT):
#   conda env create --file=ccldev.yml
#   conda activate ccldev
# Usage: bash build_ccl_pr.sh <folder>   (then CCL_PR=<folder>)
set -e
F="${1:?folder for the PR clone}"
if [ ! -d "${F}/.git" ]; then
  git clone https://github.com/DhayaaAnbajagane/CCL.git "${F}"
fi
cd "${F}"
git checkout 647ad4a401f646dec13cebcf32dd0b4d2164eb6b
PY="$(command -v python)"
cmake -H. -Bbuild -DCMAKE_BUILD_TYPE=Release -DPYTHON_EXECUTABLE="${PY}" \
  -DPYTHON_VERSION="$("${PY}" -c 'import sys; print("%d.%d.%d" % sys.version_info[:3])')"
make -Cbuild _ccllib
cp build/pyccl/_ccllib.so build/pyccl/ccllib.py pyccl/
"${PY}" -c "import sys; sys.path.insert(0, '.'); import pyccl; print('pyccl', pyccl.__version__, 'from', pyccl.__file__)"
