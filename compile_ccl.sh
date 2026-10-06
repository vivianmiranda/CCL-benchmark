#!/bin/bash
# Build the prepared PR checkout with the existing study's CMake recipe.
# GSL and FFTW come from ccldev; do not download replacement libraries.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Use: source compile_ccl.sh" >&2
  exit 1
fi

(
  cd "$(dirname "${BASH_SOURCE[0]}")" || return 1
  source ./set_installation_options.sh || return 1
  if [[ -z "${CONDA_PREFIX:-}" || ! -x .local/bin/python ||
        -n "${VIRTUAL_ENV:-}" || -n "${ROOTDIR:-}" ]]; then
    echo "Activate only ccldev and run setup_ccl.sh first." >&2
    return 1
  fi
  if [[ "$(git -C "${CCL_PATH}" rev-parse HEAD)" != \
        "${CCL_GIT_COMMIT}" ]]; then
    echo "CCL revision differs from set_installation_options.sh." >&2
    return 1
  fi

  PYTHON="$(pwd -P)/.local/bin/python"
  CCL_PATH="$(cd "${CCL_PATH}" && pwd -P)" || return 1
  export PKG_CONFIG_PATH="${CONDA_PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
  pkg-config --exists gsl fftw3 || return 1

  cmake -S "${CCL_PATH}" -B "${CCL_PATH}/build" \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="${CONDA_PREFIX}" \
    -DPYTHON_EXECUTABLE="${PYTHON}" \
    -DPYTHON_VERSION="$("${PYTHON}" -c \
      'import platform; print(platform.python_version())')" || return 1
  make -C "${CCL_PATH}/build" _ccllib || return 1
  cp "${CCL_PATH}/build/pyccl/_ccllib.so" \
    "${CCL_PATH}/build/pyccl/ccllib.py" "${CCL_PATH}/pyccl/" || return 1

  # Point only this private environment at the built PR. This replaces
  # per-command PYTHONPATH changes; stop_ccl.sh restores the old Python.
  "${PYTHON}" - "${CCL_PATH}" <<'PY' || return 1
from pathlib import Path
import sys
import sysconfig

source = Path(sys.argv[1]).resolve()
site = Path(sysconfig.get_path("purelib"))
(site / "benchmark_ccl.pth").write_text(
    f"import sys; sys.path.insert(0, {str(source)!r})\n"
)
PY
  OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 \
    "${PYTHON}" - "${CCL_PATH}" <<'PY' || return 1
from pathlib import Path
import sys
import pyccl

assert Path(pyccl.__file__).resolve().parent == Path(sys.argv[1]) / "pyccl"
print("CCL:", pyccl.__version__, pyccl.__file__)
PY
  echo "Compilation complete. Next: source start_ccl.sh"
) || return 1
return 0
