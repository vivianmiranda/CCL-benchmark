#!/bin/bash
# Prepare a private Python layer and the study's pinned CCL PR checkout.
# Source in a fresh terminal after activating the ccldev Conda base.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "Use: source setup_ccl.sh" >&2
  exit 1
fi

(
  cd "$(dirname "${BASH_SOURCE[0]}")" || return 1
  source ./set_installation_options.sh || return 1
  if [[ -z "${CONDA_PREFIX:-}" || -n "${VIRTUAL_ENV:-}" ||
        -n "${ROOTDIR:-}" ]]; then
    echo "Activate only the ccldev Conda base in a fresh shell." >&2
    return 1
  fi
  PYTHON="${CONDA_PREFIX}/bin/python"
  if [[ "$("${PYTHON}" -c \
    'import sys; print("%d.%d" % sys.version_info[:2])')" != \
    "${PYTHON_VERSION}" ]]; then
    echo "The Conda base must use Python ${PYTHON_VERSION}." >&2
    return 1
  fi

  if [[ ! -d "${CCL_PATH}" ]]; then
    mkdir -p "$(dirname "${CCL_PATH}")" || return 1
    git clone "${CCL_GIT_URL}" "${CCL_PATH}" || return 1
    git -C "${CCL_PATH}" checkout "${CCL_GIT_COMMIT}" || return 1
  fi
  if [[ "$(git -C "${CCL_PATH}" rev-parse HEAD)" != \
        "${CCL_GIT_COMMIT}" ]]; then
    echo "CCL revision differs from set_installation_options.sh." >&2
    return 1
  fi

  "${PYTHON}" -m venv --system-site-packages .local || return 1
  echo "Setup complete. Next: source compile_ccl.sh"
) || return 1
return 0
